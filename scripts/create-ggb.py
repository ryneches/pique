"""
Demonstrates how to programmatically create a Gaggle Genome Browser data file.

GGB files are SQLite databases with a schema described here:
http://gaggle.systemsbiology.net/docs/geese/genomebrowser/dev/schema/

Run it like so:
python create-ggb.py [output-file.hbgb]

Christopher Bare
Nov. 2011
"""
import sqlite3
import uuid
import argparse
import math
import re
import random
from datetime import datetime
import urllib
from contextlib import closing


# halo genome information
halo_genome = {
        "species":"Halobacterium salinarum NRC-1",
        "attributes": {
            "ucsc.db.name":"haloHalo1",
            "domain":"archaea"
        },
        "sequences":[
            {"name":"chromosome", "length":2014239, "topology":"circular", "ncbi_uid":165, "id":1},
            {"name":"PNRC200", "length":365425, "topology":"circular", "ncbi_uid":166, "id":2},
            {"name":"PNRC100", "length":191346, "topology":"circular", "ncbi_uid":13234, "id":3} ]
    }

def set_attribute(cursor, object_uuid, key, value):
    """
    Set the key/value attribute on the object identified by the uuid. For example,
    setting the 'color' attribute on a track to blue would be done like this:
    set_attribute(cursor, track_uuid, 'color', '0x0000FF')
    """
    if type(object_uuid) is not str:
        object_uuid = str(object_uuid)
    cursor.execute("""
        INSERT INTO attributes (uuid, key, value) VALUES (?, ?, ?);
    """, (object_uuid, key, value,))

def create_schema(filename):
    """
    Creates the necessary schema for a GGB data file.
    """
    conn = None
    c = None
    try:
        conn = sqlite3.connect(filename)
        with conn:
            c = conn.cursor()
            c.execute("""
                CREATE TABLE datasets (
                  uuid text PRIMARY KEY NOT NULL,
                  name text);""")
    
            c.execute("""
                CREATE TABLE datasets_sequences (
                  datasets_uuid text not null,
                  sequences_uuid text not null );""")
            c.execute("""
                CREATE TABLE sequences (
                  id integer primary key AUTOINCREMENT not null,
                  uuid text not null,
                  name text not null,
                  length integer not null,
                  topology text);""")
    
            c.execute("""
                CREATE TABLE datasets_tracks (
                  datasets_uuid text not null,
                  tracks_uuid text not null);""")
            c.execute("""
                CREATE TABLE tracks (
                  uuid text primary key not null,
                  name text not null,
                  type text not null,
                  table_name text not null);""")
    
            c.execute("""
                CREATE TABLE attributes (
                  uuid NOT NULL,
                  key text NOT NULL,
                  value);""")
    except:
        raise

def create_dataset(filename, genome):
    """
    Creates a new dataset and adds genome data into GGB data file.
    Returns the uuid of the new dataset.
    """
    conn = None
    c = None
    try:
        conn = sqlite3.connect(filename)
        with conn:
            c = conn.cursor()

            # create dataset and assign it a uuid
            dataset_uuid = uuid.uuid1()
            dataset_name = 'Test GGB dataset'
            c.execute("""
                INSERT INTO datasets (uuid, name) VALUES (?, ?);
            """, (str(dataset_uuid), dataset_name,))
            
            print "-- creating dataset uuid " + str(dataset_uuid)
        
            # add species as an attribute of the dataset
            set_attribute(c, dataset_uuid, 'species', genome['species'])
        
            # add some optional attributes of the dataset
            set_attribute(c, dataset_uuid, 'created-by', 'script create-gg.py')
            set_attribute(c, dataset_uuid, 'created-at', datetime.now().isoformat())
        
            # add more optional attributes
            if genome['attributes']:
                for key in genome['attributes']:
                    set_attribute(c, dataset_uuid, key, genome['attributes'][key])
        
            for sequence in genome['sequences']:
                sequence['uuid'] = uuid.uuid1()
                c.execute("""
                    INSERT INTO sequences (uuid, name, length, topology) VALUES (?, ?, ?, ?);
                """, (str(sequence['uuid']), sequence['name'], sequence['length'], sequence['topology'],))
                sequence['id'] = c.lastrowid
                c.execute("""
                    INSERT INTO datasets_sequences (datasets_uuid, sequences_uuid) VALUES (?, ?);
                """, (str(dataset_uuid), str(sequence['uuid']),))
                print "-- inserted sequence: %s id: %d." % (sequence['name'], sequence['id'],)
        
            conn.commit()
    
    except:
        raise

    return dataset_uuid


def add_track(filename, dataset_uuid, name, type, features, table_name=None, attributes=None):
    """
    Creates a track, adding an entry into the tracks table, associating the new track
    with the dataset and creating the features table.

    Parameters
    filename: name of SQLite db (aka GGB data file)
    dataset_uuid: identifies the dataset
    name: name of track
    type: type of track (gene, quantitative.segment, quantitative.positional)
    features: list of features
    
    Features are tuples containing the data fields required by the specified
    track type:
    gene => (sequences_id, strand, start, end, name, common_name, gene_type)
    quantitative.segment => (sequences_id, strand, start, end, value)
    quantitative.positional => (sequences_id, strand, position, value)

    The steps here are:
    1) create temp table with proper schema for the type of track we're adding
    2) insert features
    3) create feature table
    4) copy features, sorting by sequence, strand, start, end or position
       Sorting is important, because the rendering code relies on the features
       being in the table in this order. This ends up working like a poor man's
       clustered index on feature location.
    5) clean up temp table
    6) add entry in tracks and datasets_tracks tables
    7) insert attributes for track

    Returns the UUID assigned to the track.
    """
    
    # create a feature table_name if none given
    if table_name is None:
        table_name = re.sub(r'_{2,}', '_', "features_" + re.sub(r'\W', '_', name))
        table_name = table_name.rstrip('_')
    
    conn = None
    c = None
    try:
        conn = sqlite3.connect(filename)
        with conn:
            c = conn.cursor()
            
            # uniquify feature table name
            c.execute("""
                SELECT count(*) FROM tracks where table_name = ?;
            """, (table_name,))
            if c.fetchone()[0] > 0:
                c.execute("""
                    SELECT table_name FROM tracks where table_name like ?;
                """, (table_name + '%',))
                max_suffix = 0
                for row in c.fetchone():
                    m = re.match(table_name+"_(\\d+)",)
                    if m:
                        max_suffix = max(int(m.group(1)), max_suffix)
                table_name = table_name + "_" + str(max_suffix+1)
            
            if type=='gene':
                c.execute("""
                    CREATE TABLE temp (
                        sequences_id integer NOT NULL,
                        strand text NOT NULL,
                        start integer NOT NULL,
                        end integer NOT NULL,
                        name text,
                        common_name text,
                        gene_type text);
                """)
                for feature in features:
                    c.execute("""
                        INSERT INTO temp (sequences_id, strand, start, end, name, common_name, gene_type)
                        VALUES (?, ?, ?, ?, ?, ?, ?);
                    """, feature)
                c.execute("""
                    CREATE TABLE %s (
                        sequences_id integer NOT NULL,
                        strand text NOT NULL,
                        start integer NOT NULL,
                        end integer NOT NULL,
                        name text,
                        common_name text,
                        gene_type text);
                """ % (table_name,))
            
                # Copy features into feature table sorting along the way.
                c.execute("""
                    INSERT INTO %s
                    SELECT * from temp order by sequences_id, strand, start, end;
                """ % (table_name,))
                
                # default track attributes
                if attributes is None:
                    attributes = {'top':0.46, 'height':0.08, 'viewer':'Gene'}
            
            elif type=='quantitative.segment':
                c.execute("""
                    CREATE TABLE temp (
                        sequences_id integer NOT NULL,
                        strand text NOT NULL,
                        start integer NOT NULL,
                        end integer NOT NULL,
                        value numeric);
                """)
                for feature in features:
                    c.execute("""
                        INSERT INTO temp (sequences_id, strand, start, end, value)
                        VALUES (?, ?, ?, ?, ?);
                    """, feature)
                c.execute("""
                    CREATE TABLE %s (
                        sequences_id integer NOT NULL,
                        strand text NOT NULL,
                        start integer NOT NULL,
                        end integer NOT NULL,
                        value numeric);
                """ % (table_name,))
        
                # Copy features into feature table sorting along the way.
                c.execute("""
                    INSERT INTO %s
                    SELECT * from temp order by sequences_id, strand, start, end;
                """ % (table_name,))
                
                # default track attributes
                if attributes is None:
                    attributes = {'top':0.10, 'height':0.10, 'color':'0x60336699', 'viewer':'Scaling',
                        'max.value':1.0, 'min.value':-1.0, 'rangeMax':1.0, 'rangeMin':-1.0}
            
            elif type=='quantitative.positional':
                c.execute("""
                    CREATE TABLE temp (
                        sequences_id integer NOT NULL,
                        strand text NOT NULL,
                        position integer NOT NULL,
                        value numeric);
                """)
                for feature in features:
                    c.execute("""
                        INSERT INTO temp (sequences_id, strand, position, value)
                        VALUES (?, ?, ?, ?);
                    """, feature)
                c.execute("""
                    CREATE TABLE %s (
                        sequences_id integer NOT NULL,
                        strand text NOT NULL,
                        position integer NOT NULL,
                        value numeric);
                """ % (table_name,))
        
                # Copy features into feature table sorting along the way.
                c.execute("""
                    INSERT INTO %s
                    SELECT * from temp order by sequences_id, strand, position;
                """ % (table_name,))
                
                # default track attributes
                if attributes is None:
                    attributes = {'top':0.10, 'height':0.10, 'color':'0x80336699', 'viewer':'Scaling'}
            
            else:
                raise Exception("Unknown track type: " + type)
       
            c.execute("""DROP TABLE temp;""")
            c.execute("""vacuum;""")
            
            # generate a UUID for the track
            track_uuid = uuid.uuid1();
            
            c.execute("""
                INSERT INTO tracks (uuid, name, type, table_name) VALUES (?,?,?,?);
            """, (str(track_uuid), name, type, table_name,))
            c.execute("""
                INSERT INTO datasets_tracks (datasets_uuid, tracks_uuid) VALUES (?,?);
            """, (str(dataset_uuid), str(track_uuid),))
            
            if attributes:
                for key in attributes:
                    set_attribute(c, track_uuid, key, attributes[key])
        
            conn.commit()
            
            print "-- created track: %s, uuid: %s, features table: %s." % (name, str(track_uuid), table_name,)
            
    except:
        raise

    return track_uuid

def fake_track_data(genome):
    """
    A generator for fake track data.
    Yields features, which are tuples of (sequence_id, strand, start, end, value)
    """
    for sequence in genome['sequences']:
        for strand in ['+','-']:
            offset = 0 if strand=='+' else 20
            for i in range(1+offset, sequence['length'], 40):
                yield ( sequence['id'], strand, i, i+59, math.sin( i/(1000.0*math.pi)), )

def random_positional_track_data(genome):
    """
    A generator for random position track data.
    Yields features, which are tuples of (sequence_id, strand, position, value)
    """
    for sequence in genome['sequences']:
        for strand in ['+','-']:
            for i in range(1,sequence['length']/5000):
                yield ( sequence['id'], strand, random.randint(1,sequence['length']), random.lognormvariate(1, 1) )

def get_genes_from_ncbi(genome):
    """
    A generator that downloads protein coding gene annotations from NCBI Entrez
    and yields feature suitable for adding to a gene track. This might work for
    other refseq organisms, but is not tested for any besides halo.
    
    Note that this function makes use of ncbi_uids stuck in the sequences of the
    genome object. These uids can be found, at least for prokaryotic organisms by
    following the refseq links at this URL
    http://www.ncbi.nlm.nih.gov/genomes/lproks.cgi
    You'll get to a page containing a table with columns for Genome Info, Features, 
    BLAST homologs, Links and Review Info. Under features, click on "Protein coding"
    and the uid will be in the URL. Why it's so awkward to get a gene list out of
    NCBI, I don't know. Maybe there's an easier way that I'm too clueless to find.
    """
    for sequence in genome['sequences']:
        print "-- fetching genes for %s from NCBI. uid=%d" % (sequence['name'], sequence['ncbi_uid'],)
        url = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=genome&cmd=file&dopt=Protein+Table&list_uids=%d" % (sequence['ncbi_uid'],)
        with closing(urllib.urlopen(url)) as f:
            description = f.next()
            column_headers = f.next()
            for line in f:
                
                # skip blank lines
                if len(line.strip()) == 0:
                    continue
                
                fields = line.rstrip("\n").split("\t")
                strand = fields[3]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[8]
                common_name = fields[7] if fields[7]!='-' else None
                gene_type = "cds" # coding sequence
                yield (sequence['id'], strand, start, end, name, common_name, gene_type,)

def main():
    parser = argparse.ArgumentParser(description='Create a Gaggle Genome Browser data file.')
    parser.add_argument('filename', help='name of GGB file to create')
    args = parser.parse_args()
    
    create_schema(args.filename)
    dataset_uuid = create_dataset(args.filename, halo_genome)
    add_track(args.filename, dataset_uuid, "Genome", "gene", get_genes_from_ncbi(halo_genome))
    add_track(args.filename, dataset_uuid, "Sine wave", "quantitative.segment", fake_track_data(halo_genome))
    add_track(args.filename, dataset_uuid, "Random markers", "quantitative.positional", random_positional_track_data(halo_genome), 
        attributes={'top':0.30, 'height':0, 'color':'0x80FF0000', 'viewer':'Triangle marker', 'triangle.size':25})

if __name__ == "__main__":
    main()
