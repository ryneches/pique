from django.db import models
from django.core.files.storage import FileSystemStorage
from django.contrib.auth.models import User

attachment_file_storage = FileSystemStorage( location='/tmp/pique_app', base_url='/attachments' )

class Experiment( models.Model ) :
    name    = models.CharField( max_length=265, editable=False )
    created = models.DateTimeField( auto_now_add=True, editable=False )
    updated = models.DateTimeField( auto_now=True, auto_now_add=True, editable=False )
    IP      = models.FileField( upload_to='attachments', storage=attachment_file_storage )
    BG      = models.FileField( upload_to='attachments', storage=attachment_file_storage )
    mapfile = models.FileField( upload_to='attachments', storage=attachment_file_storage )
    alpha   = models.IntegerField( default=32 )
    frag    = models.IntegerField( default=300 )
    owner   = models.
