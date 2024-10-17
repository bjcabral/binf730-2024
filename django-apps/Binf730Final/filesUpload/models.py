
# Model for the upload of two sequence files
import os
from django.db import models

class SequenceFile(models.Model):
    file1 = models.FileField(upload_to='sequences/')
    file2 = models.FileField(upload_to='sequences/')
