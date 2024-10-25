# File for apps label definitions
from django.apps import AppConfig

class FileUploadsConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'fileUploads'

    def ready(self):
        import fileUploads.templatetags.custom_filters