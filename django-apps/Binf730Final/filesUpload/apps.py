# File for apps label definitions
from django.apps import AppConfig

class FilesUploadConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'filesUpload'

    def ready(self):
        import filesUpload.templatetags.custom_filters