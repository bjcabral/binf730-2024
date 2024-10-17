# URL configuration for Binf730Final project.
from . import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import path
from filesUpload import views

urlpatterns = [
    # ex: /upload/
    path('upload/', views.upload_sequence_files, name='upload_sequence_files'),
    # ex: /align/
    path('align/', views.align_sequences, name='align_sequences'),
    # ex: /select_alignment_method/
    path('select_alignment_method/', views.select_alignment_method, name='select_alignment_method'),
]
