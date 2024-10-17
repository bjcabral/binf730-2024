# URL configuration for Binf730Final project.
from django.urls import path
from . import views

urlpatterns = [
    path('upload/', views.upload_sequence_files, name='upload_sequence_files'),
    path('method-and-score-scheme/', views.method_and_score_scheme, name='method_and_score_scheme'),
    path('align/', views.align_sequences, name='align_sequences'),
]