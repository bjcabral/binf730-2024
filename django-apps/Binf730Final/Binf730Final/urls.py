# URL configuration for Binf730Final project.
from django.urls import path
from filesUpload import views

urlpatterns = [
    path('upload/', views.upload_sequence_files, name='upload_sequence_files'),
    path('method-and-score/', views.method_and_score_scheme, name='method_and_score_scheme'),
    path('align/', views.align_sequences, name='align_sequences'),
    path('calculate-distance/', views.calculate_distance, name='calculate_distance'),
    path('construct-tree/', views.construct_tree, name='construct_tree'),
]

