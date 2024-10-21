# Form for the upload of the files
from django import forms
from .models import SequenceFile

class AlignmentScoreForm(forms.Form):
    match_score = forms.IntegerField(label='Match Score: ')
    mismatch_score = forms.IntegerField(label='Mismatch: ')
    gap_score = forms.IntegerField(label='Gap Score: ')

class SequenceFileForm(forms.Form):
    SEQUENCE_INPUT_CHOICES = [
        ('manual', 'Enter sequences manually'),
        ('fasta', 'Upload FASTA files'),
    ]
    sequence_input_type = forms.ChoiceField(choices=SEQUENCE_INPUT_CHOICES, widget=forms.RadioSelect)
    number_of_sequences = forms.IntegerField(min_value=1, required=False)
    number_of_fasta_files = forms.IntegerField(min_value=1, required=False)

class AlignmentMethodForm(forms.Form):
    ALIGNMENT_CHOICES = [
        ('global', 'Global'),
        ('local', 'Local'),
    ]
    alignment_method = forms.ChoiceField(choices=ALIGNMENT_CHOICES, label="Select Alignment Method")

class DistanceMatrixForm(forms.Form):
    DISTANCE_CHOICES = [
        ('identity', 'Identity'),
        ('blastn', 'BLAST Nucleic Acid'),
        ('trans', 'Transition/Transversion'),
        ('blosum62', 'BLOSUM62'),
        ('pam250', 'PAM250'),
    ]
    distance_method = forms.ChoiceField(choices=DISTANCE_CHOICES)

class TreeConstructionForm(forms.Form):
    TREE_CHOICES = [
        ('upgma', 'UPGMA'),
        ('nj', 'Neighbor-Joining'),
    ]
    tree_method = forms.ChoiceField(choices=TREE_CHOICES)