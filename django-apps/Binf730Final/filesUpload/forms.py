# Form for the upload of the files
from django import forms
from .models import SequenceFile

class AlignmentScoreForm(forms.Form):
    match_score = forms.IntegerField(label='Match Score: ')
    mismatch_score = forms.IntegerField(label='Mismatch: ')
    gap_score = forms.IntegerField(label='Gap Score: ')

from django import forms

class SequenceFileForm(forms.Form):
    SEQUENCE_INPUT_CHOICES = [
        ('manual', 'Manual Input'),
        ('fasta', 'FASTA File Upload'),
    ]
    sequence_input_type = forms.ChoiceField(choices=SEQUENCE_INPUT_CHOICES)
    number_of_sequences = forms.IntegerField(min_value=1, required=False)
    number_of_fasta_files = forms.IntegerField(min_value=1, required=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.data.get('sequence_input_type') == 'manual':
            number_of_sequences = int(self.data.get('number_of_sequences', 0))
            for i in range(number_of_sequences):
                self.fields[f'sequence_id_{i + 1}'] = forms.CharField(label=f'Sequence ID {i + 1}')
                self.fields[f'sequence_{i + 1}'] = forms.CharField(label=f'Sequence {i + 1}', widget=forms.Textarea)
        elif self.data.get('sequence_input_type') == 'fasta':
            number_of_fasta_files = int(self.data.get('number_of_fasta_files', 0))
            for i in range(number_of_fasta_files):
                self.fields[f'fasta_file_{i + 1}'] = forms.FileField(label=f'FASTA File {i + 1}')

    def clean(self):
        cleaned_data = super().clean()
        sequence_input_type = cleaned_data.get('sequence_input_type')
        if sequence_input_type == 'manual' and not cleaned_data.get('number_of_sequences'):
            raise forms.ValidationError("Please specify the number of sequences for manual input.")
        elif sequence_input_type == 'fasta' and not cleaned_data.get('number_of_fasta_files'):
            raise forms.ValidationError("Please specify the number of FASTA files to upload.")
        return cleaned_data

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