# Form for the upload of the files
from django import forms
from .models import SequenceFile

class AlignmentScoreForm(forms.Form):
    match_score = forms.IntegerField(label='Match Score: ')
    mismatch_score = forms.IntegerField(label='Mismatch: ')
    gap_score = forms.IntegerField(label='Gap Score: ')

class SequenceFileForm(forms.ModelForm):
    class Meta:
        model = SequenceFile
        fields = ('file1', 'file2')

class AlignmentMethodForm(forms.Form):
    ALIGNMENT_CHOICES = [
        ('manual', 'Manual'),
        ('clustalw', 'ClustalW'),
    ]
    alignment_method = forms.ChoiceField(choices=ALIGNMENT_CHOICES, label="Select Alignment Method")
