# View file for the upload form
from django.shortcuts import render, redirect
from .forms import SequenceFileForm, AlignmentScoreForm, AlignmentMethodForm
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage
from Bio import SeqIO
from Bio import Align
from Bio.Align import Applications
from Bio.Align import AlignInfo
from io import StringIO


def success(request):
    return HttpResponse("File uploaded successfully.")

def upload_sequence_files(request):
    if request.method == 'POST':
        form = SequenceFileForm(request.POST, request.FILES)
        if form.is_valid():
            try:
                # Access file names
                file1 = request.FILES['file1']
                file2 = request.FILES['file2']

                # Save files to a temporary location
                fs = FileSystemStorage()
                filename1 = fs.save(file1.name, file1)
                filename2 = fs.save(file2.name, file2)

                # Parse file contents
                SeqIO.parse(fs.path(filename1), "fasta")
                SeqIO.parse(fs.path(filename2), "fasta")

                # Store file names in the session
                request.session['file1_name'] = filename1
                request.session['file2_name'] = filename2

                # Redirect to align_sequences view
                return redirect('align_sequences')
            except FileNotFoundError:
                print(f"File not found: {file1.name} or {file2.name}")
            except ValueError as e:
                print(f"Error parsing FASTA file: {e}")
                return HttpResponse("Files are not in FASTA format.")
    else:
        form = SequenceFileForm()
    return render(request, 'upload.html', {'form': form})

def align_sequences(request):
    file1_name = request.session.get('file1_name')
    file2_name = request.session.get('file2_name')

    if request.method == 'POST':
        alignment_score_form = AlignmentScoreForm(request.POST)
        if alignment_score_form.is_valid():
            match_score = alignment_score_form.cleaned_data['match_score']
            mismatch_score = alignment_score_form.cleaned_data['mismatch_score']
            gap_score = alignment_score_form.cleaned_data['gap_score']

            # Store alignment scores in the session
            request.session['match_score'] = match_score
            request.session['mismatch_score'] = mismatch_score
            request.session['gap_score'] = gap_score

            # Render the alignment method form
            alignment_method_form = AlignmentMethodForm()
            return render(request, 'alignment_method.html', {'form': alignment_method_form, 'file1_name': file1_name, 'file2_name': file2_name})
        else:
            return render(request, 'alignment_scores.html', {'form': alignment_score_form, 'file1_name': file1_name, 'file2_name': file2_name})
    else:
        # Handle GET requests by showing the alignment scores form
        alignment_score_form = AlignmentScoreForm()
        return render(request, 'alignment_scores.html', {'form': alignment_score_form, 'file1_name': file1_name, 'file2_name': file2_name})

def select_alignment_method(request):
    file1_name = request.session.get('file1_name')
    file2_name = request.session.get('file2_name')

    if request.method == 'POST':
        alignment_method_form = AlignmentMethodForm(request.POST)
        if alignment_method_form.is_valid():
            print("WE MADE IT HERE!")
            alignment_method = alignment_method_form.cleaned_data['alignment_method']
            print(f"The alignment method is {alignment_method}")
            # Load sequences for alignment
            fs = FileSystemStorage()
            seq1 = SeqIO.read(fs.path(file1_name), 'fasta')
            seq2 = SeqIO.read(fs.path(file2_name), 'fasta')

            # Get alignment scores from the session
            match_score = request.session.get('match_score')
            mismatch_score = request.session.get('mismatch_score')
            gap_score = request.session.get('gap_score')
            print(f"The alignment scores are: {match_score}, {mismatch_score}, {gap_score}")
            
            try:
            # Instantiate PairwiseAligner() based on the selected alignment method
                aligner = Align.PairwiseAligner()

            # Pass parameters to aligner
                aligner.match_score = match_score
                aligner.mismatch_score = mismatch_score
                aligner.query_gap_score = gap_score
                aligner.target_gap_score = gap_score

            # Perform the alignment
                alignments = aligner.align(seq1.seq, seq2.seq)
            except Exceptino as e:
                print("Error occured during alignment: {e}")
                return HttpResponse(f"Error occured during alignment: {e}")
            # Get the best alignment
            best_alignment = alignments

            # Prepare response
            alignment_result = str(best_alignment)
            return HttpResponse(alignment_result)
        else:
            return render(request, 'alignment_method.html', {'form': alignment_method_form, 'file1_name': file1_name, 'file2_name': file2_name})
    else:
        # Handle GET requests by showing the alignment method form
        alignment_method_form = AlignmentMethodForm()
        return render(request, 'alignment_method.html', {'form': alignment_method_form, 'file1_name': file1_name, 'file2_name': file2_name})
