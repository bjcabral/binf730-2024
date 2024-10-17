# View file for the upload form
from django.shortcuts import render, redirect
from .forms import SequenceFileForm, AlignmentScoreForm
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

                form.save()
                # return redirect('success')
                return redirect('align_sequences', file1_name=filename1, file2_name=filename2)
            except FileNotFoundError:
                print(f"File not found: {file1.name} or {file2.name}")
            except ValueError as e:
                print(f"Error parsing FASTA file: {e}")
                # files are not FASTA files
                return HttpResponse("Files are not in FASTA format.")
    else:
        form = SequenceFileForm()
    return render(request, 'upload.html', {'form': form})

def align_sequences(request, file1_name, file2_name):
    if request.method == 'POST':
        form = AlignmentScoreForm(request.POST)
        if form.is_valid():
            match_score = form.cleaned_data['match_score']
            mismatch_score = form.cleaned_data['mismatch_score']
            gap_score = form.cleaned_data['gap_score']

            # Load sequences for alignment
            fs = FileSystemStorage()
            seq1 = SeqIO.read(fs.path(file1_name), "fasta")
            seq2 = SeqIO.read(fs.path(file2_name), "fasta")

            # Instantiate PairwiseAligner()
            aligner = Align.PairwiseAligner()

            # Pass parameters to aligner
            aligner.match_score = match_score
            aligner.mismatch_score = mismatch_score
            aligner.query_gap_score = gap_score
            aligner.target_gap_score = gap_score

            alignments = aligner.align(seq1.seq, seq2.seq)

            # Print the best alignment
            best_alignment = alignments[0]
            return HttpResponse(str(best_alignment))
    else:
        # Handle GET request
        form = AlignmentScoreForm()
        return render(request, 'align.html', {'form': form, 'file1_name': file1_name, 'file2_name': file2_name})
        

