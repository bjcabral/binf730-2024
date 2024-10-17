# View file for the upload form
from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage
from Bio import SeqIO
from Bio import Align
from .forms import SequenceFileForm, AlignmentScoreForm, AlignmentMethodForm


def upload_sequence_files(request):
    if request.method == 'POST':
        form = SequenceFileForm(request.POST, request.FILES)
        if form.is_valid():
            try:
                file1 = request.FILES['file1']
                file2 = request.FILES['file2']
                fs = FileSystemStorage()
                filename1 = fs.save(file1.name, file1)
                filename2 = fs.save(file2.name, file2)
                SeqIO.parse(fs.path(filename1), "fasta")
                SeqIO.parse(fs.path(filename2), "fasta")
                request.session['file1_name'] = filename1
                request.session['file2_name'] = filename2
                return redirect('method_and_score_scheme')
            except Exception as e:
                return HttpResponse(f"Error: {e}")
    else:
        form = SequenceFileForm()
    return render(request, 'upload.html', {'form': form})


def method_and_score_scheme(request):
    file1_name = request.session.get('file1_name')
    file2_name = request.session.get('file2_name')

    if not file1_name or not file2_name:
        return HttpResponse("Files not found. Please upload again.")

    if request.method == 'POST':
        alignment_score_form = AlignmentScoreForm(request.POST)
        alignment_method_form = AlignmentMethodForm(request.POST)
        if alignment_score_form.is_valid() and alignment_method_form.is_valid():
            match_score = alignment_score_form.cleaned_data['match_score']
            mismatch_score = alignment_score_form.cleaned_data['mismatch_score']
            gap_score = alignment_score_form.cleaned_data['gap_score']
            alignment_method = alignment_method_form.cleaned_data['alignment_method']

            request.session['match_score'] = match_score
            request.session['mismatch_score'] = mismatch_score
            request.session['gap_score'] = gap_score
            request.session['alignment_method'] = alignment_method

            return redirect('align_sequences')
    else:
        alignment_score_form = AlignmentScoreForm()
        alignment_method_form = AlignmentMethodForm()
    return render(request, 'method_and_score_scheme.html',
                  {'score_form': alignment_score_form, 'method_form': alignment_method_form})


def align_sequences(request):
    file1_name = request.session.get('file1_name')
    file2_name = request.session.get('file2_name')
    match_score = request.session.get('match_score')
    mismatch_score = request.session.get('mismatch_score')
    gap_score = request.session.get('gap_score')
    # alignment_method = request.session.get('alignment_method')
    alignment_method = 'global'

    if not file1_name or not file2_name:
        return HttpResponse("Files not found. Please upload again.")

    fs = FileSystemStorage()
    seq1 = SeqIO.read(fs.path(file1_name), 'fasta')
    seq2 = SeqIO.read(fs.path(file2_name), 'fasta')

    try:
        aligner = Align.PairwiseAligner()
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.query_gap_score = gap_score
        aligner.target_gap_score = gap_score
        aligner.mode = alignment_method  # Use the alignment method from the session
        # Perform the alignment
        alignments = aligner.align(seq1.seq, seq2.seq)

        # Get the best alignment (first one in list)
        best_alignment = alignments[0]

        # Prepare response with aligned sequences formatted as HTML
        aligned_seq_html = f"<pre>{best_alignment}</pre>"

        return HttpResponse(aligned_seq_html)

    except Exception as e:
        return HttpResponse(f"Error occurred during alignment: {e}")