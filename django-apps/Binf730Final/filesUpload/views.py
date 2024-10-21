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
            sequence_input_type = request.POST.get('sequence_input_type')
            sequences = []

            if sequence_input_type == 'manual':
                number_of_sequences = int(request.POST.get('number_of_sequences'))
                for i in range(number_of_sequences):
                    sequence = request.POST.get(f'sequence_{i + 1}')
                    sequences.append(sequence)

            elif sequence_input_type == 'fasta':
                number_of_fasta_files = int(request.POST.get('number_of_fasta_files'))
                for i in range(number_of_fasta_files):
                    fasta_file = request.FILES.get(f'fasta_file_{i + 1}')
                    fs = FileSystemStorage()
                    filename = fs.save(fasta_file.name, fasta_file)
                    sequences.append(SeqIO.parse(fs.path(filename), "fasta"))

            request.session['sequences'] = sequences
            return redirect('method_and_score_scheme')
        else:
            return HttpResponse('Form is not valid.')
    else:
        form = SequenceFileForm()
    return render(request, 'upload.html', {'form': form})

def method_and_score_scheme(request):
    sequences = request.session.get('sequences')

    if not sequences:
        return HttpResponse("Sequences not found. Please upload again.")

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

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os

def align_sequences(request):
    sequences = request.session.get('sequences')
    match_score = request.session.get('match_score')
    mismatch_score = request.session.get('mismatch_score')
    gap_score = request.session.get('gap_score')
    alignment_method = request.session.get('alignment_method')

    if not sequences:
        return HttpResponse("Sequences not found. Please upload again.")

    try:
        aligner = Align.PairwiseAligner()
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.query_gap_score = gap_score
        aligner.target_gap_score = gap_score
        aligner.mode = alignment_method

        aligned_seq_html = "<h2>Pairwise Alignments:</h2>"
        all_alignments = []

        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                seq1 = sequences[i]
                seq2 = sequences[j]

                # Perform the alignment
                alignments = aligner.align(seq1, seq2)

                # Get the best alignment (first one in list)
                best_alignment = alignments[0]

                # Add the alignment to the HTML response
                aligned_seq_html += f"<h3>Alignment between Sequence {i+1} and Sequence {j+1}:</h3>"
                aligned_seq_html += f"<pre>{best_alignment}</pre>"

                # Store the alignment for writing to file
                all_alignments.append((f"Seq{i+1}", str(best_alignment[0])))
                all_alignments.append((f"Seq{j+1}", str(best_alignment[1])))

        # Create a MultipleSeqAlignment object
        msa = MultipleSeqAlignment([
            SeqRecord(Seq(seq), id=name) for name, seq in all_alignments
        ])

        # Specify the directory and filename
        directory = '/home/ubuntu/binf-730/django-apps/Binf730Final/media'
        output_file = f'{directory}/aligned.aln'

        # Ensure the directory exists
        os.makedirs(directory, exist_ok=True)

        # Write alignment to Clustal format
        AlignIO.write(msa, output_file, 'clustal')

        # Check if file was created and add to the response
        if os.path.exists(output_file):
            aligned_seq_html += f"<p>Alignment file has been saved to: {output_file}</p>"
        else:
            aligned_seq_html += "<p>Failed to save alignment file.</p>"

        return HttpResponse(aligned_seq_html)

    except Exception as e:
        return HttpResponse(f"Error occurred during alignment: {e}")