# View file for the upload form
import os
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage
from Bio import SeqIO, AlignIO, Phylo
from Bio import Align
from .forms import SequenceFileForm, AlignmentScoreForm, AlignmentMethodForm, DistanceMatrixForm, TreeConstructionForm
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO


# Specify the directory and filename
output_directory = '/home/ubuntu/binf-730/django-apps/Binf730Final/media'
aligned_file = "{output_directory}/aligned.aln"
# Save the files in fasta format in one single file when enter manually or when fasta files are uploaded
output_fasta_file = "{output_directory}/uploaded_sequences.fasta"
os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)

# Append the entered sequences, manual or fasta file into one single file for the aligner.
def append_to_fasta_file(sequence_id, sequence):
    with open(output_fasta_file, 'a') as fasta_file:
        fasta_file.write(f'>{sequence_id}\n{sequence}\n')

# function use to get the input data (sequences) from the users, either manual sequences, or fasta formatted files
# the function prompts the user for how they want to input the data. Manual means they will enter short sequences
# manually, while fasta means they will input the path to a faste formatted sequence file. First, the user will
# indicate how many sequences they will enter and the code will iterate accordingly.
def upload_sequence_files(request):
    if request.method == 'POST':
        form = SequenceFileForm(request.POST, request.FILES)
        if form.is_valid():
            sequence_input_type = request.POST.get('sequence_input_type')

            # Clear the existing file content
            open(output_fasta_file, 'w').close()

            if sequence_input_type == 'manual':
                number_of_sequences = int(request.POST.get('number_of_sequences'))
                for i in range(number_of_sequences):
                    sequence_id = request.POST.get(f'sequence_id_{i + 1}')
                    sequence = request.POST.get(f'sequence_{i + 1}')
                    append_to_fasta_file(sequence_id, sequence)

            elif sequence_input_type == 'fasta':
                number_of_fasta_files = int(request.POST.get('number_of_fasta_files'))
                fs = FileSystemStorage()
                for i in range(number_of_fasta_files):
                    fasta_file = request.FILES.get(f'fasta_file_{i + 1}')
                    filename = fs.save(fasta_file.name, fasta_file)
                    for record in SeqIO.parse(fs.path(filename), "fasta"):
                        append_to_fasta_file(record.id, str(record.seq))
                    fs.delete(filename)  # Delete the temporary file

            request.session['fasta_file'] = output_fasta_file
            return redirect('method_and_score_scheme')
        else:
            return HttpResponse('Form is not valid.')
    else:
        form = SequenceFileForm()
    return render(request, 'upload.html', {'form': form})

# This function prompts the user for the the alignment method (global or local) and the score matrix.
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

# These methods uses the Align() function from BioPython to execute the sequence alignment and save the aligned
# sequence file for later use with the cal_distance() function.
def align_sequences(request):
    aligned_file_name = '/home/ubuntu/binf-730/media/aligned.aln'
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

        # Write alignment to Clustal format
        AlignIO.write(msa, aligned_file_name, 'clustal')

        # Check if file was created and add to the response
        if os.path.exists(aligned_file_name):
            aligned_seq_html += f"<p>Alignment file has been saved to: {aligned_file_name}</p>"
        else:
            aligned_seq_html += "<p>Failed to save alignment file.</p>"

        return redirect('calculate_distance')

    except Exception as e:
        return HttpResponse(f"Error occurred during alignment: {e}")

# Function to calculate the distance matrix using the user selected method. BLOSUM250, identity, etc.
# this function collects the input from the user and uses read_aligned() and cal_dist() to read
# the aligned files in clustalo format and calculate the distance matrix.
def calculate_distance(request):
    if request.method == 'POST':
        form = DistanceMatrixForm(request.POST)
        if form.is_valid():
            distance_method = form.cleaned_data['distance_method']
            aligned_seqs = read_aligned(aligned_file_name)
            distance_matrix = cal_dist(aligned_seqs, distance_method)
            request.session['distance_matrix'] = distance_matrix.tolist()  # Convert to list for JSON serialization
            return redirect('construct_tree')
    else:
        form = DistanceMatrixForm()
    return render(request, 'calculate_distance.html', {'form': form})

# Function to construct the tree from the distance matrix using the algorithm selected by the user
# such as Neighbor Joining (nj) or UPGMA.
def construct_tree(request):
    if request.method == 'POST':
        form = TreeConstructionForm(request.POST)
        if form.is_valid():
            tree_method = form.cleaned_data['tree_method']
            distance_matrix = request.session.get('distance_matrix')
            if not distance_matrix:
                return HttpResponse("Distance matrix not found. Please calculate distances first.")
            tree = construct_tree_from_matrix(distance_matrix, tree_method)
            tree_image = draw_tree(tree)
            return render(request, 'display_tree.html', {'tree_image': tree_image})
    else:
        form = TreeConstructionForm()
    return render(request, 'construct_tree.html', {'form': form})

# I use this function to read the aligned files from disk after it is saved by the align_sequences() function
def read_aligned(aligned_seqs_file):
    return AlignIO.read(aligned_seqs_file, "clustal")

# Function to calculate the distance matrix based on the method
def cal_dist(aligned_seqs, method):
    calculator = DistanceCalculator(method)
    return calculator.get_distance(aligned_seqs)

# Function to construct the phylogenetic tree, given the distance matrix and the method: upgma or nj.
def construct_tree_from_matrix(distance_matrix, method):
    constructor = DistanceTreeConstructor()
    if method == 'upgma':
        return constructor.upgma(distance_matrix)
    elif method == 'nj':
        return constructor.nj(distance_matrix)

# Finally, we draw a nice phylogenetic tree
def draw_tree(tree):
    tree_io = StringIO()
    Phylo.draw(tree, do_show=False, write_to=tree_io, format='svg')
    tree_image = tree_io.getvalue()
    return tree_image





