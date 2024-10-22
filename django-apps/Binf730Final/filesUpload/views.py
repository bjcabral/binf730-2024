# View file for the upload form
import os
from io import StringIO

#from Bio.Align import MultipleSeqAlignment
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage
from Bio import SeqIO, AlignIO, Phylo
#from Bio import Align
from .forms import SequenceFileForm, AlignmentScoreForm, AlignmentMethodForm, DistanceMatrixForm, TreeConstructionForm, \
    SubstitutionMatrixForm
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
#from io import
from Bio.Align.Applications import ClustalOmegaCommandline
from . import config
import uuid


# Specify the directory and filename
output_directory = config.DATA_DIR
aligned_file = f"{output_directory}/{config.ALIGNED_FILE_NAME}"
# Save the files in fasta format in one single file when enter manually or when fasta files are uploaded
sequences_file = f"{output_directory}/{config.INPUT_FILE_NAME}"
os.makedirs(os.path.dirname(aligned_file), exist_ok=True)

# Append the entered sequences, manual or fasta file into one single file for the aligner.
def append_to_fasta_file(sequence_id, sequence):
    with open(sequences_file, 'a') as fasta_file:
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
            open(sequences_file, 'w').close()

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

                    request.session['fasta_file'] = sequences_file
            return redirect('method_and_score_scheme')
        else:
            return HttpResponse('Form is not valid.')
    else:
        form = SequenceFileForm()
    return render(request, 'upload.html', {'form': form})

# This function prompts the user for the alignment method (global or local) and the score matrix.
def method_and_score_scheme(request):
    fasta_file = request.session.get('fasta_file')

    if not fasta_file:
        return HttpResponse("Sequences not found. Please upload again.")

    if request.method == 'POST':
        alignment_method_form = AlignmentMethodForm(request.POST)
        substitution_matrix_form = SubstitutionMatrixForm(request.POST)

        if alignment_method_form.is_valid() and substitution_matrix_form.is_valid():
            alignment_method = alignment_method_form.cleaned_data['alignment_method']
            substitution_matrix = substitution_matrix_form.cleaned_data['substitution_matrix']

            request.session['alignment_method'] = alignment_method
            request.session['substitution_matrix'] = substitution_matrix

            if substitution_matrix == 'manual':
                match_score = substitution_matrix_form.cleaned_data['match_score']
                mismatch_score = substitution_matrix_form.cleaned_data['mismatch_score']
                gap_score = substitution_matrix_form.cleaned_data['gap_score']

                request.session['match_score'] = match_score
                request.session['mismatch_score'] = mismatch_score
                request.session['gap_score'] = gap_score

            return redirect('align_sequences')
    else:
        alignment_method_form = AlignmentMethodForm()
        substitution_matrix_form = SubstitutionMatrixForm()

    return render(request, 'method_and_score_scheme.html', {
        'method_form': alignment_method_form,
        'matrix_form': substitution_matrix_form
    })

# These methods uses the Align() function from BioPython to execute the sequence alignment and save the aligned
# sequence file for later use with the cal_distance() function.

def align_sequences(request):
    sequences = request.session.get('sequences')
    if not sequences:
        return HttpResponse("Sequences not found. Please upload again.")

    try:
        # Create processed_data directory if it doesn't exist
        processed_data_dir = os.path.join(config.DATA_DIR, 'processed_data')
        os.makedirs(processed_data_dir, exist_ok=True)

        # Write sequences to a temporary file
        input_file = os.path.join(processed_data_dir, "temp_input.fasta")
        with open(input_file, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">Seq{i+1}\n{seq}\n")

        # Generate a unique name for the aligned file
        aligned_file_name = f"aligned_{uuid.uuid4().hex[:8]}.fasta"
        aligned_file_path = os.path.join(processed_data_dir, aligned_file_name)

        # Set up Clustal Omega command
        clustalomega_cline = ClustalOmegaCommandline(
            infile=input_file,
            outfile=aligned_file_path,
            verbose=True,
            auto=True
        )

        # Run Clustal Omega
        stdout, stderr = clustalomega_cline()

        # Read the aligned sequences
        aligned_seqs = AlignIO.read(aligned_file_path, "clustal")

        # Generate HTML response
        aligned_seq_html = "<h2>Multiple Sequence Alignment:</h2>"
        aligned_seq_html += f"<pre>{aligned_seqs}</pre>"

        # Store the aligned file name in the session
        request.session['aligned_file_name'] = aligned_file_path

        # Clean up temporary input file
        os.remove(input_file)

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
            aligned_file_name = request.session.get('aligned_file_name')
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





