# BINF-730 Final Project
# Fall 2024 Semester
# Student Name: Braulio J. Cabral
# Professor: Aman Ullah
# George Mason University
# Bioinformatics and Computations Biology Graduate Certificate

# View file for the upload form
import os
from . import config
from io import StringIO
from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage
from Bio import SeqIO, AlignIO, Phylo
from .forms import SequenceFileForm, AlignmentScoreForm, AlignmentMethodForm, DistanceMatrixForm, TreeConstructionForm, \
    SubstitutionMatrixForm
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, DistanceMatrix
from Bio.Align.Applications import ClustalOmegaCommandline
import uuid
import matplotlib.pyplot as plt
from io import BytesIO
import base64

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
    fasta_file = request.session.get('fasta_file')
    if not fasta_file:
        return HttpResponse("Sequences not found. Please upload again.")

    try:
        # Create processed_data directory if it doesn't exist
        processed_data_dir = os.path.join(config.DATA_DIR, 'processed_data')
        os.makedirs(processed_data_dir, exist_ok=True)

        # Generate a unique name for the aligned file
        import uuid
        aligned_file_name = f"aligned_{uuid.uuid4().hex[:8]}.aln"
        aligned_file_path = os.path.join(processed_data_dir, aligned_file_name)

        # Set up Clustal Omega command
        clustalomega_cline = ClustalOmegaCommandline(
            infile=fasta_file,
            outfile=aligned_file_path,
            verbose=True,
            auto=True,
            outfmt="clustal"
        )

        # Run Clustal Omega
        stdout, stderr = clustalomega_cline()

        # Verify the output file exists and is not empty
        if not os.path.exists(aligned_file_path) or os.path.getsize(aligned_file_path) == 0:
            raise Exception("Alignment failed to produce output file")

        # Read the aligned sequences to verify format
        try:
            aligned_seqs = AlignIO.read(aligned_file_path, "clustal")
        except ValueError as e:
            # If reading as CLUSTAL fails, try to read as FASTA and convert
            seqs = list(SeqIO.parse(fasta_file, "fasta"))
            if not seqs:
                raise Exception("No sequences found in the input file")
            AlignIO.write(seqs, aligned_file_path, "clustal")
            aligned_seqs = AlignIO.read(aligned_file_path, "clustal")

        # Store the aligned file path in the session
        request.session['aligned_file_name'] = aligned_file_path

        return redirect('calculate_distance')

    except Exception as e:
        return HttpResponse(f"Error occurred during alignment: {str(e)}")

# Function to read the aligned sequences file
def read_aligned(aligned_seqs_file):
    return AlignIO.read(aligned_seqs_file, "clustal")

# Function to calculate the distance matrix using the user selected method. BLOSUM250, identity, etc.
# this function collects the input from the user and uses read_aligned() and cal_dist() to read
# the aligned files in clustalo format and calculate the distance matrix.
def calculate_distance(request):
    if request.method == 'POST':
        form = DistanceMatrixForm(request.POST)
        if form.is_valid():
            distance_method = form.cleaned_data['distance_method']
            aligned_file_name = request.session.get('aligned_file_name')

            if not aligned_file_name:
                return HttpResponse("Aligned sequences not found. Please align sequences first.")

            try:
                aligned_seqs = read_aligned(aligned_file_name)
                calculator = DistanceCalculator(distance_method)
                distance_matrix = calculator.get_distance(aligned_seqs)

                # Convert DistanceMatrix to a list of lists
                matrix_list = [list(row) for row in distance_matrix]

                # Store both the matrix and the names in the session
                request.session['distance_matrix'] = matrix_list
                request.session['matrix_names'] = distance_matrix.names

                # Create a formatted string representation of the matrix
                matrix_str = ""
                for i, row in enumerate(matrix_list):
                    matrix_str += f"{distance_matrix.names[i]}: {' '.join([f'{x:.4f}' for x in row])}\n"

                # Render the template with the matrix data
                return render(request, 'calculate_distance.html', {
                    'form': form,
                    'matrix': matrix_list,
                    'matrix_names': distance_matrix.names,
                    'matrix_str': matrix_str,
                    'distance_method': distance_method,
                    'show_next_step': True
                })
            except Exception as e:
                return HttpResponse(f"Error calculating distance matrix: {str(e)}")
    else:
        form = DistanceMatrixForm()
    return render(request, 'calculate_distance.html', {'form': form})
             
# Function to construct the tree from the distance matrix using the algorithm selected by the user
# such as Neighbor Joining (nj) or UPGMA.
# uses matplotlib directly.
def construct_tree(request):
    if request.method == 'POST':
        form = TreeConstructionForm(request.POST)
        if form.is_valid():
            tree_method = form.cleaned_data['tree_method']
            matrix_list = request.session.get('distance_matrix')
            matrix_names = request.session.get('matrix_names')

            if not matrix_list or not matrix_names:
                return HttpResponse("Distance matrix not found. Please calculate distances first.")

            # Reconstruct the DistanceMatrix object
            distance_matrix = DistanceMatrix(names=matrix_names)
            for i, row in enumerate(matrix_list):
                for j, value in enumerate(row[:i]):
                    distance_matrix[i, j] = value

            tree = construct_tree_from_matrix(distance_matrix, tree_method)
            tree_image = draw_tree(tree)
            return render(request, 'display_tree.html', {'tree_image': tree_image})
    else:
        form = TreeConstructionForm()
    return render(request, 'construct_tree.html', {'form': form})

# Function to construct the phylogenetic tree, given the distance matrix and the method: upgma or nj.
def construct_tree_from_matrix(distance_matrix, method):
    constructor = DistanceTreeConstructor()
    if method == 'upgma':
        return constructor.upgma(distance_matrix)
    elif method == 'nj':
        return constructor.nj(distance_matrix)

# Finally, we draw a nice phylogenetic tree using the function below called by the construct_tree(request)
# function above. Could not use this function without getting io related error. Problem with Phylo.draw. will
# use matplotlib to draw the tree instead
#def draw_tree(tree):
#    tree_io = StringIO()
#    Phylo.draw(tree, do_show=False, write_to=tree_io, format='svg')
#    tree_image = tree_io.getvalue()
#    tree_io.close()
#    return tree_image

def draw_tree(tree):
    fig, ax = plt.subplots(figsize=(10, 8))
    Phylo.draw(tree, axes=ax, do_show=False)

    # Save the plot to a BytesIO object
    buf = BytesIO()
    plt.savefig(buf, format='png')
    plt.close(fig)

    # Encode the image to base64
    image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')

    # Create the HTML img tag
    img_tag = f'<img src="data:image/png;base64,{image_base64}" alt="Phylogenetic Tree">'

    return img_tag


