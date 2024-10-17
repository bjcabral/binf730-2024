# Name: Braulio J. Cabral
# BINF 730 Biological Sequence Analysis
# Homework 2, question 2: Question 5. Write code to implement, UPGMA, WPGMA, and neighbor-joining algorithms.
# Your code should produce output from which you will draw a phylogenetic tree for some sample data.
# Please use the example from the lecture and the primate data from the lecture.
# Also, pick your own set of short sequences.

# This program combines all three methods UPGMA, WPGMA, Neighbor-joining. The program prompts the user for
# a set of sequences file or to enter short sequences manually then select the method to use.
# The output includes a distance matrix and displays a Phylogenetic Tree.

# Imports
import os
import numpy as np
import pandas as pd
from tabulate import tabulate
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

# Function to read fasta formatted sequence files
def read_fasta_sequence(input_file_path, debug = True):
    output_file_path = "/home/ubuntu/binf-730/homework2-code/concadenated_fasta_file.fa"
    if debug:
        if os.path.exists(input_file_path):
            print(f"File exist in: {input_file_path}")
        else:
            print(f"Could not find the file with path: {input_file_path}")

    records = list(SeqIO.parse(input_file_path, "fasta"))
#    sequences = [(record.id, record.description, str(record.seq)) for record in records]
    sequences = [record.format("fasta") for record in records]

    # write all records to a single file
    with open(output_file_path, 'a') as outfile:
        for sequence in sequences:
            outfile.write(sequence)

    return sequences 

# the function below calls an EBI API endpoint to use clustalw2 for multiple sequence alignment and
# returns the aligned sequences.
# Well, this was giving me a hard time due to data formatting so I used a local installation of Clustalw on Linux Ubuntu
# the function takes the path to the concadenated files for alignment and returns the path to the aligned sequences
# call the function with first the fasta file to aligned and the path to the output file

def Seq_Alignment(file_path, outfile):
    # input file
    sequences_to_align = file_path 
    # output file path
    aligned_seqs = outfile 
    # local Clustalw2 path
    clustalw_path = "/usr/bin/clustalw"

    # Clustalw command line
    clustalw_cline = ClustalwCommandline(clustalw_path, infile=sequences_to_align, align=True, outorder="ALIGNED", outfile=aligned_seqs, output="clustal")

    # Make sure clustaw exists
    assert os.path.isfile(clustalw_path), "ClustalW executable not found. Please, check path or install clustalw"
    # Time to execute the command
    try:
       stdout, stderr = clustalw_cline()
       print(f"Alignment completed, output file is {outfile}")
    except Exception as e:
                         print(f"Error executing ClustalW command: {e}")


# A function for calculating the distance matrix from the aligned sequences

def read_aligned(aligned_seqs_file):
    aligned_seqs = AlignIO.read(aligned_seqs_file, "clustal")
    return aligned_seqs 

def hamming_distance(seq1, seq2):
    #Calculate distance for a pair of sequences. This will be used in a loop to iterate over multiple sequences
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

# read the aligned files and parse the data in a format that the distance can be calculated
def cal_dist_matrix(alignments):
    sequences = [] 
    sequence_ids = []
    distances = []
    #Read all sequence records and put in array
    for record in alignments:
        sequences.append(record.seq)
        sequence_ids.append(record.id)
        #print(f"ID: {record.id}, Sequence: {record.seq}")

    # use the array to iterate over sequences and calculate distance
    for i in range(len(sequences)):
       for j in range(i + 1, len(sequences)):
           distances.append(hamming_distance(sequences[i], sequences[j]))
           print(f"Hamming distance between seuence {i + 1} and sequence {j+1}: {distances}")
    print(f"Here are the distances: {distances}")
    
    # Construct a matrix out of the distance values
    matrix_size = len(sequence_ids)

    # Create a placeholder matrix filled with zeros
    matrix = np.zeros((matrix_size, matrix_size))

    # Populate the new matrix with the distances
    index = 0
    for i in range(matrix_size):
        for j in range(i + 1, matrix_size):
            matrix[i, j] = distances[index]
            matrix[j, i] = distances[index]
            index += 1

    # Using Panda to create a dataframe so we can display the matrix beautifully :-)
    matrix_df = pd.DataFrame(matrix, columns=sequence_ids, index=sequence_ids)

    #Let's see how this looks
    print(tabulate(matrix_df, headers=sequence_ids, tablefmt="fancy_grid"))


    

# def process_upgma():


#   distance_matrix =
#  return distance_matrix

# def process_wpgma():

# def process_neighbor_joining():

# ef create_phylo():


# Main function

def main():
    concadenated_file_path = "/home/ubuntu/binf-730/homework2-code/concadenated_fasta_file.fa"
    aligned_file = "/home/ubuntu/binf-730/homework2-code/aligned_seqs.fa"

    input_method = input("Enter sequence input method (manual or fasta): ")
    sequences = []

    if input_method.lower() == 'manual':
        num_sequences = int(input("Enter the number of sequences to input manually: "))
        for i in range(num_sequences):
            seq_id = input(f"Enter ID for sequence {i + 1}: ")
            seq = input(f"Enter sequence {i + 1}: ")
            sequences.append((seq_id, seq))

        with open(concadenated_file_path, "a") as concat:
            print(f"number of sequences = {len(sequences)}")
            for id, seq in sequences:
                concat.write(f">{id}\n{seq}\n")

    elif input_method.lower() == 'fasta':
        num_files = int(input("Enter the number of FASTA files to input: "))
        for i in range(num_files):
            file_path = input(f"Enter the path to FASTA file {i + 1}: ")
            sequences.extend(read_fasta_sequence(file_path))
    else:
        raise ValueError('Invalid input method selected')


# input and output files paths
    # Test the above code
    print("sequences")
#    for_alignment = []
#    for i, seq in enumerate(sequences):
#        # print(f"Sequence {i+1}: {seq}")
#        for_alignment.append(seq)
#        # print(sequences)
#    for i, seq in enumerate(for_alignment):
#        # print(f"Alignment {i + 1}: {seq}")
#      print(seq)

    # algorithm = input("Enter the phylogenetic method to use (upgma, wpgma, or neighbor): ")
    Seq_Alignment(concadenated_file_path, aligned_file)
    alignment = read_aligned(aligned_file) 
    cal_dist_matrix(alignment)

if __name__ == "__main__":
    main()

