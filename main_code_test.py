
from Bio import SeqIO
from Bio import AlignIO
from fileinput import filename

def load_fasta(filename):
#Reads a FASTA file and returns a list of tuples (id, sequence)
with open(filename, 'r') as f:
    # Initialize variables
    records = []
    id = None
    sequence = ''
    # Iterate through lines in the file
    for line in f:
        # Remove whitespace and newline characters
        line = line.strip()
        # Check if the line is a header (starts with '>')
        if line.startswith('>'):
            # If we have an existing record, save it before starting a new one
            if id is not None:
                records.append((id, sequence))
            id = line[1:]  # Remove '>'
            sequence = ''
        # If it's not a header, it's part of the sequence
        else:
            # Append the line to the current sequence
            sequence += line
    # After the loop, save the last record if it exists
    if id is not None:
        records.append((id, sequence))
    # Return the list of records
    return records 

def global_alignment(seq_a, seq_b):

    return alignment

def percent_difference(alignment):
    return mismatches / aligned_length

def identify_breed(test_seq, ref_records):
   return best_id, best_diff, all_diffs

def p_value_best_match(test_seq, ref_records, n=1000):
    return p_value, null_dist

def build_distance_matrix(ref_records, include_test=None):
    return matrix + labels

def build_tree(distance_matrix, labels, method='upgma'):
    return tree

