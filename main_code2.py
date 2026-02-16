from Bio import SeqIO
# Define function for loading the sequence database from a FASTA file
def load_sequences(fasta_file):
    sequences = []
    # Iterate through each record in the FASTA file and append it to the list using parse function from SeqIO
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record))
    return sequences

from Bio import pairwise2
#Define function for performing global alignment between two sequences and returning the best alignment
def alignment_score(seq1, seq2):
    # Perform global alignment using the pairwise2 module from Biopython
    alignments = pairwise2.align.globalxx(seq1, seq2)
    # Get the best alignment score (the first one in the list)
    best_alignment = alignments[0]

    return best_alignment.score

def find_closest_breed(test_seq, database):
    best_score = -1
    best_record = None

    for record in database:
        score = alignment_score(test_seq, record.seq)

        if score > best_score:
            best_score = score
            best_record = record
    return best_record.id, best_score

import numpy as np
# Define function to compute p-value (z-score) for the best alignment score compared to a distribution of scores
def p_value(best_score, all_scores):
    mean = np.mean(all_scores)
    std = np.std(all_scores)
# Avoid division by zero in case all scores are the same
    if std == 0:
        return 0
# Compute z-score for the best score
    z_score = (best_score - mean) / std
    return z_score

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def build_tree(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")

    calc = DistanceCalculator('identity')
    dm = calc.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    return tree

def build_alignment_from_sequences(records, alignment_file):
    max_len = max(len(record.seq) for record in records)
    aligned_records = []

    for record in records:
        seq = str(record.seq)
        if len(seq) < max_len:
            # Simple right-pad with gaps so all sequences share length
            seq = seq + "-" * (max_len - len(seq))
        aligned_records.append(
            SeqRecord(Seq(seq), id=record.id, description="")
        )

    alignment = MultipleSeqAlignment(aligned_records)
    AlignIO.write(alignment, alignment_file, "fasta")

def main():

    # Load database and test sequence
    database = load_sequences("Data/toy_refs.fasta")
    test_seq_record = load_sequences("toy_test.fasta")[0]

    test_seq = test_seq_record.seq

    # Find closest breed
    best_score = -1
    best_id = None
    all_scores = []

    for record in database:
        score = alignment_score(test_seq, record.seq)
        all_scores.append(score)

        if score > best_score:
            best_score = score
            best_id = record.id

    print("Closest breed:", best_id)
    print("Best alignment score:", best_score)

    # Compute p-value (z-score)
    z = p_value(best_score, all_scores)
    print("Z-score:", z)

    #  phylogeny
    alignment_path = "alignment.fasta"
    if not os.path.exists(alignment_path):
        build_alignment_from_sequences(database + [test_seq_record], alignment_path)

    tree = build_tree(alignment_path)
    Phylo.draw(tree)

if __name__ == "__main__":
    main()



