# ---------------- IMPORTS ----------------
import random
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt


# ---------------- I/O FUNCTIONS ----------------
def read_fasta(filepath):
    """
    Reads a FASTA file and returns a dictionary
    with header as key and sequence as value.
    """
    sequences = {}

    with open(filepath) as f:
        name = None
        seq = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)

        if name:
            sequences[name] = "".join(seq)

    return sequences


# ---------------- SEQUENCE ANALYSIS ----------------
def percent_identity(seq1, seq2):
    """
    Calculates percentage identity between two sequences.
    Assumes sequences are same length.
    """
    matches = 0

    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            matches += 1

    return (matches / len(seq1)) * 100


def percent_difference(seq1, seq2):
    """
    Returns percent difference between sequences.
    """
    return 100 - percent_identity(seq1, seq2)


def find_closest(mystery_seq, database):
    """
    Finds closest breed match.
    Returns breed name and similarity.
    """
    best_breed = None
    best_score = -1

    for breed, seq in database.items():
        score = percent_identity(mystery_seq, seq)

        if score > best_score:
            best_score = score
            best_breed = breed

    return best_breed, best_score


# ---------------- STATISTICS ----------------
def p_value(mystery_seq, best_seq, database, n=1000):
    """
    Calculates p-value using random sampling.
    """
    observed = percent_identity(mystery_seq, best_seq)

    count = 0
    db_seqs = list(database.values())

    for _ in range(n):
        rand_seq = random.choice(db_seqs)
        score = percent_identity(mystery_seq, rand_seq)

        if score >= observed:
            count += 1

    return count / n


# ---------------- PHYLOGENY ----------------
def build_distance_matrix(sequences):
    """
    Builds a distance matrix for phylogeny.
    """
    names = list(sequences.keys())
    matrix = []

    for name1 in names:
        row = []
        for name2 in names:
            dist = percent_difference(
                sequences[name1], sequences[name2]
            ) / 100
            row.append(dist)
        matrix.append(row)

    return DistanceMatrix(names, matrix)


def build_tree(distance_matrix):
    """
    Constructs Neighbor Joining tree.
    """
    constructor = DistanceTreeConstructor()
    return constructor.nj(distance_matrix)


def plot_tree(tree):
    """
    Displays phylogenetic tree.
    """
    Phylo.draw(tree)
    plt.show()


# ---------------- MAIN PROGRAM ----------------
if __name__ == "__main__":

    print("Loading sequences...")

    database = read_fasta("dog_breeds.fa")
    mystery = read_fasta("mystery.fa")

    mystery_seq = list(mystery.values())[0]

    # ---- Closest match ----
    print("\nFinding closest breed...")
    breed, similarity = find_closest(mystery_seq, database)

    print("Closest breed:", breed)
    print("Similarity: {:.2f}%".format(similarity))
