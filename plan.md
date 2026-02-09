Plan for coursework:

Goal
- Input: reference database of dog DNA sequences (FASTA) + one test sequence
- Output: closest breed match + percent difference
- Stretch: p-value for match; reconstructed phylogeny

Assumptions (update once dataset is released)
- Reference sequences are FASTA with breed name in header
- Test sequence is FASTA or raw sequence string
- All sequences are from the same region/gene

Core pipeline (high level)
1) Load reference sequences from FASTA
2) For each reference, align test sequence to reference
3) Compute percent difference from alignment
4) Pick the lowest percent difference as the closest match
5) Report breed name + percent difference

Stretch 1: p-value for match
- Build null distribution by permuting test sequence N times
- For each permuted sequence, compute best match distance
- p-value = proportion of permuted best distances <= observed best distance

Stretch 2: phylogeny
- Compute pairwise distances among references (and optionally test)
- Build tree using UPGMA or Neighbor Joining
- Output Newick string or simple tree diagram

Function layout (draft)
def load_fasta(path):
    # returns list of (id, sequence)


def global_align(seq_a, seq_b):
    # returns alignment result for two sequences


def percent_difference(alignment):
    # returns mismatches / aligned length


def identify_closest(test_seq, ref_records):
    # returns best_id, best_diff, all_diffs


def p_value_best_match(test_seq, ref_records, n=1000):
    # returns p_value, null_dist


def build_distance_matrix(ref_records, include_test=None):
    # returns matrix + labels


def build_tree(distance_matrix, labels, method="upgma"):
    # returns Newick string or tree object

Notes
- Use Biopython for FASTA parsing, pairwise alignment, and phylogeny
- Decide how to handle gaps (count as differences or ignore) and be consistent
