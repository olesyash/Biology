from Bio.Align import PairwiseAligner
from Bio import SeqIO
from itertools import combinations
from Bio import Entrez, SeqIO
from io import StringIO
from Bio import pairwise2

Entrez.email = "jcecomputationalbiology@gmail.com"  # Provide an email address

def count_matches(seq1, seq2):
    """Count exact matches between two aligned sequences."""
    return sum(1 for i in range(len(seq1)) if seq1[i] == seq2[i])


def create_comparison_dict(alignments_lst):
    """
    Analyze differences between sequences for all alignments in the list.

    Args:
        alignments_lst: List of tuples (seq1_desc, seq2_desc, alignment)

    Returns:
        list: List of dictionaries containing comparison data and sequence descriptions
    """
    results = []

    for seq1_desc, seq2_desc, alignment in alignments_lst:
        seq1 = str(alignment[0])
        seq2 = str(alignment[1])

        matches = 0
        transitions = 0
        transversions = 0
        gap_1 = 0
        gap_2 = 0

        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                matches += 1
            elif seq1[i] == '-':
                gap_1 += 1
            elif seq2[i] == '-':
                gap_2 += 1
            elif (seq1[i].upper(), seq2[i].upper()) in [('C', 'T'), ('T', 'C'), ('A', 'G'), ('G', 'A')]:
                transitions += 1
            else:
                transversions += 1

        comparison_data = {
            'seq1_desc': seq1_desc,
            'seq2_desc': seq2_desc,
            'matches': matches,
            'transitions': transitions,
            'transversions': transversions,
            'gap_1': gap_1,
            'gap_2': gap_2,
            'sequence_length': len(seq1)
        }
        results.append(comparison_data)

    return results


def create_alignments(org_lst,  mode: str = "global", match: int = 1,
                     mismatch: int = 0, open_gap: int = 0, extend_gap: int = 0):
    """Create pairwise alignments for all sequence combinations."""

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    
    return [(seq1.description, seq2.description, aligner.align(seq1, seq2)[0])
            for seq1, seq2 in combinations(org_lst, 2)]

def create_alignments_with_globalms(org_lst, mode: str = "global", match: int = 1,  mismatch: int = 0,
                                    open_gap: int = 0, extend_gap: int = 0):
    """Create pairwise alignments for all sequence combinations."""
    alignments = []
    for seq1, seq2 in combinations(org_lst, 2):
        alignment = pairwise2.align.globalms(seq1.seq, seq2.seq, match, mismatch, open_gap, extend_gap)
        alignments.append((seq1.description, seq2.description, alignment[0]))
    return alignments

def analyze_alignmented_sequences(alignments_lst):
    """Main function to analyze sequences from a list of organisms."""
    # Calculate scores
    results = []
    for seq1_desc, seq2_desc, alignment in alignments_lst:
        matches = count_matches(str(alignment[0]), str(alignment[1]))
        normalized_score = round((matches / len(alignment[0])) * 100, 2)
        results.append((seq1_desc, seq2_desc, normalized_score))

    return results


if __name__ == "__main__":
    fasta_file = "ex2_sequences_a.fasta"
    org_lst = list(SeqIO.parse(fasta_file, "fasta"))

    # Use the globalms function
    alignments2 = create_alignments_with_globalms(org_lst)
    for aligment in alignments2:
        score = aligment[2][2]
        print(f"\nComparison between {aligment[0]}. and {aligment[1]} using globalms, score is {score}")

    alignments = create_alignments(org_lst)

    # Get both scores and detailed comparison data
    scores = analyze_alignmented_sequences(alignments)
    comparisons = create_comparison_dict(alignments)

    # Print both types of results
    for score_data in scores:
        print(f"{score_data[0]} VS {score_data[1]}: {score_data[2]}%")

    for comp in comparisons:
        print(f"\nComparison between {comp['seq1_desc']} and {comp['seq2_desc']}:")
        print(comp)

    with Entrez.efetch(db="nucleotide", id=["AF451972", "AF176731", "X90314"],
                       rettype="fasta") as handle:
        fasta_data = handle.read()
        org_lst = list(SeqIO.parse(StringIO(fasta_data), "fasta"))
        alignments = create_alignments(org_lst)
        scores = analyze_alignmented_sequences(alignments)
        comparisons = create_comparison_dict(alignments)

    print("\n")
    print("*" * 50)
    print("\n")
    # Print both types of results
    for score_data in scores:
        print(f"{score_data[0]} VS {score_data[1]}: {score_data[2]}%")

    for comp in comparisons:
        print(f"\nComparison between {comp['seq1_desc']} and {comp['seq2_desc']}:")
        print(comp)

    test = [('org1', 'org2', ('CCGT----GGATCT', 'ATGTAAAAG-GGGT'))]
    test_comparisons = create_comparison_dict(test)
    for comp in test_comparisons:
        print(f"\nComparison between {comp['seq1_desc']} and {comp['seq2_desc']}:")
        print(comp)

