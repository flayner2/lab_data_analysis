import re

from Bio.SeqIO import SeqRecord
from helpers import masked_seqs_stats


def test_seq_without_xgroups() -> None:
    seq = SeqRecord(
        seq=(
            "ACCTATAGGTTGTCGTCGACAAAGAAATGAATCAACTTCCTCTGGTGGTT"
            "CATGGCAAATGATATCTGGAACTGGTAGTTTACGTGGTTCAACAACAGCC"
            "CACACATCTATTACAGAGGGATCTAATTCTTCTGGCTCGACTAGCAAAGG"
            "TTTATTTGAAAATTTTTTACATCAAGCTCATGGATCTAGTAAAGCAATAT"
            "TGGAAGATGACGAATCCGTATCACAAGTACCTGCCCGGGCGGCCGCTCGA"
            "AAGCCG"
        ),
        id="DT319104.1",
        name="DT319104.1",
    )
    seq_class, seq_xgroups = masked_seqs_stats.find_x_regions_and_calculate_stats(seq)

    # Check if it has no xgroups
    assert len(seq_xgroups) == 0

    # Check if the seq class is invalid (i.e., 0)
    assert seq_class == 0


def test_seq_with_xgroups() -> None:

    seq = SeqRecord(
        seq=(
            "AAGCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCGA"
            "GACGGCCGCCCGGGCAGGTACACCCAAGGATTTAATCGTCAAACCATGACGGGTCTCGAAAATCGAAA"
            "CGGACAACATGACAAGGAAATGGGCCCGATGATGAACGAAGTCACCAGACCGAGATACATCAGGGACG"
            "ATAAGAATGCCAAAATTATCGACACATCGGTGGAAAC"
        ),
        id="DT319107.1",
        name="DT319107.1",
    )
    seq_class, seq_xgroups = masked_seqs_stats.find_x_regions_and_calculate_stats(seq)

    # Check the ammount of xgroups
    assert len(seq_xgroups) == 1

    # Check the length of the xgroup
    pattern = re.compile(r"X+")
    substring = pattern.search(str(seq.seq))

    assert len(substring.group(0)) == seq_xgroups[0].xgroup_len

    # Check the xgroup distance to 5'
    first = seq.seq.find("X")

    assert first == seq_xgroups[0].dist_from_5

    # Check the xgroup distance to 3'
    last = seq.seq.rfind("X")
    dist = len(seq.seq) - last - 1

    assert dist == seq_xgroups[0].dist_from_3

    # Check the seq class
    assert seq_class == 3
