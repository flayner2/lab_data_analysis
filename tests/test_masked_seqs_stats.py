import sys
from Bio.SeqIO import SeqRecord

# Configuring local imports
sys.path.insert(1, "/home/maycon/Documents/LAB/lab_data_analysis/helpers")


def test_find_x_regions_and_calculate_stats() -> None:
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
