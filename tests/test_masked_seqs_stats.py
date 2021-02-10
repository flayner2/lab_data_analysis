from Bio.SeqIO import SeqRecord
from helpers import masked_seqs_stats

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
taxon = "Polistes_canadensis"
features_list = masked_seqs_stats.find_x_regions_and_calculate_stats(seq, taxon)
features = features_list[0]

normal_seq = SeqRecord(
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
no_xgroups = masked_seqs_stats.find_x_regions_and_calculate_stats(normal_seq, taxon)


def test_num_of_xgroups() -> None:
    assert len(features_list) == 1


def test_seq_len() -> None:
    assert len(seq.seq) == features["seq_len"]


def test_seq_id() -> None:
    assert seq.id == features["seq_id"]


def test_xgroup_id() -> None:
    assert features["xgroup_id"] == 0


def test_xgroup_len() -> None:
    substring = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    assert len(substring) == features["xgroup_len"]


def test_distance_to_3prime() -> None:
    first = seq.seq.find("X")
    assert first == features["dist_from_3"]


def test_distance_to_5prime() -> None:
    last = seq.seq.rfind("X")
    dist = len(seq.seq) - last - 1
    assert dist == features["dist_from_5"]


def test_normalseq_noxgroups() -> None:
    assert len(no_xgroups) == 0


def test_get_seq_class() -> None:
    seq_class = masked_seqs_stats.get_seq_class(features_list)
    assert seq_class == 5


def test_seq_class_attribution() -> None:
    assert features["seq_class"] == 5
