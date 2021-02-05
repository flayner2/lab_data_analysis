from Bio.SeqIO import SeqRecord
from typing import Union
import re


def find_x_regions_and_calculate_stats(
    sequence: SeqRecord, taxon: str
) -> list:
    seq_features_list = []
    x_group_counter = 0

    pattern = re.compile(r"X+")
    search_res = pattern.finditer(sequence.seq)

    seq_features = {
        "seq_id": sequence.id,
        "taxon": taxon,
        "seq_len": len(sequence.seq),
    }

    for match in search_res:
        start, end = match.span()
        xgroup_len = len(match.group())

        seq_features["xgroup_id"] = x_group_counter
        seq_features["xgroup_len"] = xgroup_len
        seq_features["dist_from_3"] = start
        seq_features["dist_from_5"] = seq_features["seq_len"] - end

        seq_features_list.append(seq_features)
        x_group_counter += 1

    return seq_features_list


def get_seq_class():
    pass
