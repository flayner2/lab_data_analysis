from Bio.SeqIO import SeqRecord
import re


class XGroup:
    """
    An XGroup represents a vector-masked subsequence of a particular sequence. It is
    called XGroup because the subsequence is represented as a string o "X"s.
    """

    def __init__(self, xgroup_len: int, dist_from_5: int, dist_from_3: int) -> None:

        # Length of the sequence of "X"s
        self.xgroup_len = xgroup_len

        # The distances from 5' and 3'
        self.dist_from_5 = dist_from_5
        self.dist_from_3 = dist_from_3

    def __repr__(self) -> str:
        return f"XGroup<{self.xgroup_len}|{self.dist_from_5}|{self.dist_from_3}>"


def find_x_regions_and_calculate_stats(sequence: SeqRecord, taxon: str) -> list[dict]:
    """Takes a nucleotide sequence with masked regions represented by Xs and calculates
    the number o X regions, the length of each X region, the distances of each X region
    to the 3' and 5' ends of the sequence, the length of the sequence and the sequence
    class based on X region features

    Args:
        sequence (SeqRecord): a SeqRecord object representing a nucleotide sequence
        taxon (str): the name of the taxon the sequence belongs to

    Returns:
        list[dict]: a list of dictionaries for each X group in the sequence
    """
    seq_features_list = []
    x_group_counter = 1

    pattern = re.compile(r"X+")
    search_res = pattern.finditer(str(sequence.seq))

    seq_features = {
        "seq_id": sequence.id,
        "taxon": taxon,
        "seq_len": len(sequence.seq),
    }

    for match in search_res:
        start, end = match.span()
        xgroup_len = len(match.group())

        seq_features["seq_xgroup_count"] = x_group_counter  # TODO: move this down
        seq_features["xgroup_len"] = xgroup_len
        seq_features["dist_from_5"] = start
        seq_features["dist_from_3"] = seq_features["seq_len"] - end

        seq_features_list.append(seq_features)
        x_group_counter += 1

    if len(seq_features_list) > 0:
        seq_class = get_seq_class(seq_features_list)

        for each_dict in seq_features_list:
            each_dict["seq_class"] = seq_class

    return seq_features_list


def get_seq_class(x_groups_list: list[dict]) -> int:
    """Infers the sequence class based on its X groups characteristics

    Args:
        x_groups_list (list[dict]): a list of dictionaries containing information about
                                    the X groups of a sequence

    Returns:
        int: the sequence group code
    """
    num_of_xgroups = len(x_groups_list)

    if num_of_xgroups == 2:
        return 1
    elif num_of_xgroups > 2:
        return 2
    else:
        x_group = x_groups_list[0]
        xgrop_len = x_group["xgroup_len"]
        dist_from_5 = x_group["dist_from_5"]

        if xgrop_len <= 300 and dist_from_5 < 50:
            return 3
        elif xgrop_len > 300 and dist_from_5 < 50:
            return 4
        elif xgrop_len <= 300 and (50 < dist_from_5 <= 300):
            return 5
        elif xgrop_len > 300 and (50 < dist_from_5 <= 300):
            return 6
        else:
            return 7
