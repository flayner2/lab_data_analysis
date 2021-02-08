from Bio.SeqIO import SeqRecord
import re


def find_x_regions_and_calculate_stats(
    sequence: SeqRecord, taxon: str
) -> list:
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

        seq_features["seq_xgroup_count"] = x_group_counter
        seq_features["xgroup_len"] = xgroup_len
        seq_features["dist_from_3"] = start
        seq_features["dist_from_5"] = seq_features["seq_len"] - end

        seq_features_list.append(seq_features)
        x_group_counter += 1

    if len(seq_features_list) > 0:
        seq_class = get_seq_class(seq_features_list)

        for each_dict in seq_features_list:
            each_dict["seq_class"] = seq_class

    return seq_features_list


def get_seq_class(x_groups_list: list[dict]) -> int:
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
