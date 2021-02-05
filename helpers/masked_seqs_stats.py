from Bio.SeqIO import SeqRecord
from typing import Union


def find_x_regions_and_calculate_stats(
    sequence: SeqRecord,
) -> Union[dict[str, Union[str, int]], None]:
    pass


def get_dist_to_3prime(sequence: SeqRecord) -> Union[int, None]:
    first_x_position = sequence.seq.find("X")

    if first_x_position != -1:
        return first_x_position
    else:
        return None


def get_dist_to_5prime(sequence: SeqRecord) -> Union[int, None]:
    last_x_position = sequence.seq.rfind("X")

    if last_x_position != -1:
        return len(sequence.seq) - last_x_position - 1
    else:
        return None


def get_xregion_len():
    pass


def get_seq_class():
    pass