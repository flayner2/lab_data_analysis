import re
from Bio.Seq import Seq


def process_seq_by_class(seq: Seq, seq_class: int) -> Seq:
    assert type(seq) == Seq, "Argument {seq} must be of type Seq"

    if seq_class not in [1, 3, 6, 7]:
        return seq

    if seq_class == 1:
        return trim_x_only_keep_in_between(seq)


def trim_x_only_keep_in_between(seq: Seq) -> Seq:
    pattern = re.compile(r"X+")
    search_res = pattern.finditer(str(seq))

    x1 = next(search_res)
    x2 = next(search_res)

    _, x1_end = x1.span()
    x2_start, _ = x2.span()

    return seq[x1_end:x2_start]


def trim_x_only(seq: Seq) -> Seq:
    pass


def trim_x_and_5prime(seq: Seq) -> Seq:
    pass


def trim_x_and_3prime(seq: Seq) -> Seq:
    pass
