from Bio.Seq import Seq


def process_seq_by_class(seq: Seq, seq_class: int) -> Seq:
    if seq_class not in [1, 3, 6, 7]:
        return seq


def trim_x_only_keep_in_between(seq: Seq) -> Seq:
    pass


def trim_x_only(seq: Seq) -> Seq:
    pass


def trim_x_and_5prime(seq: Seq) -> Seq:
    pass


def trim_x_and_3prime(seq: Seq) -> Seq:
    pass
