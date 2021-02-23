"""
This script removes certain subsequences from a sequence. For example, vector-masked
subsequences represented by "X"s or polynucleotide subsequences such as poly-A and
poly-T tails.
"""
# Third-party imports
from Bio.Seq import Seq

# Own module imports
from masked_seqs_stats import XGroup
from process_est_sequences import ESTSeq


def trim_x_only_keep_in_between(seq: Seq, xgroups: list[XGroup]) -> Seq:
    """Trims both XGroups from a sequence, keeping only the sequence in between them.

    Args:
        seq (Seq): the sequence containing vector-masked "XGroups".
        xgroups (list[XGroup]): a list with all XGroups for a sequence. Should be of
        length 2.

    Returns:
        Seq: the sequence in between both XGroups.
    """

    # This function only takes care of cases when there are exactly two XGroups in the
    # sequence, so we check that.
    assert len(xgroups) == 2, (
        "This function should only be called when there are exactly two XGroups for"
        " the sequence"
    )

    # We can do this since this function is only called on sequences that have
    # exactly 2 XGroups.
    x1 = xgroups[0]
    x2 = xgroups[1]

    # Here we find the lowest end, which is where the leftmost (closer to 5') XGroup
    # should end, and the highest start, which is where the rightmost (closer to 3')
    # XGroup should start.
    start = min(x1.end, x2.end)
    end = max(x1.start, x2.start)

    # We return a slice of the sequence that only contains what is in between the end
    # of the first XGroup and the start of the last XGroup. This way we remove both
    # XGroups and the old 5' and 3' ends, leaving only the sequence in between those.
    return seq[start:end]


def trim_x_and_5prime(seq: Seq, xgroups: list[XGroup]) -> Seq:
    """Removes the leading XGroup subsequence plus the 5' end o sequence.

    Args:
        seq (Seq): the sequence containing vector-masked "XGroups".
        xgroups (list[XGroup]): a list with all XGroups for a sequence. Should be of
        length 1.

    Returns:
        Seq: the sequence without its 5' end and the leading subsequence of "X"s.
    """

    # This function only takes care of cases when there's exactly one XGroups in the
    # sequence, so we check that.
    assert len(xgroups) == 1, (
        "This function should only be called when there's only one XGroup for the"
        " sequence."
    )

    # Position of the last "X" character in the sequence.
    last_x = xgroups[0].end

    return seq[last_x:]


def trim_x_and_3prime(seq: Seq, xgroups: list[XGroup]) -> Seq:
    """Removes the lagging XGroup subsequence plus the 3' end o sequence.

    Args:
        seq (Seq): the sequence containing vector-masked "XGroups".
        xgroups (list[XGroup]): a list with all XGroups for a sequence. Should be of
        length 1.

    Returns:
        Seq: the sequence without its 3' end and the laggin subsequence of "X"s.
    """

    # This function only takes care of cases when there's exactly one XGroups in the
    # sequence, so we check that.
    assert len(xgroups) == 1, (
        "This function should only be called when there's only one XGroup for the"
        " sequence."
    )

    # Position of the first "X" character in the sequence.
    first_x = xgroups[0].start

    return seq[:first_x]


def remove_xgroup_by_seq_class(seq: Seq, seq_class: int, xgroups: list[XGroup]) -> Seq:
    """Removes XGroups from a sequence based on the sequence class.

    Args:
        seq (Seq): the sequence containing the "X" masked subsequences.
        seq_class (int): the class for the sequence based on the presence, number and
        length of its XGroups.
        xgroups (list[XGroup]): a list with all XGroups for a sequence.

    Returns:
        Seq: the sequence with trimmed subsequences based on its class.
    """

    # Sanitization check
    assert type(seq) == Seq, "Argument {seq} must be of type Seq"

    # Another sanitization check. If somehow the sequence has a nonvalid class
    # just return the untouched sequence
    if seq_class not in [1, 3, 6, 7]:
        return seq

    if seq_class == 1:
        return trim_x_only_keep_in_between(seq, xgroups=xgroups)
    elif seq_class == 3:
        return trim_x_and_5prime(seq, xgroups=xgroups)
    else:
        return trim_x_and_3prime(seq, xgroups=xgroups)


def trim_polynucleotide_subsequence(seq: Seq, al_positions: tuple[int]) -> Seq:
    start, end = al_positions

    filtered_seq = seq[: start - 1] + seq[end:]

    return filtered_seq


def trim_polynucleotides_by_dist_to_xgroups(estseq: ESTSeq) -> Seq:
    # To check if a polynucleotide subsequence is to be removed from seq
    # we need to check if it is at most 10 bases away from a XGroup.
    # MAYBE: check for the alignment z-score before trimming it
    if estseq.xgroup_list:
        for alignment in estseq.al_list:
            al_start, al_end = alignment.al_positions

            for xgroup in estseq.xgroup_list:
                x_start, x_end = xgroup.start, xgroup.end
                # This check is needed for us to see if the XGroup is
                # located before or after the alignment.
                lowest = min(al_start, x_start)

                # We then calculate the distance while aware of the
                # positions of the XGroup and the alignment relative
                # to each other.
                if lowest == al_start:
                    dist = x_start - al_end
                else:
                    dist = al_start - x_end

                # We will only trim the polynucleotide sequences if they are
                # at most 10 bases away from the XGroup.
                if dist <= 10:
                    # This check is needed because this could be the second
                    # time we try to trim the same sequence, and we want to do
                    # that on the updated sequence to keep it updated.
                    if estseq.processed_seq:
                        filtered_sequence = trim_polynucleotide_subsequence(
                            estseq.processed_seq, (al_start, al_end)
                        )
                    else:
                        filtered_sequence = trim_polynucleotide_subsequence(
                            estseq.masked_seq, (al_start, al_end)
                        )

                    estseq.set_processed_sequence(filtered_sequence)

                    # If we found a valid alignment and removed it from the sequence
                    # we don't need to check distances for that alignment anymore.
                    break
