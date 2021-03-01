"""
This script removes certain subsequences from a sequence. For example, vector-masked
subsequences represented by "X"s or polynucleotide subsequences such as poly-A and
poly-T tails.
"""

# Standard lib imports
from copy import deepcopy
from typing import Optional

# Third-party imports
from Bio.Seq import Seq

# Own module imports
from .masked_seqs_stats import XGroup
from .estseq import ESTSeq


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


def trim_subsequence(seq: Seq, positions: tuple[int, int]) -> Seq:
    """Trims a subsequence from a sequence, returning a new Seq object containing
    the sequence without the trimmed region.

    Args:
        seq (Seq): the original sequence to be trimmed.
        positions (tuple[int, int]): a 2-tuple with (start, end) positions of the
        subsequence to be removed. Note that it is assumed that the positions are
        1-indexed instead of 0-indexed.

    Returns:
        Seq: a new Seq object containing the original sequence without the trimmed
        region.
    """

    start, end = positions

    filtered_seq = seq[: start - 1] + seq[end:]

    return filtered_seq


def trim_polynucleotides_by_dist_to_xgroups(
    estseq: ESTSeq, max_dist: int = 10, z_cutoff: float = 8.0, inplace: bool = False
) -> Optional[ESTSeq]:
    """Wrapper to allow the trimming of subsequences aligned to a sequence based
    on their distance to XGroup regions and their z_score.

    Args:
        estseq (ESTSeq): the ESTSeq object with the original sequences, XGroup and
        alignment information.
        max_dist (int, optional): the maximum allowed distance, in nucleotides, from
        the alignment region to any XGroup region. Defaults to 10.
        z_cutoff (float, optional): the minimum z-score value for the alignment to be
        eligible for trimming. Defaults to 8.0.
        inplace (bool, optional): whether the changes should happen in place or a new
        ESTSeq should be returned. Defaults to False.

    Returns:
        Optional[ESTSeq]: either a new ESTSeq object with a `processed_seq` property
        without its trimmed alignment regions, or None.
    """

    # If we don't want the changes to happen inplace, copy the object.
    if not inplace:
        estseq = deepcopy(estseq)

    # To check if a polynucleotide subsequence is to be removed from seq
    # we need to check if it is at most `max_dist` bases away from a XGroup.
    if estseq.xgroup_list:
        for alignment in estseq.al_list:
            # Do a quick check to only allow trimming alignments that have a z-score
            # of at least `z_cutoff`.
            if alignment.z_score < z_cutoff:
                continue

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
                # at most `max_dist` bases away from the XGroup.
                if dist <= max_dist:
                    # This check is needed because this could be the second
                    # time we try to trim the same sequence, and we want to do
                    # that on the updated sequence to keep it updated.
                    if estseq.processed_seq:
                        filtered_sequence = trim_subsequence(
                            estseq.processed_seq, (al_start, al_end)
                        )
                    else:
                        filtered_sequence = trim_subsequence(
                            estseq.masked_seq, (al_start, al_end)
                        )

                    estseq.set_processed_seq(filtered_sequence)

                    # If we found a valid alignment and removed it from the sequence
                    # we don't need to check distances for that alignment anymore.
                    break

    # If we don't want the changes to happen inplace, we need
    # to return the new version of the estseq.
    if not inplace:
        return estseq


# TODO: Document this
def trim_polynucleotides_by_dist_to_ends(
    estseq: ESTSeq, max_dist: int = 20, z_cutoff: float = 30.0, inplace: bool = False
) -> Optional[ESTSeq]:

    # If we don't want the changes to happen inplace, copy the object.
    if not inplace:
        estseq = deepcopy(estseq)

    for alignment in estseq.al_list:
        # If the z-score is too low or if there are no alignment positions
        # we don't care about this alignment.
        if alignment.z_score < z_cutoff or not alignment.al_positions:
            continue

        for position in alignment.al_positions:
            start, end = position

            # The sequence length might be different between the clean
            # and the processed seq. We need to calculate the distance based
            # on the length of the sequence we're operating on.
            if estseq.processed_seq:
                seq_len = len(estseq.processed_seq)
            else:
                seq_len = estseq.seq_len

            dist_to_5 = start
            dist_to_3 = seq_len - end

            if dist_to_5 < max_dist or dist_to_3 < max_dist:
                if estseq.processed_seq:
                    filtered_sequence = trim_subsequence(
                        estseq.processed_seq, (start, end)
                    )
                else:
                    filtered_sequence = trim_subsequence(estseq.clean_seq, (start, end))

                estseq.set_processed_seq(filtered_sequence)

    # If we don't want the changes to happen inplace, we need
    # to return the new version of the estseq.
    if not inplace:
        return estseq
