"""
This is a script that is able to find vector-masked subsequences (represented by "X"s
in place of the masked nucleotides) of a particular sequence and produce a collection
of those subsequences as XGroup objects with their corresponding length and distances
to the sequence ends, as well as a number representing a class for that sequence based
on the number, length and distance to ends of its XGroups, if there are any.
"""
# Standard lib imports
import re

# Third party imports
from Bio.SeqIO import SeqRecord


class XGroup:
    """
    An XGroup represents a vector-masked subsequence of a particular sequence. It is
    called XGroup because the subsequence is represented as a string o "X"s.
    """

    def __init__(self, xgroup_len: int, dist_from_5: int, dist_from_3: int) -> None:
        """Creates a XGroup object that corresponds to a continuous subsequence of "X"s
        belonging to a sequence.

        Args:
            xgroup_len (int): the length of the subsequence.
            dist_from_5 (int): the distance from the beginning of the sequence to the
            first "X" character of the XGroup.
            dist_from_3 (int): the distance from the last "X" character of the XGroup
            to the end of the sequence.
        """

        # Length of the sequence of "X"s
        self.xgroup_len = xgroup_len

        # The distances from 5' and 3'
        self.dist_from_5 = dist_from_5
        self.dist_from_3 = dist_from_3

    def __repr__(self) -> str:
        """A formal representation of an XGroup object.

        Returns:
            str: a more formal representation of the XGroup. For example:

                XGroup<50|4|10>
                XGroup<221|0|424>
        """
        return f"XGroup<{self.xgroup_len}|{self.dist_from_5}|{self.dist_from_3}>"


def find_x_regions_and_calculate_stats(sequence: SeqRecord) -> tuple[int, list]:
    """Takes a nucleotide sequence with masked regions represented by "X"s and
    calculates the length of each X region, the distances of each X region to the 3'
    and 5' ends of the sequence and the sequence class based on X region features.

    Args:
        sequence (SeqRecord): a SeqRecord object representing a nucleotide sequence.

    Returns:
        tuple[int, list]: a tuple containing the sequence class as an integer from 0 to
        7 (by default, returns 0) and a list XGroup objects for each vector-masked
        subsequence in the sequence if there are any, otherwise returns an empty list.
    """

    # Search the string version of the sequence for any continous occurences of "X"s
    # of any length. With finditer() we get all the possible occurences.
    pattern = re.compile(r"X+")
    search_res = pattern.finditer(str(sequence.seq))

    seq_class, xgroups = (0, [])

    # Iterates through each subsequence of "X"s based on their positions relative
    # to the original sequence.
    for match in search_res:
        # span() returns a tuple with the start and end position of the match
        # relative to the original string.
        start, end = match.span()
        # Length of the "X" subsequence.
        xgroup_len = len(match.group())
        # The distance to 5' is just start, which is the index of the first "X".
        # In turn, we need to calculate the distance from 3' by subtracting the length
        # of the original sequence by the end position of the "X" subsequence, which is
        # the index of the last "X".
        dist_from_3 = len(sequence.seq) - end

        # Create a new XGroup object with the information collected above
        new_xgroup = XGroup(
            xgroup_len=xgroup_len, dist_from_5=start, dist_from_3=dist_from_3
        )
        # Add it to the list of XGroups
        xgroups.append(new_xgroup)

    # Only calculate the sequence class if there is any XGroup at all. If there
    # are no XGroups, leave the sequence class as 0, which will be treated as an
    # invalid class and will be ignored later.
    if len(xgroups) > 0:
        seq_class = get_seq_class(xgroups)

    return (seq_class, xgroups)


def get_seq_class(x_groups_list: list[XGroup]) -> int:
    """Infers the sequence class based on the number, length and position of its
    XGroups, if there are any.

    Args:
        x_groups_list (list[XGroup]): a list containing all XGroups from a particular
        sequence.

    Returns:
        int: the sequence class.
    """

    # Sanity check
    assert len(x_groups_list) > 0, "The list contains no XGroups"

    # Calculate how many XGroups there are for this sequence.
    num_of_xgroups = len(x_groups_list)

    # Infer the class for the sequence based on its XGroups.
    # Classes 1 and 2 only care about the number of XGroups in a sequence.
    if num_of_xgroups == 2:
        return 1
    elif num_of_xgroups > 2:
        return 2
    else:
        # Classes 3 to 7 also depend on the length of the XGropus and their distance
        # to the ends of the sequence. They all assume there is only one XGroup in
        # the sequence.
        x_group = x_groups_list[0]
        xgrop_len = x_group.xgroup_len
        dist_from_5 = x_group.dist_from_5

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
