"""
This is a simple script used to filter a set of EST sequences based on the presence
of masked vector sequences and aligned polyA and polyT sequences
"""
# Standard lib imports
import os
from typing import Optional
from copy import deepcopy
from collections import defaultdict

# Third-party imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Own modules
from helpers import filter_seqs
from helpers import masked_seqs_stats
from helpers import swat_parser
from helpers.estseq import ESTSeq


def load_seqs(taxon: str, path: str, extension: str = ".fasta") -> list[SeqRecord]:
    """Loads a set o SeqRecords for a taxon from a file in path that contains the taxon
    name and returns it as a list of SeqRecords.

    Args:
        taxon (str): the name of the taxon to look for.
        path (str): the path to the folder where the sequence file is.
        extension (str): the expected file extension for the file to load.

    Returns:
        list[SeqRecord]: a list of all SeqRecords in the file for a particular taxon.
    """
    # Find the seq file for that taxon
    for seq_file in os.listdir(path):
        if taxon in seq_file and seq_file.endswith(extension):
            # Generate the full path to the seq file
            seq_file_path = os.path.join(path, seq_file)
            # Return the seq list
            return list(SeqIO.parse(seq_file_path, "fasta"))


def load_alignments(
    taxon: str, path: str, subjects: list[str]
) -> list[swat_parser.AlignmentRecord]:
    """Generates a list of AlignmentRecord objects for all possible alignments against
    all given `subjects` for each sequence belonging to a `taxon`.

    Args:
        taxon (str): the name of the taxon to look for.
        path (str): the path to all alignment and allscores files.
        subjects (list[str]): all the subject sequences with which each sequence was
        aligned against.

    Returns:
        list[swat_parser.AlignmentRecord]: a list containing an AlignmentRecord object
        for each sequence for the `taxon`.
    """
    # A list to hold all AlignmentRecord objects for each sequence
    # of a particular taxon.
    alignments_list = []

    # Find the file for each subject sequence
    for subject in subjects:
        # Find the alignment files for a taxon
        for al_file in os.listdir(path):
            # Check if the file belongs to that taxon and that subject
            if taxon in al_file and subject in al_file:
                # Get the path to the allscores and alignments files for
                # each taxon, for each subject
                if al_file.endswith(".allscores"):
                    all_scores_file_path = os.path.join(path, al_file)
                else:
                    alignments_file_path = os.path.join(path, al_file)

        # Update the list of AlignmentRecord objects
        alignments_list.extend(
            swat_parser.SwatParser.parse_swat_allscores(all_scores_file_path, subject)
        )

        # Generate a dict of alignment positions for each sequence
        positions_dict = swat_parser.SwatParser.parse_swat_alignment_output(
            alignments_file_path
        )

        # Set the alignment positions for each AlignmentRecord
        for alignment in alignments_list:
            for positions in positions_dict[alignment.id]:
                *pos, score, z = positions

                # Check if we're actually dealing with the correct alignment
                # by also comparing the scores, to help us remove ambiguity.
                if alignment.raw_score == score and alignment.z_score == z:
                    alignment.set_alignment_positions(tuple(pos))

    return alignments_list


def create_estseq_list(taxon: str, clean_seqs: list[SeqRecord]) -> list[ESTSeq]:
    """Builds a list of ESTSeqs for each EST sequence in `clean_seqs`.

    Args:
        taxon (str): the taxon to which the ESTs belong.
        clean_seqs (list[SeqRecord]): a list of non-vector-masked EST sequences.

    Returns:
        list[ESTSeq]: a list containing an ESTSeq for each sequence in `clean_seqs`.
    """

    # Create a new ESTSeq object for each clean SeqRecord object
    estseq_list = [
        ESTSeq(seq_id=clean_seq_record.id, taxon=taxon, clean_seq=clean_seq_record.seq)
        for clean_seq_record in clean_seqs
    ]

    return estseq_list


def set_masked_seqs_for_ESTSeqs(
    seqs_list: list[ESTSeq], masked_seqs: list[SeqRecord], inplace: bool = False
) -> Optional[list[ESTSeq]]:
    """Finds the matching masked sequence for an ESTSeq object, if there's any, and
    set the objects masked_seq attribute to it.

    Args:
        seqs_list (list[ESTSeq]): a list of ESTSeq objects.
        masked_seqs (list[SeqRecord]): a list of vector-masked SeqRecords.
        inplace (bool, optional): `True` means the function will mutate the original
        list passed to `seqs_list`. `False` means it returns a new list. Defaults to
        False.

    Returns:
        Optional[list[ESTSeq]]: if `inplace` is `False`, returns the new list of ESTSeq
        objects with information about their corresponding masked sequences. If it is
        `True`, returns `None`.
    """

    assert len(masked_seqs) > 0, "No masked sequences for this taxon. Check your data."

    # If we don't wanna mutate the original list and objects,
    # we operate on a copy of it.
    if not inplace:
        seqs_list = deepcopy(seqs_list)

    # Make a dictionary from the list of masked_seqs,
    # where each id is the key and the sequence itself is the value.
    masked_seqs_map = {
        masked_seq_record.id: masked_seq_record.seq for masked_seq_record in masked_seqs
    }

    # Set the corresponding masked sequence for each ESTSeq object, if there is one.
    for estseq in seqs_list:
        estseq.set_masked_seq(masked_seqs_map.get(estseq.seq_id, None))

    # If we're not mutating the list inplace, we need to return the new list
    if not inplace:
        return seqs_list


def set_alignments_for_ESTSeqs(
    seqs_list: list[ESTSeq],
    alignments: list[swat_parser.AlignmentRecord],
    inplace: bool = False,
) -> Optional[list[ESTSeq]]:

    assert len(alignments) > 0, "No alignments found for this taxon. Check your data."

    # If we don't wanna mutate the original list and objects,
    # we operate on a copy of it.
    if not inplace:
        seqs_list = deepcopy(seqs_list)

    # Make a dictionary from the list of AlignmentRecord objects,
    # where each id is the key and the object itself is the value.
    # The best way to account for multiple entries for the same key is to
    # use a defaultdict of lists. That's why we don't use a comprehension here.
    alignments_map = defaultdict(list)
    for alignment in alignments:
        alignments_map[alignment.id].append(alignment)

    for estseq in seqs_list:
        estseq.set_alignments(alignments_map.get(estseq.seq_id, None))

    # If we're not mutating the list inplace, we need to return the new list
    if not inplace:
        return seqs_list


def set_xgroups_for_ESTSeqs(
    seqs_list: list[ESTSeq], inplace: bool = False
) -> Optional[list[ESTSeq]]:
    """Finds the XGroups belonging to a vector-masked sequence and append them to its
    list of XGroups, while also calculating and updating its sequence class.

    Args:
        seqs_list (list[ESTSeq]): a list of ESTSeq objects.
        inplace (bool, optional): `True` means the function will mutate the original
        list passed to `seqs_list`. `False` means it returns a new list. Defaults to
        False.

    Returns:
        Optional[list[ESTSeq]]: if `inplace` is `False`, returns the new list of ESTSeq
        objects with information about their corresponding masked sequences. If it is
        `True`, returns `None`.
    """

    # If we don't wanna mutate the original list and objects,
    # we operate on a copy of it.
    if not inplace:
        seqs_list = deepcopy(seqs_list)

    for estseq in seqs_list:
        if estseq.masked_seq:
            # Get the class and list of XGroups for the ESTSeq
            seq_class, xgroups = masked_seqs_stats.find_x_regions_and_calculate_stats(
                estseq.masked_seq
            )

            estseq.set_seq_class(seq_class)

            # Add each of the XGroups to the ESTSeq's list of XGroups
            for xgroup in xgroups:
                estseq.set_xgroups(xgroup)

    # If we're not mutating the list inplace, we need to return the new list
    if not inplace:
        return seqs_list


def main() -> None:
    # The path to the common directory, from which we access
    # the subdirectories containing the sequences or alignments
    # MAYBE: make this cross-OS compatible with os.path.join
    common_root_dir = "/home/maycon/Documents/LAB/eusociality/local_data"

    # The subdirs
    alignments_dir = os.path.join(common_root_dir, "alignments", "2021-02-17")
    clean_seqs_dir = os.path.join(common_root_dir, "blast_removal_results")
    masked_seqs_dir = os.path.join(common_root_dir, "crossmatch_filtered")

    # List our taxa
    taxa_list = ["Apis_mellifera", "Polistes_canadensis", "Solenopsis_invicta"]

    # Subject sequences we align against
    subjects = ["polyA", "polyT"]

    # Loop over each taxon
    for taxon in taxa_list:
        # Load the clean and masked sequences
        clean_seqs = load_seqs(taxon, clean_seqs_dir)
        masked_seqs = load_seqs(taxon, masked_seqs_dir, extension=".screen")

        # Filter the masked seqs to make sure all sequences contain "X"s
        masked_seqs = [
            masked_seq for masked_seq in masked_seqs if "X" in masked_seq.seq
        ]

        # Create a list of ESTSeq objects
        estseq_list = create_estseq_list(taxon, clean_seqs)

        # Add the vector-masked sequences for each ESTSeq that have one
        set_masked_seqs_for_ESTSeqs(
            seqs_list=estseq_list, masked_seqs=masked_seqs, inplace=True
        )

        # Add the information about the XGroups for the ESTSeqs
        set_xgroups_for_ESTSeqs(seqs_list=estseq_list, inplace=True)

        # If the taxon is not "Polistes_canadensis", for which we don't have good
        # alignment information, we don't try to load or assign its alignments.
        # MAYBE: find another way around this.
        if taxon != taxa_list[1]:
            # Load all alignments for the taxon
            alignments = load_alignments(taxon, alignments_dir, subjects)
            # Set the alignments for each ESTSeq
            set_alignments_for_ESTSeqs(estseq_list, alignments=alignments, inplace=True)

        # TODO: filter the list of ESTSeqs to remove class 4 ESTSeqs


if __name__ == "__main__":
    main()
