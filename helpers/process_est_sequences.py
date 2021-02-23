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
import filter_seqs
import masked_seqs_stats
import swat_parser


class ESTSeq:
    """
    A class to hold information about a particular EST sequence, with properties
    regarding to its masked vector sequences, sequence length, to which taxon it
    belongs, alignment information, and many other relevant stats.
    """

    def __init__(
        self,
        seq_id: str,
        taxon: str,
        clean_seq: Seq,
    ) -> None:
        """Constructs an ESTSeq object with its sequence ID, the taxon it belongs to
        and its non-vector-masked sequence. Calculates the sequence length and
        initializes Null or empty values for the rest of the properties.

        Args:
            seq_id (str): the sequence identifier.
            taxon (str): the taxon of origin for a particular sequence.
            clean_seq (Seq): the non-masked sequence.
            masked_seq (Union[Seq, None], optional): The vector masked sequence if it
            exists, else None. Defaults to None.
        """

        # Setting the attributes according to the parameters
        self.seq_id = seq_id
        self.taxon = taxon
        self.clean_seq = clean_seq

        # The vector-masked sequence, should be set later with its setter
        self.masked_seq = None

        # The processed sequence, should be set later with its setter
        # This sequence is gonna change a lot over time because we need to
        # perform mutiple filtering and trimming steps on the original sequence
        self.processed_seq = None

        # The length of the original sequence
        self.seq_len = len(self.clean_seq)

        # Sequence class based on its XGroups
        self.seq_class = 0

        # List of XGroups, initalized as an empty list
        self.xgroup_list = []
        self.xgroup_count = 0

        # List of AlignmentRecords, initalized as an empty list
        self.al_list = []
        self.al_count = 0

    def set_xgroups(self, xgroup: masked_seqs_stats.XGroup) -> None:
        """A setter to update the list of XGroups for the sequence. Also updates
        the count of XGroups.

        Args:
            xgroup (masked_seqs_stats.XGroup): a XGroup object belonging to this
            sequence.
        """
        if xgroup:
            self.xgroup_list.append(xgroup)
            self.xgroup_count = len(self.xgroup_list)

    def set_alignments(self, alignment: swat_parser.AlignmentRecord) -> None:
        """A setter to update the list of AlignmentRecords for the sequence. Also
        updates the count of AlignmentRecords.

        Args:
            alignment (swat_parser.AlignmentRecord): an AlignmentRecord object
            belonging to this sequence.
        """
        if alignment:
            # The interface for this method gives back a list,
            # that's why we extend instead of appending.
            self.al_list.extend(alignment)
            self.al_count = len(self.al_list)

    def set_masked_seq(self, masked_seq: Seq) -> None:
        """A setter to update the vector-masked sequence for the ESTSeq object.

        Args:
            masked_seq (Seq): a vector-masked version of the sequence with "X"s in
            place of the vector sequence.
        """
        if masked_seq:
            self.masked_seq = masked_seq

    def set_processed_seq(self, processed_seq: Seq) -> None:
        """A setter to update the processed sequence for the ESTSeq object.

        Args:
            processed_seq (Seq): a processed version of the sequence with "X"s in
            place of the vector sequence.
        """
        if processed_seq:
            self.processed_seq = processed_seq

    def set_seq_class(self, seq_class: int) -> None:
        """A setter to update the sequence class for the ESTSeq object.

        Args:
            seq_class (int): an integer representing the sequence class.
        """
        if seq_class:
            self.seq_class = seq_class

    def __str__(self) -> str:
        """Defines a string representation of an ESTSeq object.

        Returns:
            str: the string representation of the ESTSeq. For example:

                Sequence ABC001 for Homo sapiens
                Sequence length: 100bp
                Sequence class: 1
                X-groups count: 2
                Alignments count: 2
                Alignments: [AlignmentRecord<polyA|ABC001|8.1|3>,
                AlignmentRecord<polyT|ABC001|9.3|1>]
        """
        return (
            f"Sequence {self.seq_id} for {self.taxon}\nSequence length: "
            f"{self.seq_len}bp\nSequence class: {self.seq_class}\nX-groups count: "
            f"{self.xgroup_count}\nAlignments count: {self.al_count}\nAlignments: "
            f"{self.al_list}"
        )

    def __repr__(self) -> str:
        """Defines a formal representation of an ESTSeq object.

        Returns:
            str: a more formal representation of the ESTSeq. For example:

                ESTSeq<Homo_sapiens|ABC001>
        """
        return f"ESTSeq<{self.taxon}|{self.seq_id}>"


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
                alignment.set_alignment_positions(positions)

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
