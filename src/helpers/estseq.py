"""
A class to represent an EST sequence with features like its vector-masked sequences
and alignments.
"""


from Bio.Seq import Seq

from . import masked_seqs_stats
from . import swat_parser


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

    def set_alignments(self, alignment: list[swat_parser.AlignmentRecord]) -> None:
        """A setter to update the list of AlignmentRecords for the sequence. Also
        updates the count of AlignmentRecords.

        Args:
            alignment (list[swat_parser.AlignmentRecord]): a list of AlignmentRecord
            objects belonging to this sequence.
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

    def clear_alignments(self) -> None:
        """A helper method to reset the alignments of an ESTSeq object."""

        self.al_list = []
        self.al_count = 0

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
