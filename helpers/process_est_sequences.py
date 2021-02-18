"""
This is a simple script used to filter a set of EST sequences based on the presence
of masked vector sequences and aligned polyA and polyT sequences
"""
# Standard lib imports
import os
from typing import Union, Optional

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
        masked_seq: Union[Seq, None] = None,
    ) -> None:
        """Constructs an ESTSeq object with its sequence ID, the taxon it belongs to,
        its non-vector-masked sequence and possibly its vector-masked sequence if it
        exists. Calculates the sequence length and initializes empty values for the
        rest of the properties.

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
        self.masked_seq = masked_seq
        self.clean_seq = clean_seq

        # Calculated attributes
        self.seq_len = len(self.clean_seq)

        # Initialize default and empty values
        # Xgroups
        self.xgroup_list = []
        self.xgroup_count = len(self.xgroup_list)

        # Seq class
        self.seq_class = 0

        # Alignments
        self.al_list = []


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

        # Generate the list of AlignmentRecord objects
        alignments_list = swat_parser.SwatParser.parse_swat_allscores(
            all_scores_file_path, subject
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


# NOTE: maybe try to optimze this function
def create_estseq_list(
    taxon: str,
    clean_seqs: list[SeqRecord],
    masked_seqs: list[SeqRecord],
    alignments: list[Optional[swat_parser.AlignmentRecord]] = [],
) -> list[ESTSeq]:
    """Builds a list of ESTSeqs for each EST sequence in `clean_seqs` with informations
    about the vector-masked version of that EST (if it exists) and its alignments (if
    it has any).

    Args:
        taxon (str): the taxon to which the ESTs belong.
        clean_seqs (list[SeqRecord]): a list of non-vector-masked EST sequences.
        masked_seqs (list[SeqRecord]): a list of vector-masked EST sequences.
        alignments (list[Optional[swat_parser.AlignmentRecord]], optional): either a
        list of alignments of an EST sequence agains a set of subject sequences or an
        empty list if there are no alignments. Defaults to [].

    Returns:
        list[ESTSeq]: [description]
    """

    # A list to hold all built ESTSeq objects
    estseq_list = []

    # Create a new ESTSeq object for each clean SeqRecord object
    for clean_seq_record in clean_seqs:
        new_estseq = ESTSeq(
            seq_id=clean_seq_record.id, taxon=taxon, clean_seq=clean_seq_record.seq
        )

        for masked_seq_record in masked_seqs:
            # Find the masked SeqRecord that corresponds to the ESTSeq's id AND has an "X" in it
            if masked_seq_record.id == new_estseq.seq_id:
                # If the Seq of the masked SeqRecord contains an "X" in it
                if "X" in masked_seq_record:
                    # Then set the new ESTSeq's masked_seq to be that
                    new_estseq.masked_seq = masked_seq_record.seq

                # We found our masked seq and there's no need to search further
                # so we break
                break

        # Set the alignments list for the ESTSeq object
        if alignments:
            for alignment in alignments:
                # Find an AlignmentRecord that corresponds to our ESTSeq's id
                if alignment.id == new_estseq.seq_id:
                    # Append that AlignmentRecord to the ESTSeq's alignments list
                    new_estseq.al_list.append(alignment)

        # Finally, append that ESTSeq to the final list of ESTSeqs
        estseq_list.append(new_estseq)

    return estseq_list


# TODO: implement this
def update_xgroup_info(taxon: str, estseq_list: list[ESTSeq]) -> list[ESTSeq]:
    pass


def main():
    # The path to the common directory, from which we access
    # the subdirectories containing the sequences or alignments
    # NOTE: maybe make this cross-OS compatible with os.path.join
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

        # If the taxon is not "Polistes_canadensis", for which we don't have good
        # alignment information
        if taxon != taxa_list[1]:
            # Load all alignments for the taxon
            alignments = load_alignments(taxon, alignments_dir, subjects)
            # Create a list of ESTSeq objects
            raw_estseq_list = create_estseq_list(
                taxon, clean_seqs, masked_seqs, alignments
            )
        # Else if it is PolistesCanadensis, we create a list
        else:
            # Create a list of ESTSeq objects without alignment information
            raw_estseq_list = create_estseq_list(taxon, clean_seqs, masked_seqs)

        # TODO: implement this
        finished_estseq_list = update_xgroup_info(taxon, raw_estseq_list)


if __name__ == "__main__":
    main()
