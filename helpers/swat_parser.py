"""
A Python script to parse the results of an alignment produced by running `swat`.
It searches the given folder for a .alignment and a .allscores file. A .alignment
file can be produced by catching `swat`'s printed results to the STDOUT to a file.
A score value may be provided if you also want to use this to find the starting and
ending positions of an alignment.
"""


class Record:
    """Base Class for Record type objects"""

    def __init__(self) -> None:
        raise NotImplementedError


class Parser:
    """Base Class for Parser type objects"""

    def __init__(self) -> None:
        raise NotImplementedError


class AlignmentRecord(Record):
    """A Record containing information about its alignment against another subject
    sequence"""

    def __init__(self, id: str, raw_score: int, z_score: float, subject: str) -> None:
        """Instantiates an AlignmentRecord object with metrics describing its alignment
        score

        Args:
            id (str): the sequence identifier for the particular sequence that was
            aligned
            raw_score (int): the raw alignment score for the sequence
            z_score (float): the corrected alignment z-score for the sequence
            subject (str): the subject which this sequence was aligned against
        """
        self.id = id
        self.raw_score = raw_score
        self.z_score = z_score
        self.subject = subject

    def set_alignment_positions(self, start: int, end: int) -> None:
        """Retrieves the indexes corresponding to the starting and ending positions
        where the subject aligned with the sequence

        Args:
            start (int): the starting position of the alignment
            end (int): the ending position of the alignment
        """
        self.al_start = start
        self.al_end = end


class SwatParser(Parser):
    """A Parser for swat alignment output files, which can also calculate metrics such
    as the start and end position of a particular sequence's alignment
    """

    @staticmethod
    def parse_swat_results(
        allscores_file: str, subject: str = "None"
    ) -> list[AlignmentRecord]:
        result_records = []

        with open(allscores_file, "r") as scores_file:
            raw_records = [line.strip() for line in scores_file.readlines()]

            for record in raw_records:
                record = record.split("\t")
                result_records.append(
                    AlignmentRecord(
                        id=record[0],
                        raw_score=record[2],
                        z_score=record[3],
                        subject=subject,
                    )
                )

        return result_records
