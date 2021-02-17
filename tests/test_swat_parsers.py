import os
from helpers import swat_parser


# Constants
TEST_FOLDER = os.path.dirname(__file__)


def test_alignment_record_construction():
    """Test the construction and position assignments of an AlignmentRecord object"""

    # Initializing the object with example values from a real alignment
    record = swat_parser.AlignmentRecord(
        id="HX342487.1", raw_score=134, z_score=17.58, subject="polyT"
    )

    # Testing the objects initial properties
    assert record.id == "HX342487.1"
    assert record.raw_score == 134
    assert record.z_score == 17.58
    assert record.subject == "polyT"

    # Setting the alignment starting and ending positions
    record.set_alignment_positions(start=1, end=2)

    # Testing if we are able to set the alignment positions
    assert record.al_start == 1
    assert record.al_end == 2


def test_parsing_allscores_file():
    """Test the parsing and object creation from a mock `swat` .allscores file"""

    # Path to the input file
    file_path = f"{TEST_FOLDER}/swat_test.allscores"

    # Parse the file and get the list of AlignmentRecord objects
    alignments = swat_parser.SwatParser.parse_swat_allscores(
        allscores_file=file_path,
        subject="polyA",
    )

    # Checking if the expected number of records was produced
    assert len(alignments) == 2

    # Getting both records out of the list
    record1 = alignments[0]
    record2 = alignments[1]

    # Checking if the records are instances of the AlignmentRecord class
    assert isinstance(record1, swat_parser.AlignmentRecord)
    assert isinstance(record2, swat_parser.AlignmentRecord)

    # Testing each record's properties
    assert record1.id == "HX317466.1"
    assert record2.id == "HX287430.1"

    assert record1.raw_score == 148
    assert record2.raw_score == 132

    assert record1.z_score == 30.81
    assert record2.z_score == 27.24

    assert record1.subject == "polyA"
    assert record2.subject == "polyA"


def test_parsing_alignments_file():
    pass