from helpers import swat_parser


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
