from helpers import swat_parser


def test_alignment_record_construction():
    record = swat_parser.AlignmentRecord(
        id="HX342487.1", raw_score=134, z_score=17.58, subject="polyT"
    )

    assert record.id == "HX342487.1"
    assert record.raw_score == 134
    assert record.z_score == 17.58
    assert record.subject == "polyT"

    record.set_alignment_positions(1, 2)

    assert record.al_start == 1
    assert record.al_end == 2
