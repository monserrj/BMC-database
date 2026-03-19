"""Tests of the protein table.

These tests should be run from the repository root using pytest -v
"""

from pathlib import Path

import pytest

from scripts import db

@pytest.fixture
def dbpath(request):
    """Path to temporary database location
    
    Remember to unlink when done (if not needed for further actions)
    """
    return request.path.parent / Path("fixtures") / "tmp" / "new.sqlite"

def test_make_new_db(dbpath: Path) -> None:
    """Confirm can create (and unlink) a new empty database."""
    assert not dbpath.is_file()
    db.create_db(dbpath)
    assert dbpath.is_file()
    with dbpath.open("rb") as handle:
        magic = handle.read(16)
        assert magic == b"SQLite format 3\0"
    dbpath.unlink()



