"""Tests of the protein table.

These tests should be run from the repository root using pytest -v
"""

from pathlib import Path

import pytest
import subprocess

from scripts import db, file_and_data, readfile

@pytest.fixture
def dbpath(request) -> Path:
    """Path to temporary database location
    
    Remember to unlink when done (if not needed for further actions)
    """
    return request.path.parent / "fixtures" / "tmp" / "new.sqlite"

@pytest.fixture
def prot_minimal(request) -> Path:
    """Path to minimal protein data"""
    return request.path.parent / "fixtures" / "protein" / "prot_data_minimal_correct.csv"

@pytest.fixture
def db_info(request) -> Path:
    """Path to external database information"""
    return request.path.parent / "fixtures" / "ext_db" / "db_info_ext.csv"

@pytest.fixture
def output_add_minimal_protein_db(request) -> str:
    """Path to minimal protein data database contents"""
    return (request.path.parent / "fixtures" / "outputs" / "protein" / "add_minimal_protein_data.sql").open().read()


def test_make_new_db(dbpath: Path) -> None:
    """Confirm can create (and unlink) a new empty database."""
    assert not dbpath.is_file()
    db.create_db(dbpath)
    assert dbpath.is_file()
    with dbpath.open("rb") as handle:
        magic = handle.read(16)
        assert magic == b"SQLite format 3\0"
    dbpath.unlink()

def test_add_minimal_protein_data(dbpath: Path, prot_minimal: Path, db_info: Path, output_add_minimal_protein_db) -> None:
    """Add protein data to empty database"""
    # Create database and attach session
    assert not dbpath.is_file()
    db.create_db(dbpath)
    assert dbpath.is_file()
    session = db.get_session(dbpath)

    # Read and add protein data
    data = readfile.read_file(prot_minimal, verbose=False)
    file_and_data.link_db_csv(data, session, db_info)

    # Dump database and capture output
    result = subprocess.run(["sqlite3", str(dbpath), ".dump"], capture_output=True)
    assert result.stdout.decode("utf-8") == output_add_minimal_protein_db
   
    # Unlink database
    dbpath.unlink()