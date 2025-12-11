#!/usr/bin/env python

"""microcomp_db.py

This script is the main entry point for the microcomp_db project. It is responsible for creating the database, reading data from a CSV file, and linking the data to the database.

We'll use click to manage the command line interface and SQLAlchemy to manage the database connection and operations.
"""

# This script allows for the creation of the database with data using db.py an data_addition.py. Instructions
# followed and explanations kept to help following through
import sys
from pathlib import Path

import click

# Import the database file, csv and the data addition file to create
# and populate database
from db import create_db, get_session
from readfile import read_file
from file_and_data import link_db_csv


@click.command()
@click.option(
    "--dbpath",
    default=Path("bmc.db"),
    type=click.Path(exists=False, path_type=Path),
)
@click.option(
    "--csvpath",
    default=Path("../data/raw/prot_info/prot_data_minimal_correct_UP.csv"),
    type=click.Path(exists=False, path_type=Path),
)
@click.option(
    "--force/--no-force",
    default=False,
)
@click.option(
    "--verbose/--no-verbose",
    default=False,
)
def main(dbpath: Path, csvpath: Path, force: bool, verbose: bool):
    print(dbpath)
    """Main function to run the script."""
    # Create the database if it doesn't exist, or we're forcing overwrite
    if force is True and dbpath.exists():  # overwrite database
        print(f"Overwriting {dbpath}")
    elif dbpath.exists():
        print(f"Not overwriting {dbpath} (exiting)")
        sys.exit(0)
    create_db(dbpath)
    session = get_session(dbpath)

    # Read CSV file to import data
    # print(f"{csvpath=}")
    data = read_file(csvpath, verbose)

    # Link the data from the CSV file to the database
    link_db_csv(data, session)


# The syntax here is "special"
# Every file in Python has a __name__ attribute
# If the file is being run as the main program, __name__ is set to "__main__"
# This happens if the file is run from the command line as
# python scripts/monolith.py
# If the file is being imported, __name__ is set to the name of the file
# This is a common Python idiom
if __name__ == "__main__":
    # Call the main function
    main()
