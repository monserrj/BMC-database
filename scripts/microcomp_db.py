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
from readfile import cli_open_csvfile
from file_and_data import link_db_csv
from data_addition import add_data


@click.command()
@click.option(
    "--dbpath",
    default=Path("bmc.db"),
    type=click.Path(exists=False, path_type=Path),
)
@click.option(
    "--force/--no-force",
    default=False,
)
def main(dbpath: Path, force: bool):
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
    print(session)
    sys.exit()

    # Add data from CSV file into bmc.db
    add_data()

    # Read CSV file to import data
    cli_open_csvfile()
    # Link the data from the CSV file to the database
    link_db_csv()


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
