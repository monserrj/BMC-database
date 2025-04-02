#!/usr/bin/env python

# This script allows for the creation of the database with data using db.py an data_addition.py. Instructions
# followed and explanations kept to help following through
from pathlib import Path

# Import the database file, csv and the data addition file to create
# and populate database
from db import create_db
from readfile import cli_open_csvfile
from file_and_data import link_db_csv
from data_addition import add_data

# Temporary code to delete local database while we debug
Path("bmc.db").unlink()

# The syntax here is "special"
# Every file in Python has a __name__ attribute
# If the file is being run as the main program, __name__ is set to "__main__"
# This happens if the file is run from the command line as
# python scripts/monolith.py
# If the file is being imported, __name__ is set to the name of the file
# This is a common Python idiom
if __name__ == "__main__":
    # Function to create database
    create_db()
    # Function to add data from csv file into bmc.db
    add_data()
    # Function to read csv file to import data
    cli_open_csvfile()
    # Function to link the data from the csv file to the database
    link_db_csv()
