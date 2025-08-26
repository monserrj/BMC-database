#!/usr/bin/env python

'''This script allows for the reading of excel file for manual data addition to the database. Instructions
followed and explanations kept to help following through
'''

# Open .csv to be able to work with .csv files
import csv

# Import click to create special functions to interact with database
import click

# Import path to be able to open .csv files in different folders
from pathlib import Path


def read_file(filepath: Path, verbose: bool):
    mydata = []
    # print(f"{filepath=}")
    try:
        with open(filepath, newline="", encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            next(reader)  # Skip header row

            for row in reader:
                # Convert empty strings to Null
                row = [None if val == "" else val for val in row]
                if (
                    row[0] is not None
                ):  # LP: we could do with more sanity checking of data here
                    mydata.append(tuple(row))
                else:
                    click.echo("No sequence data, discarding")
    except Exception as e:
        click.echo(f"Error reading file: {e}")
        raise click.Abort  # Aborts with error

    # Show first row for verification:
    if verbose:
        click.echo(str(mydata[0]) + "\n")
        click.echo("Data read correctly")
    return mydata


# Create a click function
@click.command()  # Declare function as a click command
# @click.version_option("0.1.0", prog_name="select_csvfile")  # Defines the name and version of the click function
@click.argument(
    "filepath",
    type=click.Path(  # Define type of arg",
        exists=True,  # makes sure the file exists
        file_okay=True,  # makes sure the input path points to a file
        readable=True,  # Make sure the content is readable
        path_type=Path,  # Return the input into a path object
    ),
)  # To define path/file as an argument and make Click treat any input as a path object.
@click.option("--verbose", "-v", is_flag=True, help="Provide more informative output.")
# Define the click function
def cli_open_csvfile(filepath: Path, verbose: bool):
    """Prompt to enter the path and filename of csv file and then read
    the data"""
    return read_file(filepath, verbose)


if __name__ == "__main__":
    cli_open_csvfile()
