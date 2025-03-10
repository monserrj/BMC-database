#!/usr/bin/env python

# This script allows for the reading of excel file for manual data addition to the database. Instructions
# followed and explanations kept to help following through

# Open .csv to be able to work with .csv files
import csv
# Import click to create special functions to interact with database
import click
# Import path to be able to open .csv files in different folders
from pathlib import Path

# Script below here:

# Create a click function
@click.command() # Declare function as a click command
# @click.version_option("0.1.0", prog_name="select_csvfile") # Defines the name and version of the click function
@click.argument("filepath",
    type=click.Path( # Define type of arg",
        exists=True, # makes sure the file exists
        file_okay=True, # makes sure the input path points to a file
        readable=True, # Make sure the content is readable
        path_type=Path, # Return the input into a path object
        ),
    ) # To define path/file as an argument and make Click treat any input as a path object.


# Define the click function
def cli_open_csvfile (filepath):
    """ Prompt to enter the path and filename of csv file and then read
    the data"""
    mydata = []
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
        # filepath = click.prompt("Enter the filepath to add data to your database:", type=click.File)
        # I do not need this as I want a more direct approach readfile.py --infile path/to/my_file.csv
        click.echo("Data read correctly")
    except Exception as e:
        click.echo(f"Error reading file: {e}")
    
    # Show first row for verification:
    click.echo(mydata[1])
    click.echo() # Add an extra line to the end of the output

# Provisional: Check the file selection and reading process


if __name__ == "__main__":
    cli_open_csvfile()