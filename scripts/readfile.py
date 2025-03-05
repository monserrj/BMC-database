#!/usr/bin/env python

# This script allows for the reading of excel file for manual data addition to the database. Instructions
# followed and explanations kept to help following through

# Open .csv to be able to work with .csv files
import csv

# Import path to be able to open .csv files in different folders
from pathlib import Path

# For exiting system for trouble shooting import sys
import sys

# Script below here:
def select_csvfile ():
    """ Prompt to enter the path and filename of csv file"""
    file_path = input("Enter the full path of your csv file: ").strip()
    
    # Convert the string path to a Path object
    file_path = Path(file_path)
    # Check if the file exists
    if not file_path.exists():
        print(f"Error: The file '{file_path}' does not exist. Please check the path and try again.")
        return None
    return file_path

def read_csvfile (file_path):
    """Read data from the specified CSV file"""
    mydata = []
    try:
        with open(file_path, newline="", encoding="utf-8") as csvfile:
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
                    print("No sequence data, discarding")
        return mydata
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found")
        return None
    except Exception as e:
        print(f"Error {e} while reading the file")
        return None

# Provisional: Check the file selection and reading process
def main ():
    read_csvfile()
    
file_path = select_csvfile()
if file_path:
    data = read_csvfile(file_path)
    # Show first row for verification
    print("File successfully read. Data preview:", data[:1])