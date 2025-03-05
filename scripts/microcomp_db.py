#!/usr/bin/env python

# This script allows for the creation of the database with data using db.py an data_addition.py. Instructions
# followed and explanations kept to help following through

# Open .csv to be able to work with .csv files
import csv

# Import path to be able to open .csv files in different folders
from pathlib import Path

# Import a sessionmaker to create a session object
from sqlalchemy.orm import sessionmaker

# To be able to make exceptions in code (try/except):
# from sqlalchemy.exc import IntegrityError, PendingRollbackError
# For exiting system for trouble shooting import sys
import sys

# Import the database file and the data addition file to create
# and populate database
from db import create_db

from readfile import select_csvfile, read_csvfile


def main ():
    create_db()
    select_csvfile()
    read_csvfile()
    load_data()


# The syntax here is "special"
# Every file in Python has a __name__ attribute
# If the file is being run as the main program, __name__ is set to "__main__"
# This happens if the file is run from the command line as
# python scripts/monolith.py
# If the file is being imported, __name__ is set to the name of the file
# This is a common Python idiom
if __name__ == "__main__":
    main()


    # Add the data to the database
    for idx, (
        protseq,
        NCBIid,
        uniprot,
        struct,
        genename,
        namerank,
        dnaseq,
        taxref,
        taxdb,
        spec,
        genu,
        fam,
        order,
        phyl,
        classt,
        stra,
        pdbid,
        pdb_1,
        pdb_2,
        pdb_3,
        pathid,
        KOid,
    ) in enumerate(mydata):
        
        # Check what data is available:
        # print(f"This is before adding session.query {protseq=}, {NCBIid=},{uniprot=}, {struct=}")
        print(f"\nStarting next loop ({idx=}) with {uniprot=}, {KOid=}")

        try:  # attempt to add to db, if we raise an error, we roll back
            # Add protein data
            protein = protein_addition(
                session,
                protseq=protseq,
                NCBIid=NCBIid,
                uniprot=uniprot,
                struct=struct,
                dnaseq=dnaseq,
            )
            print(f"\nProtein record returned: {protein}")
        
            # Add name data
            name = name_addition(
                session,
                genename=genename,
                namerank=namerank,
                protein = protein
            )
        
            print(f"\nProtein record returned: {protein}")
            print(f"Name record returned: {name}")
        
            # Add taxonomy data
            # We expect a single tax_id per protein, so if the protein
            # already is linked to a taxonomy, we skip adding the taxonomy
            if not len(protein.taxonomies):
                taxonomy = taxonomy_addition(
                    session,
                    taxref=taxref,
                    taxdb = taxdb,
                    spec = spec,
                    genu = genu,
                    fam = fam,
                    order = order,
                    phyl = phyl,
                    classt = classt,
                    stra = stra,
                    protein = protein
                )
                print(f"\nProtein record returned: {protein}")
                print(f"Name record returned: {name}")
                print(f"Taxonomy record returned: {taxonomy}")
            session.commit()

        except Exception as exc:
            print(f"Error adding data: {exc}")
            print("Rolling back changes and skipping to next entry")
            session.rollback()