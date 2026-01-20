#!/usr/bin/env python

"""This script links the addition of data manually to the BMC using csv from readfile to database described in db.
Following the data addition structure from data_addition. Instructions
followed and explanations kept to help following through

To be able to make exceptions in code (try/except):
from sqlalchemy.exc import IntegrityError, PendingRollbackError
For exiting system for trouble shooting import sys
"""

from db import (
    protein_addition,
    #name_addition,
    #cds_addition,
    #xref_addition,
    #xdatabase_addition,
)

# To make boolean values from the csv true booleans:
def parse_bool(value):
    """Convert CSV boolean-like values to Python bool."""
    if isinstance(value, bool):
        return value  # already boolean
    # Accept TRUE/FALSE (case-insensitive) or yes/no
    return str(value).strip().lower() in ("true", "yes")

# Making a function for the addition of data to the different tables created with db.py following data_instructions:
def link_db_csv(mydata, session):
    """Function for each data type addition created
    Add the data to the database in a loop, but we'll have to
    check if the data entered already exist
    and update the corresponding tables accordingly.
    We need to check if the sequence, structure, and accession already exist,
    and update the corresponding tables accordingly if they do not.
    We can then update the linker tables by adding the corresponding items.
    """

    # Can I make this bit more efficient with a loop?
    for idx, row in enumerate(mydata):
        (
            protseq,
            protacc,
            struct,
            canonical,
            protname,
            # external dbs come next in sets of 4 (dbname, url, acc, type)
            *db_refs,
        ) = row
        # Check what data is available:
        print(f"\nStarting next loop ({idx=}) with {protacc=}")
        # Convert canonical to boolean
        canonical = parse_bool(canonical)

        try:  # attempt to add to db, if we raise an error, we roll back
            # Add protein data
            protein = protein_addition(
                session,
                protseq=protseq,
                protacc=protacc,
                struct=struct,
                canonical=canonical,
            )
            print(f"\nProtein record returned: {protein}")
            
            # gene = cds_addition(
            #     session,
            #     cdstype=cdstype,
            #     dnaseq=dnaseq,
            #     protein=protein,
            # )
            # print(f"\nProtein record returned: {protein}")
            # print(f"\nGene record returned: {gene}")

            # Add name data
            # name = name_addition(session, protname=protname, protein=protein)

            # print(f"\nProtein record returned: {protein}")
            # print(f"Name record returned: {name}")

            # Add external reference data
            # Handle external DB refs
            # They come in groups of 4: db, url, acc, type
        #     for i in range(0, len(db_refs), 4):
        #         dbname, dburl, accessionid, dbtype = db_refs[i : i + 4]

        #         # skip empty slots if some rows don’t have all 7 filled
        #         if not any([dbname, dburl, accessionid, dbtype]):
        #             continue

        #         xref = xref_addition(
        #             session,
        #             accessionid=accessionid,
        #             dbname=dbname,
        #             dburl=dburl,
        #             protein=protein,
        #         )
        #     print(f"\nProtein record returned: {protein}")
        #     print(f"Name record returned: {xref}")
            
            # Commit the transaction after successful addition
            print("Committing the session")
            session.commit()

        except Exception as exc:
            print(f"Error adding data: {exc}")
            print("Rolling back changes and skipping to next entry")
            session.rollback()
