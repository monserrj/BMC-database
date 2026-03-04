#!/usr/bin/env python

"""This script links the addition of data manually to the BMC using csv from readfile to database described in db.
Following the data addition structure from data_addition. Instructions
followed and explanations kept to help following through

To be able to make exceptions in code (try/except):
from sqlalchemy.exc import IntegrityError, PendingRollbackError
For exiting system for trouble shooting import sys
"""
from more_itertools import chunked
import logging
from db import (
    add_protein,
    #name_addition,
    #cds_addition,
    add_xref,
    add_xdatabase,
)
from enums import ExDatabase, enum_from_str

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
            struct,
            canonical,
            name,
            # external dbs come next in sets of 2 (dbname, acc)
            *db_refs,
            dnaseq,
        ) = row
        
        # Check what data is available:
        logging.info(f"\nStarting next loop ({idx=}) with data: ({protseq[5]})")
        # Convert canonical to boolean
        canonical = parse_bool(canonical)

        try:  # attempt to add to db, if we raise an error, we roll back
            # Add protein data
            # accession values are now auto-generated; ignore CSV column
            protein = add_protein(
                session,
                protseq=protseq,
                struct=struct,
                canonical=canonical,
            )
            logging.info(f"\nProtein record returned: {protein}")
            if protein is None:
                logging.warning("Skipping row due to failed protein insertion: %5s", protseq)
                # nothing further to link, move to next entry
                session.rollback()
                continue
            
            # gene = cds_addition(
            #     session,
            #     cdstype=cdstype,
            #     dnaseq=dnaseq,
            #     protein=protein,
            # )
            # logging.info(f"\nProtein record returned: {protein}")
            # logging.info(f"Gene record returned: {gene}")

            # Add name data
            # name = name_addition(session, protname=protname, protein=protein)

            # logging.info(f"\nProtein record returned: {protein}")
            # logging.info(f"Name record returned: {name}")

            # Add external reference data
            # Handle external DB refs
            for db_name, db_acc in chunked(db_refs, 2):
                if not db_name or not db_acc:
                    continue  # skip empty db refs
                
                # Convert db_name string to ExDatabase enum with explicit typo catching
                try:
                    db_enum = enum_from_str(ExDatabase, db_name)
                except ValueError as e:
                    logging.error(f"Invalid database name '{db_name}': {e}. Skipping this xref.")
                    continue  # skip this xref but continue processing others
                
                # Known enum value: auto-add it without confirmation prompt
                xdb = add_xdatabase(session, db_enum, confirm=False)
                
                # Add xref if database was successfully retrieved/added
                if xdb:
                    xref = add_xref(
                        session,
                        xdb=xdb,
                        protein=protein,
                        xrefacc=db_acc,
                    )
                    logging.info(f"\nExternal reference record returned: {xref}")
            
            # Commit the transaction after successful addition
            logging.info("Committing the session")
            session.commit()
            # detach all objects to avoid cross-iteration state issues
            # ask Leighton unsure
            session.expunge_all()

        except Exception as exc:
            logging.error(f"Error adding data: {exc}")
            logging.info("Rolling back changes and skipping to next entry")
            session.rollback()
