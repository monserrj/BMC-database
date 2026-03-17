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
from pathlib import Path

from db import (
    Xdatabase,
    add_protein,
    # name_addition,
    # cds_addition,
    add_xref,
    add_xdatabase,
    Xref,
    CdsXref,
    DATABASE_TARGETS,
)
from readfile import read_file

# To make boolean values from the csv true booleans:
def parse_bool(value):
    """Convert CSV boolean-like values to Python bool."""
    if isinstance(value, bool):
        return value  # already boolean
    # Accept TRUE/FALSE (case-insensitive) or yes/no
    return str(value).strip().lower() in ("true", "yes")

# Making a function for the addition of data to the different tables created with db.py following data_instructions:
def load_databases_from_csv(session, dbinfo_path: Path, verbose: bool = False):
    """Load external database info from a CSV and add to the Xdatabase table.

    The CSV is expected to have a header row, then rows like:
        name,type,url

    Databases are mapped to CDS or PROTEIN targets using DATABASE_TARGETS during import
    (mapping is only used at import time, not stored in database).
    Any invalid rows are logged and skipped.
    """

    dbinfo = read_file(dbinfo_path, verbose)
    for idx, row in enumerate(dbinfo):
        try:
            xname, xtype, xurl = row[0], row[1], row[2]
        except (ValueError, IndexError):
            logging.error("Invalid database info row %s: %s", idx, row)
            continue

        if not xname:
            logging.warning("Skipping empty database name at row %s", idx)
            continue

        xdb = add_xdatabase(
            session,
            xname=xname,
            xurl=xurl,
            xtype=xtype,
            require_password=False,
        )
        if not xdb:
            logging.warning("Failed to add or find database %s", xname)
    session.commit()


def link_db_csv(mydata, session, dbinfo_path: Path | None = None):
    """Import protein data and link references from a CSV.

    If dbinfo_path is provided, it is first loaded into the database.
    """

    # Optionally load external database info first
    if dbinfo_path:
        logging.info(f"Loading databases from CSV: {dbinfo_path}")
        load_databases_from_csv(session, dbinfo_path)
        logging.info("Database loading complete")

    # Process protein rows
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
        logging.info(f"\nStarting next loop ({idx=}) with protein prefix: {protseq[:10]}...")
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

                # Lookup the database record by name
                xdb = session.query(Xdatabase).filter_by(xref_db_name=db_name).first()

                if not xdb:
                    # Try to add the database (will prompt for password and info if not in approved list)
                    xdb = add_xdatabase(
                        session,
                        xname=db_name,
                        require_password=True,
                    )
                    if not xdb:
                        logging.error(
                            "Failed to add database '%s'; skipping this xref.",
                            db_name,
                        )
                        continue

                # Link xref based on database target (CDS or PROTEIN)
                # Use DATABASE_TARGETS mapping to determine linking behavior
                is_cds_db = DATABASE_TARGETS.get(db_name, False)  # Default to PROTEIN if not in mapping
                
                if is_cds_db:
                    # This database links to CDS (genes)
                    logging.debug(f"Database '{db_name}' is for CDS linking")
                    # TODO: Link to CDS when CDS data is available
                    # For now, this will be implemented when CDS is uncommented
                else:
                    # This database links to PROTEIN
                    logging.debug(f"Database '{db_name}' is for PROTEIN linking")
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

