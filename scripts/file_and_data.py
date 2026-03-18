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
    add_cds,
    Cds,
    # name_addition,
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
        # Parse CSV row with optional originseq column
        # Format: protseq, struct, canonical, name, *db_refs, dnaseq, [originseq]
        # originseq is optional and comes after dnaseq when present

        try:
            protseq = row[0]
            struct = row[1]
            canonical = row[2]
            name = row[3]

            rest = row[4:]

            def _looks_like_dna(val: str) -> bool:
                return (
                    isinstance(val, str)
                    and len(val.strip()) > 30
                    and set(val.strip().upper()) <= set("ACGTN")
                )

            # Expect the last two columns (if present) to be:
            #   - dna_seq (required)
            #   - originseq (optional; may be empty)
            # All preceding columns are treated as pairs of db refs.
            if len(rest) >= 2:
                dnaseq = rest[-2]
                originseq = rest[-1] if rest[-1] not in (None, "") else None
                db_refs = rest[:-2]
            elif len(rest) == 1:
                dnaseq = rest[0]
                originseq = None
                db_refs = []
            else:
                dnaseq = None
                originseq = None
                db_refs = []

            # If the originseq isn't a valid DNA string, ignore it (don't create a new CDS)
            if originseq is not None and not _looks_like_dna(originseq):
                logging.warning(
                    "Row %s: originseq does not look like DNA; ignoring value '%s'",
                    idx,
                    originseq,
                )
                originseq = None

            # Validate dna sequence as DNA too
            if dnaseq is not None and not _looks_like_dna(dnaseq):
                logging.warning(
                    "Row %s: dna_seq does not look like DNA; skipping row (value='%s')",
                    idx,
                    dnaseq,
                )
                dnaseq = None
        except (ValueError, IndexError) as e:
            logging.error(f"Row {idx}: Error parsing row: {e}")
            continue
        
        # Check what data is available:
        logging.info(f"\nStarting next loop ({idx=}) with protein prefix: {protseq[:10]}...")
        # Convert canonical to boolean
        canonical = parse_bool(canonical)

        try:  # attempt to add to db, if we raise an error, we roll back
            # Add protein data first
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

            # Now handle CDS if dna_seq is valid
            cds = None
            if dnaseq and dnaseq.strip():
                # Check if this CDS sequence already exists
                existing_cds = session.query(Cds).filter(Cds.cds_seq == dnaseq).first()
                if existing_cds:
                    logging.warning(
                        f"CDS sequence already exists (accession {existing_cds.cds_accession}); "
                        f"skipping CDS addition to maintain 1:1 CDS-Protein relationship"
                    )
                    # Don't rollback, just skip CDS
                else:
                    # CDS is new, add it
                    cds = add_cds(
                        session,
                        cdsseq=dnaseq,
                        originseq=originseq,  # Optional origin sequence for engineered/modified genes
                        protein=protein,
                    )
                    logging.info(f"\nCDS record returned: {cds}")
            else:
                logging.debug("No valid DNA sequence provided for this row; skipping CDS addition")
                # Don't rollback, just skip CDS

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
                    if cds:
                        logging.debug(f"Database '{db_name}' is for CDS linking")
                        # Create xref for this CDS database reference
                        xref = session.query(Xref).filter(Xref.xref_acc_ext == db_acc).first()
                        if not xref:
                            xref = Xref(xref_acc_ext=db_acc, xref_db=xdb)
                            session.add(xref)
                            session.flush()
                            logging.info(f"External reference {db_acc=} from database {db_name} added")
                        
                        # Link to CDS
                        link_cds = session.query(CdsXref).filter_by(cds_id=cds.cds_id, xref_id=xref.xref_id).first()
                        if not link_cds:
                            link_cds = CdsXref(cds_id=cds.cds_id, xref_id=xref.xref_id)
                            session.add(link_cds)
                            session.flush()
                            logging.info(f"Linked CDS {cds.cds_accession} <-> Xref {db_acc}")
                        logging.info(f"\nExternal reference record returned for CDS: {xref}")
                    else:
                        logging.warning(f"Database '{db_name}' is for CDS linking but no CDS data available; skipping xref")
                else:
                    # This database links to PROTEIN
                    logging.debug(f"Database '{db_name}' is for PROTEIN linking")
                    xref = add_xref(
                        session,
                        xdb=xdb,
                        protein=protein,
                        xrefacc=db_acc,
                    )
                    logging.info(f"\nExternal reference record returned for PROTEIN: {xref}")
            
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

