#!/usr/bin/env python

"""This script is the development of the BMC db. Instructions
followed and explanations needed kept to help following through
the process."""

#######################################################
# (1) Imports
#######################################################

import logging  # We'll use this to report program state and other things to the user
import sys

from pathlib import Path  # for type hints
from typing import Optional


from enums import DatabaseType, StructProtType, ModificationType, enum_from_str

# Import  SQLAlchemy classes needed with a declarative approach.
from sqlalchemy import (
    Boolean,
    Column,
    Enum,
    Integer,
    String,
    ForeignKey,
    PrimaryKeyConstraint,
    UniqueConstraint,
    or_,
    create_engine,
    func,
)
from sqlalchemy.orm import (
    DeclarativeBase,
    relationship,
    Mapped,
    mapped_column,
    sessionmaker,
)
from sqlalchemy.exc import IntegrityError
import os
from getpass import getpass

#######################################################
# (2) Generic SQLAlchemy configuration
#######################################################


class Base(DeclarativeBase):
    """Declarative base class"""

    pass


# Move logger here, so it is consistent across modules and avoid accidental root logger use
logger = logging.getLogger(__name__)
Session = sessionmaker()  # We also need a Session object for database connections


#######################################################
# (3) Define the database structure
#######################################################


class Protein(Base):
    """Each row describes a unique protein.

    The table stores

    - protein sequence,
    - unique accession ID
    - type of protein structure
    - canonical or isoform
    """

    __tablename__ = "protein"  # this is the name that will be used in SQLite

    # Define Protein table columns using Declarative:
    prot_id: Mapped[int] = mapped_column(
        primary_key=True, autoincrement=True
    )  # Autopopulated ID for local table
    prot_accession: Mapped[str] = mapped_column(
        nullable=False, unique=True
    )  # A unique accession number, must be present
    prot_seq: Mapped[str] = mapped_column(nullable=False, unique=True)
    is_canonical: Mapped[Optional[bool]] = mapped_column(
        default=True
    )  # True/False flag for whether this protein is canonical
    struct_prot_type: Mapped[Optional[Enum]] = Column(Enum(StructProtType))

    # Define relationships to other tables in Declarative
    references: Mapped["ProteinXref"] = relationship(back_populates="protein")
    cds: Mapped["Cds"] = relationship(back_populates="protein")
    names: Mapped["ProteinName"] = relationship(back_populates="protein")
    isoforms: Mapped["Isoforms"] = relationship(
        back_populates="isoform", foreign_keys="[Isoforms.isoform_prot_id]"
    )
    canonicals: Mapped["Isoforms"] = relationship(
        back_populates="canonical", foreign_keys="[Isoforms.canonical_prot_id]"
    )
    complexes: Mapped["ProteinComplex"] = relationship(back_populates="protein")
    interaction_1: Mapped["Ppi"] = relationship(
        back_populates="protein_1",
        foreign_keys="[Ppi.prot_id_1]",
    )
    interaction_2: Mapped["Ppi"] = relationship(
        back_populates="protein_2",
        foreign_keys="[Ppi.prot_id_2]",
    )

    # Define string representation of protein row
    def __str__(self):
        outstr = [
            f"Protein ID: {self.prot_id}",
            f"Protein accession: {self.prot_id}, Protein sequence: {self.prot_seq}",
            f"Protein structure type: {self.struct_prot_type}",
            f"Protein is canonical: {self.is_canonical}",
        ]
        return "\n".join(outstr)


class Xdatabase(Base):
    """Each row describes an external database"""

    __tablename__ = "xdatabase"

    # Define database row using Declarative
    xref_db_id: Mapped[int] = mapped_column(
        primary_key=True, autoincrement=True
    )  # Local id for database
    xref_db_name: Mapped[str] = mapped_column(nullable=False, unique=True)  # Database name
    xref_db_type: Mapped[Enum] = Column(Enum(DatabaseType))  # Database type (sequence, structure, function, taxonomy)
    xref_db_url: Mapped[Optional[str]] = mapped_column(nullable=False, unique=True)  # Database URL
    
    # Define relationships using Declarative
    xref: Mapped[list["Xref"]] = relationship(back_populates="xref_db")




class Xref(Base):  # All external references
    """Each row describes an external database cross-reference"""

    __tablename__ = "xref"

    # Define table content:
    xref_id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    xref_acc_ext: Mapped[str] = mapped_column(
        nullable=False, unique=True
    )  # External DB accession
    xref_db_id: Mapped[int] = mapped_column(
        ForeignKey("xdatabase.xref_db_id"),
        nullable=False,
    )

    # Define relationships using Declarative
    xref_db: Mapped["Xdatabase"] = relationship(back_populates="xref")
    proteins: Mapped["ProteinXref"] = relationship(back_populates="xref")
    cdss: Mapped["CdsXref"] = relationship(back_populates="x_ref")


class ProteinXref(Base):
    """Linker table, protein to database cross-reference"""

    __tablename__ = "protein_xref"
    # Require a unique combination of protein ID and crossreference ID
    __table_args__ = (PrimaryKeyConstraint("prot_id", "xref_id"),)

    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),  # primary_key=True
    )
    xref_id: Mapped[int] = mapped_column(
        ForeignKey("xref.xref_id"),
    )

    # Define relationships using Declarative
    protein: Mapped["Protein"] = relationship(back_populates="references")
    xref: Mapped["Xref"] = relationship(back_populates="proteins")


class Cds(Base):
    """Each row describes a unique CDS.

    The table stores

    - gene ID,
    - DNA sequence,
    - origin sequences if engineered,
    - accession number
    - and the corresponding protein.
    """

    __tablename__ = "cds"

    # Define table content  using declarative:
    cds_id: Mapped[int] = mapped_column(
        primary_key=True, autoincrement=True
    )  # Autopopulated ID for local table
    cds_seq: Mapped[str] = mapped_column(
        nullable=False, unique=True
    )  # Gene DNA sequence
    cds_accession: Mapped[str] = mapped_column(
        nullable=False, unique=True
    )  # Unique accession number for CD
    prot_id: Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), nullable=False)

    # Introduce all relationship between tables: Check this?
    protein: Mapped["Protein"] = relationship(back_populates="cds")
    origins: Mapped["Origin_cds"] = relationship(
        back_populates="origin", foreign_keys="[Origin_cds.origin_id]"
    )
    modifieds: Mapped["Origin_cds"] = relationship(
        back_populates="modified", foreign_keys="[Origin_cds.modified_id]"
    )
    references: Mapped["CdsXref"] = relationship(back_populates="cds")
    modifications: Mapped["CdsModification"] = relationship(back_populates="cds")


class Origin_cds(Base):
    """Each row describes the link between original and modified CDS.

    The table stores

    - original DNA sequence id,
    - the corresponding modified CDS id.
    """

    __tablename__ = "origin_cds"
    __table_args__ = (PrimaryKeyConstraint("origin_id", "modified_id"),)

    # Define table content:
    origin_id: Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id"), nullable=False
    )  # Foreign key, original CDS sequence (cds_id key)
    modified_id: Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id"), nullable=False
    )  # Foreign key, modified CDS sequence (cds_id key)

    # Introduce all relationship between tables:
    modified: Mapped["Cds"] = relationship(
        "Cds", back_populates="modifieds", foreign_keys=[modified_id]
    )
    origin: Mapped["Cds"] = relationship(
        "Cds", back_populates="origins", foreign_keys=[origin_id]
    )


class CdsXref(Base):
    """Linker table, CDS to database cross-reference"""

    __tablename__ = "cds_xref"
    __table_args__ = (PrimaryKeyConstraint("cds_id", "xref_id"),)

    # Define table content:
    cds_id: Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id"),
    )
    xref_id: Mapped[int] = mapped_column(
        ForeignKey("xref.xref_id"),
    )

    # Introduce all relationship between tables:
    cds: Mapped["Cds"] = relationship(back_populates="references")
    x_ref: Mapped["Xref"] = relationship(back_populates="cdss")


class Modification(Base):
    """Each row represents a unique modification applied to a CDS.


    The table stores

    - modification ID,
    - modification type (enum): truncated, extended, fusion, synthetic, mutated, domesticated,
    - description of the modification.
    """

    __tablename__ = "modification"
    # __table_args__ = (UniqueConstraint("modification_description", "modification_type"),)

    # Define table content using declarative:
    modification_id: Mapped[int] = mapped_column(
        primary_key=True, autoincrement=True
    )  # Autopopulated ID for local table
    modification_description: Mapped[str] = mapped_column(nullable=False)

    modification_type: Mapped[Enum] = Column(Enum(ModificationType))

    # Introduce all relationship between tables:
    modifications: Mapped["CdsModification"] = relationship(
        back_populates="modification"
    )


class CdsModification(Base):
    """Each row represents a unique modification applied to a CDS.

    The table stores

    - modification ID,
    - CDS ID.
    """

    __tablename__ = "cds_modification"
    __table_args__ = (PrimaryKeyConstraint("cds_id", "modification_id"),)

    # Define table content:
    cds_id: Mapped[int] = mapped_column(ForeignKey("cds.cds_id"))
    modification_id: Mapped[int] = mapped_column(
        ForeignKey("modification.modification_id")
    )

    # Introduce all relationship between tables:
    modification: Mapped["Modification"] = relationship(back_populates="modifications")
    cds: Mapped["Cds"] = relationship(back_populates="modifications")


class Name(Base):
    """Each row describes a protein name.

    The table stores

    - name ID,
    - protein name.


    """

    __tablename__ = "name"

    # Define table content in Declarative:
    name_id: Mapped[int] = mapped_column(
        primary_key=True, autoincrement=True
    )  # primary key column
    prot_name: Mapped[str] = mapped_column(nullable=False)

    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinName"]] = relationship(back_populates="name")


class ProteinName(Base):
    """Each row describes a protein name associated to a protein.

    The table stores

    - protein ID,
    - name ID.
    Several proteins can share the same name
    A name can be shared by several proteins
    """

    __tablename__ = "protein_name"
    __table_args__ = (PrimaryKeyConstraint("prot_id", "name_id"),)

    # Define table content:
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    name_id: Mapped[int] = mapped_column(
        ForeignKey("name.name_id"),
    )

    # Introduce all relationship between tables:
    name: Mapped["Name"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="names")


class Isoforms(Base):
    """Each row describes a protein isoform relationship with its canonical protein.

    The table stores

    - canonical protein ID,
    - isoform protein ID.

    Several isoforms can be associated to a canonical protein.
    An isoform can only be associated to one canonical protein.
    """

    __tablename__ = "isoforms"
    __table_args__ = (PrimaryKeyConstraint("canonical_prot_id", "isoform_prot_id"),)

    # Define table content:
    canonical_prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    isoform_prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )

    # Introduce all relationship between tables:
    isoform: Mapped["Protein"] = relationship(
        "Protein", back_populates="isoforms", foreign_keys=[isoform_prot_id]
    )
    canonical: Mapped["Protein"] = relationship(
        "Protein", back_populates="canonicals", foreign_keys=[canonical_prot_id]
    )


class Complex(Base):
    """Each row describes a unique BMC-like complex.

    The table stores

    - complex ID,
    - complex name,
    - complex type (PDU, EUT...),
    - complex source (native, engineered, predicted, theoretical),
    - whether it is active,
    - whether it has been experimentally tested.
    """

    __tablename__ = "complex"

    # Define table content in Declarative:
    complex_id: Mapped[int] = mapped_column(
        primary_key=True, autoincrement=True
    )  # Autopopulated ID for local table
    complex_accession: Mapped[str] = mapped_column(nullable=False, unique=True)
    complex_type: Mapped[str] = mapped_column(
        nullable=False
    )  # Classification undecided (pdu,eut,grm..)
    is_active: Mapped[Optional[bool]] = mapped_column(nullable=True)  # Active/Inactive.
    is_exp_tested: Mapped[Optional[bool]] = mapped_column(nullable=True)  # Y/N.
    #    complex_source = Column(SQLEnum(ComplexSource, name="complex_source_enum"), nullable=False)  # Not currently defined

    # Introduce all relationship between tables:
    proteins: Mapped["ProteinComplex"] = relationship(
        "ProteinComplex",
        back_populates="complex",
    )
    interactions: Mapped["Ppi_complex"] = relationship(
        "Ppi_complex",
        back_populates="complex",
    )


class ProteinComplex(Base):
    """Each row describes a unique protein-complex relationship.

    The table stores

    - protein ID,
    - complex ID,
    - whether the protein is essential for the complex assembly,
    - copy number of the protein in the complex.
    """

    __tablename__ = "protein_complex"
    __table_args__ = (PrimaryKeyConstraint("prot_id", "complex_id"),)

    # Define table content:
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    complex_id: Mapped[int] = mapped_column(
        ForeignKey("complex.complex_id"),
    )

    is_essential: Mapped[Optional[bool]] = mapped_column(
        default=False, nullable=True
    )  # Whether the protein is essential for the complex assembly
    stoichiometry: Mapped[Optional[int]] = mapped_column(
        nullable=True
    )  # Copy number of the protein in the complex

    # Introduce all relationship between tables:
    complex: Mapped["Complex"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="complexes")


class Interaction(Base):
    """Each row describes a unique protein-protein interaction.

    The table stores

    - interaction ID,
    - interaction type,
    - interaction description.
    """

    __tablename__ = "interaction"
    __table_args__ = (UniqueConstraint("interact_type", "interact_description"),)

    # Define table content in declarative:
    interact_id: Mapped[int] = mapped_column(
        primary_key=True, autoincrement=True
    )  # Autopopulated ID for local table
    interact_type: Mapped[str] = mapped_column(
        nullable=False
    )  # Type of interaction (e.g: electrostatic, hydrophobic, etc)

    interact_description: Mapped[str] = mapped_column(
        nullable=True
    )  # Description of the interaction

    # Introduce all relationship between tables:
    ppis: Mapped["Ppi"] = relationship(back_populates="interaction")


class Ppi(Base):
    """Each row describes a unique protein-protein interaction between two  specific proteins.

    The table stores

    - ppiID (to assign a unique identifier to each interaction),
    - interaction ID,
    - protein ID 1,
    - protein ID 2.
    """

    __tablename__ = "Ppi"

    # Define table content:
    ppi_id: Mapped[int] = mapped_column(
        primary_key=True, unique=True, autoincrement=True
    )
    interact_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("interaction.interact_id"),
    )
    prot_id_1: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    prot_id_2: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )

    # Introduce all relationship between tables:
    complex: Mapped["Ppi_complex"] = relationship(back_populates="ppi")
    protein_1: Mapped["Protein"] = relationship(
        "Protein", back_populates="interaction_1", foreign_keys=[prot_id_1]
    )
    protein_2: Mapped["Protein"] = relationship(
        "Protein", back_populates="interaction_2", foreign_keys=[prot_id_2]
    )
    interaction: Mapped["Interaction"] = relationship(back_populates="ppis")


class Ppi_complex(Base):
    """Each row describes a unique protein-protein interaction associated to a complex.

    The table stores

    - ppiID,
    - complex ID.
    """

    __tablename__ = "ppi_complex"
    __table_args__ = (PrimaryKeyConstraint("ppi_id", "complex_id"),)

    # Define table content:
    ppi_id: Mapped[int] = mapped_column(
        ForeignKey("Ppi.ppi_id"),
    )
    complex_id: Mapped[int] = mapped_column(
        ForeignKey("complex.complex_id"),
    )

    # Introduce all relationship between tables:
    ppi: Mapped["Ppi"] = relationship(back_populates="complex")
    complex: Mapped["Complex"] = relationship(back_populates="interactions")


""" Functions to add data to the database"""


# # Function to add protein data
def add_protein(session, protseq, struct, canonical):
    """
    Args:
        session: SQLAlchemy session
        protseq (str): Protein sequence.
        struct (str): Protein structure type.
        canonical (boolean): Canonical if true, false if isoform

    Accessions are auto-generated in the format PRXXXX where XXXX is a
    zero-padded number derived from the protein's primary key. This ensures
    sequential identifiers starting at PR0001.

    Explanation on how the code works:
    1. If the protein already exists, store it in the `protein` variable
    2. If the protein does not exist, create a new protein object and add it to
    the session
    3. Commit the session to the database

    NOTE: SQLAlchemy will autoflush the session when we query the database, so we
        do not need to manually flush the session before committing
        (https://docs.sqlalchemy.org/en/20/orm/session_basics.html#session-flushing)
    This will also automatically commit the session if no exceptions are
        raised but, if errors are raised, we will need to rollback the session
        and continue to the next entry.
    """
    logger.debug(f"\nNow in {add_protein.__name__}")

    logger.debug(f"Before query, {protseq[:10]=}..., {struct=}, {canonical=}")

    # Create a new protein object
    try:
        # first check if the sequence already exists
        protein = session.query(Protein).filter(Protein.prot_seq == protseq).first()
        if protein:
            logger.info("Protein sequence already exists, returning existing accession %s", protein.prot_accession)
            return protein
        else:
            logger.debug(f"After query, {protein=}")
            
            # Enforce protein structure type
            try:
                struct = StructProtType[struct.upper()]
            except KeyError:
                logger.warning(
                    "Skipping invalid struct type: %s", struct,
                )
                return None

            # determine next accession number by looking at current max prot_id
            max_id = session.query(func.coalesce(func.max(Protein.prot_id), 0)).scalar()
            next_id = max_id + 1
            accession = f"PR{next_id:04d}"

            # create new protein with generated accession
            protein = Protein(
                prot_seq=protseq,
                prot_accession=accession,
                struct_prot_type=struct,
                is_canonical=canonical,
            )
            session.add(protein)
            session.flush()  # assign prot_id
            logger.info("Generated accession %s for protein ID %s", protein.prot_accession, protein.prot_id)

            return protein
    except Exception as exc:
        logger.exception("Failed to add protein accession=%s", accession)
        logger.exception(exc)
        session.rollback()
        raise


# Mapping of known databases to whether they link to CDS (True) or Protein (False)
# CDS (genes): NCBI, NCBITAX, GO (when linked to genes)
# Protein: Uniprot, KO (KEGG Orthology), PDB
DATABASE_TARGETS = {
    "NCBI": True,      # CDS
    "NCBITAX": True,   # CDS
    "GO": True,        # CDS
    "Uniprot": False,  # PROTEIN
    "KO": False,       # PROTEIN
    "PDB": False,      # PROTEIN
}

# Function to add Xdatabase data 
def add_xdatabase(session, xname, xurl=None, xtype=None, require_password=True):

#     ''' Args:
#     session: SQLAlchemy session to the database
#     xname (str): Database name.
#     xurl (str): Database URL.
#     xtype (str): Database type (sequence, structure, function, taxonomy)
#     confirm (bool): Whether to ask the user to confirm adding a new database type if xtype is not recognised

#     Explanation on how the code works:
#     1. Check if the database already exists,  store it in the `xdb` variable
#     2. If the database does not exist, check that this new database typews added (Y/N) questions
#     3. If yes, create a new db object and add it to the session
#     4. Commit the session to the database
#     '''

    logger.debug(f"\nNow in {add_xdatabase.__name__}")
    
    logger.debug(f"Before query, {xname=}")

    # Check existence by name first, then by URL
    xdb = (
        session.query(Xdatabase)
        .filter_by(xref_db_name=xname)
        .first()
    )
    if not xdb and xurl:
        xdb = (
            session.query(Xdatabase)
            .filter_by(xref_db_url=xurl)
            .first()
        )

    if xdb:
        logger.info("Database %s already exists", xname)
        return xdb

    # Validate/normalize xtype to DatabaseType enum
    if isinstance(xtype, str):
        try:
            xtype = DatabaseType[xtype.strip().upper()]
        except KeyError:
            valid = [e.name for e in DatabaseType]
            logger.warning(
                "Invalid database type '%s' (must be one of %s)",
                xtype,
                valid,
            )
            return None

    if not require_password:
        # Add database directly without password (for trusted CSV sources)
        try:
            xdb = Xdatabase(
                xref_db_name=xname,
                xref_db_url=xurl,
                xref_db_type=xtype,
            )
            session.add(xdb)
            session.flush()
            logger.info("Database '%s' added from trusted source", xname)
            return xdb
        except Exception as exc:
            logger.exception("Failed to add database '%s' — rolling back", xname)
            session.rollback()
            raise

    # Insert new database (requires password and may need additional info)
    logger.info("This external database %s does not exist, to add to database provide password", xname)
    password = getpass("Password: ")
    expected_password = os.getenv("DB_ADD_PASSWORD")

    if expected_password is None:
        # how to set password: export DB_ADD_PASSWORD="..."
        logger.error("Environment variable DB_ADD_PASSWORD is not set")
        return None

    if password == expected_password:
        # Password correct, now check if we need additional information
        if not xurl or not xtype:
            logger.info("Database '%s' approved. Please provide additional information:", xname)
            if not xurl:
                xurl = input("URL: ").strip()
            if not xtype:
                xtype_str = input("Type (Sequence/Structure/Function/Taxonomy): ").strip()
                try:
                    xtype = DatabaseType[xtype_str.upper()]
                except KeyError:
                    valid = [e.name for e in DatabaseType]
                    logger.error("Invalid database type '%s' (must be one of %s)", xtype_str, valid)
                    return None
            target_str = input("Link to CDS (gene) or Protein? [CDS/Protein]: ").strip().upper()
            # Note: This is only used at import time, not stored in database

        try:
            xdb = Xdatabase(
                xref_db_name=xname,
                xref_db_url=xurl,
                xref_db_type=xtype,
            )
            session.add(xdb)
            session.flush()

            logger.info("Database '%s' successfully added", xname)
            return xdb

        except Exception as exc:
            logger.exception("Failed to add database '%s' — rolling back", xname)
            logger.exception(exc)
            session.rollback()
            raise

    else:
        logger.warning("Incorrect password provided, database '%s' not added", xname)
        return None

# Function to add xref data and linking to cds and protein
def add_xref(session, xdb, protein, xrefacc, cds=None):

    ''' Args:
        xdb: Xdatabase ORM object (external DB already added)
        protein: Protein ORM object
        xrefacc (str): accession string for this external reference
        cds: Cds ORM object (optional, currently unused)

    Explanation on how the code works:
    1. Check if the xref already exists,  store it in the `xref` variable
    2. If the xref does not exist, create a new xref object and add it to the
    session

    Need to think on how to do that, condition depending on type of database?
    3. Check if the xref is already associated with the cds_xref, and if not
    create a new `cds_xref` object and add it to the session
    4. Check if the xref is already associated with the prot_xref, and if not
    create a new `prot_xref` object and add it to the session

    4. Commit the session to the database
    '''

    logger.debug(f"\nNow in {add_xref.__name__} with {xrefacc=}")

    logger.debug(f"Before query, {xdb=}, {xrefacc=}, {protein=}")

    # if no protein provided, skip linking steps later
    if protein is None:
        logger.warning("No protein object provided to add_xref; xref will not be linked to protein.")

    # Create a new xref object
    try:
        xref = (
            session.query(Xref)
            .filter(Xref.xref_acc_ext == xrefacc).first()
        )
        logger.debug(f"After query, {xref=}")

        # Add xref if it is not already present
        if not xref:
            # use the relationship attribute 'xref_db' to link to Xdatabase
            xref = Xref(xref_acc_ext=xrefacc, xref_db=xdb)
            session.add(xref)
            session.flush()
            logger.info(f"External reference {xrefacc=} from database {xdb=} added")
        else:
            logger.info(f"This external reference {xrefacc=} from database {xdb=} has already being added")
        logger.info(f"External reference row returned: External reference {xrefacc=}")

        # Associate the reference and protein information in the prot_xref
        logger.debug(f"{xref.proteins=}, {type(xref.proteins)}")

        # 2. Link to Protein (ProteinXref) if we have a protein object
        if protein is not None:
            # use explicit primary key values to avoid ORM dependency issues
            link_prot = (
                session.query(ProteinXref)
                .filter_by(prot_id=protein.prot_id, xref_id=xref.xref_id)
                .first()
            )

            if not link_prot:
                link_prot = ProteinXref(prot_id=protein.prot_id, xref_id=xref.xref_id)
                session.add(link_prot)
                session.flush()
                logger.info(f"Linked Protein {protein} <-> Xref {xref} (flush assigned link_prot={link_prot})")
            else:
                logger.info(f"Protein {protein} already linked to Xref {xref}")
        else:
            logger.debug("Skipping protein linkage because no protein provided.")

        # # 3. Link to CDS (CdsXref)
        # link_cds = (
        #     session.query(CdsXref)
        #     # Unsure about this?
        #     .filter_by(cds==cds and xref==xref)
        #     .first()
        # )

        # if not link_cds:
        #     link_cds = Cds_xref(cds=cds, xref=xref)
        #     session.add(link_cds)
        #     print(f"Linked CDS {cds} <-> Xref {xref}")
        # else:
        #     print(f"CDS {cds} already linked to Xref {xref}")
        return xref  # Return the reference row we just added to the db/otherwise dealt with
    except Exception as exc:
        logger.exception("Failed to add xref %s for database %s", xrefacc, xdb)
        logger.exception(exc)
        session.rollback()
        raise



# Function to add CDS data
def add_cds(session, cdsseq, originseq, protein):
    """Add a CDS (Coding DNA Sequence) record to the database.

    Args:
        session: SQLAlchemy session
        cdsseq (str): Gene DNA sequence (potentially modified)
        originseq (str): Original CDS sequence if engineered/modified (optional)
        protein: Protein ORM object this CDS encodes

    Accessions are auto-generated in the format CDXXXX where XXXX is a
    zero-padded number derived from the CDS's primary key.

    If originseq is provided and differs from cdsseq:
    - Creates/retrieves both original and modified CDS records
    - Links them via Origin_cds table
    If originseq is None or same as cdsseq:
    - Creates only the current CDS record

    Explanation on how the code works:
    1. Check if the current CDS sequence already exists
    2. If it exists, return existing CDS
    3. If originseq provided and different from cdsseq:
       a. Find or create the origin CDS
       b. Create the modified CDS
       c. Link them via Origin_cds table
    4. Else: Create only the current CDS
    5. Return the CDS row
    """
    logger.debug(f"\nNow in {add_cds.__name__}")
    logger.debug(f"Before query, {cdsseq[:10]=}..., {originseq[:10] if originseq else None}..., {protein=}")

    try:
        # First check if the current sequence already exists
        cds = session.query(Cds).filter(Cds.cds_seq == cdsseq).first()
        if cds:
            logger.info("CDS sequence already exists, returning existing accession %s", cds.cds_accession)
            return cds

        logger.debug(f"After query, {cds=}")

        # Handle origin sequence if provided and different from current sequence
        origin_cds = None
        if originseq and originseq.strip() and originseq != cdsseq:
            logger.debug("CDS has origin sequence, processing origin")
            # Find or create the origin CDS
            origin_cds = session.query(Cds).filter(Cds.cds_seq == originseq).first()
            if not origin_cds:
                # Create origin CDS with same protein
                max_id = session.query(func.coalesce(func.max(Cds.cds_id), 0)).scalar()
                origin_accession = f"CD{max_id + 1:04d}"
                origin_cds = Cds(
                    cds_seq=originseq,
                    cds_accession=origin_accession,
                    prot_id=protein.prot_id,
                )
                session.add(origin_cds)
                session.flush()
                logger.info("Created origin CDS with accession %s", origin_cds.cds_accession)
            else:
                logger.info("Origin CDS already exists with accession %s", origin_cds.cds_accession)

        # Determine next accession number for the modified CDS
        max_id = session.query(func.coalesce(func.max(Cds.cds_id), 0)).scalar()
        next_id = max_id + 1
        accession = f"CD{next_id:04d}"

        # Create new CDS with generated accession
        cds = Cds(
            cds_seq=cdsseq,
            cds_accession=accession,
            prot_id=protein.prot_id,
        )
        session.add(cds)
        session.flush()  # assign cds_id

        logger.info("Generated accession %s for CDS ID %s", cds.cds_accession, cds.cds_id)

        # Link origin to modified if origin exists
        if origin_cds:
            try:
                origin_link = Origin_cds(
                    origin_id=origin_cds.cds_id,
                    modified_id=cds.cds_id,
                )
                session.add(origin_link)
                session.flush()
                logger.info(
                    "Linked origin CDS %s (ID %s) to modified CDS %s (ID %s)",
                    origin_cds.cds_accession,
                    origin_cds.cds_id,
                    cds.cds_accession,
                    cds.cds_id,
                )
            except IntegrityError as e:
                logger.warning("Origin-modified CDS link already exists: %s", e)
                session.rollback()

        logger.debug(f"CDS row returned: {cds}")
        return cds
    except Exception as exc:
        logger.exception("Failed to add CDS with sequence=%s", cdsseq[:10])
        logger.exception(exc)
        session.rollback()
        raise
# Add the add_modification function
def add_modification(session, modif_descrip, modif_type, cds):
    """
    Adds or retrieves a modification and associates it with a CDS.
    
    Args:
        session: SQLAlchemy session
        modif_descrip (str): Description of the modification.
        modif_type (str): Type of modification (e.g., truncated, extended, fusion, synthetic, mutated, domesticated).
        cds: CDS object to associate with the modification.
        
    Returns:
        CdsModification object representing the link between CDS and modification.
        
    Explanation on how the code works:
    1. Check if a modification with the same description and type already exists.
    2. If it exists, retrieve its ID; if not, create a new modification record.
    3. Check if the link between the CDS and modification already exists in cds_modification.
    4. If the link does not exist, create a new record to associate the CDS with the modification.
    5. Flush and return the link.
    """
    logger.debug(f"\nNow in {add_modification.__name__}")
    logger.debug(f"Before query, {modif_descrip=}, {modif_type=}, {cds=}")
    
    try:
        # Skip if no description provided
        if not modif_descrip or not modif_type:
            logger.debug("Skipping modification: description or type is missing")
            return None
        
        # Enforce modification type enum conversion first
        try:
            modif_type_enum = ModificationType[modif_type.strip().upper()]
        except KeyError:
            valid = [e.name for e in ModificationType]
            logger.warning(
                "Invalid modification type '%s' (must be one of %s)",
                modif_type,
                valid,
            )
            return None
        
        # Check if modification already exists
        modification = session.query(Modification).filter(
            Modification.modification_description == modif_descrip,
        ).first()
        
        if modification:
            logger.info(f"Modification already exists with ID {modification.modification_id}")
        else:
            # Create new modification record
            modification = Modification(
                modification_description=modif_descrip,
                modification_type=modif_type_enum
            )
            session.add(modification)
            session.flush()
            logger.info(f"Added new modification with ID {modification.modification_id}")
        
        # Now link modification to CDS in cds_modification
        cds_mod_link = session.query(CdsModification).filter(
            CdsModification.cds_id == cds.cds_id,
            CdsModification.modification_id == modification.modification_id
        ).first()
        
        if cds_mod_link:
            logger.info(f"CDS-modification link already exists")
        else:
            try:
                cds_mod_link = CdsModification(
                    cds_id=cds.cds_id,
                    modification_id=modification.modification_id
                )
                session.add(cds_mod_link)
                session.flush()
                logger.info(f"Linked CDS {cds.cds_accession} to modification '{modif_descrip}'")
            except IntegrityError as e:
                logger.warning("CDS-modification link already exists: %s", e)
                session.rollback()
        
        return cds_mod_link
        
    except Exception as exc:
        logger.exception("Failed to add modification with description=%s", modif_descrip)
        logger.exception(exc)
        session.rollback()
        raise

# # Function to add name data 
# def add_name(session, protname, protein):

#     ''' Args:
#     protname (str): Protein name.
#     protein: Protein information added with add_protein function

#     Explanation on how the code works:
#     1. Check if the name already exists,  store it in the `name` variable
#     2. If the name does not exist, create a new name object and add it to the
#        session
#     3. Check if the name is already associated with the protein, and if not
#        create a new `proteinname` object and add it to the session
#     4. Commit the session to the database
#     '''

#     print(f"\nNow in {add_name.__name__}")

#     with session.no_autoflush:
#         print(f"Before query, {protname=}")

#         # Create a new name object
#         name = (
#             session.query(Name)
#             .filter(Name.prot_name == protname)
#             .first()
#         )
#         print(f"After query, {name=}")

#         # Add name if it is not already present
#         if not name:
#             name = Name(prot_name=protname)
#             session.add(name)
#             session.flush()
#             print(f"Name {protname=} added")
#         else:
#             print(f"This name {protname} has already being added")
#         print(f"Name row returned: {name}")

#         print(f"{name.proteins=}, {type(name.proteins)}")

#         if name not in name.proteins:
#             proteinname = ProteinName()
#             proteinname.name = name
#             protein.names.append(proteinname)
#             print(f"{name.proteins=}")
#             print(f"\nLinked name {name.name_id} to Protein {protein.prot_id}")
#             print(f"{proteinname=}")
#         else:
#             print(
#                 f"Name {name.name_id} is already linked to Protein {protein.prot_id}"
#             )
#         print(f"{name}")
#         print(f"Linked name from protein: {proteingene.gene}")

#         return name  # Return the gene row we just added to the db/otherwise dealt with


# Add modification support
def add_modification(session, modif_descrip, modif_type, cds):
    """Add or retrieve a modification and link it to a CDS.

    Args:
        session: SQLAlchemy session
        modif_descrip: Human-readable modification description.
        modif_type: Modification type string (must match ModificationType enum).
        cds: CDS ORM object to link the modification to.

    Returns:
        CdsModification ORM object (link), or None if it could not be created.
    """

    logger.debug(f"\nNow in {add_modification.__name__}")

    try:
        if not modif_descrip or not modif_type or not cds:
            logger.debug("Skipping modification: missing description/type/CDS")
            return None

        try:
            modif_type_enum = ModificationType[modif_type.strip().upper()]
        except KeyError:
            valid = [e.name for e in ModificationType]
            logger.warning(
                "Invalid modification type '%s' (must be one of %s)",
                modif_type,
                valid,
            )
            return None

        modification = (
            session.query(Modification)
            .filter(
                Modification.modification_description == modif_descrip,
                Modification.modification_type == modif_type_enum,
            )
            .first()
        )

        if modification:
            logger.info("Modification already exists with ID %s", modification.modification_id)
        else:
            modification = Modification(
                modification_description=modif_descrip,
                modification_type=modif_type_enum,
            )
            session.add(modification)
            session.flush()
            logger.info("Added new modification with ID %s", modification.modification_id)

        cds_mod_link = (
            session.query(CdsModification)
            .filter(
                CdsModification.cds_id == cds.cds_id,
                CdsModification.modification_id == modification.modification_id,
            )
            .first()
        )

        if cds_mod_link:
            logger.info("CDS-modification link already exists")
            return cds_mod_link

        cds_mod_link = CdsModification(
            cds_id=cds.cds_id,
            modification_id=modification.modification_id,
        )
        session.add(cds_mod_link)
        session.flush()
        logger.info(
            "Linked CDS %s to modification '%s'",
            cds.cds_accession,
            modif_descrip,
        )
        return cds_mod_link

    except Exception as exc:
        logger.exception("Failed to add modification with description=%s", modif_descrip)
        logger.exception(exc)
        session.rollback()
        raise


#######################################################
# (4) Helper functions that use the definitions to
#     construct or interact with the database
#######################################################


# Function to create the database and tables
def create_db(dbpath):
    """Function to create all the tables from the database"""
    # Create a database engine to connect to the database.
    # This creates a new empty database file called bmc.db in the current directory.

    logger.info("Initialising database")

    dbpath = Path(dbpath).resolve()  # ensure it is a Path object
    dbpath.parent.mkdir(parents=True, exist_ok=True)  # make sure folder exists

    db_URL = f"sqlite:///{dbpath.as_posix()}"

    logger.debug("Binding to session at %s", db_URL)
    engine = create_engine(db_URL)

    try:
        Base.metadata.create_all(bind=engine)
        logger.info("Database and tables created successfully")

    except Exception as exc:
        logger.error("Error creating database: %s", exc)
        raise


# Function to get a live session to the database
def get_session(dbpath):
    """Returns live session to database."""
    logger.debug("Opening database session for %s", dbpath)
    dbpath = Path(dbpath).resolve()  # ensure it is a Path object
    db_URL = f"sqlite:///{dbpath.as_posix()}"
    engine = create_engine(db_URL)
    Session.configure(bind=engine)
    return Session()


#######################################################
# (5) Code that runs when this file is called as
#     a script
#######################################################


if __name__ == "__main__":
    # Set up logging
    logger = logging.getLogger()
    logformatter = logging.Formatter(
        "[%(levelname)s] %(funcName)s %(asctime)s %(message)s"
    )
    logger.setLevel(logging.DEBUG)

    # Add a handler that writes to the console
    termhandler = logging.StreamHandler(sys.stdout)
    termhandler.setFormatter(logformatter)
    logger.addHandler(termhandler)

    # Path to output database
    outdbpath = Path("db.sqlite3")

    # Create database
    logger.info("Creating database at %s", outdbpath)
    create_db(outdbpath)
    logger.info("Created database at %s", outdbpath)

    # Render database as an ER diagram
    #from eralchemy import render_er

    # render_er(Base, "er_diagram.pdf")

## How to populate parent child relationships for CDS
# Suppose we have CDS with the following relationships
#
# seq        acc     parent
# ATG...TAC  CDS_1   NULL
# ATG...TAC  CDS_2   CDS_1
# ATG...TAC  CDS_3   CDS_2
#
# We'd parse/populate (in a function) as follows:
#
# line 1
#  add CDS acc, CDS seq, and because parent is NULL, origin_cds is NULL
# line 2
#  as parent is CDS_1, query Cds table to get cds_id corresponding to CDS_1 (call this ID*)
#  add CDS acc, CDS seq, and ID* (as parent)
# line 3
#  as parent is CDS_2, query Cds table to get cds_id corresponding to CDS_2 (call this ID**)
#  add CDS acc, CDS seq, and ID** (as parent)
