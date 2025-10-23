#!/usr/bin/env python

"""This script is the development of the BMC db. Instructions
followed and explanations needed kept to help following through
the process."""

import logging  # We'll use this to report program state and other things to the user
import sys

from pathlib import Path  # for type hints
from .enums import DatabaseType,   StructProtType, ModificationType, enum_from_str
# Import  SQLAlchemy classes needed with a declarative approach.
from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy import (
    Boolean,
    Column,
    Integer,
    String,
    ForeignKey,
    PrimaryKeyConstraint,
    UniqueConstraint,
    String,
    or_,
)
from sqlalchemy.orm import relationship, Mapped, mapped_column

# Import a sessionmaker to create a session object
from sqlalchemy.orm import sessionmaker

# Create_engine function to create an engine object
from sqlalchemy import create_engine

from typing import Optional

# To create enumerations for some of the vocabulary restricted fields
from enum import Enum

# Database set up:

# Create a base class to inherit from.
Base = declarative_base()

# Session set up:
Session = sessionmaker()  # we also need a session object


# Declare tables
class Protein(Base):
    """Table representing a protein This table will store the protein sequence,
    unique accession number, type of protein structure and whether is
    canonical and isoform.
    """

    __tablename__ = "protein"  # this is the name that will be used in SQLite

    # Define table content:
    prot_id = Column(Integer, primary_key=True, autoincrement=True)
    prot_accession = Column(
        Integer, nullable=False, unique=True
    )  # Unique accession number for proteins ( I want to automate it, any suggestions on resources for this?)
    prot_seq = Column(String, nullable=False, unique=True)
    struct_prot_type = Column(SQLEnum(StructProtType, name="struct_prot_type_enum"))
    is_canonical = Column(Boolean, default=True)  # Only false if it is an isoform

    # Define type of output for protein
    def __str__(self):
        outstr = [
            f"Protein ID: {self.prot_id}",
            f"Protein accession: {self.prot_id}Protein sequence: {self.prot_seq}",
            f"Protein structure type: {self.struct_prot_type}",
            f"Protein is canonical: {self.is_canonical}",
        ]
        return "\n".join(outstr)

    # Unsure if needed, but to keep the code consistent with the rest of the db
    # Introduce all relationship between tables:
    references: Mapped["ProteinXref"] = relationship(back_populates="protein")
    cds: Mapped["Cds"] = relationship(back_populates="protein")
    names: Mapped["ProteinName"] = relationship(back_populates="protein")
    isoforms: Mapped["Isoforms"] = relationship(back_populates="isoform")
    canonicals: Mapped["Isoforms"] = relationship(back_populates="canonical")
    complexes: Mapped["ProteinComplex"] = relationship(back_populates="protein")
    interaction_1: Mapped["Prot_prot_interact"] = relationship(
        back_populates="protein_1"
    )
    interaction_2: Mapped["Prot_prot_interact"] = relationship(
        back_populates="protein_2"
    )


class Xdatabase(Base):
    """Descriptions of external databases (linked by Xref accessions)"""

    __tablename__ = "xdatabase"

    # Define table content:
    xref_db_id = Column(Integer, primary_key=True)
    xref_db_name = Column(String, nullable=False, unique=True)
    xref_href = Column(String, nullable=False, unique=True)  # URL for database access
    xref_type = Colun(SQLEnum(DatabaseType, name ="database_type_enum", nullable=False))

    # Introduce all relationship between tables:
    xref: Mapped["Xref"] = relationship(back_populates="xref_db_id")
    databases: Mapped[list["Xref"]] = relationship()


class Xref(Base):  # All external references
    """Table representing a singe database cross-reference"""

    __tablename__ = "xref"

    # Define table content:
    xref_id = Column(Integer, primary_key=True, autoincrement=True)
    xref_acc_ext = Column(
        String, unique=True, nullable=False
    )  # Accession from external db
    xref_db_id: Mapped[int] = mapped_column(
        ForeignKey("xdatabase.xref_db_id"),
        nullable=False,
    )

    # Introduce all relationships between tables
    xref_db: Mapped["Xdatabase"] = relationship(back_populates="databases")
    proteins: Mapped["ProteinXref"] = relationship(back_populates="xref")
    cdss: Mapped["CdsXref"] = relationship(back_populates="cds")


class ProteinXref(Base):
    """Linker table, protein to database cross-reference"""

    __tablename__ = "protein_xref"
    # For linker tables, we need to set a PrimaryKeyConstraint,
    # so SQLite knows what the primary key is for the table
    __table_args__ = (PrimaryKeyConstraint("prot_id", "xref_id"),)

    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),  # primary_key=True
    )
    xref_id: Mapped[int] = mapped_column(
        ForeignKey("xref.xref_id"),
    )

    # Introduce all relationship between tables:
    protein: Mapped["Protein"] = relationship(back_populates="references")
    xref: Mapped["Xref"] = relationship(back_populates="proteins")


class Cds(Base):
    """Table representing the gene details linked to a protein
    This table will store the gene ID, DNA sequence,
    origin sequences if engineered, accession number and the corresponding
    protein
    """

    __tablename__ = "cds"

    # Define table content:
    cds_id = Column(Integer, primary_key=True)
    cds_seq = Column(String, nullable=False, unique=True)
    cds_accession = Column(
        String, nullable=False, unique=True
    )  # Unique accession number for CDS

    prot_id: Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), nullable=False)

    # Introduce all relationship between tables: Check this?
    protein: Mapped["Protein"] = relationship(back_populates="cds")
    origins: Mapped["Origin"] = relationship(back_populates="origin")
    cdss: Mapped["Origin"] = relationship(back_populates="cds")
    references: Mapped["CdsXref"] = relationship(back_populates="cds_xref")
    modifications: Mapped["CdsModification"] = relationship(back_populates="cds")


class Origin(Base):
    """Table representing the original DNA sequence of a modified CDS
    This table will store the original DNA sequence and the corresponding
    modified CDS sequence id
    """

    __tablename__ = "origin"
    __table_args__ = (PrimaryKeyConstraint("origin_id", "cds_id"),)

    # Define table content:
    # Should I change names here?
    origin_id: Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id"), nullable=False
    )  # Original CDS sequence (cds_id key)
    cds_id: Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id"), nullable=False
    )  # Foreign key, modified CDS sequence (cds_id key)

    # Introduce all relationship between tables:
    cds: Mapped["Cds"] = relationship(back_populates="cds")
    origin: Mapped["Cds"] = relationship(back_populates="origins")


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
    """Table representing the modification of the engineered CDS sequences
    This modifications will be vocabulary restricted:
    truncated, fusion, synthetic, mutated, domesticated"""

    __tablename__ = "modification"

    # Define table content:
    modification_id = Column(Integer, primary_key=True, autoincrement=True)
    modification_type = Column(SQLEnum(ModificationType, name="modification_type_enum"), nullable=False)
    modification_description = Column(String, nullable=False)

    UniqueConstraint(modification_type, modification_description)

    # Introduce all relationship between tables:
    modifications: Mapped["CdsModification"] = relationship(
        back_populates="modification"
    )


class CdsModification(Base):
    """Linker table, CDS to modification cross-reference"""

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
    """Table representing a protein name
    Each name_ID represents a protein name. Several name strings
    given to an unique protein, and several proteins sharing same name
    """

    __tablename__ = "name"

    # Define table content:
    name_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column
    prot_name = Column(String, nullable=False)

    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinName"]] = relationship(back_populates="name")


class ProteinName(Base):
    """Linker table, protein to name cross-reference
    This table will store the protein ID and the name ID"""

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
    """Table representing the isoforms of a protein
    This table will store the canonical protein ID and the isoform protein ID.
    The canonical protein ID is the one that is considered the main protein,
    while the isoform protein ID is a variant of the canonical protein.
    """

    __tablename__ = "isoforms"
    __table_args__ = (PrimaryKeyConstraint("canonical__prot_id", "isoform__prot_id"),)

    # Define table content:
    canonical__prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    isoform__prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )

    # Introduce all relationship between tables:
    isoform: Mapped["Protein"] = relationship(back_populates="isoforms")
    canonical: Mapped["Protein"] = relationship(back_populates="canonicals")


class Complex(Base):
    """Table representing the complex that can be form by the interaction
    between several proteins, including native BMC or engineered ones
    This table will store the complex features, including the type (e.g: pdu, eut),
    whether is has enzymatic activty or not, if it has been experimentally tested
    whether it assembles or not, the origin of the complex in a restricted vocabulary (meaning whether is it
    a native complex, engineered or created with a theorical or bioinformatic approach)
    """

    __tablename__ = "complex"

    # Define table content:
    complex_id = Column(Integer, primary_key=True, autoincrement=True)
    complex_name = Column(String, nullable=False)
    complex_type = Column(
        String, nullable=False
    )  # Classification undecided (pdu,eut,grm..)
    is_active = Column(String)  # Active/Inactive. Boolean?
    is_exp_tested = Column(String)  # Y/N. Boolean?
    complex_source = Column(SQLEnum(ComplexSource, name="complex_source_enum"), nullable=False)

    UniqueConstraint(
        "complex_name", "complex_type", "is_active", "is_exp_tested", "complex_source"
    )

    # Introduce all relationship between tables:
    proteins: Mapped["ProteinComplex"] = relationship(
        "ProteinComplex",
        secondary="protein_complex",
        back_populates="complexes",
        lazy="dynamic",
    )


class ProteinComplex(Base):
    """Linker table, protein to complex cross-reference
    This table will store the protein ID and the complex ID"""

    __tablename__ = "protein_complex"
    __table_args__ = (PrimaryKeyConstraint("prot_id", "complex_id"),)

    # Define table content:
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    complex_id: Mapped[int] = mapped_column(
        ForeignKey("complex.complex_id"),
    )
    is_essential = Column(Boolean, default=False, nullable=True)
    # Whether the protein is essential for the complex assembly
    copy_number = Column(Integer, nullable=True)
    # Copy number of the protein in the complex

    # Introduce all relationship between tables:
    complex: Mapped["Complex"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="complexes")


class Interaction(Base):
    """Table representing the different interactions
    between several proteins.
    """

    __tablename__ = "interaction"

    # Define table content:
    interact_id = Column(Integer, primary_key=True, autoincrement=True)
    interact_type = Column(
        String, nullable=False
    )  # Type of interaction (e.g: electrostatic, hydrophobic, etc)
    interact_description = Column(String)

    UniqueConstraint(interact_type, interact_description)

    # Introduce all relationship between tables:
    interactions: Mapped["Prot_prot_interact"] = relationship(
        back_populates="interaction"
    )


class Prot_prot_interact(Base):
    """Table representing the different specific interactions
    between several proteins inside a complex and during assembly
    This table will store the different proteins known to interact
    specifically with each other.
    """

    __tablename__ = "prot_prot_interact"
    __table_args__ = (PrimaryKeyConstraint("interact_id", "prot_id_1", "prot_id_2"),)

    # Define table content:
    interact_id: Mapped[int] = mapped_column(
        ForeignKey("interaction.interact_id"),
    )
    prot_id_1: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    prot_id_2: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )

    # To enforce unique no repeated protein to protein interactions are created
    UniqueConstraint("prot_id_1", "prot_id_2", "interaction_id")

    # Introduce all relationship between tables:
    protein_1: Mapped["Protein"] = relationship(back_populates="interaction_1")
    protein_2: Mapped["Protein"] = relationship(back_populates="interaction_2")
    interaction: Mapped["Interaction"] = relationship(back_populates="interactions")


""" Functions to add data to the database"""

# # Function for protein data addition
# def protein_addition(session, protaccession, protseq, struct, canonical):

#     '''
#     Args:
#     protseq (str): Protein sequence.
#     protaccession (str): Protein unique accession number
#     struct (str): Protein structure type.
#     canonical (boolean): Canonical if true, false if isoform

#     Explanation on how the code works:
#     1. If the protein already exists, store it in the `protein` variable
#     2. If the protein does not exist, create a new protein object and add it to
#        the session
#     3. Commit the session to the database

#     NOTE: SQLAlchemy will autoflush the session when we query the database, so we
#           do not need to manually flush the session before committing
#           (https://docs.sqlalchemy.org/en/20/orm/session_basics.html#session-flushing)
#            This will also automatically commit the session if no exceptions are
#           raised but, if errors are raised, we will need to rollback the session
#           and continue to the next entry.
#     '''

#     print(f"\nNow in {protein_addition.__name__}")

#     print(
#         f"Before query, {protseq[:10]=}..., {protaccession=}, {struct=}, {canonical=}"
#     )

    # Check and convert xtype
    # if isinstance(xtype, str):
    #     try:
    #         xtype = enum_from_str(StructProtType, xtype)
    #     except ValueError as e:
    #         print(f"Invalid protein type: {xtype!r} — {e}")
    #         return None

#     # Create a new protein object
#     protein = (
#         session.query(Protein)
#         .filter(Protein.prot_seq == protseq)
#         .filter(Protein.prot_accession == protaccession) # Do I need this?
#         .first()
#     )
#     print(f"After query, {protein=}")

#     # Add protein if it is not already present
#     if not protein:
#         protein = Protein(
#             prot_seq=protseq,
#             prot_accession=protaccession,
#             struct_prot_type=struct,
#             is_canonical=canonical,
#         )
#         session.add(protein)
#         session.flush()  # This sends the changes to the database, so prot_id is assigned
#     else:
#         print(f"Protein with accession number {protaccession} already exists")

#     print(f"Protein row returned: {protein}")

#     # except Exception as exc:
#     #     print(f"Error committing protein/gene combination: {exc}")
#     #     print("Rolling back changes and skipping to next entry\n")
#     #     session.rollback()
#     # Rollback makes it that when there is a "fail",
#     # like not unique uniprot reference, its "forgets" the error and keeps going.
#     return (
#         protein  # Return the protein row we just added to the db
#     )

# # Function for Xdatabase data addition
# def xdatabase_addition(session, xname, href, xtype):

#     ''' Args:
#     xdb_name (str): Database name.
#     href (str): Database URL.
#     xdb_type (str): Database type (sequence, structure, function, taxonomy)

#     Explanation on how the code works:
#     1. Check if the database already exists,  store it in the `xdb` variable
#     2. If the database does not exist, create a new db object and add it to the
#        session
#     3. Commit the session to the database
#     '''

#     print(f"\nNow in {xdatabase_addition.__name__}")
#     print(f"Before query, {xname=}, {href=}, {xtype=}")

    # Check and convert xtype
    # if isinstance(xtype, str):
    #     try:
    #         xtype = enum_from_str(DatabaseType, xtype)
    #     except ValueError as e:
    #         print(f"Invalid database type: {xtype!r} — {e}")
    #         return None

    # elif not isinstance(xtype, DatabaseType):
    #     print(f"xtype must be a string or DatabaseType, got {type(xtype)}")
    #     return None
    # print(f"Converted xtype to Enum: {xtype}")
    
#     # Create a new db object
#     xdb = (
#         session.query(Xdatabase)
#         .filter(Xdatabase.xref_db_name == xname)
#         .filter(Xdatabase.xref_href == href)
#         .first()
#     )
#     print(f"After query, {xdb=}")

#     # Add db if it is not already present
#     if not xdb:
#         xdb = Xdatabase(xref_db_name=xname, xref_href=href, xref_type=xtype)
#         session.add(xdb)
#         session.flush()
#         print(f"Database {xname} added")
#     else:
#         print(f"This database {xname} has already being added")
#     print(f"Database row returned: {xdb}")

#     return xdb  # Return the db row we just added to the db/otherwise dealt with

# # Function for CDS data addition
# def cds_addition(session, cdsseq, cdsaccession, cdsorigin, protein):

#     ''' Args:
#     cdsseq (str): Gene DNA sequence.
#     cdsaccession: (str) Unique accession number CDS
#     cds_origin (str): Original CDS sequence if the sequence is modified.
#     protein: Protein information added with protein_addition function

#     Explanation on how the code works:
#     1. Check if the cds already exists,  store it in the `cds` variable
#     2. If the cds does not exist, create a new cds object and add it to the
#        session
#     3. Check if the cds is already associated with the protein, and if not
#        create a new `proteingene` object and add it to the session
#     4. Commit the session to the database
#     '''

#     print(f"\nNow in {cds_addition.__name__}")

#     print(f"Before query, {cdsseq[:10]=}...,{cdsaccession=} {cdsorigin=}")

#     # Create a new cds object
#     cds = (
#         session.query(Cds)
#         .filter(Cds.cds_seq == cdsseq)
#         .first()
#     )
#     print(f"After query, {cds=}")

#     # Add cds if it is not already present
#     if not cds:
#         print(f"Before adding seq and accession {cds=}")
#         cds = Cds(
#             cds_seq=cdsseq,
#             cds_accession=cdsaccession,
#             protein=protein, # Use relationship not FK directly
#             origin=cdsorigin) # Use relationship not FK directly
#         print(f" After adding seq and accession{cds=}")
#         print(f"{protein.prot_id=}, {cds.cds_id=}, {cds.prot_id}")
#         print(f"{protein.cdss=}")
#         session.add(cds)
#         session.flush()
#         print(f"Cds {cdsseq[:10]=} with accession {cdsaccession=} added")
#     else:
#         print(f"This CDS {cdsseq[:10]=} with accession {cdsaccession=} has already being added")

#     print(f"Cds row returned: {cds}")

#     return cds  # Return the cds row we just added to the db/otherwise dealt with

# # Function for xref data addition and linking to cds and protein
# def xref_addition(session, xdb, cds, protein, xrefacc):

#     ''' Args:
#         xdb: Xdatabase ORM object (external DB already added)
#         cds: Cds ORM object
#         protein: Protein ORM object
#         xrefacc (str): accession string for this external reference

#     Explanation on how the code works:
#     1. Check if the xref already exists,  store it in the `xref` variable
#     2. If the xref does not exist, create a new xref object and add it to the
#        session

#     Need to think on how to do that, condition depending on type of database?
#     3. Check if the xref is already associated with the cds_xref, and if not
#        create a new `cds_xref` object and add it to the session
#     4. Check if the xref is already associated with the prot_xref, and if not
#        create a new `prot_xref` object and add it to the session

#     4. Commit the session to the database
#     '''

#     print(f"\nNow in {xref_addition.__name__} with {xrefacc=}")

#     print(f"Before query, {xdb=}, {xrefacc=}")

#     # Create a new xref object
#     xref = (
#         session.query(Xref)
#         .filter(Xref.xref_acc_ext == xrefacc)
#         .first()
#     )
#     print(f"After query, {xref=}")

#     # Add xref if it is not already present
#     if not xref:
#         xref = Xref(xref_acc_ext=xrefacc, XDatabase=xdb)
#         session.add(xref)
#         session.flush()
#         print(f"External reference {xrefacc=} from database {xdb=} added")
#     else:
#         print(f"This external reference {xrefacc=} from database {xdb=} has already being added")
#     print(f"External reference row returned: External reference {xrefacc=}")

#     # Associate the reference and protein information in the prot_xref
#     print(f"{xref.proteins=}, {type(xref.proteins)}")

#     # 2. Link to Protein (ProteinXref)
#     link_prot = (
#         session.query(ProteinXref)
#         # Unsure about this?
#         .filter_by(protein==protein and xref==xref)
#         .first()
#     )

#     if not link_prot:
#         link_prot = ProteinXref(protein=protein, xref=xref)
#         session.add(link_prot)
#         print(f"Linked Protein {protein} <-> Xref {xref}")
#     else:
#         print(f"Protein {protein} already linked to Xref {xref}")

#     # 3. Link to CDS (CdsXref)
#     link_cds = (
#         session.query(CdsXref)
#         # Unsure about this?
#         .filter_by(cds==cds and xref==xref)
#         .first()
#     )

#     if not link_cds:
#         link_cds = Cds_xref(cds=cds, xref=xref)
#         session.add(link_cds)
#         print(f"Linked CDS {cds} <-> Xref {xref}")
#     else:
#         print(f"CDS {cds} already linked to Xref {xref}")


#     return xref  # Return the renference row we just added to the db/otherwise dealt with

# # Function for name data addition
# def name_addition(session, protname, protein):

#     ''' Args:
#     protname (str): Protein name.
#     protein: Protein information added with protein_addition function

#     Explanation on how the code works:
#     1. Check if the name already exists,  store it in the `name` variable
#     2. If the name does not exist, create a new name object and add it to the
#        session
#     3. Check if the name is already associated with the protein, and if not
#        create a new `proteinname` object and add it to the session
#     4. Commit the session to the database
#     '''

#     print(f"\nNow in {name_addition.__name__}")

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


# Function to create the database and tables
def create_db(dbpath: Path):
    """Function to create all the tables from the database"""
    # Create a database engine to connect to the database.
    # This creates a new empty database file called bmc.db in the current directory.
    db_URL = f"sqlite:///{dbpath}"
    engine = create_engine(db_URL)
    Base.metadata.create_all(bind=engine)
    print("Database and tables created successfully")


# Function to get a live session to the database
def get_session(dbpath: Path):
    """Returns live session to database."""
    db_URL = f"sqlite:///{dbpath}"
    engine = create_engine(db_URL)
    Session.configure(bind=engine)
    return Session()


if __name__ == "__main__":
    # Set up logging
    logger = logging.getLogger()
    logformatter = logging.Formatter("%(asctime)s %(message)s")
    logger.setLevel(logging.DEBUG)

    # Add a handler that writes to the console
    termhandler = logging.StreamHandler(sys.stdout)
    termhandler.setFormatter(logformatter)
    logger.addHandler(termhandler)

    # Path to output database
    outdbpath = Path("db.sqlite3")

    # Create database
    logger.info("Creating and populating database at %s", outdbpath)
    create_db(outdbpath)


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
