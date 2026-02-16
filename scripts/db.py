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

from enums import DatabaseType, StructProtType, enum_from_str

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
    String,
    or_,
    create_engine,
)
from sqlalchemy.orm import (
    DeclarativeBase,
    relationship,
    Mapped,
    mapped_column,
    sessionmaker,
)

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
    xref_db_name: Mapped[str] = mapped_column(
        nullable=False, unique=True
    )  # Database with crossref
    xref_href: Mapped[str] = mapped_column(
        nullable=False, unique=True
    )  # URL to database
    #    xref_type = Column(SQLEnum(DatabaseType, name ="database_type_enum", nullable=False))   # Not currently defined

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

    #    modification_type = Column(SQLEnum(ModificationType, name="modification_type_enum"), nullable=False)   # Not currently defined

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


# # Function for protein data addition
def protein_addition(session, protacc, protseq, struct, canonical):
    """
    Args:
    protseq (str): Protein sequence.
    protacc (str): Protein unique accession number
    struct (str): Protein structure type.
    canonical (boolean): Canonical if true, false if isoform

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
    logger.debug(f"\nNow in {protein_addition.__name__}")

    print(f"Before query, {protseq[:10]=}..., {protacc=}, {struct=}, {canonical=}")

    # Create a new protein object
    try:
        protein = session.query(Protein).filter(Protein.prot_seq == protseq).first()
        print(f"After query, {protein=}")

        # Enforce protein structure type
        if struct.upper() not in StructProtType:  # invalid protein structure type
            logger.warning(
                "Skipping invalid struct type: %s (accession: %s)", struct, protacc
            )
        else:
            struct = StructProtType[struct.upper()]

        if protein:
            logger.info(
                "Protein already exists: accession=%s id=%s",
                protacc,
                protein.prot_id,
            )
        # Add protein if it is not already present
        if not protein:
            protein = Protein(
                prot_seq=protseq,
                prot_accession=protacc,
                struct_prot_type=struct,
                is_canonical=canonical,
            )
            session.add(protein)
            session.flush()  # This sends the changes to the database, so prot_id is assigned
            logger.info(
                "Protein %s added with accession %s",
                protein.prot_id,
                protacc,
            )
            return protein
    except Exception as exc:
        logger.exception("Failed to add protein accession=%s", protacc)
        logger.exception(exc)
        session.rollback()
        raise


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
    from eralchemy import render_er

    render_er(Base, "er_diagram.pdf")

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
