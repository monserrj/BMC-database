#!/usr/bin/env python

# This script is the development of the BMC db. Instructions
# followed and explanations needed kept to help following through
from pathlib import Path  # for type hints

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

    # Introduce all relationship between tables:
    # names: Mapped[list["ProteinName"]] = relationship()
    # cds: Mapped[list["Cds"]] = relationship()
    # xref_prots: Mapped[list["ProteinXref"]] = relationship()

    # I am sure this is wrong, fix when TS
    # isoforms: Mapped[list["Isoforms"]] = relationship()
    # canonicals: Mapped[list["Isoforms"]] = relationship()
    # references: Mapped[list["ProteinXref"]] = relationship()

    # Define table content:
    prot_id = Column(Integer, primary_key=True, autoincrement=True)
    prot_accession = Column(
        Integer, nullable=False, unique=True
    )  # Unique accession number for proteins
    prot_seq = Column(String, nullable=False, unique=True)
    struct_prot_type = Column(
        Integer, nullable=True
    )  # BMC type H/T/P or non structural
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


class Xdatabase(Base):
    """Descriptions of external databases (linked by Xref accessions)"""

    __tablename__ = "xdatabase"

    databases: Mapped[list["Xref"]] = relationship()

    xref_db_id = Column(Integer, primary_key=True)
    xref_db = Column(String, nullable=False, unique=True)  # Database name
    xref_href = Column(String, nullable=False, unique=True)  # URL for database access
    xref_type = Column(
        String, nullable=False
    )  # types of database: seq, struct, funct, tax'

    xref: Mapped["Xref"] = relationship(back_populates="xref_db_id")


class Xref(Base):  # All external references
    """Table representing a singe database cross-reference"""

    __tablename__ = "xref"

    # Check if this can be done
    proteins: Mapped[list["ProteinXref"]] = relationship()
    # cdss: Mapped[list["CdsXref"]] = relationship()

    # Define table content:
    xref_id = Column(Integer, primary_key=True, autoincrement=True)
    xref_acc = Column(String, unique=True, nullable=False)  # Accession from external db
    xref_db_id: Mapped[int] = mapped_column(
        ForeignKey("xdatabase.xref_db_id"),
        nullable=False,
    )  # External reference db key

    xref_db: Mapped["Xdatabase"] = relationship(back_populates="databases")


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

    protein: Mapped["Protein"] = relationship(back_populates="references")
    xref: Mapped["Xref"] = relationship(back_populates="proteins")


def create_db(dbpath: Path):
    """Function to create all the tables from the database"""
    # Create a database engine to connect to the database.
    # This creates a new empty database file called bmc.db in the current directory.
    db_URL = f"sqlite:///{dbpath}"
    engine = create_engine(db_URL)
    Base.metadata.create_all(bind=engine)
    print("Database and tables created successfully")


def get_session(dbpath: Path):
    """Returns live session to database."""
    db_URL = f"sqlite:///{dbpath}"
    engine = create_engine(db_URL)
    Session.configure(bind=engine)
    return Session()


if __name__ == "__main__":
    create_db(Path("db-lp.sqlite3"))


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
