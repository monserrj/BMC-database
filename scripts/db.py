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

    # Unsure if needed, but to keep the code consistent with the rest of the db
    # Introduce all relationship between tables:
    references: Mapped["ProteinXref"] = relationship(back_populates="protein")
    cds: Mapped["Cds"] = relationship(back_populates="protein")
    names: Mapped["ProteinName"] = relationship(back_populates="protein")
    isoforms: Mapped["Isoforms"] = relationship(back_populates="isoform")
    canonicals: Mapped["Isoforms"] = relationship(back_populates="canonical")
    complexes: Mapped["ProteinComplex"] = relationship(back_populates="protein")
    interaction_1: Mapped["Prot_prot_interact"] = relationship(back_populates="protein_1")
    interaction_2: Mapped["Prot_prot_interact"] = relationship(back_populates="protein_2")


class Xdatabase(Base):
    """Descriptions of external databases (linked by Xref accessions)"""

    __tablename__ = "xdatabase"
    
    # Define table content:
    xref_db_id = Column(Integer, primary_key=True)
    xref_db = Column(String, nullable=False, unique=True)  # Database name
    xref_href = Column(String, nullable=False, unique=True)  # URL for database access
    xref_type = Column(
        String, nullable=False
    )  # types of database: seq, struct, funct, tax'

    # Introduce all relationship between tables:
    xref: Mapped["Xref"] = relationship(back_populates="xref_db_id")
    databases: Mapped[list["Xref"]] = relationship()

class Xref(Base):  # All external references
    """Table representing a singe database cross-reference"""

    __tablename__ = "xref"

    # Define table content:
    xref_id = Column(Integer, primary_key=True, autoincrement=True)
    xref_acc = Column(String, unique=True, nullable=False)  # Accession from external db
    xref_db_id: Mapped[int] = mapped_column(
        ForeignKey("xdatabase.xref_db_id"),
        nullable=False,
    )  # External reference db key

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

    #Define table content:
    cds_id = Column(Integer, primary_key=True)
    cds_seq = Column(String, nullable=False, unique=True)
    cds_accession = Column(
        String, nullable=False, unique=True
    )  # Unique accession number for CDS
    
    prot_id : Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), nullable=False
        ) 

    cds_origin : Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id"), nullable=True
    ) # Foreign key, origin DNA sequence (cds_id key)
    
    # Introduce all relationship between tables:
    protein: Mapped["Protein"] = relationship(back_populates="cds")
    origin : Mapped["Cds"] = relationship(back_populates="cds", remote_side=[cds_id], post_update=True)
    # It says to add this as self-referential foreign key
    # https://docs.sqlalchemy.org/en/20/orm/self_referential.html
    references: Mapped["CdsXref"] = relationship(back_populates="cds_xref")
    modifications: Mapped["CdsModif"] = relationship(back_populates="cds")

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
    
class Modif(Base):
    """Table representing the modification of the engineered CDS sequences
    This modifications will be vocabulary restricted: 
    truncated, fusion, synthetic, mutated, domesticated """
    
    __tablename__ = "modif"
       
    # Define table content:
    modif_id = Column(Integer, primary_key=True, autoincrement=True)
    modif_type = Column(String, nullable=False) # vocabulary restricted
    modif_description = Column(String, nullable=False) # Description of the specific modification
    
    UniqueConstraint (modif_type, modif_description)

    # Introduce all relationship between tables:
    modifications: Mapped["CdsModif"] = relationship(back_populates="modification")

class CdsModif(Base):
    """Linker table, CDS to modification cross-reference"""

    __tablename__ = "cds_modif"
    __table_args__ = (PrimaryKeyConstraint("cds_id", "modif_id"),)

    # Define table content:
    cds_id: Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id")
    )
    modif_id: Mapped[int] = mapped_column(
        ForeignKey("modif.modif_id")
    )
    
    # Introduce all relationship between tables:
    modification: Mapped["Modif"] = relationship(back_populates="modifications")
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

class Isoforms (Base):
    """Table representing the isoforms of a protein
    This table will store the canonical protein ID and the isoform protein ID.
    The canonical protein ID is the one that is considered the main protein,
    while the isoform protein ID is a variant of the canonical protein.
    """

    __tablename__ = "isoforms"
    __table_args__ = (PrimaryKeyConstraint("canonical__prot_id", "isoform__prot_id"),)

    # Define table content:
    canonical__prot_id : Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    isoform__prot_id : Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
        )
    
    # Introduce all relationship between tables:
    isoform : Mapped["Protein"] = relationship(back_populates="isoforms")
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
    complex_name = Column(String, nullable=False)  # Name of the complex
    complex_type = Column(String, nullable=False) # Classification undecided (pdu,eut,grm..)
    is_active = Column(String) # Active/Inactive. Boolean?
    is_exp_tested = Column(String) #Y/N. Boolean?
    complex_source = Column(String, nullable=False) 
    # Native/engineered/theoretical...
    
    UniqueConstraint("complex_name", "complex_type", "is_active",
                     "is_exp_tested", "complex_source")
    
    # Introduce all relationship between tables:
    proteins: Mapped["ProteinComplex"] = relationship(
        "ProteinComplex", secondary="protein_complex", back_populates="complexes", lazy="dynamic"
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

class Prot_prot_interact(Base):
    """Table representing the different specific interactions
    between several proteins inside a complex and during assembly
    This table will store the different proteins known to interact
    specifically with each other.
    """

    __tablename__ = "prot_prot_interact"
    
    # Define table content:
    ppi_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    prot_id_1: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )
    prot_id_2: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"),
    )

    # To enforce unique no repeated protein to protein interactions are created
    UniqueConstraint("prot_id_1", "prot_id_2")
    
    # Introduce all relationship between tables:
    protein_1: Mapped["Protein"] = relationship(back_populates="interaction_1")
    protein_2: Mapped["Protein"] = relationship(back_populates="interaction_2")


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
    create_db(Path("db.sqlite3"))











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
