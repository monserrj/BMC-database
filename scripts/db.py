#!/usr/bin/env python

# This script is the development of the BMC db. Instructions
# followed and explanations needed kept to help following through

# Import  SQLAlchemy classes needed with a declarative approach.
from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy import (
    Column,
    Integer,
    String,
    ForeignKey,
    UniqueConstraint,
    String,
)
from sqlalchemy.orm import relationship, Mapped, mapped_column

# Create_engine function to create an engine object
from sqlalchemy import create_engine

from typing import Optional

# Database set up:

# Create a base class to inherit from.
Base = declarative_base()

# Create a database engine to connect to the database.
# This creates a new empty database file called bmc.db in the current directory.
db_URL = "sqlite:///bmc.db"
engine = create_engine(db_URL)

# # Session set up:
# Session = sessionmaker()  # we also need a session object
# Session.configure(bind=engine)
# session = Session()

# Database creation:

# Create the tables in the database.
# Tables with one-to-many and many-to-many relationships must be created
# before creating other tables, to satisfy the logic of the code.
class ProteinName(Base):
    __tablename__ = "protein_name"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), primary_key=True
    )
    name_id: Mapped[int] = mapped_column(
        ForeignKey("name.name_id"), primary_key=True
    )  # Not primary key as can be repeated
    name_rank: Mapped[Optional[int]]
    name: Mapped["Name"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="names")
    # UniqueConstraint("prot_id", "name_id", "name_rank")


class ProteinTaxonomy(Base):
    __tablename__ = "protein_taxonomy"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), primary_key=True
    )
    tax_id: Mapped[int] = mapped_column(ForeignKey("taxonomy.tax_id"), primary_key=True)
    taxonomy: Mapped["Taxonomy"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="taxonomies")
    # UniqueConstraint("prot_id", "tax_id")


class ProteinPdb(Base):
    __tablename__ = "protein_pdb"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), primary_key=True
    )
    pdb_id: Mapped[int] = mapped_column(
        ForeignKey("pdb.pdb_id")
    )  # Recheck if this is primary key too
    pdb: Mapped["Pdb"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="pdbs")
    UniqueConstraint("prot_id", "pdb_id")


class ProteinDomain(Base):
    __tablename__ = "protein_domain"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), primary_key=True
    )
    dom_id: Mapped[int] = mapped_column(ForeignKey("domain.dom_id"))
    domain: Mapped["Domain"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="domains")
    UniqueConstraint("prot_id", "dom_id")


# class Isoforms (Base):
#     __tablename__ = "isoforms"
#     canonical__prot_id : Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), primary_key=True)
#     isoform__prot_id : Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), primary_key=True)
#     isoform : Mapped["Protein"] = relationship(back_populates="isoforms")
#     canonical: Mapped["Protein"] = relationship(back_populates="canonicals") # Recheck this


class ProteinFunction(Base):
    __tablename__ = "protein_function"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), primary_key=True
    )
    go_id: Mapped[int] = mapped_column(ForeignKey("function.go_id"))
    function: Mapped["Function"] = relationship(back_populates="proteins")
    UniqueConstraint("prot_id", "go_id")


class ProteinPath(Base):
    __tablename__ = "protein_path"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), primary_key=True
    )
    path_id: Mapped[int] = mapped_column(ForeignKey("enzyme_path.path_id"))
    path: Mapped["EnzymePath"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="paths")
    UniqueConstraint("prot_id", "path_id")


# proteincomplex = Table(
#      "protein_complex",
#      Base.metadata,
#     Column("prot_id", Integer, ForeignKey("protein.prot_id")),
#     Column("complex_id", Integer, ForeignKey("complex.complex_id")),
#     Column("prot_essential_assembly", Integer),
#     Column("ppi_ID", Integer, ForeignKey("protein.prot_id")),
#     Column("copy_number", Integer),
# )


# Create a class corresponding to each database table. No classes
# are instantiated until we create an instance of the Base class -
# this gets done after the tables are created.
# Each class inherits from the Base class we created earlier.
class Protein(Base):
    """Table representing a protein This table will store the protein sequence, links to accessions
    in other databases, and links to the type of protein structure.
    The protein sequence is taken to be unique, so we'll add a
    UniqueConstraint to the table to enforce this.
    """

    __tablename__ = "protein"  # this is the name that will be used in SQLite
    # Introduce all relationship between tables:
    names: Mapped[list["ProteinName"]] = relationship()
    taxonomies: Mapped[list["ProteinTaxonomy"]] = relationship()
    pdbs: Mapped[list["ProteinPdb"]] = relationship()
    domains: Mapped[list["ProteinDomain"]] = relationship()
    # isoforms: Mapped[list["Isoforms"]] = relationship()
    # canonicals: Mapped[list["Isoforms"]] = relationship() # I am sure this is wrong, fix when TS
    functions: Mapped[list["ProteinFunction"]] = relationship()
    paths: Mapped[list["ProteinPath"]] = relationship()
    # Define table content:
    prot_id = Column(Integer, primary_key=True, autoincrement=True)
    prot_seq = Column(String, nullable=False, unique=True)
    locus_NCBI_id = Column(String, unique=True, nullable=True)
    uniprot_id = Column(String, unique=True, nullable=True)
    struct_prot_type = Column(Integer, nullable=True)
    dna_seq = Column(
        String, nullable=False, unique=True
    )  # New addition for testing Name instead of gen table
    # Introduce back_populates so when a relationship between different tables is
    # introduced, they information will be back-populated to be consistent across
    # all tables. Relationships must be introduced in both related tables (e.g.:
    # Gene and Protein, with relationship based on table ProteinGene)
    # complexes = relationship(
    #     "Complex", secondary=proteincomplex, back_populates="proteins", lazy="dynamic"
    #     )
    # interacts = relationship(
    #     "Protprotinteract",  secondary=proteincomplex, back_populates="proteins", lazy="dynamic"
    #     )
    # Define type of output for protein
    def __str__(self):
        outstr = [
            f"Protein ID: {self.prot_id}",
            f"Protein sequence: {self.prot_seq}",
            f"NCBI ID: {self.locus_NCBI_id}",
            f"Uniprot ID: {self.uniprot_id}",
            f"Protein structure type: {self.struct_prot_type}",
            f"DNA sequence: {self.dna_seq}",
        ]
        return "\n".join(outstr)


class Name(Base):
    """Table representing a gene name
    Each gene_ID represents a gene name. Several name strings
    given to an unique protein, and several proteins sharing same name
    """

    __tablename__ = "name"
    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinName"]] = relationship()
    # Define table content:
    name_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column
    gene_name = Column(String, nullable=False)
    # dna_seq = Column(String, nullable=False) removed as it will be part of the protein identity


class Taxonomy(Base):
    """Table representing the taxon accession of a protein
    This table will store the taxon origin of the protein sequence, e.g.:
    specie, genus, family... and a accession number to a database with details
    about that organism
    """

    __tablename__ = "taxonomy"
    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinTaxonomy"]] = relationship()
    # Define table content:
    tax_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    tax_ref = Column(String, unique=True, nullable=False)  # accession number in db
    tax_db = Column(String, nullable=False)  # Name of database used (e.g.: NCBI, GTDB)
    species = Column(String, nullable=False)
    genus = Column(String)
    family = Column(String)
    order_tax = Column(String)
    phylum = Column(String)
    class_tax = Column(String)
    strain = Column(String)
    # To enforce unique taxonomy references
    __table_args__ = (UniqueConstraint("species", "strain"),)


class Pdb(Base):
    """Table representing the Pdb accession of a protein
    This table will store the different pdb accession number that represent the
    structure of a protein
    """

    __tablename__ = "pdb"
    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinPdb"]] = relationship()
    # Define table content:
    pdb_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    pdb_acc_1 = Column(String, unique=True)  # primary accession number in pdb
    pdb_acc_2 = Column(String, nullable=False)  # accession number
    pdb_acc_3 = Column(String, nullable=False)  # accession number


class Domain(Base):
    """Table representing the conserved domain family of a protein
    This table will store the different conserved domain accession number that represent the
    structure of a protein, including the reference database where the accession number
    was taken
    """

    __tablename__ = "domain"
    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinDomain"]] = relationship()
    # Define table content:
    dom_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    dom_ref = Column(
        String, unique=True, nullable=False
    )  # domain accession in external db
    dom_db = Column(Integer, nullable=False)  # external database name e.g. pfam, CDD
    # To enforce unique domain family references
    __table_args__ = (UniqueConstraint("dom_id", "dom_ref"),)


class Function(Base):
    """Table representing a function of a protein
    This table will store the different Gene Ontology accession numbers, type and
    description of the function.
    The types can be MF (molecular function), CC (cellular compartment),
    or BP (biological process)
    """

    __tablename__ = "function"
    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinFunction"]] = relationship()
    # Define table content:
    go_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    go_ref = Column(String, unique=True, nullable=False)  # accession number in GO
    go_type = Column(String, nullable=False)  # GO type (MF,CC,BP)
    go_description = Column(String, nullable=False)  # text description of function
    # To enforce unique function references
    __table_args__ = (UniqueConstraint("go_id", "go_ref"),)


class EnzymePath(Base):
    """Table representing the enzymatic reaction in which the protein
    participates
    This table will store the different Kegg Ontology accession numbers
    related wit an specific function of a protein. Thus, a protein can have
    more than one reference
    """

    __tablename__ = "enzyme_path"
    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinPath"]] = relationship()
    # Define table content:
    path_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column
    KO_ref = Column(String, unique=True, nullable=False)  # accession number in KO
    # To enforce unique enzymatic pathway references
    __table_args__ = (UniqueConstraint("path_id", "KO_ref"),)


# class Complex(Base):
#     """Table representing the complex that can be form by the interaction
#     between several proteins, including native BMC or engineered ones
#     This table will store the complex features, including the type (e.g: pdu, eut),
#     whether is has enzymatic activty or not, if it has been experimentally tested
#     whether it assembles or not, the origin of the complex (meaning whether is it
#     a native complex, engineered or created with a theorical or bioinformatic approach)
#     """
#     __tablename__ = "complex"
#     complex_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
#     complex_type = Column(String) # Classification undecided (pdueut,grm..)
#     complex_activity = Column(String, nullable=False) # Active/Inactive
#     assembly_exp_tested = Column(String, nullable=False) #Y/N. If Y reference paper?
#     complex_source = Column(String, nullable=False) #Native/engineered/theoretical...
#     # A many-to-many relationship between Protein and enzymatic activity
#     proteins = relationship(
#         "Protein", secondary=proteincomplex, back_populates="complexes", lazy="dynamic")
# # To enforce unique no repeated complexes are created
#     __table_args__ = (UniqueConstraint("complex_id", "complex_type", "complex_activity", "assembly_exp_tested", "complex_source"),)
# class Prot_port_interact(Base):
#     """Table representing the different specific interactions
#     between several proteins inside a complex and during assembly
#     This table will store the different proteins known to interact specifically.
#     Prot_ID_1 will be the main protein, with the rest of the proteins being the
#     ones interacting with it
#     """
#     __tablename__ = "prot_prot_interact"
#     ppi_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
#     prot_id_1 = Column(Integer, ForeignKey("Protein.prot_id"))
#     prot_id_2 = Column(Integer, ForeignKey("Protein.prot_id"))
#     prot_id_3 = Column(Integer, ForeignKey("Protein.prot_id"))
#     prot_id_4 = Column(Integer, ForeignKey("Protein.prot_id"))
#     prot_id_5 = Column(Integer, ForeignKey("Protein.prot_id"))
#     prot_id_6 = Column(Integer, ForeignKey("Protein.prot_id"))
#     prot_id_7 = Column(Integer, ForeignKey("Protein.prot_id"))
#     # A many-to-many relationship between Protein and enzymatic activity
#     proteins = relationship(
#         "Protein", secondary=proteincomplex, back_populates="interacts", lazy="dynamic")
# # To enforce unique no repeated protein to protein interactions are created
#     __table_args__ = (UniqueConstraint("prot_id_1", "prot_id_2", "prot_id_3", "prot_id_4", "prot_id_5", "prot_id_6", "prot_id_7"),)


def create_db():
    """Function to create all the tables from the database"""
    Base.metadata.create_all(bind=engine)
    print("Database and tables created successfully")

if __name__ =="__main__":
    create_db()