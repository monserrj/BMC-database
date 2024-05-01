#!/usr/bin/env python 

# This script is the development of the BMC db following the 
# BMC-database-skeleton.sql file. Instructions from create_Db.py file
# followed and explanations needed kept to help following through

# Import  SQLAlchemy classes needed with a declarative approach.

from sqlalchemy.orm import declarative_base
from sqlalchemy import Column, Integer, String, Table, ForeignKey, UniqueConstraint, VARCHAR
from sqlalchemy.orm import relationship, backref

# Create_engine function to create an engine object
from sqlalchemy import create_engine

# Import a sessionmaker to create a session object
from sqlalchemy.orm import sessionmaker


# Database creation:

# Create a base class to inherit from.
Base = declarative_base()

# Create a database engine to connect to the database.
# This creates a new empty database file called bmc.db in the current directory.
engine = create_engine("sqlite:///bmc.db")


# Create the tables in the database.
# Tables with one-to-many and many-to-many relationships must be created
# before creating other tables, to satisfy the logic of the code.


# need to create an unique constraint so protein_gene table field combinations
# are always unique (only 1 priority name)
proteingene = Table(
    "protein_gene",  # This name will be used in SQLite
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("gen_id", Integer, ForeignKey("gene.gen_id")),
    Column("name_rank", Integer), 
    # To enforce unique combinations of protein, gene ID and rank
    # to ensure several names from a protein/gene are
    # not made principal
    __table_args__ = (UniqueConstraint("prot_id", "gene_id", "name_rank"))
)

proteintaxon = Table(
    "protein_taxon",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column ("tax_id", Integer, ForeignKey("taxon.tax_id")),
)            
            
proteinpdb = Table(
    "protein_pdb",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("pdb_id", Integer, ForeignKey("pdb.pdb_id")),
)

proteindomain = Table(
    "protein_domain",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("dom_id", Integer, ForeignKey("domain.dom_id")),
)

isoforms = Table(
    "isoforms",
    Base.metadata,
    Column("canonical_prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("isoform_prot_id", Integer, ForeignKey("protein.prot_id")),
)

proteinGO = Table(
    "protein_GO",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("go_id", Integer, ForeignKey("function_GO.go_id")),
)

proteinpath = Table(
    "protein_path",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("path_id", Integer, ForeignKey("enzyme_path.path_id")),
)

proteincomplex = Table(
    "protein_complex",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("complex_id", Integer, ForeignKey("complex.complex_id")),
    Column("prot_essential_assembly", Integer),
    Column("interact_prot_ID", Integer, ForeignKey("protein.prot_id")),
    Column("copy_number", Integer),
)


# Create a class corresponding to each database table. No classes
# are instantiated until we create an instance of the Base class -
# this gets done after the tables are created.
# Each class inherits from the Base class we created earlier.
class Protein(Base):
    """Table representing a protein

    This table will store the protein sequence, links to accessions
    in other databases, and links to the type of protein structure.

    The protein sequence is taken to be unique, so we'll add a
    UniqueConstraint to the table to enforce this.
    """

    __tablename__ = "protein"  # this is the name that will be used in SQLite

    prot_id = Column(Integer, primary_key=True)  # primary key column
    prot_seq = Column(String, nullable=False, unique=True)  # sequence string
    locus_NCBI_ID = Column(VARCHAR)
    uniprot_ID = Column (VARCHAR)
    struct_prot_type = Column (Integer, nullable = True)
  
class Gene(Base):
    """Table representing a gene name and DNA sequence

    Each gene_ID represents a gene name and sequences. Several name strings
    given to an unique sequence
    """

    __tablename__ = "gene"
    
    gene_id = Column(Integer, primary_key=True)  # primary key column
    gene_name = Column(Integer, nullable=False)
    dna_seq = Column(String, nullable=False)

  # Define relationships after defining columns
    # A one-to-many relationship between Protein and Gene
    protein_gene = relationship("gene", secondary=proteingene)

class Taxon(Base):
    """Table representing the taxon accession of a protein

    This table will store the taxon origin of the protein sequence, e.g.:
    specie, genus, family... and a accession number to a database with details
    about that organism
    """

    __tablename__ = "taxon"
    tax_id = Column(Integer, primary_key=True)  # primary key column
    tax_ref = Column(VARCHAR, unique=True, nullable=False) # accession number in db
    tax_db = Column(String, nullable=False) # Name of database used (e.g.: NCBI, GTDB)
    specie = Column(String, nullable=False)
    genus = Column(String)
    family = Column(String)
    order_tax = Column(String)
    phylum = Column(String)
    class_tax = Column(String)
    strain = Column(String) 

# A many-to-one relationship between Protein and Taxonomy
    protein_taxon = relationship("taxon", secondary=proteintaxon)

# To enforce unique taxon references
    __table_args__ = (UniqueConstraint("tax_id", "tax_ref"))
        
class Pdb(Base):
    """Table representing the Pdb accession of a protein

    This table will store the different pdb accession number that represent the
    structure of a protein
    """

    __tablename__ = "pdb"
    pdb_id = Column(Integer, primary_key=True)  # primary key column
    pdb_acc_1 = Column(VARCHAR, unique=True) # primary accession number in pdb
    pdb_acc_2 = Column(VARCHAR, nullable=False) # accession number
    pdb_acc_3 = Column(VARCHAR, nullable=False) # accession number 

# A one-to-many relationship between Protein and Pdb structure
    protein_pdb = relationship("pdb", secondary=proteinpdb)  
    
class Domain(Base):
    """Table representing the conserved domain family of a protein

    This table will store the different conserved domain accession number that represent the
    structure of a protein, including the reference database where the accession number
    was taken
    """

    __tablename__ = "domain"
    dom_id = Column(Integer, primary_key=True)  # primary key column
    dom_ref = Column(VARCHAR, unique=True, nullable=False) # domain accession in external db
    dom_db = Column(Integer, nullable=False) # external database name e.g. pfam, CDD

# A many-to-many relationship between Protein and domain family
    protein_domain = relationship("domain", secondary=proteindomain)  
# To enforce unique domain family references
    __table_args__ = (UniqueConstraint("dom_id", "dom_ref"))
    
class Function(Base):
    """Table representing a function of a protein

    This table will store the different Gene Ontology accession numbers, type and
    description of the function.
    The types can be MF (molecular function), CC (cellular compartment),
    or BP (biological process)
    """

    __tablename__ = "Function"
    go_id = Column(Integer, primary_key=True)  # primary key column
    go_ref = Column(VARCHAR, unique=True, nullable=False) # accession number in GO
    go_type = Column(String, nullable=False) # GO type (MF,CC,BP)
    go_description = Column(String, nullable=False) # text description of function

# A many-to-many relationship between Protein and function
    protein_GO = relationship("function", secondary=proteinGO)  
# To enforce unique function references
    __table_args__ = (UniqueConstraint("go_id", "go_ref"))
    
class Enzyme_path(Base):
    """Table representing the enzymatic reaction in which the protein
    participates

    This table will store the different Kegg Ontlogy accession numbers
    related wit an specific function of a protein. Thus, a protein can have
    more than one reference
    """

    __tablename__ = "enzyme_path"
    path_id = Column(Integer, primary_key=True)  # primary key column
    KO_ref = Column(VARCHAR, unique=True, nullable=False) # accession number in KO

# A many-to-many relationship between Protein and enzymatic activity
    proteinpath = relationship("enzyme_path", secondary=proteinpath)  
# To enforce unique enzymatic pathway references
    __table_args__ = (UniqueConstraint("path_id", "KO_ref"))
    
 class Complex(Base):
    """Table representing the complex that can be form by the interaction
    between several proteins, including native BMC or engineered ones

    This table will store the complex features, including the type (e.g: pdu, eut),
    whether is has enzymatic activty or not, if it has been experimentally tested
    whether it assembles or not, the origin of the complex (meaning whether is it
    a native complex, engineered or created with a theorical or bioinformatic approach)
    """

    __tablename__ = "complex"
    complex_id = Column(Integer, primary_key=True)  # primary key column
    complex_type = Column(String) # Classification undecided (pdueut,grm..)
    complex_activity = Column(String, nullable=False) # Active/Inactive
    assembly_exp_tested = Column(String, nullable=False) #Y/N. If Y reference paper?
    complex_source = Column(String, nullable=False) #Native/engineered/theoretical...

# A many-to-many relationship between Protein and enzymatic activity
    proteincomplex = relationship("complex", secondary=proteincomplex)
# To enforce unique no repeated complexes are created
    __table_args__ = (UniqueConstraint("complex_id", "complex_type", "complex_activity", "assembly_exp_tested", "complex_source"))  