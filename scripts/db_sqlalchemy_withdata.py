#!/usr/bin/env python

# This script is the development of the BMC db following the 
# BMC-database-skeleton.sql file. Instructions from create_Db.py file
# followed and explanations needed kept to help following through

# Open .csv to be able to work with .csv files
import csv
# Import path to be able to open .csv files in different folders
from pathlib import Path

# Import  SQLAlchemy classes needed with a declarative approach.
from sqlalchemy.orm import declarative_base
from sqlalchemy import Column, Integer, String, Table, ForeignKey, UniqueConstraint, String
from sqlalchemy.orm import relationship

# Create_engine function to create an engine object
from sqlalchemy import create_engine

# Import a sessionmaker to create a session object
from sqlalchemy.orm import sessionmaker

# To be able to make exceptions in code (try/except): 
# from sqlalchemy.exc import IntegrityError, PendingRollbackError
# For exitying system for trouble shooting import sys
import sys

# Database creation:

# Create a base class to inherit from.
Base = declarative_base()

# Create a database engine to connect to the database.
# This creates a new empty database file called bmc.db in the current directory.
engine = create_engine("sqlite:///bmc.db")


# Create the tables in the database.
# Tables with one-to-many and many-to-many relationships must be created
# before creating other tables, to satisfy the logic of the code.


# Need to create an unique constraint so protein_gene table field combinations
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
   UniqueConstraint("prot_id", "gen_id", "name_rank"), # Remove table_args if not class
)

proteintaxon = Table(
    "protein_taxon",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column ("tax_id", Integer, ForeignKey("taxon.tax_id")),
    UniqueConstraint("prot_id", "tax_id"),
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
    Column("go_id", Integer, ForeignKey("function.go_id")),
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
    Column("interact_prot_id", Integer, ForeignKey("protein.prot_id")),
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
    locus_NCBI_id = Column(String)
    uniprot_id = Column (String)
    struct_prot_type = Column (Integer, nullable = True)
    # Introduce back_populates so when a relationship between different tables is
    # introduced, they information will be backpopulated to be consistant accross
    # all tables. Relationships must be introduced in both related tables (e.g.: 
    # Gene and Protein, with relationship based on table Proteingene)
    gene = relationship("Gene", secondary=proteingene,
                                back_populates="protein",
                                lazy="dynamic")
    taxon = relationship("Taxonomy", secondary=proteintaxon,
                                back_populates="protein",
                                lazy="dynamic")
    pdb = relationship("Pdb", secondary=proteinpdb,
                                back_populates="protein",
                                lazy="dynamic")
    domain = relationship("Domain", secondary=proteindomain,
                                back_populates="protein",
                                lazy="dynamic")
    function = relationship("Function", secondary=proteinGO,
                                back_populates="protein",
                                lazy="dynamic")
    path = relationship("Enzyme_path", secondary=proteinpath,
                                back_populates="protein",
                                lazy="dynamic")
    complexmain = relationship("Complex", foreign_keys=['prot_id'],
                                back_populates="protein",
                                lazy="dynamic")
    complexinteract = relationship("Complex", foreign_keys=['interact_prot_id'],
                                back_populates="protein",
                                lazy="dynamic")
class Gene(Base):
    """Table representing a gene name and DNA sequence

    Each gene_ID represents a gene name and sequences. Several name strings
    given to an unique sequence
    """

    __tablename__ = "gene"
    
    gen_id = Column(Integer, primary_key=True)  # primary key column
    gen_name = Column(Integer, nullable=False)
    dna_seq = Column(String, nullable=False)

    # Define relationships after defining columns
    # A one-to-many relationship between Protein and Gene
    protein = relationship("Protein", secondary=proteingene,
                                back_populates="gene",
                                lazy="dynamic")

class Taxonomy(Base):
    """Table representing the taxon accession of a protein

    This table will store the taxon origin of the protein sequence, e.g.:
    specie, genus, family... and a accession number to a database with details
    about that organism
    """

    __tablename__ = "taxon"
    tax_id = Column(Integer, primary_key=True)  # primary key column
    tax_ref = Column(String, unique=True, nullable=False) # accession number in db
    tax_db = Column(String, nullable=False) # Name of database used (e.g.: NCBI, GTDB)
    species = Column(String, nullable=False)
    genus = Column(String)
    family = Column(String)
    order_tax = Column(String)
    phylum = Column(String)
    class_tax = Column(String)
    strain = Column(String) 

    # A many-to-one relationship between Protein and Taxonomy
    protein = relationship("Protein", secondary=proteintaxon,
                                back_populates="taxon",
                                lazy="dynamic")
    
    # To enforce unique taxon references
    __table_args__ = (UniqueConstraint("tax_id", "tax_ref"),)
        
class Pdb(Base):
    """Table representing the Pdb accession of a protein

    This table will store the different pdb accession number that represent the
    structure of a protein
    """

    __tablename__ = "pdb"
    pdb_id = Column(Integer, primary_key=True)  # primary key column
    pdb_acc_1 = Column(String, unique=True) # primary accession number in pdb
    pdb_acc_2 = Column(String, nullable=False) # accession number
    pdb_acc_3 = Column(String, nullable=False) # accession number 

    # A one-to-many relationship between Protein and Pdb structure
    protein = relationship("Protein", secondary=proteinpdb,
                                back_populates="pdb",
                                lazy="dynamic") 
    
class Domain(Base):
    """Table representing the conserved domain family of a protein

    This table will store the different conserved domain accession number that represent the
    structure of a protein, including the reference database where the accession number
    was taken
    """

    __tablename__ = "domain"
    dom_id = Column(Integer, primary_key=True)  # primary key column
    dom_ref = Column(String, unique=True, nullable=False) # domain accession in external db
    dom_db = Column(Integer, nullable=False) # external database name e.g. pfam, CDD

    # A many-to-many relationship between Protein and domain family
    protein = relationship("Protein", secondary=proteindomain,
                                back_populates="domain",
                                lazy="dynamic")

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
    go_id = Column(Integer, primary_key=True)  # primary key column
    go_ref = Column(String, unique=True, nullable=False) # accession number in GO
    go_type = Column(String, nullable=False) # GO type (MF,CC,BP)
    go_description = Column(String, nullable=False) # text description of function

    # A many-to-many relationship between Protein and function
    protein = relationship("Protein", secondary=proteinGO,
                                back_populates="function",
                                lazy="dynamic")
    # To enforce unique function references
    __table_args__ = (UniqueConstraint("go_id", "go_ref"),)
    
class Enzyme_path(Base):
    """Table representing the enzymatic reaction in which the protein
    participates

    This table will store the different Kegg Ontlogy accession numbers
    related wit an specific function of a protein. Thus, a protein can have
    more than one reference
    """

    __tablename__ = "enzyme_path"
    path_id = Column(Integer, primary_key=True)  # primary key column
    KO_ref = Column(String, unique=True, nullable=False) # accession number in KO

    # A many-to-many relationship between Protein and enzymatic activity
    protein = relationship("Protein", secondary=proteinpath,
                                back_populates="enzyme_path",
                                lazy="dynamic")
    # To enforce unique enzymatic pathway references
    __table_args__ = (UniqueConstraint("path_id", "KO_ref"),)
    
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
    proteins = relationship("Protein", secondary=proteincomplex,
                                back_populates="complex",
                                lazy="dynamic")

    # To enforce unique no repeated complexes are created
    __table_args__ = (UniqueConstraint("complex_id", "complex_type", "complex_activity", "assembly_exp_tested", "complex_source"),)

# Now that we have defined the tables, we can create the tables in the
# database.
Base.metadata.create_all(engine)

# Add some data to populate the database.
# Using csv files for data.
# Add the data to the database in a loop, but we'll have to
# check if the data entered already exist
# and update the corresponding tables accordingly.


# Open the csv files:
# All files are stored in the same directory
# datadir = Path("../data/raw")

# Open CSV file prot_data
# Define path for all data files
raw_dir = Path(__file__).resolve().parent.parent / 'data'/'raw'
#Define path to the prot_data file
prot_data_file = raw_dir / 'prot_data.csv'

mydata = []
with open(prot_data_file, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip header row
    for row in reader:
        mydata.append(tuple(row))

# Start the session
Session = sessionmaker()  # we also need a session object
Session.configure(bind=engine)
session = Session()

# Add the data to the database
# We need to check if the sequence, structure, and accession already exist,
# and update the corresponding tables accordingly if they do not.
# We can then update the linker tables by adding the corresponding items.
for (prot, protseq, NCBIid, uniprot, struct, gen, name, dnaseq, 
     tax, taxref, taxdb, spec, genu, fam, order, phyl, 
     classt, stra, pdbid, pdb_1, pdb_2, pdb_3, path, KOid) in mydata:
    
    # Check what data is available:
    # print(prot,protseq,NCBIid)

#         # Create a new protein object
#     new_prot = (
#         session.query(Protein)
#         .filter(Protein.prot_id == prot)
#         .filter(Protein.prot_seq == protseq)
#         .filter(Protein.locus_NCBI_id == NCBIid)
#         .filter(Protein.uniprot_id == uniprot)
#         .filter(Protein.struct_prot_type == struct)
#         .first()
#         )
#     if not isinstance(new_prot, Protein):
#         new_prot = Protein(prot_id=prot,
#                         prot_seq = protseq,
#                         locus_NCBI_id = NCBIid,
#                         uniprot_id = uniprot,
#                         struct_prot_type = struct)
#         session.add(new_prot)
#         session.commit()

#     # Create a new gene object
#     new_gene = (
#         session.query(Gene)
#         .filter(Gene.gen_id == gen)
#         .filter(Gene.gen_name == name)
#         .filter(Gene.dna_seq == dnaseq)
#         .first()
#     )
#     if not isinstance(new_gene, Gene):
#         new_gene = Gene(gen_id=gen, gen_name = name, dna_seq = dnaseq)
#         session.add(new_gene)
#         session.commit()

#     # Create new tax
#     # new_tax = (
#     #     session.query(Taxonomy)
#     #     .filter(Taxonomy.tax_id == tax)
#     #     .filter(Taxonomy.tax_ref == taxref)
#     #     .filter(Taxonomy.tax_db == taxdb)
#     #     .filter(Taxonomy.species == spec)
#     #     .filter(Taxonomy.genus == genu)
#     #     .filter(Taxonomy.family == fam)
#     #     .filter(Taxonomy.order_tax == order)
#     #     .filter(Taxonomy.phylum == phyl)
#     #     .filter(Taxonomy.class_tax == classt)
#     #     .filter(Taxonomy.strain == stra)
#     #     .first
#     # )
#     # if not isinstance(new_tax, Taxonomy):
#     #     new_tax = Taxonomy(tax_id = tax,
#     #                     tax_ref = taxref,
#     #                     ltax_db = taxdb,
#     #                     species = spec,
#     #                     genus = genu,
#     #                     family = fam,
#     #                     order_tax = order,
#     #                     phylum = phyl,
#     #                     class_tax = classt,
#     #                     strain = stra)
#         # session.add(new_tax)
#         # session.commit()

#     new_pdb = (
#         session.query(Pdb)
#         .filter(Pdb.pdb_id == pdbid)
#         .filter(Pdb.pdb_acc_1 == pdb_1)
#         .filter(Pdb.pdb_acc_2 == pdb_2)
#         .filter(Pdb.pdb_acc_3 == pdb_3)
#         .first
#     )
#     if not isinstance(new_pdb, Pdb):
#         new_pdb = Pdb(pdb_id = pdbid, pdb_acc_1 = pdb_1,
#                     pdb_acc_2 = pdb_2, pdb_acc_3 = pdb_3)
#         session.add(new_pdb)
#         session.commit()

#     new_path = (
#         session.query(Enzyme_path)
#         .filter(Enzyme_path.path_id == path)
#         .filter(Enzyme_path.KO_ref == KOid)
#         .first
#     )
#     if not isinstance(new_path, Enzyme_path):
#         new_path = Enzyme_path(path_id = path, KO_ref = KOid)
#         session.add(new_path)
#         session.commit()

#     # Add the protein to the Protein_gene
#     new_prot.Protein_gene.append(new_prot)
#     session.commit()
#     # Should I add all of them?



# # # Now we can query the database to see if the data has been added correctly
# # # Unsure whether I understand this correctly
# # for protein in session.query(Protein):
# #     print(f"\nPROTEIN: {protein.prot_id}")
# #     print("PROTEIN AND GENES:")
# #     for gen in protein.protein_gene:
# #         print(f"\t{gen.gen_id}, {gen.prot_id}")
# #     print("GENES:")
# #     for gene in protein.gene:
# #         print(f"\t{gene.gen_seq}, {gene.gen_id}, {gene.gen_name}")