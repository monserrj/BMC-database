#!/usr/bin/env python

# This script is the development of the BMC db following the
# BMC-database-skeleton.sql file. Instructions from create_Db.py file
# followed and explanations needed kept to help following through

# Open .csv to be able to work with .csv files
import csv

# Import path to be able to open .csv files in different folders
from pathlib import Path

from typing import Optional

# Import  SQLAlchemy classes needed with a declarative approach.
from sqlalchemy.orm import declarative_base
from sqlalchemy import (
    Column,
    Integer,
    String,
    Table,
    ForeignKey,
    UniqueConstraint,
    String,
)
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

# Temporary code to delete local database while we debug
Path("bmc.db").unlink()

# Create a database engine to connect to the database.
# This creates a new empty database file called bmc.db in the current directory.
engine = create_engine("sqlite:///bmc.db")

# Create the tables in the database.
# Tables with one-to-many and many-to-many relationships must be created
# before creating other tables, to satisfy the logic of the code.

# Need to create an unique constraint so protein_gene table field combinations
# are always unique (only 1 priority name)
# proteingene = Table(
#     "protein_gene",  # This name will be used in SQLite
#     Base.metadata,
#     Column("prot_id", Integer, ForeignKey("protein.prot_id")),
#     Column("gene_id", Integer, ForeignKey("gene.gene_id")),
#     Column("name_rank", Integer),
#     # To enforce unique combinations of protein, gene ID and rank
#     # to ensure several names from a protein/gene are
#     # not made principal
#     UniqueConstraint(
#         "prot_id", "gene_id", "name_rank"
#     ),  # Remove table_args if not class
# )


from sqlalchemy.orm import Mapped, mapped_column


class ProteinGene(Base):
    __tablename__ = "protein_gene"
    protid: Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), primary_key=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("gene.gene_id"), primary_key=True)
    name_rank: Mapped[Optional[int]]
    gene: Mapped["Gene"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="genes")


proteintaxon = Table(
    "protein_taxon",
    Base.metadata,
    Column("prot_id", Integer, ForeignKey("protein.prot_id")),
    Column("tax_id", Integer, ForeignKey("taxon.tax_id")),
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
    Column("path_id", Integer, ForeignKey("enzymepath.path_id")),
)

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
    """Table representing a protein

    This table will store the protein sequence, links to accessions
    in other databases, and links to the type of protein structure.

    The protein sequence is taken to be unique, so we'll add a
    UniqueConstraint to the table to enforce this.
    """

    __tablename__ = "protein"  # this is the name that will be used in SQLite
    genes: Mapped[list["ProteinGene"]] = relationship()

    prot_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column. Added autoincrement
    prot_seq = Column(String, nullable=False, unique=True)  # sequence string
    locus_NCBI_id = Column(String)
    uniprot_id = Column(String, unique=True)
    struct_prot_type = Column(Integer, nullable=True)
    # Introduce back_populates so when a relationship between different tables is
    # introduced, they information will be backpopulated to be consistant accross
    # all tables. Relationships must be introduced in both related tables (e.g.:
    # Gene and Protein, with relationship based on table Proteingene)
    # genes = relationship(
    #     "Gene", secondary=proteingene, back_populates="proteins", lazy="dynamic"
    # )
    taxons = relationship(
        "Taxonomy", secondary=proteintaxon, back_populates="proteins", lazy="dynamic"
    )
    pdbs = relationship(
        "Pdb", secondary=proteinpdb, back_populates="proteins", lazy="dynamic"
    )
    domains = relationship(
        "Domain", secondary=proteindomain, back_populates="proteins", lazy="dynamic"
    )
    functions = relationship(
        "Function", secondary=proteinGO, back_populates="proteins", lazy="dynamic"
    )
    paths = relationship(
        "Enzymepath", secondary=proteinpath, back_populates="proteins", lazy="dynamic"
    )
    # complexes = relationship(
    #     "Complex", secondary=proteincomplex, back_populates="proteins", lazy="dynamic"
    #     )
    # interacts = relationship(
    #     "Protprotinteract",  secondary=proteincomplex, back_populates="proteins", lazy="dynamic"
    #     )

    def __str__(self):
        outstr = [
            f"Protein ID: {self.prot_id}",
            f"Protein sequence: {self.prot_seq}",
            f"NCBI ID: {self.locus_NCBI_id}",
            f"Uniprot ID: {self.uniprot_id}",
            f"Protein structure type: {self.struct_prot_type}",
        ]
        return "\n".join(outstr)


class Gene(Base):
    """Table representing a gene name and DNA sequence

    Each gene_ID represents a gene name and sequences. Several name strings
    given to an unique sequence
    """

    __tablename__ = "gene"
    proteins: Mapped[list["ProteinGene"]] = relationship()

    gene_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column
    gene_name = Column(String, nullable=False)
    dna_seq = Column(String, nullable=False)

    # Define relationships after defining columns
    # A one-to-many relationship between Protein and Gene
    # proteins = relationship(
    #     "Protein", secondary=proteingene, back_populates="genes", lazy="dynamic"
    # )


class Taxonomy(Base):
    """Table representing the taxon accession of a protein

    This table will store the taxon origin of the protein sequence, e.g.:
    specie, genus, family... and a accession number to a database with details
    about that organism
    """

    __tablename__ = "taxon"
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

    # A many-to-one relationship between Protein and Taxonomy
    proteins = relationship(
        "Protein", secondary=proteintaxon, back_populates="taxons", lazy="dynamic"
    )

    # To enforce unique taxon references
    __table_args__ = (UniqueConstraint("tax_id", "tax_ref"),)


class Pdb(Base):
    """Table representing the Pdb accession of a protein

    This table will store the different pdb accession number that represent the
    structure of a protein
    """

    __tablename__ = "pdb"
    pdb_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    pdb_acc_1 = Column(String, unique=True)  # primary accession number in pdb
    pdb_acc_2 = Column(String, nullable=False)  # accession number
    pdb_acc_3 = Column(String, nullable=False)  # accession number

    # A one-to-many relationship between Protein and Pdb structure
    proteins = relationship(
        "Protein", secondary=proteinpdb, back_populates="pdbs", lazy="dynamic"
    )


class Domain(Base):
    """Table representing the conserved domain family of a protein

    This table will store the different conserved domain accession number that represent the
    structure of a protein, including the reference database where the accession number
    was taken
    """

    __tablename__ = "domain"
    dom_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    dom_ref = Column(
        String, unique=True, nullable=False
    )  # domain accession in external db
    dom_db = Column(Integer, nullable=False)  # external database name e.g. pfam, CDD

    # A many-to-many relationship between Protein and domain family
    proteins = relationship(
        "Protein", secondary=proteindomain, back_populates="domains", lazy="dynamic"
    )

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
    go_id = Column(Integer, primary_key=True, autoincrement=True)  # primary key column
    go_ref = Column(String, unique=True, nullable=False)  # accession number in GO
    go_type = Column(String, nullable=False)  # GO type (MF,CC,BP)
    go_description = Column(String, nullable=False)  # text description of function

    # A many-to-many relationship between Protein and function
    proteins = relationship(
        "Protein", secondary=proteinGO, back_populates="functions", lazy="dynamic"
    )
    # To enforce unique function references
    __table_args__ = (UniqueConstraint("go_id", "go_ref"),)


class Enzymepath(Base):
    """Table representing the enzymatic reaction in which the protein
    participates

    This table will store the different Kegg Ontlogy accession numbers
    related wit an specific function of a protein. Thus, a protein can have
    more than one reference
    """

    __tablename__ = "enzymepath"
    path_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column
    KO_ref = Column(String, unique=True, nullable=False)  # accession number in KO

    # A many-to-many relationship between Protein and enzymatic activity
    proteins = relationship(
        "Protein", secondary=proteinpath, back_populates="paths", lazy="dynamic"
    )
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
raw_dir = Path(__file__).resolve().parent.parent / "data" / "raw" / "prot_info"
# Define path to the prot_data file
prot_data_file = raw_dir / "prot_data_repeated_uniprot_stop.csv"

mydata = []
with open(prot_data_file, newline="") as csvfile:
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
for (
    protseq,
    NCBIid,
    uniprot,
    struct,
    geneid,
    name,
    namerank,
    dnaseq,
    taxid,
    taxref,
    taxdb,
    spec,
    genu,
    fam,
    order,
    phyl,
    classt,
    stra,
    pdbid,
    pdb_1,
    pdb_2,
    pdb_3,
    pathid,
    KOid,
) in mydata:
    # Check what data is available:
    print(
        f"This is before adding session.query {protseq=}, {NCBIid=},{uniprot=}, {struct=}"
    )

    # 1. If the protein already exists, we store it in the `protein` variable
    # 2. If the protein does not exist, we create a new protein object and add it to
    #    the session
    # 3. We then check if the gene already exists, and if so we store it in the `gene`
    #    variable
    # 4. If the gene does not exist, we create a new gene object and add it to the
    #    session
    # 5. We then check if the gene is already associated with the protein, and if not
    #    we create a new `proteingene` object and add it to the session
    # 6. We then commit the session to the database
    ## NOTE: SQLAlchemy will autoflush the session when we query the database, so we
    ##       do not need to manually flush the session before committing
    ##       (https://docs.sqlalchemy.org/en/20/orm/session_basics.html#session-flushing)
    ##       This will also automatically commit the session if no exceptions are
    ##       raised but, if errors are raised, we will need to rollback the session
    ##       and continue to the next entry.

    try:
        print(f"Before query, {protseq[:10]=}..., {NCBIid=}, {uniprot=}, {struct=}")
        protein = (
            session.query(Protein)
            # .filter(Protein.prot_id == protid) So it is added automatically
            .filter(Protein.prot_seq == protseq)
            .filter(Protein.locus_NCBI_id == NCBIid)
            .filter(Protein.uniprot_id == uniprot)
            # .filter(Protein.struct_prot_type == struct) It can be the same to others (is just the type: Hexamer, pentamer...)
            .first()
        )
        print(f"After query, {protein=}")
        if not protein:
            protein = Protein(
                prot_seq=protseq,
                locus_NCBI_id=NCBIid,
                uniprot_id=uniprot,
                struct_prot_type=struct,
            )
            session.add(protein)
        else:
            print(f"Protein with prot id XXXX and NCBI_id {NCBIid} already exist")

        print(f"{protseq=}, {NCBIid=}, {uniprot=}, {struct=}")

        # Create a new gene object
        print(f"Before query, {geneid=}, {name=}, {dnaseq=}")
        gene = (
            session.query(Gene)
            # .filter(Gene.gene_id == geneid) automatically assigned
            .filter(Gene.gene_name == name)
            # .filter(Gene.dna_seq == dnaseq) Do not matter if the dna seq is the same (same prot can have different names)
            .first()
        )
        # Ensure gene is not None before proceeding
        if not gene:
            gene = Gene(gene_name=name, dna_seq=dnaseq)
            session.add(gene)
            print(f"Gene {name=} added")
        else:
            print(
                f"This gene name {name} has already being added to this gene ID {geneid}"
            )
        print(f"{protein.genes=}, {type(protein.genes)}")
        # print(f"{protein.genes.first()=}, {protein.genes.all()=}")

        genes_associated_with_protein = protein.genes

        if gene not in genes_associated_with_protein:
            print(f"{protein.prot_id=}, {gene.gene_id=}, {namerank=}")
            # protein_gene = proteingene.insert().values(
            #     gene_id=gene.gene_id,
            #     prot_id=protein.prot_id,
            #     name_rank=namerank
            # )
            proteingene = ProteinGene(name_rank=int(namerank))
            proteingene.gene = gene
            print(f"{proteingene=}")
            # # protein.genes.append(gene)  # Add the gene to the protein's genes collection
            # # rank = proteingene(name_rank=namerank)
            # session.add(protein_gene)
            # session.commit()  # Commit the changes
            protein.genes.append(proteingene)
            print(f"Linked Gene {gene.gene_id} to Protein {protein.prot_id}")
            print(f"{proteingene=}")
        else:
            print(f"Gene {gene.gene_id} is already linked to Protein {protein.prot_id}")
            print(geneid, name, dnaseq)

        # Try to commit our changes
        print("Committing changes")
        session.commit()
    except Exception as exc:
        print(f"Error committing protein/gene combination: {exc}")
        print("Rolling back changes and skipping to next entry")
        session.rollback()
        # sys.exit()

    continue  # Temporary skip of adding taxonomy data while we debug

    # Create new tax
    existing_tax = (
        session.query(Taxonomy)
        # .filter(Taxonomy.tax_id == taxid) ID automatically assigned
        .filter(Taxonomy.tax_ref == taxref)
        # .filter(Taxonomy.tax_db == taxdb) Same database for several references
        # .filter(Taxonomy.species == spec)
        # .filter(Taxonomy.genus == genu)
        # .filter(Taxonomy.family == fam)
        # .filter(Taxonomy.order_tax == order)
        # .filter(Taxonomy.phylum == phyl)
        # .filter(Taxonomy.class_tax == classt)
        .filter(Taxonomy.strain == stra)
        .first  # Strain should be unique no?
    )
    if not existing_tax:
        new_tax = Taxonomy(
            tax_ref=taxref,
            ltax_db=taxdb,
            species=spec,
            genus=genu,
            family=fam,
            order_tax=order,
            phylum=phyl,
            class_tax=classt,
            strain=stra,
        )
        session.add(new_tax)
        session.commit()
        if new_tax not in prot.taxon:
            prot.taxon.append(new_tax)
            session.commit()
            print(f"Added Taxon {taxid} to Protein {protid}")
        else:
            print(f"Taxon {taxid} already associated with Protein {protid}")
    else:
        print(f"This taxon has already being added {taxid, taxref}")

    print(taxid, taxref, taxdb, spec, genu, fam, order, phyl, classt, stra)

    # Create new pdb
    existing_pdb = (
        session.query(Pdb)
        # .filter(Pdb.pdb_id == pdbid) Added automatically
        .filter(Pdb.pdb_acc_1 == pdb_1)
        .filter(Pdb.pdb_acc_2 == pdb_2)
        .filter(Pdb.pdb_acc_3 == pdb_3)
        .first
    )
    if not existing_pdb:
        new_pdb = Pdb(pdb_acc_1=pdb_1, pdb_acc_2=pdb_2, pdb_acc_3=pdb_3)
        session.add(new_pdb)
        session.commit()
        if new_pdb not in prot.pdb:
            prot.pdb.append(new_pdb)
            session.commit()
            print(f"Added Pdb structure {pdbid} to Protein {protid}")
        else:
            print(f"Pdb structure {pdbid} already associated with Protein {protid}")

    else:
        print(f" This pdb entry already exist {pdbid}")

    print(pdbid, pdb_1, pdb_2, pdb_3)

    # Create new enzyme path
    existing_path = (
        session.query(Enzymepath)
        # .filter(Enzymepath.path_id == pathid) Added automatically
        .filter(Enzymepath.KO_ref == KOid)
        .first
    )
    if not existing_path:
        new_path = Enzymepath(KO_ref=KOid)
        session.add(new_path)
        session.commit()
        if new_path not in prot.path:
            prot.path.append(new_path)
            session.commit()
            print(f"Added EnzymePath {pathid} to Protein {protid}")
        else:
            print(f"EnzymePath {pathid} already associated with Protein {protid}")
    else:
        print(f"Existing gene ontology reference {KOid}")

    print(pathid, KOid)

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
