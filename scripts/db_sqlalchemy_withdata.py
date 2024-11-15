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
from sqlalchemy.orm import Mapped, mapped_column


class ProteinGene(Base):
    __tablename__ = "protein_gene"
    prot_id: Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), primary_key=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("gene.gene_id"), primary_key=True)
    name_rank: Mapped[Optional[int]]
    gene: Mapped["Gene"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="genes")

 
class ProteinTaxonomy (Base):
    __tablename__ = "protein_taxonomy"
    prot_id : Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), primary_key=True)
    tax_id : Mapped[int] = mapped_column(ForeignKey("taxonomy.tax_id"))
    taxonomy: Mapped["Taxonomy"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="taxonomies")
    UniqueConstraint("prot_id", "tax_id"),

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
    taxonomies: Mapped[list["ProteinTaxonomy"]] = relationship()

    prot_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column. Added autoincrement
    prot_seq = Column(String, nullable=False, unique=True)  # sequence string
    locus_NCBI_id = Column(String,unique=True, nullable=True)
    uniprot_id = Column(String, unique=True, nullable=True)
    struct_prot_type = Column(Integer, nullable=True)
    
    # Introduce back_populates so when a relationship between different tables is
    # introduced, they information will be backpopulated to be consistant accross
    # all tables. Relationships must be introduced in both related tables (e.g.:
    # Gene and Protein, with relationship based on table Proteingene)
    # genes = relationship(
    #     "Gene", secondary=proteingene, back_populates="proteins", lazy="dynamic"
    # )
    # taxons = relationship(
    #     "Taxonomy", secondary=proteintaxon, back_populates="proteins", lazy="dynamic"
    # )
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
    
class Taxonomy(Base):
    """Table representing the taxon accession of a protein
    This table will store the taxon origin of the protein sequence, e.g.:
    specie, genus, family... and a accession number to a database with details
    about that organism
    """
    
    __tablename__ = "taxonomy"
    proteins: Mapped[list["ProteinTaxonomy"]] = relationship()
    
    tax_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column
    tax_ref = Column(
        String, unique=True, nullable=False
    )  # accession number in db
    tax_db = Column(
        String, nullable=False
    )  # Name of database used (e.g.: NCBI, GTDB)
    species = Column(String, nullable=False)
    genus = Column(String)
    family = Column(String)
    order_tax = Column(String)
    phylum = Column(String)
    class_tax = Column(String)
    strain = Column(String)
    # To enforce unique taxonomy references
    __table_args__ = (UniqueConstraint("tax_id", "tax_ref"),)
    __table_args__ = (UniqueConstraint("species", "strain"),)


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

# Making a function for the addition of data to protein, gene and protein_gene tables:
def protein_gene_data_addition (session, protseq, NCBIid, uniprot, struct, name, namerank, dnaseq): 
    # Args:
    # protseq (str): Protein sequence.
    # NCBIid (str): NCBI locus ID.
    # uniprot (str): UniProt ID.
    # struct (str): Protein structure type.
    # name (str): Gene name.
    # dnaseq (str): Gene DNA sequence.
    # namerank (int): Rank of the gene name associated with the protein.
    
    # Explanation on how the code works:
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
    print(f"\nNow in {protein_gene_data_addition.__name__}")

    ## (2) LP: This could be refactored into a function, called from here
    try:
        print(f"Before query, {protseq[:10]=}..., {NCBIid=}, {uniprot=}, {struct=}")
        # Create a new protein object
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
        # Add protein if it is not already present
        if not protein:
            protein = Protein(
                prot_seq=protseq,
                locus_NCBI_id=NCBIid,
                uniprot_id=uniprot,
                struct_prot_type=struct,
            )
            session.add(protein)
            session.flush () # This sends the changes to the database, so prot_id is assigned
        else:
            print(f"Protein with prot id XXXX and NCBI_id {NCBIid} already exists")
            
        print(f"Protein row returned: {protein}")
        #print(f"{protseq[:10]=}, {NCBIid=}, {uniprot=}, {struct=}")
        ## (2)
        
        ## (3) LP: This could be refactored into a function, called from here
        # Create a new gene object
        print(f"Before query, {name=}, {dnaseq[:10]=}")
        gene = (
            session.query(Gene)
            # .filter(Gene.gene_id == geneid) automatically assigned
            .filter(Gene.gene_name == name)
            .filter(Gene.dna_seq == dnaseq)
            .all() # I want both conditions to be meet for it not to added, I don't care if the gene name is the same if the seuqence is different and viceversa
        )
        # Add gene if it is not already present
        if not gene:
            gene = Gene(gene_name=name, dna_seq=dnaseq)
            session.add(gene)
            session.flush ()
            print(f"Gene {name=} added")
        else:
            print(
                f"This gene name {name} has already being added"
            )
        print(f"Gene row returned: {gene}")
        ## (3)

        ## (4) LP: Could refactor this to a function, called from here
        # Associate the gene and protein information in the protein_gene table
        print(f"{protein.genes=}, {type(protein.genes)}")
        # print(f"{protein.genes.first()=}, {protein.genes.all()=}")
        
        gene_associated_with_protein = protein.genes  # LP: could avoid assigning this
        
        if gene not in gene_associated_with_protein:
            print(f"{protein.prot_id=}, {gene.gene_id=}, {namerank=}")
            proteingene = ProteinGene(name_rank=int(namerank))
            proteingene.gene = gene
            # print(f"{proteingene=}")
            print(f"{protein.genes=}")
            protein.genes.append(proteingene)
            print(f"{protein.genes=}")
            print(f"Linked Gene {gene.gene_id} to Protein {protein.prot_id}")
            print(f"{proteingene=}")
            
        else:
            print(f"Gene {gene.gene_id} is already linked to Protein {protein.prot_id}")
            print(f"{name}, {dnaseq [:10]=}")
        print(f"Linked gene from protein: {proteingene.gene}")
        ## (4)

        # Try to commit our changes
        print("Committing changes")
        session.commit()
        
    except Exception as exc:
        print(f"Error committing protein/gene combination: {exc}")
        print("Rolling back changes and skipping to next entry")
        session.rollback()
        # Rollback makes it that when there is a "fail", 
        # like not unique uniprot reference, its "forgets" the error and keeps going.

    return protein  # Return the protein row we just added to the db/otherwise dealt with

# Making a function for the addition of data to tax and prot_tax tables:
def taxonomy_data_addition (session, taxref, taxdb, spec, genu, fam, order, phyl, classt, stra, protein):
    # Args:
    # taxref (str): Taxonomy reference ID
    # taxdb (str): Taxonomy database used for reference ID
    # spec (str): Specie (taxonomy classification)
    # genu (str): Genus (taxonomy classification)
    # fam (str): Family (taxonomy classification)
    # order (str): Order (taxonomy classification)
    # phyl (str): Phylogeny (taxonomy classification)
    # classt (str): Class (taxonomy classification)
    # stra (str): Strain (taxonomy classification)
    
    # Explanation on how the code works:
    # 1. We check if the taxonomy already exists, and if so we store it in the `tax`
    #    variable
    # 2. If the taxonomy does not exist, we create a new gene object and add it to the
    #    session
    # 3. We then check if the taxonomy is already associated with the protein being added, and if not
    #    we create a new `proteintax` object and add it to the session
    # 4. We then commit the session to the database
    print(f"\nNow in {taxonomy_data_addition.__name__}")

    # LP: Had to turn off autoflushing to suppress an error here
    # See https://github.com/sqlalchemy/sqlalchemy/discussions/12049=
    with session.no_autoflush:

        try:
            ## (5) LP: Could refactor this to a function, called from here
            print(f"Before query, {taxref=}, {taxdb=}, {spec=}, {genu=}, {fam=}, {order=}, {phyl=}, {classt=}, {stra=}")
            # Create a new taxonomy object
            taxonomy = (
                session.query(Taxonomy)
                .filter(Taxonomy.tax_ref == taxref)
                .filter(Taxonomy.strain == stra) # Strain should be unique no?
                .first()
            )
            print(f"After query, {taxonomy=}")
            # Add taxonomy if it is not already present
            if not taxonomy:
                taxonomy = Taxonomy(
                    tax_ref=taxref,
                    tax_db=taxdb,
                    species=spec,
                    genus=genu,
                    family=fam,
                    order_tax=order,
                    phylum=phyl,
                    class_tax=classt,
                    strain=stra,
                )
                session.add(taxonomy)
                session.flush()
            else:
                print(f"This taxonomy {taxref} already exists")
                
            print(f"{taxref=}, {taxdb=}, {spec=}, {genu=}, {fam=}, {order=}, {phyl=}, {classt=}, {stra=}")
            print(f"Taxonomy row returned: {taxonomy}")
            ## (5)
            
            # print(f"{taxonomy.proteins=}")

            # Associate the taxonomy with the protein information in the protein_tax table
            print(f"{taxonomy.proteins=}, {type(taxonomy.proteins)}")

            taxonomy_associated_with_protein = taxonomy.proteins  # LP don't need to assign
            if taxonomy not in taxonomy_associated_with_protein:
                print(f"{protein.prot_id=}, {taxonomy.tax_id=}")
                proteintaxonomy = ProteinTaxonomy()
                print(f"{taxonomy.proteins=}")
                proteintaxonomy.taxonomy = taxonomy
                # proteintaxonomy.prot_id = protein
                print(f"{proteintaxonomy=}")
                print(f"{taxonomy.proteins=}")
                #taxonomy.proteins.append(proteintaxonomy)
                protein.taxonomies.append(proteintaxonomy)
                print(f"{taxonomy.proteins=}")
                print(f"Linked Taxonomy {taxonomy.tax_id} to Protein {protein.prot_id}")
                print(f"{proteintaxonomy=}")
            else:
                print(f"Taxonomy {taxonomy.tax_id} is already linked to Protein {protein.prot_id}")
                print(f"{taxref}, {spec}, {stra}")
            
            # Try to commit our changes;
            print("Committing changes")
            session.commit()
            
        except Exception as exc:
            print(f"Error committing protein/gene combination: {exc}")
            print("Rolling back changes and skipping to next entry")
            session.rollback()

# Script below here
if __name__ == "__main__":

    # Now that we have defined the tables, we can create the tables in the
    # database.
    Base.metadata.create_all(engine)

    # Add some data to populate the database.
    # Function for each data type addition created
    # Add the data to the database in a loop, but we'll have to
    # check if the data entered already exist
    # and update the corresponding tables accordingly.
    # We need to check if the sequence, structure, and accession already exist,
    # and update the corresponding tables accordingly if they do not.
    # We can then update the linker tables by adding the corresponding items.

    # Start the session
    Session = sessionmaker()  # we also need a session object
    Session.configure(bind=engine)
    session = Session()

    # Open the csv file
    # Define path for data file directory
    raw_dir = Path(__file__).resolve().parent.parent / "data" / "raw" / "prot_info"
    # Define path to the data file
    prot_data_file = raw_dir / "prot_data_minimal_correct.csv"

    ## (1) LP: this bit could be refactored into a function
    mydata = []
    with open(prot_data_file, newline="") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row
        for row in reader:
            # Convert empty strings to Null
            row = [None if val == "" else val for val in row]
            mydata.append(tuple(row))
    ## (1)

    # Add the data to the database
    for (
        protseq,
        NCBIid,
        uniprot,
        struct,
        name,
        namerank,
        dnaseq,
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
        print(f"This is before adding session.query {protseq=}, {NCBIid=},{uniprot=}, {struct=}")
        
        # Add protein-gene data
        protein = protein_gene_data_addition(
            session,
            protseq=protseq,
            NCBIid=NCBIid,
            uniprot=uniprot,
            struct=struct,
            name=name,
            dnaseq=dnaseq,
            namerank=namerank
        )
        print(f"Protein record returned: {protein}")
        
        # Add taxonomy data
        taxonomy_data_addition(
            session,
            taxref=taxref,
            taxdb = taxdb,
            spec = spec,
            genu = genu,
            fam = fam,
            order = order,
            phyl = phyl,
            classt = classt,
            stra = stra,
            protein = protein
        )


    # # Create new pdb
    # existing_pdb = (
    #     session.query(Pdb)
    #     # .filter(Pdb.pdb_id == pdbid) Added automatically
    #     .filter(Pdb.pdb_acc_1 == pdb_1)
    #     .filter(Pdb.pdb_acc_2 == pdb_2)
    #     .filter(Pdb.pdb_acc_3 == pdb_3)
    #     .first
    # )
    # if not existing_pdb:
    #     new_pdb = Pdb(pdb_acc_1=pdb_1, pdb_acc_2=pdb_2, pdb_acc_3=pdb_3)
    #     session.add(new_pdb)
    #     session.commit()
    #     if new_pdb not in prot.pdb:
    #         prot.pdb.append(new_pdb)
    #         session.commit()
    #         print(f"Added Pdb structure {pdbid} to Protein {protid}")
    #     else:
    #         print(f"Pdb structure {pdbid} already associated with Protein {protid}")

    # else:
    #     print(f" This pdb entry already exist {pdbid}")

    # print(pdbid, pdb_1, pdb_2, pdb_3)

    # # Create new enzyme path
    # existing_path = (
    #     session.query(Enzymepath)
    #     # .filter(Enzymepath.path_id == pathid) Added automatically
    #     .filter(Enzymepath.KO_ref == KOid)
    #     .first
    # )
    # if not existing_path:
    #     new_path = Enzymepath(KO_ref=KOid)
    #     session.add(new_path)
    #     session.commit()
    #     if new_path not in prot.path:
    #         prot.path.append(new_path)
    #         session.commit()
    #         print(f"Added EnzymePath {pathid} to Protein {protid}")
    #     else:
    #         print(f"EnzymePath {pathid} already associated with Protein {protid}")
    # else:
    #     print(f"Existing gene ontology reference {KOid}")

    # print(pathid, KOid)

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
