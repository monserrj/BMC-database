#!/usr/bin/env python

# This script is the development of the BMC db. Instructions
# followed and explanations needed kept to help following through
from pathlib import Path  # for type hints

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


# Function for protein data addition
def protein_addition(session, protseq, NCBIid, uniprot, struct, dnaseq):
    # Args:
    # protseq (str): Protein sequence.
    # NCBIid (str): NCBI locus ID.
    # uniprot (str): UniProt ID.
    # struct (str): Protein structure type.
    # dnaseq (str): Gene DNA sequence.

    # Explanation on how the code works:
    # 1. If the protein already exists, store it in the `protein` variable
    # 2. If the protein does not exist, create a new protein object and add it to
    #    the session
    # 3. Commit the session to the database

    ## NOTE: SQLAlchemy will autoflush the session when we query the database, so we
    ##       do not need to manually flush the session before committing
    ##       (https://docs.sqlalchemy.org/en/20/orm/session_basics.html#session-flushing)
    ##       This will also automatically commit the session if no exceptions are
    ##       raised but, if errors are raised, we will need to rollback the session
    ##       and continue to the next entry.

    print(f"\nNow in {protein_addition.__name__}")

    print(
        f"Before query, {protseq[:10]=}..., {NCBIid=}, {uniprot=}, {struct=}, {dnaseq[:10]=}"
    )
    # Create a new protein object
    protein = (
        session.query(Protein)
        # Prot_id added automatically
        .filter(Protein.prot_seq == protseq)
        .filter(Protein.locus_NCBI_id == NCBIid)
        .filter(Protein.uniprot_id == uniprot)
        # Protein struct is just the type (e.g.: hexamer), will be repeated
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
            dna_seq=dnaseq,
        )
        session.add(protein)
        session.flush()  # This sends the changes to the database, so prot_id is assigned
    else:
        print(f"Protein with prot id XXXX and NCBI_id {NCBIid} already exists")

    print(f"Protein row returned: {protein}")

    # Try to commit our changes
    # print("Committing changes")
    # session.commit()

    # except Exception as exc:
    #     print(f"Error committing protein/gene combination: {exc}")
    #     print("Rolling back changes and skipping to next entry\n")
    #     session.rollback()
    # Rollback makes it that when there is a "fail",
    # like not unique uniprot reference, its "forgets" the error and keeps going.
    return (
        protein  # Return the protein row we just added to the db/otherwise dealt with
    )


# Function for name data addition
def name_addition(session, genename, namerank, protein):
    # Args:
    # name (str): Gene name.
    # namerank (int): Rank of the gene name associated with the protein.
    # protein: Protein information added with protein_addition function

    # Explanation on how the code works:
    # 1. Check if the name already exists,  store it in the `name` variable
    # 2. If the name does not exist, create a new name object and add it to the
    #    session
    # 3. Check if the name is already associated with the protein, and if not
    #    create a new `proteinname` object and add it to the session
    # 4. Commit the session to the database

    print(f"\nNow in {name_addition.__name__}")

    with session.no_autoflush:
        # try:
        print(f"Before query, {genename=}")

        # Create a new name object
        name = (
            session.query(Name)
            # name_id automatically assigned
            .filter(Name.gene_name == genename)
            .first()
        )
        print(f"After query, {name=}")

        # Add name if it is not already present
        if not name:
            name = Name(gene_name=genename)
            session.add(name)
            session.flush()
            print(f"Name {genename=} added")
        else:
            print(f"This gene name {genename} has already being added")
        print(f"Name row returned: {name}")

        # Associate the gene name and protein information in the protein_gene table
        # Leighton suggested to make a specific function only for merged tables,
        # I think I want them linked because if I add a gene I want the information
        # instantly available for my protein and linked. need to double check with LP
        print(f"{name.proteins=}, {type(name.proteins)}")

        if name not in name.proteins:
            print(f"{protein.prot_id=}, {name.name_id=}, {namerank=}")
            proteinname = ProteinName(name_rank=int(namerank))
            # Name rank in there as is a new addition to the proteinname table
            # not in name or protein
            print(f"{name.proteins=}")
            proteinname.name = name
            print(f"{proteinname=}")
            print(f"{name.proteins=}")
            protein.names.append(proteinname)
            print(f"{name.proteins=}")
            print(f"\nLinked Gene name {name.name_id} to Protein {protein.prot_id}")
            print(f"{proteinname=}")
        else:
            print(
                f"Gene name {name.name_id} is already linked to Protein {protein.prot_id}"
            )
        print(f"{name}")
        # print(f"Linked gene from protein: {proteingene.gene}")

        # Try to commit our changes
        # print("Committing changes")
        # session.commit()

        # except Exception as exc:
        #     print(f"Error committing protein/gene combination: {exc}")
        #     print("Rolling back changes and skipping to next entry")
        #     session.rollback()
        # Rollback makes it that when there is a "fail",
        # like not unique uniprot reference, its "forgets" the error and keeps going.
        return name  # Return the gene row we just added to the db/otherwise dealt with


# Function for taxonomy data addition:
def taxonomy_addition(
    session, taxref, taxdb, spec, genu, fam, order, phyl, classt, stra, protein
):
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
    # protein: Protein information added with protein_addition function

    # Explanation on how the code works:
    # 1. Check if the taxonomy already exists, and if so store it in the `tax`
    #    variable
    # 2. If the taxonomy does not exist, create a new gene object and add it to the
    #    session
    # 3. Then check if the taxonomy is already associated with the protein being added, and if not
    #    create a new `proteintax` object and add it to the session
    # 4. Commit the session to the database

    print(f"\nNow in {taxonomy_addition.__name__}")

    # LP: Had to turn off autoflushing to suppress an error here
    # See https://github.com/sqlalchemy/sqlalchemy/discussions/12049=
    with session.no_autoflush:
        print(
            f"Before query, {taxref=}, {taxdb=}, {spec=}, {genu=}, {fam=}, {order=}, {phyl=}, {classt=}, {stra=}"
        )
        # Create a new taxonomy object
        taxonomy = session.query(Taxonomy).filter(Taxonomy.tax_ref == taxref).first()
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
            print(f"Taxonomy {taxref=} added")
        else:
            print(f"This taxonomy {taxref} already exists")

        print(f"Taxonomy row returned: {taxonomy}")

        # Associate the taxonomy with the protein information in the protein_tax table
        print(f"{taxonomy.proteins=}, {type(taxonomy.proteins)}")

        if taxonomy not in taxonomy.proteins:
            print(f"{protein.prot_id=}, {taxonomy.tax_id=}")
            proteintaxonomy = ProteinTaxonomy()
            print(f"{taxonomy.proteins=}")
            proteintaxonomy.taxonomy = taxonomy
            print(f"{proteintaxonomy=}")
            print(f"{taxonomy.proteins=}")
            protein.taxonomies.append(proteintaxonomy)
            print(f"{taxonomy.proteins=}")
            print(f"Linked Taxonomy {taxonomy.tax_id} to Protein {protein.prot_id}")
            print(f"{proteintaxonomy=}")
            session.flush()
        else:
            print(
                f"Taxonomy {taxonomy.tax_id} is already linked to Protein {protein.prot_id}"
            )
        print(f"{taxref}, {spec}, {stra}")

        #     # Try to commit our changes;
        #     print("Committing changes")
        #     session.commit()

        # except Exception as exc:
        #     print(f"Error committing protein/gene combination: {exc}")
        #     print("Rolling back changes and skipping to next entry")
        #     session.rollback()
        return taxonomy  # Return the taxonomy row we just added to the db/otherwise dealt with


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
    create_db()
