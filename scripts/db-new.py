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

# Database creation:

# Create the tables in the database.
# Tables with one-to-many and many-to-many relationships must be created
# before creating other tables, to satisfy the logic of the code.
class ProteinName(Base):
    __tablename__ = "protein_name"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), # primary_key=True
    )
    name_id: Mapped[int] = mapped_column(
        ForeignKey("name.name_id"),
    )  # Not primary key as can be repeated
    name: Mapped["Name"] = relationship(back_populates="proteins")
    protein: Mapped["Protein"] = relationship(back_populates="names")
    UniqueConstraint("prot_id", "name_id")

class Isoforms (Base):
    __tablename__ = "isoforms"
    canonical__prot_id : Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), primary_key=True)
    isoform__prot_id : Mapped[int] = mapped_column(ForeignKey("protein.prot_id"), primary_key=True)
    
    #Recheck how this works as right now is wrong
    isoform : Mapped["Protein"] = relationship(back_populates="isoforms")
    canonical: Mapped["Protein"] = relationship(back_populates="canonicals") # Recheck this

class Xref(Base): # All external references
    """Table representing the protein and CDS details from different db
    """
    
    __tablename__ = "xref"
    
    # Check if this can be done
    proteins: Mapped[list["ProteinXref"]] = relationship()
    cdss: Mapped[list["CdsXref"]] = relationship()
    
    # Define table content:
    xref_id = Column(Integer,primary_key=True, autoincrement=True)
    xref_acc = Column(String, unique=True, nullable=False) # Accession from external db
    xref_db_id = Column(Integer, nullable=False, ForeignKey="XrefIndex.xref_db_id") # External reference db key
    
    xref_db: Mapped["XrefIndex"] = relationship(back_populates="indexes")

class ProteinXref(Base):
    __tablename__ = "protein_xref"
    prot_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), # primary_key=True
    )
    xref_id: Mapped[int] = mapped_column(
        ForeignKey("xref.xref_id"),
    )
    
    protein: Mapped["Protein"] = relationship(back_populates="references")
    xref: Mapped["Xref"] = relationship(back_populates="proteins")
    UniqueConstraint("prot_id", "xref_id")


# Unsure whther class Cds goes here or below
class Cds(Base):
    """Table representing the gene details linked to a protein
    This table will store the gene ID, DNA sequence, 
    origin sequences if engineered, accession number and the corresponding
    protein
    """
    __tablename__ = "cds"
    
    references: Mapped[list["ProteinXref"]] = relationship()
    modifications: Mapped[list["CdsXref"]] = relationship()
    
    
    cds_id = Column(Integer, primary_key=True)
    cds_seq = Column(String, nullable=False, unique=True)
    # Figure out how this works:
    cds_accession = Column(String, nullable=False, Unique=True)
    
    prot_id = Column(Integer, ForeignKey("Protein.prot_id"))
    #prot_id: Mapped[int] = mapped_column(
    #    ForeignKey("protein.prot_id"), # primary_key=True
    #)
    cds_origin = Column(Integer, ForeignKey("Cds.cds_id")) # Foreign key, origin DNA sequence
    
    protein: Mapped["Protein"] = relationship(back_populates="cds")
    origin : Mapped["Cds"] = relationship(back_populates="cds", remote_side=[cds_id], post_update=True)
    # It says to add this as self-referential foreign key
    # https://docs.sqlalchemy.org/en/20/orm/self_referential.html

class CdsXref(Base):
    __tablename__ = "cds_xref"
    # cds_xref_id: Mapped[int] = mapped_column(
        
    # ) No need to assign and id for this i think
    cds_id: Mapped[int] = mapped_column(
        ForeignKey("cds.cds_id"), primary_key=True
    )
    xref_id: Mapped[int] = mapped_column(ForeignKey("xref.xref_id"), # primary_key=True
    )
    x_ref: Mapped["Xref"] = relationship(back_populates="cdss")
    cds: Mapped["Cds"] = relationship(back_populates="references")

class CdsModif(Base):
    __tablename__ = "cds_modif"
    cds_id: Mapped[int] = mapped_column(
        ForeignKey("protein.prot_id"), primary_key=True
    )
    modif_id: Mapped[int] = mapped_column(ForeignKey("modif.modif_id"), primary_key=True)
    
    modification: Mapped["Modif"] = relationship(back_populates="cdss")
    cds: Mapped["Protein"] = relationship(back_populates="modifications")




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
    """Table representing a protein This table will store the protein sequence, 
    unique accession number, type of protein structure and whether is
    canonical and isoform.
    """

    __tablename__ = "protein"  # this is the name that will be used in SQLite
    
    # Introduce all relationship between tables:
    names: Mapped[list["ProteinName"]] = relationship()
    cds: Mapped[list["Cds"]] = relationship()
    xref_prots: Mapped[list["ProteinXref"]] = relationship()
    
    # I am sure this is wrong, fix when TS
    isoforms: Mapped[list["Isoforms"]] = relationship()
    canonicals: Mapped[list["Isoforms"]] = relationship()
    
    references: Mapped[list["ProteinXref"]] = relationship()
    
    # Define table content:
    prot_id = Column(Integer, primary_key=True, autoincrement=True)
    prot_accession =  Column(Integer, pnullable=False, unique=True) # Unique accession number for proteins
    prot_seq = Column(String, nullable=False, unique=True)
    struct_prot_type = Column(Integer, nullable=True) # BMC type H/T/P or non structural
    is_canonical = Column(bool, default=True) # Only false if it is an isoform
    
    # Define type of output for protein
    def __str__(self):
        outstr = [
            f"Protein ID: {self.prot_id}",
            f"Protein accession: {self.prot_id}"
            f"Protein sequence: {self.prot_seq}",
            f"Protein structure type: {self.struct_prot_type}",
            f"Protein is canonical: {self.is_canonical}",
        ]
        return "\n".join(outstr)

class Name(Base):
    """Table representing a gene name
    Each gene_ID represents a protein name. Several name strings
    given to an unique protein, and several proteins sharing same name
    """

    __tablename__ = "name"
    # Introduce all relationship between tables:
    proteins: Mapped[list["ProteinName"]] = relationship()
    # Define table content:
    name_id = Column(
        Integer, primary_key=True, autoincrement=True
    )  # primary key column
    prot_name = Column(String, nullable=False)

class XrefIndex(Base):
    """Index table for all external databases"""
    
    __tablename__ = "xref_index"
    
    indexes: Mapped[list["Xref"]] =relationship()
    
    xref_db_id = Column(Integer, primary_key=True)
    xref_db = Column(String, nullable=False, unique=True) # Database name
    xref_href = Column(String, nullable=False, unique=True) # URL for database access
    xref_type = Column (String, nullable=False) # types of database: seq, struct, funct, tax'

class Modif(Base):
    """Table representing the modification of the engineered CDS sequences
    """
    
    __tablename__ = "modif"
    # Introduce all relationship between tables:
    cdss: Mapped[list["CdsModif"]] = relationship()
    
    # Define table content:
    modif_id = Column(Integer, primary_key=True, autoincrement=True)
    modif_type = Column(String, nullable=False) # vocabulary restricted: truncated, fusion, synthetic, mutated, domesticated
    modif_description = Column(Integer, nullable=False) # Description of the specific modification
    UniqueConstraint (modif_type,modif_description)

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
#     # A many-to-many relationship between Protein and enzymatic activity
#     proteins = relationship(
#         "Protein", secondary=proteincomplex, back_populates="interacts", lazy="dynamic")
# # To enforce unique no repeated protein to protein interactions are created
#     __table_args__ = (UniqueConstraint("prot_id_1", "prot_id_2", "prot_id_3", "prot_id_4", "prot_id_5", "prot_id_6", "prot_id_7"),)




# Function for protein data addition
def protein_addition(session, protaccession, protseq, struct, canonical):
    # Args:
    # protseq (str): Protein sequence.
    # protaccession (str): Protein unique accession number
    # struct (str): Protein structure type.
    # canonical (boolean): Canonical if true, false if isoform

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
        f"Before query, {protseq[:10]=}..., {protaccession=}, {struct=}, {canonical=}"
    )
    # Create a new protein object
    protein = (
        session.query(Protein)
        # Prot_id added automatically
        .filter(Protein.prot_seq == protseq)
        .filter(Protein.prot_accession == protaccession) # Do I need this?
        .first()
    )
    print(f"After query, {protein=}")

    # Add protein if it is not already present
    if not protein:
        protein = Protein(
            prot_seq=protseq,
            prot_accession=protaccession,
            struct_prot_type=struct,
            is_canonical=canonical,
        )
        session.add(protein)
        session.flush()  # This sends the changes to the database, so prot_id is assigned
    else:
        print(f"Protein with accession number {protaccession} already exists")

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
def name_addition(session, protname, protein):
    # Args:
    # protname (str): Protein name.
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
        print(f"Before query, {protname=}")

        # Create a new name object
        name = (
            session.query(Name)
            .filter(Name.prot_name == protname)
            .first()
        )
        print(f"After query, {name=}")

        # Add name if it is not already present
        if not name:
            name = Name(gene_name=protname)
            session.add(name)
            session.flush()
            print(f"Name {protname=} added")
        else:
            print(f"This gene name {protname} has already being added")
        print(f"Name row returned: {name}")

        print(f"{name.proteins=}, {type(name.proteins)}")

        if name not in name.proteins:
            print(f"{protein.prot_id=}, {name.name_id=}")
            proteinname = ProteinName()
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
        print(f"Linked gene from protein: {proteingene.gene}")

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

# Isoform_addition



# Function for external reference index of databases addition
def xrefindex_addition(session, xrefdb, href, xreftype):
    # Args:
    # xrefdb (str): Database name.
    # href (str): Database URL.
    # xreftype (str): Database type (sequence, structure, function, taxonomy)

    # Explanation on how the code works:
    # 1. Check if the database already exists,  store it in the `db` variable
    # 2. If the database does not exist, create a new db object and add it to the
    #    session
    # 3. Commit the session to the database

    print(f"\nNow in {xrefindex_addition.__name__}")

    with session.no_autoflush:
        print(f"Before query, {xrefdb=}, {href[10:]=},{xreftype=} ")

        # Create a new db object
        db = (
            session.query(XrefIndex)
            .filter(XrefIndex.xref_db == xrefdb)
            .filter(XrefIndex.xref_href == href)
            .first()
        )
        print(f"After query, {xrefdb=}, {href[10:]=},{xreftype=} ")

        # Add db if it is not already present
        if not db:
            db = XrefIndex(xref_db=xrefdb, xref_href=href, xref_type=xreftype)
            session.add(db)
            session.flush()
            print(f"Database renferece {xrefdb=} added")
        else:
            print(f"This database reference {xrefdb=} has already being added")
        print(f"Name row returned: {db}")

        return db  # Return the gene row we just added to the db/otherwise dealt with


def cds_addition(session, cdsseq, cdsaccession, cdsorigin, protein):
    # Args:
    # cdsseq (str): Gene DNA sequence.
    # cdsaccession: (str) Unique accession number CDS
    # cds_origin (str): Original CDS sequence if the sequence is modified.
    # N/A if native
    # protein: Protein information added with protein_addition function

    # Explanation on how the code works:
    # 1. Check if the cds already exists,  store it in the `cds` variable
    # 2. If the cds does not exist, create a new cds object and add it to the
    #    session
    # 3. Check if the cds is already associated with the protein, and if not
    #    create a new `proteingene` object and add it to the session
    # 4. Commit the session to the database

    print(f"\nNow in {cds_addition.__name__}")

    #with session.no_autoflush:
    print(f"Before query, {cdsseq[:10]=}...,{cdsaccession=} {cdsorigin=}")

    # Create a new cds object
    cds = (
        session.query(Cds)
        # cds_id and accession automatically assigned
        .filter(Cds.cds_seq == cdsseq)
        .first()
    )
    print(f"After query, {cds=}")

    # Add cds if it is not already present
    # I need to merge cds_origin to a protein ID while adding cds info
    # Unsure this is correct ( i got quite confuse come back to it)
    if not cds:
        print(f"Before adding seq and accession {cds=}")
        cds = Cds(
            cds_seq=cdsseq,
            cds_acccession=cdsaccession,
            protein=protein,
            cds_origin=cdsorigin)
        print(f" After adding seq and accession{cds=}")
        print(f"{protein.prot_id=}, {cds.cds_id=}")
        print(f"{protein.cdss=}")
        # I dont need to append as  it is not a many to many relationship apparently
        session.add(cds)
        session.flush()
        print(f"Cds {cdsseq[:10]=} with accession {cdsaccession=} added")
    else:
        print(f"This CDS {cdsseq[:10]=} with accession {cdsaccession=} has already being added")
    
    print(f"Cds row returned: {cds}")

    return cds  # Return the cds row we just added to the db/otherwise dealt with

def xref_addition(session, db, cds, protein, xrefacc):
    # Args:
    # Db: External database details added with xrefindex_addition function.
    # Cds: CDS information added with cds_addition function.
    # protein: Protein information added with protein_addition function
    # xrefacc (str): Unique accession number of that external reference
    
    # Explanation on how the code works:
    # 1. Check if the xref already exists,  store it in the `xref` variable
    # 2. If the xref does not exist, create a new xref object and add it to the
    #    session
    
    ''' Neeed to think on how to do that, condition depending on type of database?'''
    # 3. Check if the xref is already associated with the cds_xref, and if not
    #    create a new `cds_xref` object and add it to the session
    # 4. Check if the xref is already associated with the prot_xref, and if not
    #    create a new `prot_xref` object and add it to the session
    
    # 4. Commit the session to the database

    print(f"\nNow in {xref_addition.__name__}")

    print(f"Before query, {db=}, {xref_acc=}")

    # Create a new xref object
    xref = (
        session.query(Xref)
        # xref_gene_id automatically assigned
        .filter(Xref.xref_acc == xrefacc)
        .first()
    )
    print(f"After query, {xref=}")

    # Add xref if it is not already present
    if not xref:
        xref = Xref(xref_acc=xrefacc, XrefIndex=db)
        session.add(xref)
        session.flush()
        print(f"External reference {xrefacc=} from database {db=} added")
    else:
        print(f"This external reference {xrefacc=} from database {db=} has already being added")
    print(f"External reference row returned: External reference {xrefacc=}")

    # Associate the reference and protein information in the prot_xref
    print(f"{xref.proteins=}, {type(xref.proteins)}")

    if xref not in xref.proteins and xref_type
        """"Continue from here"""
    else:
        print(
            f"Gene locus reference {gref.gref_id} is already linked to Gene {cds.cds_id}"
        )
    print(f"{gene_xref}")

    return xref  # Return the renference row we just added to the db/otherwise dealt with


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
