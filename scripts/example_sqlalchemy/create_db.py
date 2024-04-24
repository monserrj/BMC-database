#!/usr/bin/env python
#
# This script is an example of using the SQLAlchemy ORM to create
# a database, following the schema in the BMC-database-skeleton.sql
# file - though not exactly. This will be a reduced schema (time
# constraints) and will have a different structure where I think
# that's either pedagogically useful, or where I think it's an
# approach you might want to consider.

# We need to import the SQLAlchemy classes we need.

# We'll be using a declarative approach to define our schema, so we
# import declarative_base and some helper classes.
from sqlalchemy.orm import declarative_base

# We'll need the create_engine function to create an engine object
from sqlalchemy import create_engine

# The columns in each table will inherit from SQLalchemy classes
# corresponding to the datatype of the column, so we need to import
# these classes.
from sqlalchemy import Column, Integer, String

# There will be at least one many-to-many relationship in our schema,
# so we need to import Table.
from sqlalchemy import Table

# We'll need to import ForeignKey to create foreign key relationships
from sqlalchemy import ForeignKey

# We'll need to import sessionmaker to create a session object
from sqlalchemy.orm import sessionmaker

# We'll need to import the relationship function to create the
# many-to-many and one-to-many relationship between Protein and
# the other tables
from sqlalchemy.orm import relationship

# We'll need to import the backref function to create the one-to-many
# relationship between Protein and ProteinAccession
from sqlalchemy.orm import backref

# We'll need to import UniqueConstraint to enforce a unique constraint
# on tables where we want to enforce uniqueness in a column
from sqlalchemy import UniqueConstraint

# Let's start to create the database.
# As we're using the declarative system, we need to create a base
# class to inherit from.
Base = declarative_base()

# We need to create a database engine to connect to the database.
# We'll use the create_engine function to create an engine object.
# This creates a new empty database file called bmc.db in the
# current directory.
engine = create_engine("sqlite:///bmc.db")


# Now we need to create the tables in the database.
# There is a limited amount of information in this schema, but we will have
# a one-to-many relationship between Protein and ProteinAccession, and
# a many-to-many relationship between Protein and ProteinStructure.
# We'll create tables to represent the one-to-many and many-to-many
# relationships before creating these tables, to satisfy the logic of the code.
proteinstructure = Table(  # linker table for many-to-many relationship
    "protein_structure",  # This name will be used in SQLite
    Base.metadata,
    Column("protein_id", Integer, ForeignKey("protein.protein_id")),
    Column("structure_id", Integer, ForeignKey("structure.structure_id")),
)

proteinaccession = Table(  # linker table for one-to-many relationship
    "protein_accession",
    Base.metadata,
    Column("protein_id", Integer, ForeignKey("protein.protein_id")),
    Column("accession_id", Integer, ForeignKey("accession.accession_id")),
)


# We create a class corresponding to each database table. No classes
# are instantiated until we create an instance of the Base class -
# this gets done after the tables are created.
# Each class inherits from the Base class we created earlier.
# I'm making this different from BMC-database-skeleton.sql
# Note that the one-to-many relationship between Protein and
# ProteinAccession is represented by the ForeignKey column in the
# Protein class (accession).
# The many-to-many relationship between Protein and ProteinStructure
# will be represented by the protein_structure table.
class Protein(Base):
    """Table representing a protein

    This table will store the protein sequence, links to accessions
    in other databases, and links to the type of protein structure.

    The protein sequence is taken to be unique, so we'll add a
    UniqueConstraint to the table to enforce this.
    """

    __tablename__ = "protein"  # this is the name that will be used in SQLite

    protein_id = Column(Integer, primary_key=True)  # primary key column
    seq = Column(String, nullable=False, unique=True)  # sequence string

    # We define relationships after defining columns
    # A one-to-many relationship between Protein and
    # ProteinAccession
    protein_accession = relationship("Accession", secondary=proteinaccession)
    # A many-to-many relationship between Protein and ProteinStructure
    # To define the relationship here, we needed to create a table
    # to represent the many-to-many relationship first
    protein_structure = relationship("Structure", secondary=proteinstructure)


class Accession(Base):
    """Table representing a protein accession

    Each accession is an accession ID/value, and the database it
    derives from
    """

    __tablename__ = "accession"
    # We need to enforce unique combinations of accession and database
    __table_args__ = (UniqueConstraint("accession", "database"),)

    accession_id = Column(Integer, primary_key=True)  # primary key column
    accession = Column(String, nullable=False)
    database = Column(String, nullable=False)


class Structure(Base):
    """Table representing a kind of protein structure

    This table will store the kind of protein structure, e.g.
    hexamer, pentamer, trimer, or not structural
    """

    __tablename__ = "structure"
    structure_id = Column(Integer, primary_key=True)  # primary key column
    structure = Column(String, unique=True, nullable=False)


# Now that we have defined the tables, we can create the tables in the
# database.
Base.metadata.create_all(engine)

# At this point we could check that the tables have been created by
# looking at the database file with sqlite3, or by using the
# SQLAlchemy inspector.

# Let's create some fake data to populate the database.
# Note that we've duplicated structures, sequences, and accessions
# We'll add the data to the database in a loop, but we'll have to
# check if the sequence, structure, and accession already exist
# and update the corresponding tables accordingly
mydata = [
    ("FCLEPPY", "accession1", "database1", "hexamer"),
    ("DEFGH", "accession2", "database2", "pentamer"),
    ("IKLMN", "accession3", "database2", "trimer"),
    ("PQRSTVWY", "accession4", "database1", "pentamer"),
    ("PQRSTVWY", "accession5", "database1", "trimer"),
    ("DEFGH", "accession4", "database2", "pentamer"),
]

# Start the session
Session = sessionmaker()  # we also need a session object
Session.configure(bind=engine)
session = Session()

# Add the data to the database
# We need to check if the sequence, structure, and accession already exist,
# and update the corresponding tables accordingly if they do not.
# We can then update the linker tables by adding the corresponding items.
for seq, acc, db, struct in mydata:
    # Create a new structure object
    new_struct = session.query(Structure).filter(Structure.structure == struct).first()
    if not isinstance(new_struct, Structure):
        new_struct = Structure(structure=struct)
        session.add(new_struct)
        session.commit()
    # Create a new accession object
    new_acc = (
        session.query(Accession)
        .filter(Accession.accession == acc)
        .filter(Accession.database == db)
        .first()
    )
    if not isinstance(new_acc, Accession):
        new_acc = Accession(accession=acc, database=db)
        session.add(new_acc)
        session.commit()
    # Create a new protein object
    new_protein = session.query(Protein).filter(Protein.seq == seq).first()
    if not isinstance(new_protein, Protein):
        new_protein = Protein(seq=seq)
        session.add(new_protein)
        session.commit()
    # Add the protein to the accession
    new_protein.protein_accession.append(new_acc)
    session.commit()
    # Add the protein to the structure
    new_protein.protein_structure.append(new_struct)
    session.commit()

# Now we can query the database to see if the data has been added correctly
# We'll query the Protein table and print out the sequence, the accessions,
# and the structures
for protein in session.query(Protein):
    print(f"\nPROTEIN: {protein.seq}")
    print("ACCESSIONS:")
    for acc in protein.protein_accession:
        print(f"\t{acc.accession}, {acc.database}")
    print("STRUCTURES:")
    for struct in protein.protein_structure:
        print(f"\t{struct.structure}")
