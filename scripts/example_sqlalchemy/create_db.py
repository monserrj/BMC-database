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
Session = sessionmaker()  # we also need a session object

# We need to create a database engine to connect to the database.
# We'll use the create_engine function to create an engine object.
# This creates a new empty database file called bmc.db in the
# current directory.
engine = create_engine("sqlite:///bmc.db")


# Now we need to create the tables in the database.
# There is a limited amount of information in this schema, but we will have
# a one-to-many relationship between Protein and ProteinAccession, and
# a many-to-many relationship between Protein and ProteinStructure.
# We'll need to create a table to represent the many-to-many relationship
# before creating these tables, to satisfy the logic of the code.
protein_structure = Table(
    "protein_structure",
    Base.metadata,
    Column("protein_id", Integer, ForeignKey("protein.protein_id")),
    Column("structure_id", Integer, ForeignKey("structure.structure_id")),
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

    __tablename__ = "protein"
    __table_args__ = (UniqueConstraint("seq"),)

    protein_id = Column(Integer, primary_key=True)  # primary key column
    seq = Column(String)
    accession = Column(Integer, ForeignKey("protein_accession.accession_id"))
    type = Column(Integer)

    # We define relationships after defining columns
    # A one-to-many relationship between Protein and
    # ProteinAccession
    protein_accession = relationship(
        "ProteinAccession", backref=backref("protein", uselist=False)
    )
    # A many-to-many relationship between Protein and ProteinStructure
    # To define the relationship here, we needed to create a table
    # to represent the many-to-many relationship first
    protein_structure = relationship(
        "ProteinStructure", secondary=protein_structure, backref="protein"
    )


class ProteinAccession(Base):
    """Table representing a protein accession

    Each accession is an accession ID/value, and the database it
    derives from
    """

    __tablename__ = "protein_accession"
    accession_id = Column(Integer, primary_key=True)  # primary key column
    accession = Column(String)
    database = Column(String)


class ProteinStructure(Base):
    """Table representing a kind of protein structure

    This table will store the kind of protein structure, e.g.
    hexamer, pentamer, trimer, or not structural
    """

    __tablename__ = "structure"
    structure_id = Column(Integer, primary_key=True)  # primary key column
    structure = Column(String)


# Now that we have defined the tables, we can create the tables in the
# database.
Base.metadata.create_all(engine)

# At this point we could check that the tables have been created by
# looking at the database file with sqlite3, or by using the
# SQLAlchemy inspector.
