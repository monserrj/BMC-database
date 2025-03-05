#!/usr/bin/env python

# This script is the addition of data manually to the BMC db. Instructions
# followed and explanations kept to help following through

# Open .csv to be able to work with .csv files
# import csv

# Import path to be able to open .csv files in different folders
from pathlib import Path

from typing import Optional

# Import  SQLAlchemy classes needed with a declarative approach.
from sqlalchemy.orm import declarative_base

# from sqlalchemy import (
#     Column,
#     Integer,
#     String,
#     Table,
#     ForeignKey,
#     UniqueConstraint,
#     String,
# )
# from sqlalchemy.orm import relationship

# Create_engine function to create an engine object
from sqlalchemy import create_engine

# Import a sessionmaker to create a session object
from sqlalchemy.orm import sessionmaker

# To be able to make exceptions in code (try/except):
# from sqlalchemy.exc import IntegrityError, PendingRollbackError
# For exiting system for trouble shooting import sys
import sys

# Create a base class to inherit from.
Base = declarative_base()

# Temporary code to delete local database while we debug
# Path("bmc.db").unlink()

# Create a database engine to connect to the database.
# This creates a new empty database file called bmc.db in the current directory.
# engine = create_engine("sqlite:///bmc.db")

# Create the tables in the database.
# Tables with one-to-many and many-to-many relationships must be created
# before creating other tables, to satisfy the logic of the code.
# from sqlalchemy.orm import Mapped, mapped_column

# Making a function for the addition of data to the different tables created with db.py:

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
    
    print(f"Before query, {protseq[:10]=}..., {NCBIid=}, {uniprot=}, {struct=}, {dnaseq[:10]=}")
        
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
        session.flush () # This sends the changes to the database, so prot_id is assigned
    else:
        print(f"Protein with prot id XXXX and NCBI_id {NCBIid} already exists")
            
    print(f"Protein row returned: {protein}")
        
    # Try to commit our changes
    print("Committing changes")
    session.commit()
        
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
        #try:
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
            session.flush ()
            print(f"Name {genename=} added")
        else:
            print(
                    f"This gene name {genename} has already being added"
                )
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
        ## (5) LP: Could refactor this to a function, called from here
        print(
            f"Before query, {taxref=}, {taxdb=}, {spec=}, {genu=}, {fam=}, {order=}, {phyl=}, {classt=}, {stra=}"
        )
        # Create a new taxonomy object
        taxonomy = (
            session.query(Taxonomy)
            .filter(Taxonomy.tax_ref == taxref)
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
            print(f"Taxonomy {taxref=} added")
        else:
            print(f"This taxonomy {taxref} already exists")
                
        print(f"Taxonomy row returned: {taxonomy}")
            
        ## (5)
            
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
        return taxonomy # Return the taxonomy row we just added to the db/otherwise dealt with



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
    raw_dir = (
        Path(__file__).resolve().parent.parent 
        / "data" 
        / "raw" 
        / "prot_info" 
        / "incorrect_name"
    )
    # Define path to the data file
    prot_data_file = raw_dir / "prot_data_2_genename.csv"

    ## (1) LP: this bit could be refactored into a function
    mydata = []
    with open(prot_data_file, newline="") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row
        for row in reader:
            # Convert empty strings to Null
            row = [None if val == "" else val for val in row]
            if (
                row[0] is not None
            ):  # LP: we could do with more sanity checking of data here
                mydata.append(tuple(row))
            else:
                print("No sequence data, discarding")
    ## (1)

    # Add the data to the database
    for idx, (
        protseq,
        NCBIid,
        uniprot,
        struct,
        genename,
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
    ) in enumerate(mydata):
        
        # Check what data is available:
        # print(f"This is before adding session.query {protseq=}, {NCBIid=},{uniprot=}, {struct=}")
        print(f"\nStarting next loop ({idx=}) with {uniprot=}, {KOid=}")

        try:  # attempt to add to db, if we raise an error, we roll back
            # Add protein data
            protein = protein_addition(
                session,
                protseq=protseq,
                NCBIid=NCBIid,
                uniprot=uniprot,
                struct=struct,
                dnaseq=dnaseq,
            )
            print(f"\nProtein record returned: {protein}")
        
            # Add name data
            name = name_addition(
                session,
                genename=genename,
                namerank=namerank,
                protein = protein
            )
        
            print(f"\nProtein record returned: {protein}")
            print(f"Name record returned: {name}")
        
            # Add taxonomy data
            # We expect a single tax_id per protein, so if the protein
            # already is linked to a taxonomy, we skip adding the taxonomy
            if not len(protein.taxonomies):
                taxonomy = taxonomy_addition(
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
                print(f"\nProtein record returned: {protein}")
                print(f"Name record returned: {name}")
                print(f"Taxonomy record returned: {taxonomy}")
            session.commit()

        except Exception as exc:
            print(f"Error adding data: {exc}")
            print("Rolling back changes and skipping to next entry")
            session.rollback()

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
