#!/usr/bin/env python

# This script describes how the addition of data for each table is done manually to the BMC db. Instructions
# followed and explanations kept to help following through

# Import path to be able to open .csv files in different folders
from pathlib import Path
# Import  SQLAlchemy classes needed with a declarative approach.
from sqlalchemy import create_engine

# Import a sessionmaker to create a session object
from sqlalchemy.orm import sessionmaker

# To be able to make exceptions in code (try/except):
# from sqlalchemy.exc import IntegrityError, PendingRollbackError
# For exiting system for trouble shooting import sys
# import sys

from db import db_URL
# Create a base class to inherit from.
engine = create_engine(db_URL)
SessionLocal = sessionmaker(bind=engine)

# Making a function to explain how to add data to the different tables created with db.py:
def add_data():
    session = SessionLocal()
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