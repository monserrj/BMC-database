#!/usr/bin/env python

# This script links the addition of data manually to the BMC using csv from readfile to database described in db.
# Following the data addition structure from data_addition. Instructions
# followed and explanations kept to help following through

# To be able to make exceptions in code (try/except):
# from sqlalchemy.exc import IntegrityError, PendingRollbackError
# For exiting system for trouble shooting import sys
from db import protein_addition, name_addition, cds_addition, gref_addition


# Making a function for the addition of data to the different tables created with db.py following data_instructions:
def link_db_csv(mydata, session):
    # Function for each data type addition created
    # Add the data to the database in a loop, but we'll have to
    # check if the data entered already exist
    # and update the corresponding tables accordingly.
    # We need to check if the sequence, structure, and accession already exist,
    # and update the corresponding tables accordingly if they do not.
    # We can then update the linker tables by adding the corresponding items.
    for idx, (
        protseq,
        accessionid,
        dbname,
        dburl,
        uniprot,
        struct,
        genename,
        cdstype,
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
        goref,
        gotype,
        godescription,
    ) in enumerate(mydata):
        # Check what data is available:
        print(f"\nStarting next loop ({idx=}) with {uniprot=}, {KOid=}")

        try:  # attempt to add to db, if we raise an error, we roll back
            # Add protein data
            protein = protein_addition(
                session,
                protseq=protseq,
                uniprot=uniprot,
                struct=struct,
            )
            print(f"\nProtein record returned: {protein}")
            
            gene = cds_addition(
                session,
                cdstype=cdstype,
                dnaseq=dnaseq,
                protein=protein,
            )
            print(f"\nProtein record returned: {protein}")
            print(f"\nGene record returned: {gene}")

            # Add name data
            name = name_addition(
                session, genename=genename, protein=protein
            )

            print(f"\nProtein record returned: {protein}")
            print(f"Name record returned: {name}")
            
            # Add locus reference data
            gref = gref_addition(
                session, accessionid=accessionid, dbname=dbname, dburl=dburl, protein=protein
            )

            print(f"\nProtein record returned: {protein}")
            print(f"Name record returned: {gref}")

            # # Add taxonomy data
            # # We expect a single tax_id per protein, so if the protein
            # # already is linked to a taxonomy, we skip adding the taxonomy
            # if not len(protein.taxonomies):
            #     taxonomy = taxonomy_addition(
            #         session,
            #         taxref=taxref,
            #         taxdb=taxdb,
            #         spec=spec,
            #         genu=genu,
            #         fam=fam,
            #         order=order,
            #         phyl=phyl,
            #         classt=classt,
            #         stra=stra,
            #         protein=protein,
            #     )
            #     print(f"\nProtein record returned: {protein}")
            #     print(f"Gene record returned: {gene}")
            #     print(f"Name record returned: {name}")
            #     print(f"Taxonomy record returned: {taxonomy}")
            # session.commit()
            
            # # Add function data
            # function = function_addition(
            #     session, goref=goref, gotype=gotype, godescription=godescription, protein=protein
            # )
            # print(f"\nProtein record returned: {protein}")
            # print(f"Function record returned: {function}")

            # # Add pdb data
            # pdb = pdb_addition(
            #     session, pdb_1=pdb_1, pdb_2=pdb_2, pdb_3=pdb_3, protein=protein
            # )
            # print(f"\nProtein record returned: {protein}")
            # print(f"Function record returned: {pdb}")
        except Exception as exc:
            print(f"Error adding data: {exc}")
            print("Rolling back changes and skipping to next entry")
            session.rollback()