#!/usr/bin/env python

# This script links the addition of data manually to the BMC using csv from readfile to database described in db.
# Following the data addition structure from data_addition. Instructions
# followed and explanations kept to help following through

# To be able to make exceptions in code (try/except):
# from sqlalchemy.exc import IntegrityError, PendingRollbackError
# For exiting system for trouble shooting import sys
from db import protein_addition, taxonomy_addition, name_addition, gene_addition, function_addition


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
            
            gene = gene_addition(
                session,
                NCBIid=NCBIid,
                dnaseq=dnaseq,
                protein=protein,
            )
            print(f"\nProtein record returned: {protein}")
            print(f"\nGene record returned: {gene}")

            # Add name data
            name = name_addition(
                session, genename=genename, namerank=namerank, protein=protein
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
                    taxdb=taxdb,
                    spec=spec,
                    genu=genu,
                    fam=fam,
                    order=order,
                    phyl=phyl,
                    classt=classt,
                    stra=stra,
                    protein=protein,
                )
                print(f"\nProtein record returned: {protein}")
                print(f"Gene record returned: {gene}")
                print(f"Name record returned: {name}")
                print(f"Taxonomy record returned: {taxonomy}")
            session.commit()
            
            # Add function data
            function = function_addition(
                session, goref=goref, gotype=gotype, godescription=godescription, protein=protein
            )
            print(f"\nProtein record returned: {protein}")
            print(f"Function record returned: {function}")

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
