# README.md `raw`

This folder would hold the original raw data as it was collected, such as:

- sequence read files (e.g. `fastq`)
  Different categories of data, or data collected on different dates or for different samples might be stored in suitably-named subfolders, e.g.by table type in the BMC db
  
- .csv with manually added data for the db to read in sql-alchemy. Csv files are inside several folder.
  Nomenclature:
    - start file is the type of data introduced
    - minimal for when minimal information for correct readability is added
    - complete for when all the information compiled from literature search and added to database_start.xlsx file is present
    - repeated_x: to state what type of information is repeated in the data (e.g.: full mean a whole line of info)
    - stop: indicates that this should be detected by the code and stop the repetition.
      
  Content:
  - domain_data:
    - domain_data_complete.csv : complete domain information available for each protein id with the corresponding prot_id reference.
    - domain_data_minimal.csv : minimal domain information available for each protein id with the corresponding prot_id reference (no repeats/incorrect info).
    - domain_data_repeated_domid : minimal domain information available for each protein id with the corresponding prot_id reference with repeated domain id.
    - domain_data_repeated_domref : minimal domain information available for each protein id with the corresponding prot_id reference with repeated domain reference.
    - domain_data_repeated_full : minimal domain information available for each protein id with the corresponding prot_id reference with repeated full doamin information.
    - domain_data_repeated_protid : minimal domain information available for each protein id with the corresponding prot_id reference with repeated protein id.
  - go_data:
    - go_data_complete.csv : complete gene ontology information available for each protein id with the corresponding prot_id reference.
    - go_data_minimal.csv : minimal gene ontology information available for each protein id with the corresponding prot_id reference (no repeats/incorrect info).
    - go_data_repeated_full.csv : minimal gene ontology information available for each protein id with the corresponding prot_id reference with repeated go entry.
    - go_data_repeated_godescription.csv : minimal gene ontology information available for each protein id with the corresponding prot_id reference with repeated go description for one go entry.
    - go_data_repeated_goid.csv : minimal gene ontology information available for each protein id with the corresponding prot_id reference with repeated go id for one go entry.
    - go_data_repeated_goref.csv : minimal gene ontology information available for each protein id with the corresponding prot_id reference with repeated go reference for one go entry.
    - go_data_repeated_protid.csv : minimal gene ontology information available for each protein id with the corresponding prot_id reference with repeated prot id for a different go entry.
  - isoform_data:
    - isoform_data_complete.csv : complete list of isoform and canonical proteins (by protein id) for all the proteins entered.
    - isoform_data_long.csv : extra long list of isoform and canonical proteins (by protein id) for different proteins entered (as existing list only includes one).
    - isoform_data_repeated_canonical.csv : complete list of isoform and canonical proteins (by protein id) with a canonical protein id repeated.
    - isoform_data_repeated_full.csv : complete list of isoform and canonical proteins (by protein id) with a full entry repeated.
    - isoform_data_repeated_isoform.csv : complete list of isoform and canonical proteins (by protein id) with a isoform protein id repeated.
    - isoform_data_samecan&iso.csv : complete list of isoform and canonical proteins (by protein id) with a the sampe protein id entered as canonical and isoform.
  - prot_info folder:
    - incorrect_gene folder: 
        - prot_data_gene_repeated_full_stop.csv : with the minimal amount of prot sequences and gene information added, with one gene details repeated.
        - prot_data_gene_repeated_geneid.csv : with the minimal amount of prot sequences and gene information added with repeated gene id for one gene reference.
        - prot_data_gene_repeated_genename_stop.csv :  with the minimal amount of prot sequences and gene information added with repeated gene name for one gene reference.
        - prot_data_gene_repeated_geneseq.csv :  with the minimal amount of prot sequences and gene information added with repeated gene sequence for one gene reference.
    - incorrect_path folder:
        - prot_data_path_repeated_full_stop.csv : with the minimal amount of prot sequences and enzymatic path information added, with one path details repeated.
        - prot_data_path_repeated_ko_stop.csv :  with the minimal amount of prot sequences and enzymatic path information added with repeated kegg orthology reference for one path.
        - prot_data_path_repeated_path.csv :  with the minimal amount of prot sequences and enzymatic path information added with repeated path id for one path.        - 
    - incorrect_pdb folder:
        - prot_data_pdb_repeated_pdb1_stop.csv : with the minimal amount of prot sequences and pdb structural information added, with the first pdb reference details repeated in another pdb id.
        - prot_data_pdb_repeated_pdb1cross2.csv : with the minimal amount of prot sequences and pdb structural information added with repeated pdb ref in 2 different pdb_acc_1 and 2 respectively for two pdb ids.
        - prot_data_pdb_repeated_pdb1cross3.csv :  with the minimal amount of prot sequences and pdb structural information added with repeated pdb ref in 2 different pdb_acc_1 and 3 respectively for two pdb ids.
        - prot_data_pdb_repeated_pdb1to2.csv :  with the minimal amount of prot sequences and pdb structural information added with repeated pdb ref in 2 different pdb_acc_1 and 2 in the same pdb_id.
        - prot_data_pdb_repeated_pdb1to3.csv :  with the minimal amount of prot sequences and pdb structural information added with repeated pdb ref in 2 different pdb_acc_1 and 3 in the same pdb_id.
        - prot_data_pdb_repeated_pdb2_stop.csv :  with the minimal amount of prot sequences and pdb structural information added, with the second pdb reference details repeated in another protein.
        - prot_data_pdb_repeated_pdb2cross3.csv : with the minimal amount of prot sequences and pdb structural information added with repeated pdb ref in 2 different pdb_acc_2 and 3 respectively for two pdb ids.
        - prot_data_pdb_repeated_pdb2to3.csv : with the minimal amount of prot sequences and pdb structural information added with repeated pdb ref in 2 different pdb_acc_2 and 3 in the same pdb_id.
        - prot_data_pdb_repeated_pdb3_stop.csv : with the minimal amount of prot sequences and pdb structural information added, with the third pdb reference details repeated in another protein.
        - prot_data_pdb_repeated_full_stop.csv : with the minimal amount of prot sequences and pdb structural information added with repeated pdb info twice.
        - prot_data_pdb_repeated_pdbid.csv : with the minimal amount of prot sequences and pdb structural information added with repeated pdb id in for different pdb details.
    - incorrect_tax folder: taxonomy database used is always the same and repeated as it should not affect it. Family, order, class and phylum are also repeated as they are common for some species
        - prot_data_tax_repeated_full_stop.csv : with the minimal amount of prot sequences and taxonomy information added with repeated taxonomy information (all details for 2 tax reference).
        - prot_data_tax_repeated_genus.csv : with the minimal amount of prot sequences and taxonomy information added with repeated genus for one taxonomy reference.
        - prot_data_tax_repeated_species.csv :  with the minimal amount of prot sequences and taxonomy information added with repeated specie for one taxonomy reference.
        - prot_data_tax_repeated_strain_stop.csv :  with the minimal amount of prot sequences and taxonomy information added with repeated strain for one taxonomy reference.
        - prot_data_tax_repeated_taxid.csv :  with the minimal amount of prot sequences and taxonomy information added with repeated tax id for one taxonomy reference.
        - prot_data_tax_repeated_taxref_stop.csv :  with the minimal amount of prot sequences and taxonomy information added with repeated taxonomy reference for one taxonomy reference.
    - prot_data_complete.csv : with all prtoein info from original excel file database_start.xlsx
    - prot_data_minimal_correct.csv : with the minimal amount of prot information for protein info, taxonomy, path, pdb and gene added with no error introduced for first testing.
    - prot_data_repeated_full_stop.csv : with the minimal amount of prot sequences with a whole repeated protein added.
    - prot_data_repeated_ncbi_stop.csv : with the minimal amount of prot sequences with a repeated ncbi reference added.
    - prot_data_repeated_protid.csv : with the minimal amount of prot sequences with a repeated protein id added.
    - prot_data_repeated_protseq_stop.csv : with the minimal amount of prot sequences with a repeated protein sequence added.
    - prot_data_repeated_prottype.csv : with the minimal amount of prot sequences with a repeated protein type (e.g.: BMC-T,non structural,etc) added.
    - prot_data_repeated_uniprot_stop.csv : with the minimal amount of prot sequences with a repeated uniprot reference added.
- complex_data.csv: Pending to be modofied file with complex information once complex tables from database are restructure.
- database_start.xlsx : Excel file with original information researched to start designing the database and its different tables.
- prot_complex_data.csv: Pending to be modofied file with complex information once complex tables from database are restructure.
