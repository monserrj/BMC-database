.mode column
.header on

CREATE TABLE Protein (
    prot_ID integer,
    prot_seq integer,
    locus_NCBI_ID varchar not null,
    uniprot_ID varchar not null,
    structural_prot_type integer, -- to define if it is an hexamer, pentamer, pseudohexamer (trimer) or not structural,
    primary key(prot_ID)
);

CREATE TABLE Gene (
    gene_name integer,
    gene_ID integer,
    dna_seq integer,-- If I have the DNA seq here I am duplicating info
    primary key(gene_ID)
);

CREATE TABLE Protein_gene (
    prot_ID integer,
    gene_ID integer,
    name_rank integer,
    -- primary gene gets rank 1; other names get 2, 3, ... Changed from bolean to be able to use UNIQUE
    -- to make a protein:gene relationship primary (promote it), then:
    -- 1. move primary relationship to name_rank of current largest number + 1
    -- 2. mover promoted gene to have name_rank 1
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(gene_ID) references Gene(gene_ID),
    UNIQUE(prot_id, gene_ID, name_rank)
    -- to make sure no protein has more than one main name
);

CREATE TABLE Taxon (
    tax_ID integer,
    tax_ref varchar not null,
    -- ref number in remote db
    tax_db varchar not null,
    -- stores whether NCBI or GTDB
    specie integer,
    genus integer,
    family integer,
    order_tax integer,
    phylum integer,
    class integer,
    strain integer,
    primary key(tax_ID)
);

CREATE TABLE Protein_taxon (
    prot_ID integer,
    tax_ID integer,
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(tax_ID) references Taxon(tax_ID)
);

CREATE TABLE Pdb (
-- if empty alphafold or other software will be used to predict structure
    pdb_ID integer,
    pdb_rscb_acc_1 varchar,
    pdb_rscb_acc_2 varchar,
    -- primary_pdb boolean. It does not matter which one is primary. Not needed
    primary key(pdb_ID)
);

CREATE TABLE Protein_pdb (
    pdb_ID integer,
    prot_ID integer,
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(pdb_ID) references Pdb(pdb_ID)
);

CREATE TABLE Domain (
    domain_ID integer,
    dom_ref_fam varchar not null,
    -- domain accession in remote db
    dom_ext_source integer not null,
    -- domain database, such as pfam, CDD...
    primary key(domain_ID)
);

CREATE TABLE Protein_domain (
    domain_ID integer,
    prot_ID integer,
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(domain_ID) references Domain(domain_ID)
);

CREATE TABLE Isoforms (
    canonical_prot_ID integer,
    isoform_prot_ID integer,
    -- Isoform ID is related to those proteins that may have to RBD
    -- 2 proteins from he same DNA seq, such as pduB and pduB'
    foreign key(canonical_prot_ID) references Protein(prot_ID),
    foreign key(isoform_prot_ID) references Protein(prot_ID)
);

CREATE TABLE Function_GO (
    go_ID integer, -- added as there was no ID assigned
    go_ext_ID varchar,
    -- ref number in GO db
    go_type varchar,
    -- MF (molecular function), CC (cellular compartment), or BP (biological process)
    go_description varchar,
    -- text description of the function
    primary key(go_ID)
);

CREATE TABLE Protein_GO (
    go_ID integer,
    prot_ID integer,
    foreign key(go_ID) references Function_GO(go_ID),
    foreign key(prot_ID) references Protein(prot_ID)
);

CREATE TABLE Enzyme_path (
    path_ID integer,
    KO varchar,
    -- ref number from Kegg Ontology db
    primary key(path_ID)
);

CREATE TABLE Protein_path (
    path_ID integer,
    prot_ID integer,
    foreign key(path_ID) references Enzyme_path(path_ID),
    foreign key(prot_ID) references Protein(prot_ID)
);

CREATE TABLE Complex (
    complex_ID integer,
    -- can be the native complex, e.g.: pdu operon, or a synthetic one, e.g.: pduABJMNK
    complex_type integer,
    -- need to decide which classification to use. Axen paper?
    complex_activity integer,
    -- Active/inactice
    assembly_exp_tested integer,
    -- Is the assembly of this compartment experimentally tested (Y/N)?
    -- Consider whether if Y answer paper references should be added. How?
    complex_source integer,
    -- Native complex, experimental complex, theorical complex
    primary key(complex_ID)
);-- Should I add classification for native complext?

CREATE TABLE Protein_complex(
    complex_ID integer,
    -- Several proteins can belong to the same complex and one protein to several complexes
    prot_ID integer,
    prot_essential_assembly integer, -- is this protein essential for the assembly of the complex (Y/N)?
    prot_interact_ID integer, -- ask Leighton, is about direct interation not complex
    copy_number integer, -- protein copy number in complex (stecheometry)
     -- structural_prot_type integer, -- move to protein table to avoid redundances
    foreign key(complex_ID) references Complex(complex_ID),
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(prot_interact_ID) references Protein(prot_ID)
);