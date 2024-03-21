.mode column.header on CREATE TABLE Protein (
    prot_ID integer,
    prot_seq integer,
    locus_NCBI_ID varchar not null,
    uniprot_ID varchar not null,
    primary key(prot_ID)
);

CREATE TABLE Gene (
    gene_name integer,
    gene_ID integer,
    dna_seq integer,
    primary key(gene_ID)
);

-- to make a protein:gene relationship primary (promote it), then:
-- 1. move primary relationship to name_rank of current largest number + 1
-- 2. mover promoted gene to have name_rank 1
CREATE TABLE Protein_gene (
    prot_ID integer,
    gene_ID integer,
    --primary_name boolean,
    name_rank integer,
    -- primary gene gets rank 1; other names get 2, 3, ...
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(gene_ID) references Gene(gene_ID),
    UNIQUE(prot_id, gene_ID, name_rank)
);

CREATE TABLE Taxon (
    tax_ID integer,
    tax_ref varchar not null,
    -- could be NCBI or GTDB
    taxdb varchar not null,
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
    pdb_ID integer,
    pdb_rscb_acc_1 varchar,
    pdb_rscb_acc_2 varchar,
    -- primary_pdb boolean,
    primary key(pdb_ID)
);

CREATE TABLE Protein_pdb (
    pdb_ID integer,
    prot_ID integer,
    -- primary_pdb boolean,
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(pdb_ID) references Pdb(pdb_ID)
);

CREATE TABLE Domain (
    domain_ID integer,
    -- domain accession in remote db
    dom_ref_fam varchar not null,
    ext_source integer not null,
    -- doman database
    primary key(domain_ID)
);

CREATE TABLE Protein_domain (
    domain_ID integer,
    prot_ID integer,
    foreign key(prot_ID) references Protein(prot_ID),
    foreign key(domain_ID) references Domain(domain_ID)
);

CREATE TABLE ISOFORMS (
    canonical_prot_ID integer,
    isoform_prot_ID integer,
    foreign key(canonical_prot_ID) references Protein(prot_ID),
    foreign key(isoform_prot_ID) references Protein(prot_ID)
);

/*
 CREATE TABLE Function_MF (
 MF_ID integer,
 GO_MF_1 varchar,
 description_MF_1 integer,
 GO_MF_2 varchar,
 description_MF_2 integer,
 primary key(MF_ID)
 );
 
 CREATE TABLE Function_CC (
 CC_ID integer,
 GO_CC_1 varchar,
 description_CC_1 integer,
 GO_CC_2 varchar,
 description_CC_2 integer,
 primary key(CC_ID)
 );
 
 CREATE TABLE Function_BP (
 BP_ID integer,
 GO_BP_1 varchar,
 description_BP_1 integer,
 GO_BP_2 varchar,
 description_BP_2 integer,
 primary key(BP_ID)
 );
 */
CREATE TABLE Function_GO (
    goID integer,
    go_type varchar,
    -- MF, CC, or BP
    go_description varchar,
    primary key(goID)
);

CREATE TABLE Protein_GO (
    goID integer,
    protein_ID integer,
    foreign key(goID) references Function_GO(goID),
    foreign key(protein_ID) references Protein(prot_ID),
);

CREATE TABLE Enzyme_path (
    path_ID integer,
    KO varchar,
    primary key(path_ID)
);

CREATE TABLE Protein_Path (
    path_ID integer,
    protein_ID integer,
    foreign key(path_ID) references Enzyme_path(path_ID),
    foreign key(protein_ID) references Protein(prot_ID),
)
/*
 CREATE TABLE Protein_function (
 BP_ID integer,
 CC_ID integer,
 MF_ID integer,
 path_ID integer,
 prot_ID integer,
 foreign key(BP_ID) references Function_BP(BP_ID),
 foreign key(CC_ID) references Function_CC(CC_ID),
 foreign key(MF_ID) references Function_MF(MF_ID),
 foreign key(path_ID) references Enzyme_path(path_ID),
 foreign key(prot_ID) references Protein(prot_ID)
 );
 */
CREATE TABLE Complex (
    complex_ID integer,
    complex_type integer,
    complex_activity,
    assembly_exp_tested integer,
    primary key(complex_ID)
);

CREATE TABLE Protein_complex(
    complex_ID integer,
    prot_ID integer,
    prot_essential_assembly integer,
    -- prot_interact_ID integer,
    copy_number integer,
    structural_prot_type integer,
    foreign key(complex_ID) references Complex(complex_ID),
    foreign key(prot_ID) references Protein(prot_ID),
    -- foreign key(prot_interact_ID) references Protein(prot_ID)
);