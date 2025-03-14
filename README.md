# R PACKAGE: veriNA3d

Paper: https://academic.oup.com/bioinformatics/article/35/24/5334/5529936
Citation details below.

This is the new repository where I will make updates on the veriNA3d R package (**stable** version), if necessary.
It can still be found at the original site: http://mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d, but the mainteinance there is discontinued.

## Description

VeriNA3d is an R package for the analysis of Nucleic Acid structural data. The software was developed on top of bio3d (Grant et al, 2006) with a higher level of abstraction. In addition of single-structure analyses, veriNA3d also implements pipelines to handle whole datasets of mmCIF/PDB structures. As far as we know, no similar software has been previously distributed, thus it aims to fill a gap in the data mining pipelines of PDB structural data analyses.

## Installation
---------------

Instructions

1- Make sure you have all the dependencies already installed in R. If not the case, open R and run:
&nbsp;

    install.packages(c("bio3d", "circlize", "jsonlite", "plot3D", "MASS", "RColorBrewer", "RANN"))

2- Install veriNA3d according with your R version:

&nbsp;

    install.packages("https://github.com/Digp/veriNA3d-stable/raw/master/veriNA3d.tar.gz", repos = NULL)

3- To start using it, just load the package!
&nbsp;

    library(veriNA3d)

## Documentation
----------------

### Dataset level

`getRNAList`: Get list of representative/non-redundant RNA structures organized in Equivalence Classes (source: Leontis & Zirbel, 2012).

`getAltRepres`: Apply filters (e.g. just protein-RNA structures) to select other representants from the members of each class.

`represAsDataFrame`: From the output of getRNAList or getAltRepres, generate a data.frame in which each row corresponds to a RNA chain, rather than an Equivalence Class.

`pipeNucData`: From a list of RNA structures/chains computes and returns structural data at the level of the nucleotide.

`pipeProtNucData`: From a list of protein-RNA structures computes and returns the interaction sites distances and atoms.

`applyToPDB`: Applies a desired function to a list of PDB IDs.

`queryEntryList`: Returns the whole list of PDB IDs in the database.

`queryObsoleteList`: Returns the list of Obsolete PDB IDs in the database.

`cleanByPucker`: From the output of pipeNucData subsets a desired subset of nucleotides in a given puckering conformation.
&nbsp;

&nbsp;


### Single-structure level

#### **Functions to query PDB data using the PDBe (EMBL-EBI) REST API or a mirror API from the MMB Lab** (All of them take a PDB ID as input)

`queryAuthors`: List of authors.

`queryReldate`: Release date.

`queryDepdate`: Deposition date.

`queryRevdate`: Revision date.

`queryDescription`: Author description.

`queryCompType`: Compound type (e.g. Nuc or Prot-nuc).

`queryChains`: Chain information.

`queryEntities`: Entitity information.

`countEntities`: In a given pdbID it counts the total number of each different kind of entity (RNA, DNA, Protein ...).

`queryFormats`: File formats for the given ID.

`queryHeader`: PDB Header.

`queryHetAtms`: HETATM entities in structure (includes modified residues, ions and ligands).

`hasHetAtm`: Checks wether a a given structure contains a particular HETATM entity. It makes use of queryHetAtms.

`queryModres`: Modified residues.

`queryLigands`: Ligands in structure.

`queryOrgLigands`: Ligands in structure (substracting ions).

`queryResol`: Resolution (if applicable).

`queryTechnique`: Experimental Technique.

`queryStatus`: Released/Obsolete and related status information.

`queryNDBId`: Cross-reference NDB ID.

`queryAPI`: Subfunction of all the previous, which can be used to make alternative queries.
&nbsp;

#### **Classify PDB structures** (PDB ID as input)

`classifyRNA`: Categorizes a structure in different RNA groups.

`classifyDNA`: Categorizes a structure in different DNA groups.
&nbsp;

#### **Input mmCIF data**

`cifDownload`: Downloads structure from Protein Data Bank.

`cifParser`: Reads the 14th common sections of all mmCIF files in the PDB and generates a CIF S4 object.

`cifAsPDB`: Wrapper of cifParser that generates a pdb object (bio3d compatible S3 object).
&nbsp;

#### **CIF accessors** (Find descriptions in mmCIF dicctionary: http://mmcif.wwpdb.org/)

`cifAtom_site`: Access the coordinates of a CIF object (read by cifParser). The resulting object is not compatible with bio3d functions, see cifAsPDB for that.

`cifAtom_sites`

`cifAtom_type`

`cifAudit_author`

`cifAudit_conform`

`cifChem_comp`

`cifDatabase_2`

`cifEntity`

`cifEntry`

`cifExptl`

`cifPdbx_database_status`

`cifStruct`

`cifStruct_asym`

`cifStruct_keywords`
&nbsp;

#### **Structure analysis**

`selectModel`: Selects the model of interest.

`findBindingSite`: Same as pipeProtNucData for a single structure.

`measureEntityDist`: Measures distances between given entities.

`measureElenoDist`: Measures distances between given atoms.

`trimSphere`: Trim a pdb object and a surrounding sphere of atoms.

`trimByID`: Same as trimSphere using the IDs and output of pipeNucData.

`checkNuc`: Checks the integrity of all the nucleotides in a given Nucleic Acid structure.

`measureNuc`: Measures a defult/desired set of distances, angles and torsional angles for a given Nucleic Acid structure.

`rVector`: Computes the rVectors between all nucleobases of a structure (source: Bottaro et al, 2014).

`eRMSD`: Compares structures with the same number of residues using the rVectors (source: Bottaro et al, 2014).

`RMSD`: Compares structures with the RMSD measure.

`dssr`: Wrapper of DSSR software (source: Lu et al, 2015), if installed.
&nbsp;

&nbsp;

### Exploratory analysis

`findHDR`: Finds High Density Regions in a 2D scatter plot

`plot2D`: Scatter plot of angles

`plot3Ddens`: 3D view of the density of 2D data.

`plotCategorical`

`plotCircularDistribution`


## Developers
-------------

Diego Gallego

Eric Matamoros

Leonardo Darré
&nbsp;

&nbsp;

*Molecular Modeling and Bioinformatics Group.*


## Citation
----------
Diego Gallego, Leonardo Darré, Pablo D Dans, Modesto Orozco, VeriNA3d: an R package for nucleic acids data mining, Bioinformatics, Volume 35, Issue 24, 15 December 2019, Pages 5334–5336, https://doi.org/10.1093/bioinformatics/btz553


## License
----------

GPL-3
