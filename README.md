
*Reference selector*

© 2017 Merly Escalona (<merlyescalona@uvigo.es>)

University of Vigo, Spain, http://darwin.uvigo.es


# SimPhy/NGSphy refselector

For simulations of targeted-sequencing experiments under a known species/gene
tree distribution, extracts the reference sequences that would have been used as
target in the probe desing.

# Assumptions

- We are working under a [SimPhy](https://github.com/adamallo/simphy) - [NGSphy](https://github.com/merlyescalona/ngsphy) simulation pipeline scenario.
Following the same hierarchical folder structure.

- Also, it is assumed that the SimPhy folder project has been compressed using [simphycompress](https://github.com/merlyescalona/simphycompress) and the
length of the concatenation N sequence used is known. To know more about the simulation pipeline scenario go to:
    - [SimPhy's repository](https://github.com/adamallo/simphy), and/or check:
        - Mallo D, de Oliveira Martins L, Posada D (2016) SimPhy: Phylogenomic Simulation of Gene, Locus and Species Trees. *Syst. Biol.* **65**(2) 334-344. doi: [http://dx.doi.org/10.1093/sysbio/syv082](http://dx.doi.org/10.1093/sysbio/syv082)
    - [INDELible's site: ](http://abacus.gene.ucl.ac.uk/software/indelible/), and/or check:
        - Fletcher, W. and Yang, Z.( 2009) INDELible: a flexible simulator of biological sequence evolution. *Mol. Biol. Evol.* **26**(8):1879-1888. doi: [https://doi.org/10.1093/molbev/msp098](https://doi.org/10.1093/molbev/msp098)
    - [NGSphy'S repository](https://github.com/merlyescalona/ngsphy), and/or check:
        - Escalona M, Rocha S, Posada D (2018) NGSphy: phylogenomic simulation of next-generation sequencing data. *Bioinformatics.* bty146. doi: [https://doi.org/10.1093/bioinformatics/bty146](https://doi.org/10.1093/bioinformatics/bty146)
    - [simphycompress repository](https://github.com/merlyescalona/simphycompress)

- Species tree replicates are filtered based on the number of loci (number of sequences/ FASTA files) existing
in each folder.

- True sequences from the SimPhy/INDELible simulation process do not contain N's.

# Input

- [SimPhy](https://github.com/adamallo/simphy) folder path
- prefix of the existing [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files
- prefix for the output files
- method indicating how to obtain the reference sequences
- (optional) length of the N sequence that will be used to separate the sequences when concatenated
- (optional) file with the description of the sequences that will be used as reference.

## Methods for reference selection

Specified method to obtain the reference sequence. Values range from `0-4` ( Default: `0`), where:

- (0): 	Considers the outgroup sequence as the reference loci.
- (1): 	Extracts a specific sequence per locus
- (2): 	Selects a random sequence from the ingroups. Same sequence throughout the loci.
- (3): 	Selects randomly a specie and generates a consensus sequence of the sequences belonging to that species.
- (4): 	Generates a consensus sequences from all the sequences involved (will need parameter `-sdf/--seq-desc-file`)

```
**NOTE:** 	The higher the method number, the longer it will take to generate the reference loci.
```

## Reference description file

Each description should be in a separate line. The order in which the descriptions are organized
is the order that will be considered for the specific replicate (i.e. line 1, replicate 1).
If there are less descriptions than species tree replicates, the remaining references will be
considered as sequence `1_0_0`.

# Output

- The output will be a directory of [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files
- There should be as many [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files as replicates have been generated for the current [SimPhy](https://github.com/adamallo/simphy) project
- Each file will contain all the selected loci, either concatenated or as a multiple alignment file


# Install

- Clone this repository

```
git clone git@github.com:merlyescalona/refselector.git
```

- Chance your current directory to the downloaded folder:

```
cd refselector
```

- Install:

```
python setup.py install --user
```

# Usage

# Documentation

Go to the [wiki](https://github.com/merlyescalona/refselector/wiki)
