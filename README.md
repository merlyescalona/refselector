
*Reference selector*

© 2017 Merly Escalona (<merlyescalona@uvigo.es>)

University of Vigo, Spain, http://darwin.uvigo.es


# SimPhy/NGSphy refselector

For simulations of targeted-sequencing experiments under a known species/gene
tree distribution, extracts the reference sequences that would have been used as
target in the probe desing.

# Assumptions

- We are working under a [SimPhy](https://github.com/adamallo/simphy) - [NGSphy](https://github.com/merlyescalona/ngsphy) simulation pipeline scenario.
Meaning, it follows hierarchical [SimPhy](https://github.com/adamallo/simphy) 's folder structure and sequence
labeling.

To know more about the simulation pipeline scenario go to:

- [SimPhy: A comprehensive simulator of gene family evolution ](https://github.com/adamallo/simphy)
- [NGSphy: phylogenomic simulation of next-generation sequencing data](https://github.com/merlyescalona/ngsphy)

# Input

- [SimPhy](https://github.com/adamallo/simphy) folder path
- prefix of the existing [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files
- prefix for the output files
- method indicating how to obtain the reference sequences
- (optional) file with the description of the sequences that will be used as reference
- (optional) length of the N sequence that will be used to separate the sequences when concatenated

# Output

- The output will be a directory of [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files
- There should be as many [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files as replicates have been generated for the current [SimPhy](https://github.com/adamallo/simphy) project
- Each file will contain all the selected loci, either concatenated or as a multiple alignment file

# Documentation

Go to the [wiki](https://github.com/merlyescalona/refselector/wiki)
