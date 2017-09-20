
# SimPhy/NGSphy refselector
*Reference selector*

© 2017 Merly Escalona (<merlyescalona@uvigo.es>)

University of Vigo, Spain, http://darwin.uvigo.es

For simulations of targeted-sequencing experiments under a known species/gene
tree distribution, extracts the reference sequences that would have been used as
target in the probe desing.

# Assumptions

- We are working under a SimPhy - NGSphy simulation pipeline scenario.
Meaning, it follows hierarchical SimPhy's folder structure and sequence
labeling.

# Input

- SimPhy folder path
- prefix of the existing FASTA files
- prefix for the output files
- method indicating how to obtain the reference sequences
- (optional) file with the description of the sequences that will be used as reference
- (optional) length of the N sequence that will be used to separate the sequences when concatenated

# Output

- The output will be a directory of FASTA files
- There should be as many FASTA files as replicates have been generated for the current SimPhy project
- Each file will contain all the selected loci, either concatenated or as a multiple alignment file

# Documentation

Go to the [wiki](https://github.com/merlyescalona/refselector/wiki)
