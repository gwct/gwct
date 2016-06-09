# GWCT
### Genome-Wide Convergence Tester

## Author
#### Gregg Thomas

## About

#### GWCT is a program to count convergenct substitutions in a set of species.

#### This is version Beta 1.0
#### The only real dependency is Python 2.7, but you must also have ancestral sequence reconstructions in a certain format (see below).

Other optional functions depend on Newick Utilities and R

## Installation

Simply download the program and run it. You may want to add the GWCT folder to your $PATH variable for ease of use.

## Usage

For input, you need to have a set of aligned protein sequences for a group of species and a phylogenetic tree.

First: Use the gwct_codeml.py script to reconstruct ancestral sequences:

`python gwct_codeml.py -i [directory containing multiple FASTA alignments] -c [PAML path] -t [tree file] -a 1`

This will create a directory which will be used as the input for GWCT. Run GWCT as follows:

`python gwct.py -i [same input directory]`
