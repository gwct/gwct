# GWCT
### Genome-Wide Convergence Tester

## Author
#### Gregg Thomas

## About

#### GWCT is a program to count convergent amino acid substitutions in a set of species.

#### This is version Beta 1.1
#### The only real dependency is Python 2.7, but you must also have ancestral sequence reconstructions in a certain format (see below).

Other optional functions need the following software:

- Tree pruning for PAML depends on [Newick Utilities](http://cegg.unige.ch/newick_utils)
- Generation of C/D plots depends on [R](https://www.r-project.org/)

## Installation

Simply download the program and run it. You may want to add the GWCT folder to your $PATH variable for ease of use.

## Basic Usage

For input, you need to have a set of aligned protein sequences for a group of species and a phylogenetic tree containing those species.

First: Use the gwct_codeml.py script to reconstruct ancestral sequences:

`python gwct_codeml.py -i [directory containing multiple FASTA alignments] -c [PAML path] -t [tree file] -o [PAML output directory]`

This will run PAML on all files and format the ancestral sequences in a way that can be parsed by the GWCT. 

Then run the GWCT as follows:

`python gwct.py -i [PAML output directory from before] -t [set of target species] -p [ancestral probability threshold] -o [GWCT output directory]`

This will produce three files: conv-sites.txt, div-sites.txt, and uniq-sites.txt that contain a list of convergent, divergent, and unique substitutions in the target species. Note: All target species must be in all alignments.

Or run it as:

`python gwct.py -i [PAML output directory from before] -w [1,2] -p [ancestral probability threshold] -o [GWCT output directory]`

This will produce counts for all pairwise comparisons of tips (`-w 1`) or all nodes (`-w 2`) along with a C/D plot (if R is installed).

## Definitions

**Convergent substitutions** are defined as those in which there is a substitution along two or more target lineages (as inferred by ancestral reconstructions) **to the same state**.

**Divergent substitutions** are defined as those in which there is a substitution along two or more target lineages (as inferred by ancestral reconstructions) **to a different state**.

**Unique substitutions** are defined when the target lineages contain one allele and all other lineages contain other alleles. Unique substitutions are only counted on the tips of the tree.

## Options

#### gwct_codeml.py

| Option | Description | 
| ------ | ----------- |
| -i | A directory containing one or more FASTA formatted alignments. Note that these files MUST have the .fa extension to be read as FASTA files. |
| -p | The FULL path to your PAML directory (ie /Users/Momo/bin/paml/) |
| -t | A file containing a single Newick formatted, rooted species tree. The tree must contain all species that are present in your alignments. |
| --prune | Set this flag and the script will compare each alignment to the species tree. If not all species in the tree are in the alignment, the tree will be pruned. Requires Newick Utilities. |
| -v | Control the amount of output printed to the screen. Print all PAML output (1) or just print a progess bar (0). Default: 1 |
| -o | Output directory name. If the directory is not present, the script will created it for you. If this option is left blank, the output directory will be determined automatically (based on the input directory name). |


#### gwct.py

| Option | Description | 
| ------ | ----------- |
| -i | An input directory containing the contents of a gwct_codeml.py run. In otherwords, the input directory for gwct.py is the output directory from gwct_codeml.py. |
| -t | A list of target species for which you want to count convergent sites. This should be a space separated list of species. For example, "spec1 spec2 spec3" will count convergent sites between those three species. To specify an internal node as a target, use a comma separated list of the species that defined that node. For example, "spec1 spec2,spec3" will count convergent sites between spec1 and the node ancestral to spec2 and spec3. |
| -w | An option to count convergent substitutions between all pairs of tips (1) or all nodes (2) in the tree and generate a C/D plot (if R is installed!). If this option is set, no targets should be given. If this option is set to count convergent sites betweeen all nodes, unique substitutions will not be counted. |
| -p | A probability threshold ranging from 0.0 to 1.0. If set, only ancestral states at or above that threshold will be used. Default: 0.0. |
| -n | The number of processes the GWCT should use. Note that one process is reserved for the main function, so entering `-n 2` is equivalent to `-n 1` and you probably won't see a speed-up. Default: 1. |
| -o | Output directory name. | 
