# **lordFAST: sensitive and Fast Alignment Search Tool for LOng noisy Read sequencing Data**

lordFAST is a sensitive tool for mapping long reads with high error rates.
lordFAST is specially designed for aligning reads from PacBio sequencing technology but provides the  user  the ability to change alignment parameters depending on the reads and application.

## How to install?

[![Build Status](https://travis-ci.org/vpc-ccg/lordfast.svg?branch=master)](https://travis-ci.org/vpc-ccg/lordfast)

### Requirements
- GCC ≥ 4.4.7
- zlib

In order to install lordFAST, please download the source code from [https://github.com/vpc-ccg/lordfast](https://github.com/vpc-ccg/lordfast) or alternatively clone the repository by running the following command:

    $ git clone https://github.com/vpc-ccg/lordfast.git
Now the code can be compiled easily by running `make` command line which builds the binary file `lordfast`.

    $ cd lordfast
    $ make

## How to run?

### SYNOPSIS

    lordfast --index FILE [OPTIONS]
    lordfast --search FILE --seq FILE [OPTIONS]

### OPTIONS
Run `lordfast -h` or `man ./HELP.man` to see available options.

#### Indexing options
    -I, --index STR
        Path to the reference genome file in FASTA format which is supposed to be indexed. [required]

#### Mapping options
    -S, --search STR
        Path to the reference genome file in FASTA format. [required]

    -s, --seq STR
        Path to the file containing read sequences in FASTA/FASTQ format. [required]

    -o, --out STR
        Write output to STR file rather than standard output. [stdout]

    -t, --threads INT
        Use INT number of CPU cores. Pass 0 to use all the available cores. [1]

#### Advanced options
    -k, --minAnchorLen INT
        Minimum required length of anchors to be considered. [14]

    -n, --numMap INT
        Perform alignment for at most INT candidates. [10]

    -l, --minReadLen INT
        Do not try to map any read shorter than INT bp and report them as unmapped. [1000]

    -c, --anchorCount INT
        Consider INT anchoring positions on the long read. [1000]

    -m, --maxRefHit INT
        Ignore anchoring positions with more than INT reference hits. [1000]

    -R, --readGroup STR
        SAM read group line in a format like '@RG\tID:foo\tSM:bar'. []

    -a, --chainAlg INT
        Chaining algorithm to use. Options are "dp-n2" and "clasp". [dp-n2]

    --noSamHeader
        Do not print sam header in the output.

#### Other options
    -h, --help
        Prints this help file.

    -v, --version
        Prints the version of software.


### EXAMPLES
#### Indexing reference genome:

    $ ./lordfast --index refgen.fasta

#### Mapping to the reference genome:

    $ ./lordfast --search refgen.fa --seq reads.fastq > map.sam
    $ ./lordfast --search refgen.fa --seq reads.fastq --threads 4 > map.sam

## Publication
*Haghshenas E., Sahinalp S.C. and Hach F.*, "lordFAST: sensitive and Fast Alignment Search Tool for LOng noisy Read sequencing Data" Bioinformatics (2018) DOI: [10.1093/bioinformatics/bty544](https://doi.org/10.1093/bioinformatics/bty544)

## Bugs
Please report the bugs through lordFAST's issues page at [https://github.com/vpc-ccg/lordfast/issues](https://github.com/vpc-ccg/lordfast/issues).

## Contact
Ehsan Haghshenas (ehaghshe AT sfu DOT ca)

## Copyright and License
This software is released under GNU General Public License (v3.0)\
Copyright (c) 2018 Simon Fraser University
- BWA (used for the BWT-based index) is developed by Heng Li and is licensed under GPL
- Edlib (used for global alignment) is developed by Martin Sosic and is licensed under MIT
- ksw (used for alignment extension) is licensed under MIT
- clasp (can be used for chaining) is developed and copyrighted by Christian Otto
