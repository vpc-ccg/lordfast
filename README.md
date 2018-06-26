# **lordFAST: sensitive and ultra-Fast Alignment Search Tool for LOng noisy Read sequencing Data**

lordFAST is a sensitive tool for mapping long reads with high errorrates.
lordFAST is specially designed for aligning reads from PacBio sequencing technology but provides the  user  the ability to change alignment parameters depending on the reads and application.

## How to install?
In order to install lordFAST, please download the source code from [https://github.com/vpc-ccg/lordfast](https://github.com/vpc-ccg/lordfast) or alternatively clone the repository by running the following command:

    $ git clone https://github.com/vpc-ccg/lordfast.git
Now the code can be compiled easily by running "make" command line which builds the binary file "lordfast".

    $ cd lordfast
    $ make

## How to run?

### SYNOPSIS

    lordfast --index FILE [OPTIONS]
    lordfast --search FILE --seq FILE [OPTIONS]

### OPTIONS
#### INDEXING OPTIONS

    -I, --index STR
        Path to the reference genome file in FASTA format which is supposed to be indexed. [required]

#### MAPPING OPTIONS
    -S, --search STR
        Path to the reference genome file in FASTA format. [required]

    -k, --minAnchorLen INT
        Minimum required length of anchors to be considered. [14]

    -o, --out STR
        Write output to STR file rather than standard output. [stdout]

    -m, --maxRefHit INT
        Ignore anchoring positions with more than INT reference matches. [1000]

    -t, --threads INT
        The number of cores for mapping the sequences.
        Pass 0 to use all the available cores in the system. [1]

#### GENERAL OPTIONS
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

## BUGS
Please report the bugs through lordfast's issues page at [https://github.com/vpc-ccg/lordfast/issues](https://github.com/vpc-ccg/lordfast/issues).

## Authors
Ehsan Haghshenas (ehaghshe AT sfu DOT ca)
</br>
Faraz Hach (fhach AT sfu DOT ca)

Copyright (c) 2018, Simon Fraser University, All rights reserved.
