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
#### GENERAL OPTIONS

	-h, --help
		Prints this help file.
		
	-v, --version
		Prints the version of software.

####INDEXING OPTIONS

    --ws window_size
	    Index the reference genome with sliding a window of size window_size
	    (default: 14)


####MAPPING OPTIONS

    --threads t
	    Use t number of cores for mapping the sequences
	    (default: 1)
	    Use 0 to use all the available cores in the system
	    
    --seq file
	    Set the input sequence to file
	    
    -o file
	    Output the mapping record into file
	    (default: output)

    --sl segment_length
        Consider  segments  of length segment_length for sampling seeds.
        (default: 500)

    --sc seed_count
        Sample seeds from first seed_count locations  of  each  segment.
        (default: 250)


###EXAMPLES
####Indexing reference genome:

    $ ./pacfast --index refgen.fasta
    $ ./pacfast --index refgen.fasta --ws 15
####Mapping to the reference genome:

    $ ./pacfast --search refgen.fa --seq reads.fastq
    $ ./pacfast --search refgen.fa --seq reads.fastq --threads 4

##Authors
Faraz Hach (fhach AT sfu DOT ca)  
Ehsan Haghshenas (ehaghshe AT sfu DOT ca)  
Iman Sarrafi (isarrafi AT sfu DOT ca)

Copyright (c) <2015-2020>, Simon Fraser University All rights reserved.
