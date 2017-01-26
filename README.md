# GNUMAP: Genomic Next-generation Universal MAPper version 4.0

## Installation

  1. Download the source code [from the repo](https://github.com/byucsl/gnumap/archive/bwt.zip) or `git clone` [https://github.com/byucsl/gnumap.git](https://github.com/byucsl/gnumap.git).
  2. In the `gnumap/` directory run `make`.
  3. _Optional:_ add the path to the GNUMAP binary (`./bin/gnumap`) to your PATH variable.

### Mac Installation

To install on a Mac, clone the repository with the following command:
  
  `git clone -b mac https://github.com/byucsl/gnumap.git`
  
## Quick Start

  1. Run `./bin/gnumap` to view a help menu explaining how to run GNUMAP.
  2. To perform an alignment using the example dataset, run the following command:
  
  `./bin/gnumap -g examples/Cel_gen.fa -o gnumap.out -a .9 examples/Cel_gen.reads.1.fq`

  * The -g parameter specifies the genome in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
  * The -o parameter gives the path and file-name in the [SAM format](https://samtools.github.io/hts-specs/).
  * The -a parameter takes a percentage (in the form of a floating point number) and specifies the minimum alignment score that will be accepted for mapped reads.
  * The last parameter is the file containing the reads needed to be mapped in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format). 
  _Note:_ one may list multiple read files and GNUMAP will map the reads from each read file.

## Common Parameters

Here are some of the common parameters used, to see a complete list of parameters refer to the file `./docs/DOCUMENTATION`.

  * -g, --genome=STRING          Genome .fa file(s)
  * -o, --output=STRING          Output file
  * -v, --verbose=INT            Verbose (default=0)
  * -c, --num_proc=INT           Number of processors to run on
  * -m, --mer_size=INT           Mer size (default=0)
  * -j, --jump=INT               The number of bases to jump in the sequence indexing
                                 (default: mer_size)
  * -k, --num_seed=INT           The total number of seed hits that must match to a
                                 location before it is considered for alignment
                                 (default: 2)
  * -h, --max_kmer=INT           Kmers in the reference genome that occur more than this
                                 will not be used in the read mapping
  * --no_nw                      This will disable the Needleman-Wunsch alignments and
                                 only use hit count as the basis for alignment. Score is
				                 calculated by summing the number of hits for a position

## Reference

The paper describing the most recent changes was presented at [BioT 2016](http://biotconf.org) and is titled *GNUMAP 4.0: Space and Time Efficient NGS Read Mapping Using the FM-Index*.

GNUMAP is made by the [Computational Sciences Laboratory](http://csl.cs.byu.edu/) at [Brigham Young University](http://byu.edu). 
