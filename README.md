# GNUMAP: Genomic Next-generation Universal MAPper version 4.0

## Installation

  1. Download the source code [from the repo](https://github.com/byucsl/gnumap/archive/bwt.zip) or `git clone` [https://github.com/byucsl/gnumap.git](https://github.com/byucsl/gnumap.git).
  2. In the gnumap/ directory run `mkdir bin && mkdir obj && make`.
  3. _Optional:_ add the path to the GNUMAP binary (`./bin/gnumap`) to your PATH variable.

## Quick Start

  1. Run `./bin/gnumap` to view a help menu explaining how to run GNUMAP.
  2. To perform an alignment using the example dataset, run the following command:
  
  `./bin/gnumap -g examples/Cel_gen.fa -o gnumap.out.sam -a .9 examples/Cel_gen.reads.1.fq`

    * The -g parameter specifies the genome in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
    * The -o parameter gives the path and file-name in the [SAM format](https://samtools.github.io/hts-specs/).
    * The -a parameter takes a percentage (in the form of a floating point number) and specifies the minimum alignment score that will be accepted for mapped reads.
    * The last parameter is the file containing the reads needed to be mapped in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).
