#!/bin/bash

#PBS -l nodes=30:ppn=1,pmem=24gb,walltime=3:00:00
#PBS -m bea
#PBS -N gnumap-sim
#PBS -M nathanlclement@gmail.com

machfile=$PBS_NODEFILE
nproc=30


PROG="/fslgroup/fslg_genome/software/gnumap_MPI/bin/gnumap-tester"
#GENOME="/fslhome/greece/genomes/hg19/chrX.fa"
GENOME="/fslgroup/fslg_genome/compute/human/hg19/gnumap/gnumap.chrX.m10.saved"
OUTPUT="/fslgroup/fslg_genome/software/gnumap_MPI/junk.out"
#SEQFILES="$(ls /fslgroup/fslg_genome/compute/SNP/sim/chrX-sim*.reads.fa)"
SEQFILES="examples/


#/usr/mpi/fsl_openmpi_gcc-1.4.2/bin/mpiexec -np $nproc -machinefile $machfile \
	$PROG -g \"$(echo $GENOME | sed -e 's/ /,/g')\" -o $OUTPUT \
	-v 1 --snp -m 14 -j 21 -h 10000 -a .95 -p -c 1 \
	--gap_penalty=-3 \
	--read=$GENOME
	\"$(echo $SEQFILES | sed -e 's/ /,/g')\"
