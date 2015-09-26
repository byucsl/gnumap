#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8048M   # memory per CPU core

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes bin/gnumap-plain -g /fslhome/masaki/compute/human_genome/GCA_000001405.18_GRCh38.p3_genomic.fna --print_all_sam -o y.sim_reads.filtered examples/1.fq
valgrind  --track-origins=yes --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes bin/gnumap-plain -g /fslhome/masaki/compute/human_genome/GCA_000001405.18_GRCh38.p3_genomic.fna --print_all_sam -o y.sim_reads.filtered examples/1.fq
