#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2048M   # memory per CPU core

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


bin/gnumap-plain -c 16 -g /fslgroup/fslg_genome/compute/human/hg19/chrX.fa --print_all_sam -o chrx.mappings.SAMPLEFILE /fslhome/masaki/fsl_groups/fslg_genome/compute/bodily_paul/ScaffScaffOnHuman/data/NA19240/SAMPLEFILE.fastq
