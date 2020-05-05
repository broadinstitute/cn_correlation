#!/bin/sh
# UGER parameters
#$ -N corrperm
#$ -b y

# set up environment
source /broad/software/scripts/useuse
reuse .matlab_2014a_mcr  #! Matlab_2012b_MCR not on RHEL6

# arguments are: $1 = executable; $2 = input_dir; $3 = output_dir; $4 = iterations
# $SGE_TASK_ID contains the chunk number

# <permute-exe> <output-dir> <ref-dir> <chunk-no> <perms-in-chunk>
if [ ! -e $3/idx_cell.chunk.$SGE_TASK_ID.mat ]
then
    echo $1 $2 $3 chunk.$SGE_TASK_ID $4
    $1 $2 $3 chunk.$SGE_TASK_ID $4
fi


