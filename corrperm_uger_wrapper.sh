#!/bin/sh
# UGER array task parameters
#$ -N corrperm
#$ -b y

#~Broad dotkit code:
#~ source /broad/software/scripts/useuse
#~ reuse .matlab_2014a_mcr

MLR=/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2014a
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MLR/bin/glnxa64:$MLR/sys/java/jre/glnxa64/jre/lib/amd64:$MLR/sys/java/jre/glnxa64/jre/lib/amd64/server:$MLR/sys/java/jre/glnxa64/jre/lib/amd64/jli:$MLR/sys/os/glnxa64:$MLR/runtime/glnxa64
export LD_LIBRARY_PATH

# wrapper arguments are: $1 = executable; $2 = input_dir; $3 = output_dir; $4 = iterations
# $SGE_TASK_ID contains the chunk number

# <permute-exe> <output-dir> <ref-dir> <chunk-no> <perms-in-chunk>
if [ ! -e $3/idx_cell.chunk.$SGE_TASK_ID.mat ]
then
    echo $1 $2 $3 chunk.$SGE_TASK_ID $4
    $1 $2 $3 chunk.$SGE_TASK_ID $4
fi


