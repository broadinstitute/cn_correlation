#!/bin/bash
# set environment variables for generating multiprocessing environment
# arguments MPE_CMD WRAPPER INPUT_DIR OUTPUT_DIR TAGOUT NITERS [SEED] ...

# MPE_CMD is path to file that will generate platform-specific command 
MPE_CMD=$1
shift

# generate environment variables for command generation
MATLAB_WRAP=$1
EXECUTABLE=$2
INPUT_DIR=$3
OUTPUT_DIR=$4
TAGOUT=$5
NITERS=$6
SEED=$7
OUTFILE=$OUTPUT_DIR/$TAGOUT.out.txt
ERRFILE=$OUTPUT_DIR/$TAGOUT.err.txt
WORK_DIR=$INPUT_DIR
source $MPE_CMD $*
