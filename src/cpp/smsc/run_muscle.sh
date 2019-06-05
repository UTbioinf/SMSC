#!/bin/bash

source $(dirname $0)/config.sh

infile=$1
outfile=$2

${MUSCLE_EXEC} -in ${infile} -out ${outfile} -maxiters ${max_iters} -quiet
