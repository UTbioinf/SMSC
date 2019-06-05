#!/bin/bash

source $(dirname $0)/config.sh

infile=$1
outfile=$2
java -Xmx${mhap_mem} -XX:MaxPermSize=${mhap_mem} -server -jar ${mhap_path} -s ${infile} > ${outfile} 2> /dev/null

