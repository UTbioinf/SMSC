#!/bin/bash

source $(dirname $0)/config.sh

exe_prenuc=${mummer_path}/aux_bin/prenuc
exe_mummer=${mummer_path}/mummer
exe_mgaps=${mummer_path}/mgaps
exe_postnuc=${mummer_path}/aux_bin/postnuc
exe_delta_filter=${mummer_path}/delta-filter

arg0=$0

run_prenuc=0
run_mummer=0
run_mgaps_post_nuc=0
run_nucmer=0
run_delta_filter=0
run_final_clean=0

is_forward=1

ref_file=
qry_file=
prefix=

min_match_length=20
min_cluster=200
max_gap=20
diagdiff=5
diagfactor=.12
breaklen=150

cnt_choices=0

function help_msg(){
    echo "Usage: ${arg0} [options]"
    echo "options:"
    echo -e "\t-P               Run prenuc"
    echo -e "\t-M               Run mummer"
    echo -e "\t-C               Run mgaps and post_nuc"
    echo -e "\t-N               Run nucmer"
    echo -e "\t-D               Run delta-filter"
    echo -e "\t-F               Final clean"
    echo -e "\t-R               Run reverse"
    echo -e "\t-r <file>        Set reference file name"
    echo -e "\t-q <file>        Set query file name"
    echo -e "\t-p <prefix>      Set prefix"
    echo -e "\t-l <len>         Set minimum match length"
    echo -e "\t-c <cluster>     Set minimum cluster size"
    echo -e "\t-g <gap>         Set maximum gap size"
    echo -e "\t-d <diagdiff>    Set maximum diagonal difference"
    echo -e "\t-f <diagfactor>  Set maximum diagonal difference factor"
    echo -e "\t-b <breaklen>    Set break length"
    echo -e "\t-h               Show help msg"
}

function inc_and_check_cnt_choice(){
    cnt_choices=$((cnt_choices + 1))
    if [ "${cnt_choices}" -gt "1" ]
    then
        echo "Too many running options. Only one of -P -M -C -N can be set" >&2
        exit 1
    fi
}

function check_set(){
    if [ "$1" == "x" ]
    then
        echo "-$2 should be set" >&2
        exit 1
    fi
}

while getopts ":PMCNDFRr:q:p:l:c:g:d:f:b:h" opt; do
    case $opt in
        P)
            run_prenuc=1
            inc_and_check_cnt_choice
            ;;
        M)
            run_mummer=1
            inc_and_check_cnt_choice
            ;;
        C)
            run_mgaps_post_nuc=1
            inc_and_check_cnt_choice
            ;;
        N)
            run_nucmer=1
            inc_and_check_cnt_choice
            ;;
        D)
            run_delta_filter=1
            inc_and_check_cnt_choice
            ;;
        F)
            run_final_clean=1
            inc_and_check_cnt_choice
            ;;
        R)
            is_forward=0
            ;;
        r)
            ref_file=$OPTARG
            ;;
        q)
            qry_file=$OPTARG
            ;;
        p)
            prefix=$OPTARG
            ;;
        l)
            min_match_length=$OPTARG
            ;;
        c)
            min_cluster=$OPTARG
            ;;
        g)
            max_gap=$OPTARG
            ;;
        d)
            diagdiff=$OPTARG
            ;;
        f)
            diagfactor=$OPTARG
            ;;
        b)
            breaklen=$OPTARG
            ;;
        h)
            help_msg
            exit 1
            ;;
        \?)
            echo "Invalid option -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

shift $((OPTIND-1))



# check_set x${ref_file}    r
# check_set x${prefix}      p
# check_set x${qry_file}    q

if [ "${run_prenuc}" -eq "1" ]
then
    check_set x${ref_file} r
    check_set x${prefix}   p
    ${exe_prenuc} ${ref_file} > ${prefix}.ntref
fi

if [ "${run_nucmer}" -eq "1" ]
then
    check_set x${ref_file}    r
    check_set x${prefix}      p
    check_set x${qry_file}    q
    if [ "${is_forward}" -eq "1" ]
    then
        ${exe_mummer} -maxmatch -L -l ${min_match_length} -n ${prefix}.ntref ${qry_file} 2> ${prefix}.err | \
        ${exe_mgaps} -l ${min_cluster} -s ${max_gap} -d ${diagdiff} -f ${diagfactor} | \
        ${exe_postnuc} -b ${breaklen} ${ref_file} ${qry_file} ${prefix}.forward
    else
        ${exe_mummer} -maxmatch -L -l ${min_match_length} -n ${prefix}.ntref ${qry_file} 2> ${prefix}.err | \
        ${exe_mgaps} -l ${min_cluster} -s ${max_gap} -d ${diagdiff} -f ${diagfactor} | \
        ${exe_postnuc} -b ${breaklen} ${ref_file} ${qry_file} ${prefix}.backward
    fi
fi

if [ "${run_mummer}" -eq "1" ]
then
    check_set x${prefix}      p
    check_set x${qry_file}    q
    ${exe_mummer} -maxmatch -L -l ${min_match_length} -n ${prefix}.ntref ${qry_file} 2> ${prefix}.err \
        > ${prefix}.mums
fi

if [ "${run_mgaps_post_nuc}" -eq "1" ]
then
    check_set x${ref_file}    r
    check_set x${prefix}      p
    check_set x${qry_file}    q
    if [ "${is_forward}" -eq "1" ]
    then
        ${exe_mgaps} -l ${min_cluster} -s ${max_gap} -d ${diagdiff} -f ${diagfactor} \
            < ${prefix}.lmums | \
        ${exe_postnuc} -b ${breaklen} ${ref_file} ${qry_file} ${prefix}.forward
        rm -f ${prefix}.lmums ${prefix}.mums
    else
        ${exe_mgaps} -l ${min_cluster} -s ${max_gap} -d ${diagdiff} -f ${diagfactor} \
            < ${prefix}.rmums | \
        ${exe_postnuc} -b ${breaklen} ${ref_file} ${qry_file} ${prefix}.backward
        rm -f ${prefix}.rmums ${prefix}.mums
    fi
fi

if [ "${run_delta_filter}" -eq "1" ]
then
    check_set x${prefix}      p
    if [ "${is_forward}" -eq "1" ]
    then
        ${exe_delta_filter} -g ${prefix}.forward.delta > ${prefix}.forward.filter
        # rm -f ${prefix}.forward.delta
    else
        ${exe_delta_filter} -g ${prefix}.backward.delta > ${prefix}.backward.filter
        # rm -f ${prefix}.backward.delta
    fi
fi

if [ "${run_final_clean}" -eq "1" ]
then
    check_set x${prefix}      p
    rm -f ${prefix}.ntref
    rm -f ${prefix}.err
fi
