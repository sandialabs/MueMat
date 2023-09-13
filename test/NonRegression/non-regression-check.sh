#!/bin/sh
#
# This script check diffs between two directories of outputs
#
# Usage:
# ./non-regression-check.sh [dir] [ref] [fast=0,1] [difftool] 
#  
# options:
# - fast: 1 == compare only convergence results
#
# default: 
# - dir = 'current-version'
# - ref = lastest directory of REF/
# - fast = 0
# - difftool = tkdiff

if [ $# -ge 1 ]; then CUR=$1;      else CUR=current-version; fi
if [ $# -ge 2 ]; then REF=$2;      else REF=REF/`ls REF/ | tail -1`; fi
if [ $# -ge 3 ]; then fast=$3;     else fast=0; fi
if [ $# -ge 4 ]; then difftool=$4; else difftool=tkdiff; fi

command -v $difftool &> /dev/null || { echo "Error: you need $difftool. Aborting."; exit 0; }

for i in `find $REF -type f -printf '%f\n'`; do 

    if [ -f $CUR/$i ]
        then
        if [ $fast -eq 1 ] 
            then
            grep ": ||r||=" $REF/$i > /tmp/muemat-nonregression-ref
            grep ": ||r||=" $CUR/$i > /tmp/muemat-nonregression-cur
        else
            cp $REF/$i /tmp/muemat-nonregression-ref
            cp $CUR/$i /tmp/muemat-nonregression-cur
        fi
        
        diff /tmp/muemat-nonregression-ref /tmp/muemat-nonregression-cur > /dev/null 2> /dev/null
        if [ $? -ne 0 ]
            then
            echo -e "[\e[1;31m!!\e[0m] $i"
            $difftool $REF/$i $CUR/$i
        else
            echo -e "[OK] $i"
        fi
    else
        echo -e "[\e[1;31m??\e[0m] $i"
    fi
    
done
