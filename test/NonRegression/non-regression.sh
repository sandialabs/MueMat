#!/bin/sh
#
# This script runs tests of config.txt and check diffs
#
# Usage:
# ./non-regression.sh [dir] [ref] [fast=0,1] [difftool] 
#  
# options:
# - fast: 1 == compare only convergence results
#
# default: 
# - dir = 'current-version'
# - ref = lastest directory of REF/
# - fast = 0
# - difftool = tkdiff
#
# how to create a new version of reference:
# mv current-version REF/`date +%F`-`git rev-parse HEAD | tail -c 10`
#

command -v parallel &> /dev/null || { echo "Error: you need GNU/parallel (http://www.gnu.org/software/parallel/). Aborting."; exit 0; }

if [ $# -ge 1 ]; then ./non-regression-run.sh $1; else ./non-regression-run.sh; fi
./non-regression-check.sh $*
