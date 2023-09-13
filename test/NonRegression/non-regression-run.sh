#!/bin/sh
#
# This script runs tests of config.txt and stores outputs in the directory $outputdir
# You need GNU/parallel to run this script
#
# Usage:
# ./non-regression-run.sh [dir] (default: dir==current-version)

if [ $# -ge 1 ]; then outputdir=$1; else outputdir=current-version; fi

command -v parallel &> /dev/null || { echo "Error: you need GNU/parallel (http://www.gnu.org/software/parallel/). Aborting."; exit 0; }

outputdir=`pwd`/$outputdir # the following script need a full path

mkdir -p $outputdir
parallel "cd ../..; echo {}; matlab -nodisplay -nojvm -r \"try, {}; catch end; quit\" | tail -n+6 > $outputdir/{}.txt; cd - > /dev/null" < config.txt


