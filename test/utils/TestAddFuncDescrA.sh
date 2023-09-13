#!/bin/bash

# UNIT TEST

rm -f *.out.m
rm -f *.err.m

for mfile in *.m
do

  echo $mfile
  perl ../../utils/AddFuncDescr.pl $mfile > $mfile.out.m 2> $mdile.err.m
  diff $mfile.out.m REF/$mfile.out.m

done

rm -f *.err.m
