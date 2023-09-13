#!/bin/sh

#Check src/*/Contents.m

# This script will only run from the MueMat home directory.
where=`basename $(pwd)`
if [ ${where} != "MueMat" ]; then
  echo "error: run this script from the MueMat home directory"
  exit
fi

cd src
# Loop over src/ directories
for dn in *; do
  if [ -d "$dn" ]; then 
      # echo "processing directory ${dn}"
      cd $dn
 
      # Check if every class of the directory is in Contents.m.
      # If a class description is not in Contents.m, a line is added automatically at the end of Contents.m.
      # Additional files (M-file of function etc.) could also be in Contents.m but it is optional.
      # Script is working even if Contents.m doesn't exist.
      for classname in `egrep -l classdef *.m`; do 
	  grep -L '%[ ]*'`basename $classname .m`'[ ]*-.*' Contents.m 2>&1 >/dev/null 
	  
	  if [ $? -ne '0' ]; then
	      # Not found -> Add the following line to Contents.m: 
	      # % classname - first comment found in the M-File (dangerous ...)
	      echo "Class $classname missing in $dn/Contents.m -> added at the end of Contents.m"
	      echo '% ' `basename $classname .m` - `cat $classname  | grep % | head -1 | cut -d '%' -f2-` >> Contents.m
	  fi
      done

      cd - >/dev/null
  fi #if [ -d "$dn" ]
done # for dn in *

cd ..
