#!/bin/sh

# Updates the files doc/source/ref_*.m so that all classes get included in the documentation.
# This script  will only run from the MueMat home directory.
# See also: makeclassdoc.m

where=`basename $(pwd)`
if [ ${where} != "MueMat" ]; then
  echo "error: run this script from the MueMat home directory"
  exit
fi
docDir="../../doc/source"
cd src
#echo `pwd`
#echo
for dn in *; do
  #echo "processing directory ${dn}"
  if [ -d "$dn" ]; then
    #echo "cwd = ${dn}"
    cd ${dn}
    helpfile="${docDir}/ref_${dn}.m"
    rm -f ${helpfile}
    echo "%% Help for ${dn}" > ${helpfile}
    echo "%" >> ${helpfile}
    for fn in *.m; do
        # Look for any line beginning with classdef, and grab the second field as the classname.
        #echo "   ==> grepping ${fn} <=="
        cn=`grep "^[ ]*classdef" ${fn} | awk '{printf "%s\n",$2}'`
        if [ "${cn}" != "" ]; then
          #echo "           % * <matlab:doc('${cn}') ${cn}>"
          for i in ${cn}; do
            echo "% * <matlab:doc('${i}') ${i}>" >> ${helpfile}
          done
        fi
    done
    cd ..
  fi
done
