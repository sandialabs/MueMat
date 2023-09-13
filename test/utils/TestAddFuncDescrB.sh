#!/bin/bash

returnvalue=0

# INVARIANT TEST

if [ -e tmp ]
    then
    echo "Please remove directory tmp"
    exit 1
fi

cp -r ../../src tmp

cd tmp
for mfile in `find . -name "*.m"`
do
  #echo "tmp/$mfile"
  perl ../../../utils/AddFuncDescr.pl $mfile > $mfile.out1.m 2> $mfile.err1.m

  # reapply AddFuncDescr on its output
  perl ../../../utils/AddFuncDescr.pl $mfile.out1.m > $mfile.out2.m 2> $mfile.err2.m

  # test
  diff $mfile.out1.m $mfile.out2.m 2>&1 > /dev/null
  if [ $? -ne 0 ]
      then
      echo "diff tmp/$mfile.out1.m tmp/$mfile.out2.m"
      returnvalue=1
  fi
done
cd -

exit $returnvalue
