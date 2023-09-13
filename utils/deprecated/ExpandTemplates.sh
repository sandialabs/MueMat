#!/bin/sh
#
#   THIS SCRIPT IS CURRENTLY UNUSED IN MUEMAT
#
#
#   A small shell script to do wimpy templates. The script looks
#   for filenames of the form *[*|*] where * is any string. It sticks
#   two new files in the directory PostTemplate for each found file.
#   The new files are identical to the original except occurences of the
#   strings TEMPLATE[aaa|bbb] are replaced by aaa in the first file and
#   bbb in the second file. For example,
#         Vec(TEMPLATE[first|count]:TEMPLATE[last|nextcount-1],:);
#   in the file DCBlkApply[|SubIn] spawns a file PostTemplate/DCBlkApply.m 
#   with
#         Vec(first:last,:);
#   and a second file PostTemplate/DCBlkApplySubIn.m with
#         Vec(count:nextcount-1,:);
#
# Note: this could probably be generalized easily to handle files of
# type   *[*|*|*]  or   *[*|*|*|*]  
#
cd PostTemplate
DIR=`pwd | sed "s/.*PostTemplate/PostTemplate/"`
ls *.m
echo
echo "Remove these files[y/n]?"
read YesOrNo
echo $YesOrNo_
if test `expr $YesOrNo`_ = 'y_' &&
   test `expr $DIR`_ = 'PostTemplate_'
then
  rm *.m
fi
echo
cd ..
for i in *\[*\|*\]
do
echo $i
FILENAME=`echo $i | sed "s/\[\([^|]*\)|\([^|]*\)\]/\1.m/g"`
cat $i | sed "s/TEMPLATE\[\([^|]*\)|\([^|]*\)\]/\1/g" > PostTemplate/$FILENAME
FILENAME=`echo $i | sed "s/\[\([^|]*\)|\([^|]*\)\]/\2.m/g"`
cat $i | sed "s/TEMPLATE\[\([^|]*\)|\([^|]*\)\]/\2/g" > PostTemplate/$FILENAME
done



#cp DmatCBlkApply.template DmatCBlkApply\[FullOut\|SubOut\]
