#!/bin/sh
#

# CONFIG
filename=('../DBlkApply.m'  '../DinvBlkApply.m')
opt=(     ''             '-D INV')
nb=$(seq 0 1)



# PREPARE THE TEMPLATE
sed 's/^[ ]*$/%-%/g' DxBlkApply.tmpl.m   > DxBlkApply.tmpl.m.1  # 1 to keep spaces on the file
sed '/^%/ s/ /=+=/g' DxBlkApply.tmpl.m.1 > DxBlkApply.tmpl.m.2  # 2 preprocess comments but prevent gcc preprocessor to remove white space in comments
#sed 's/%/\/\//'      DxBlkApply.tmpl.m.1 > DxBlkApply.tmpl.m.2 # 3 don't preprocess matlab comments
tmpl='DxBlkApply.tmpl.m.2'



# GENERATE FILES
for i in $nb
do

# echo ${filename[$i]} ${opt[$i]}
  f=${filename[$i]}

  gcc ${opt[$i]} -x c -E -P -C $tmpl > $f
  cat -s $f > $f.1; mv $f.1 $f
# sed -i 's/\/\//%/'      $f # 3
  sed -i '/^%/ s/=+=/ /g' $f # 2
  sed -i 's/%-%//g'       $f # 1
  sed -i '/./,$!d'        $f # delete all leading blank lines at top of file

done



# CLEAN
rm *.1 *.2