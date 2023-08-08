#!bin/bash

cd ~/LICHENS_Hessian/
wd=$( pwd )

for i in `seq 1 4`
do
#soft link directories for the four partitions
ln -s $wd/alignments/Telos_concat_part$i".aln" p0$i
ln -s $wd/trees/*.tre p0$i
#copy control files
cp $wd/control_file/mcmctree_LICHENS.ctl baseml_method1/p0$i
cd baseml_method1/p0$i
#set tree and alignment names in the control files
tree_name=$( ls *tre )
aln_name=$( ls *aln )
sed -i -e 's/ALN/'${aln_name}'/' *ctl 
sed -i -e 's/TREE/'${tree_name}'/' *ctl
cd $wd; done
