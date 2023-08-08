#!/bin/bash

cd ~/LICHENS_553T_Hessian/
wd=$( pwd )

for i in p*
do
cp mcmctree_LICHENS.ctl $i
cp Telos553T_IQTree_rooted_baseml.tree $i # if using IQTree
cd $i
tree_name=$( ls *tree )
aln_name=$( ls *phy )
sed -i -e 's/ALN/'${aln_name}'/' *ctl 
sed -i -e 's/TREE/'${tree_name}'/' *ctl
cd ..
done
