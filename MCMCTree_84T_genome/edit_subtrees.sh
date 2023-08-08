for i in *
do
cd $i
#reformat root age label
sed -i 's/B-0.999-1.001-/B(0.999,1.001)/' *.tre
#add root age label to taxa that dont have the caliciales outgroup
sed -i "s/;/\'B(0.999,1.001)\';/" *.tre
#add newick header to tree
cat <(head -n 1 *.tree) Telos85T_4parts_partitions12.raxml.rooted.subset.tre > Telos85T_4parts_partitions12.raxml.rooted.subset.formatted.tre
rm Telos85T_4parts_partitions12.raxml.rooted.subset.tre
cd ..
done
