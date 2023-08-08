#!/bin/bash

for i in `seq 1 32`
do
mkdir run$i
mkdir run$i/mcmctree_GBM
mkdir run$i/mcmctree_ILN
cp mcmctree.ctl run$i/mcmctree_GBM
cp mcmctree_ILN.ctl run$i/mcmctree_ILN
cp in.BV* run$i/mcmctree_GBM/in.BV
cp in.BV* run$i/mcmctree_ILN/in.BV
done
