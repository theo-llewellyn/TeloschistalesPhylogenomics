#!/bin/bash

cd OG0003067
sed -i 's/\*\.aln/partitions12_OG0003067_Telos.aln/' */*ctl
cd ../OG0001822
sed -i 's/\*\.aln/partitions12_OG0001822_Telos.aln/' */*ctl
cd ../OG0002170
sed -i 's/\*\.aln/partitions12_OG0002170_Telos.aln/' */*ctl
cd ../OG0002971
sed -i 's/\*\.aln/partitions12_OG0002971_Telos.aln/' */*ctl
