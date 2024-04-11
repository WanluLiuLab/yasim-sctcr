#!/usr/bin/env bash
set -ue
rm -fr to_jiashan
mkdir to_jiashan
cp -R sim to_jiashan
cp -R data to_jiashan
cp -R ref to_jiashan

rm -fr to_jiashan/ref/Homo_sapiens.GRCh38.cdna.all.fa
rm -fr to_jiashan/ref/Homo_sapiens.GRCh38.pep.all.fa

tar cvf to_jiashan.tar to_jiashan
xz -9 -T0 -vvv to_jiashan.tar
