#!/bin/sh
#$ -N mQTL_meta
#$ -m beas
#$ -l h_vmem=32G
#$ -pe sharedmem 2
#$ -l h_rt=10:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

module load igmm/apps/metal/2011-03-25

cd /exports/eddie/scratch/xshen33/mQTL/mQTL_result

while read mfile; do
  metal ${mfile}
done < ls.metal.script.txt