#!/bin/bash
#SBATCH -p normal
#SBATCH -J 50-1ADMIX
#SBATCH --time 480:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o out_ad.log
#SBATCH -e e_ad.log
for x in {2..20}; do
admixture -B100 --cv=10 -s time -j15 Extract50aim.bed $x > output$x 
done