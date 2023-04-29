#!/bin/bash
#SBATCH --job-name runseur
#SBATCH -p GPUv100s
#SBATCH -N 1
#SBATCH --mem 380928      # Memory Requirement (MB)
#SBATCH -t 2-0:0:00
#SBATCH -o runseur_%j.out
#SBATCH -e runseur_%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user stefano.berto@utsouthwestern.edu

module purge && module load shared slurm python/3.7.x-anaconda   
source activate stefanoconda3

pyscenic aucell sce_wt.loom Regulome_wt.csv --output AUC_wt.csv --num_workers 40

pyscenic aucell sce_het.loom Regulome_het.csv --output AUC_het.csv --num_workers 40

# END OF THE CODE