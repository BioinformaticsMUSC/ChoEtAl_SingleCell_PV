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

pyscenic ctx sce_wt_adj.csv \
mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v9-nr.mgi-m0.001-o0.0.tbl \
--expression_mtx_fname sce_wt.loom \
--output Regulome_wt.csv \
--mask_dropouts \
--num_workers 40


pyscenic ctx sce_het_adj.csv \
mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v9-nr.mgi-m0.001-o0.0.tbl \
--expression_mtx_fname sce_het.loom \
--output Regulome_het.csv \
--mask_dropouts \
--num_workers 40

# END OF THE CODE