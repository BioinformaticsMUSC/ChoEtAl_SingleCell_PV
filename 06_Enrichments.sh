# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585, WGCNA list = 6029)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir output_relabel/DGE_MAST/enrichments/

cp utils/geneset/*.RData output_relabel/DGE_MAST/enrichments/
cp utils/Enrichment.r output_relabel/DGE_MAST/enrichments/
cp output_relabel/DGE_MAST/DGE_Sign_ForEnrich_HumanID.txt output_relabel/DGE_MAST/enrichments/
cp output_relabel/DGE_MAST/DGE_Sign_ForEnrich_MouseID.txt output_relabel/DGE_MAST/enrichments/

cd output_relabel/DGE_MAST/enrichments/

mkdir STATS/

# scRNA Disorders
Rscript Enrichment.r -g DGE_Sign_ForEnrich_HumanID.txt -l ALZ_SingleCell_DEGs.RData -p -b 15585 -o STATS/ALZ_SingleCell -W 6 -H 5
Rscript Enrichment.r -g DGE_Sign_ForEnrich_HumanID.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_SingleCell -W 6 -H 5
Rscript Enrichment.r -g DGE_Sign_ForEnrich_HumanID.txt -l GeneSets_BroadTrans_ASD_SingleCell.RData -p -b 15585 -o STATS/ASD_BroadTrans_SingleCell -W 6 -H 5

# Neuropsy
Rscript Enrichment.r -g DGE_Sign_ForEnrich_HumanID.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS -W 6 -H 3
Rscript Enrichment.r -g DGE_Sign_ForEnrich_HumanID.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS -W 6 -H 5
Rscript Enrichment.r -g DGE_Sign_ForEnrich_HumanID.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/ASD_Sfari -W 6 -H 2
Rscript Enrichment.r -g DGE_Sign_ForEnrich_HumanID.txt -l GeneSets_BroadTranscript_ASD.RData -p -b 15585 -o STATS/BROAD_ASD -W 6 -H 6

# Mef2c cKO
Rscript Enrichment.r -g DGE_Sign_ForEnrich_MouseID.txt -l GeneSets_Mef2c_cKO.RData -p -b 15585 -o STATS/MEF2C_cKO -W 6  -H 3
Rscript Enrichment.r -g DGE_Sign_ForEnrich_MouseID.txt -l GeneSets_Mef2c_Het.RData -p -b 15585 -o STATS/MEF2C_Het -W 6  -H 3

# Mef2c ChIP-seq
Rscript Enrichment.r -g DGE_Sign_ForEnrich_MouseID.txt -l MEF2C_WithReelin_GeneSets.RData -p -b 15585 -o STATS/MEF2C_ChIP -W 6 -H 3

# Mef2c single cell
Rscript Enrichment.r -g DGE_Sign_ForEnrich_MouseID.txt -l GeneSets_Mef2c_PVcutrun.RData -p -b 15585 -o STATS/MEF2C_PV_CutRunFishell -W 6 -H 3
Rscript Enrichment.r -g DGE_Sign_ForEnrich_MouseID.txt -l GeneSets_Mef2c_SSTcutrun.RData -p -b 15585 -o STATS/MEF2C_SST_CutRunFishell -W 6 -H 3


rm *.RData
