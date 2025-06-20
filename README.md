# Protein–protein interactions shape trans-regulatory impact of genetic variation on protein expression and complex traits
Code used for analysis of the study, which is published on BioRxiv: https://doi.org/10.1101/2024.10.02.616321
### Code details
01_most_sig_cis_trans.R: Strongest cis- and trans-effects at the protein and mRNA level
02_h2: Scripts calculating and summarizing the cis and trans h2 of protein and mRNA expression
03_pqtl_pLI.R: pLI of pQTL target genes
04_gwas_coloc_ukb_cc.R and 04_gwas_coloc_ukb_quant.R: Colocalization between UKB pQTLs and GWAS of categorical traits and quantitative traits
05_qtl_pi1.R: pi1 values evaluating the replication of QTLs across datasets
06_trans_pqtl_ppi_tf.R: Enrichment of trans-pQTLs in TF and PPI mechanisms
07_snpEff.R and 07_snpEff.sh: SNP annotation of trans-pQTLs and trans-eQTLs
08_ppi_interface.R: Enrichment of coding trans-QTLs at PPI interfaces
09_mcl_clust.sh: Cluster protein pairs into PPI cluster using MCL
10_PPI_trans_pqtl_pipeline: Pipelines to identify PPI trans-pQTLs
11_rm_hotspot_cutoff.R: PPI trans-pQTL hotspots tend to overlap with cell-composition GWAS
12_ppi_pqtl_summ.R: Summary of PPI trans-pQTLs
13_ppi_qtl_p_in_cell_comp.R: p value inflation after removing cell-composition effect
14_rand_vs_mcl_mod.R: Comparison between random clusters and PPI clusters
15_comp_ukb_interval.R Replication of PPI trans-pQTLs in INTERVAL
16_gwas_coloc_pco_cc.R and 16_gwas_coloc_pco_quant.R: Colocalization between PPI trans-pQTLs and GWAS of categorical traits and quantitative traits
17_ppi_coloc_summ.R: Summary of colocalization for PPI trans-pQTLs
