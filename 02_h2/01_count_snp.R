library(data.table)

setwd('/project/xuanyao/jinghui/pqtl/05_h2')
### annot of all SNPs
# maf = 0.05
# baseline_annot = c()
# for (i in 1:22) {
#   annot_i = fread(paste0('00_ref/baseline/baseline.', i, '.annot.gz'))
#   frq_i = fread(paste0('00_ref/1000G_frq/1000G.mac5eur.', i, '.frq.gz'))
#   annot_i = annot_i[frq_i$FRQ < 1 - maf, ] # remove maf < 0.05
#   baseline_annot = rbind(baseline_annot, annot_i)
#   print(i)
# }
# 
# fwrite(baseline_annot, '/project/xuanyao/jinghui/pqtl/05_h2/00_ref/baseline.frq005.annot.gz', sep = '\t')

baseline_annot = fread('00_ref/baseline.frq005.annot.gz')

#### gene meta of each dataset
## GTEx blood
# gene_meta = fread('../../gtex/06_qtl_z/gtex_blood_meta.txt')
## DGN
gene_meta = fread('../../gtex/06_qtl_z/dgn_gene_meta.txt')
## INTERVAL prot
# gene_meta = fread('00_ref/Sun_2018_prot_w_coor.txt')
# colnames(gene_meta)[c(1,7)] = c('gene_id', 'start')
## UKB prot
# gene_meta = fread('00_ref/ukb_prot_w_coor.txt')
# colnames(gene_meta)[c(1,5)] = c('gene_id', 'start')
# gene_meta$chr = paste0('chr', gene_meta$chr)

cis_window = 1000000
trans_window = 5000000
trans_window2 = 1000000

args = commandArgs(trailingOnly = T)
chr_i = as.numeric(args[1])

gene_meta = gene_meta[gene_meta$chr == paste0('chr', chr_i), ]
for (i in 1:nrow(gene_meta)) {
  gene_i = gene_meta$gene_id[i]
  tss_i = gene_meta$start[i]
  
  cis_chr_snp = baseline_annot$SNP[which(baseline_annot$CHR == chr_i)]
  cis_1mb_snp = baseline_annot$SNP[which(baseline_annot$CHR == chr_i & 
                                           baseline_annot$BP > tss_i - cis_window &
                                           baseline_annot$BP < tss_i + cis_window)]
  trans_snp1 = baseline_annot$SNP[which(baseline_annot$CHR != chr_i)]
  trans_snp2 = baseline_annot$SNP[which(baseline_annot$CHR == chr_i & 
                                          (baseline_annot$BP < tss_i - trans_window |
                                             baseline_annot$BP > tss_i + trans_window))]
  trans_5mb_snp = c(trans_snp1, trans_snp2)
  
  trans_snp2 = baseline_annot$SNP[which(baseline_annot$CHR == chr_i & 
                                          (baseline_annot$BP < tss_i - trans_window2 |
                                             baseline_annot$BP > tss_i + trans_window2))]
  trans_1mb_snp = c(trans_snp1, trans_snp2)
  
  cis_chr_annot = baseline_annot[baseline_annot$SNP %in% cis_chr_snp, ]
  cis_1mb_annot = baseline_annot[baseline_annot$SNP %in% cis_1mb_snp, ]
  trans_chr_annot = baseline_annot[baseline_annot$SNP %in% trans_snp1, ]
  trans_5mb_annot = baseline_annot[baseline_annot$SNP %in% trans_5mb_snp, ]
  trans_1mb_annot = baseline_annot[baseline_annot$SNP %in% trans_1mb_snp, ]
  
  nsnp_cis_chr = colSums(cis_chr_annot[, -(1:4)])
  nsnp_cis_1mb = colSums(cis_1mb_annot[, -(1:4)])
  nsnp_trans_chr = colSums(trans_chr_annot[, -(1:4)])
  nsnp_trans_5mb = colSums(trans_5mb_annot[, -(1:4)])
  nsnp_trans_1mb = colSums(trans_1mb_annot[, -(1:4)])
  
  nsnp = data.frame(annot = names(nsnp_cis_chr),
                    nsnp_cis_chr = nsnp_cis_chr, 
                    nsnp_cis_1mb = nsnp_cis_1mb, 
                    nsnp_trans_chr = nsnp_trans_chr, 
                    nsnp_trans_5mb = nsnp_trans_5mb,
                    nsnp_trans_1mb = nsnp_trans_1mb)
  
  fwrite(nsnp, paste0('03_n_snp/ukb/', gene_i, '.txt'), sep = '\t', row.names = F)
}





