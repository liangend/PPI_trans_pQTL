library(data.table)
setwd('/project/xuanyao/jinghui')
## GTEx gene list
# gene_meta = fread('gtex/06_qtl_z/gtex_blood_meta.txt')
## DGN gene list
gene_meta = fread('gtex/06_qtl_z/dgn_gene_meta.txt')
gene_meta = gene_meta[gene_meta$chr != 'chrX', ]

## number of SNP in each annotation for each chromosome
cis_snp = c()
trans_1mb_snp = c()
trans_5mb_snp = c()
for (i in 1:nrow(gene_meta)) {
  gene_i = gene_meta$gene_id[i]
  n_snp = fread(paste0('pqtl/05_h2/03_n_snp/dgn/', gene_i, '.txt'))
  cis_snp = rbind(cis_snp, n_snp$nsnp_cis_1mb)
  trans_1mb_snp = rbind(trans_1mb_snp, n_snp$nsnp_trans_1mb)
  trans_5mb_snp = rbind(trans_5mb_snp, n_snp$nsnp_trans_5mb)
}

## calculate the average cis-coefficients (per-snp h2 for each annotation) across all proteins
cis_files = list.files('pqtl/05_h2/02_h2_results/dgn_h2/cis_1mb/', 
                       pattern = 'results')
coeff_cis = c()
for (i in cis_files) {
  h2_i = fread(paste0('pqtl/05_h2/02_h2_results/dgn_h2/cis_1mb/', i))
  coeff_i = h2_i$Coefficient
  coeff_cis = rbind(coeff_cis, coeff_i)
}
annot = h2_i$Category
colnames(coeff_cis) = annot
coeff_cis_ave = colMeans(coeff_cis) 

## calculate cis-h2
cis_h2 = cis_snp %*% coeff_cis_ave

## obtain the trans-coefficients (per-snp h2 for each annotation) across all proteins
trans_files = list.files('pqtl/05_h2/02_h2_results/dgn_h2/trans_5mb/', 
                         pattern = 'results')
coeff_trans = c()
for (i in trans_files) {
  h2_i = fread(paste0('pqtl/05_h2/02_h2_results/dgn_h2/trans_5mb/', i))
  coeff_i = h2_i$Coefficient
  coeff_trans = rbind(coeff_trans, coeff_i)
}
colnames(coeff_trans) = annot
coeff_trans_ave = colMeans(coeff_trans)

## calculate trans-h2
trans_1mb_h2 = trans_1mb_snp %*% coeff_trans_ave
trans_1mb_h2_ind = rowSums(trans_1mb_snp * coeff_trans)
trans_5mb_h2 = trans_5mb_snp %*% coeff_trans_ave
trans_5mb_h2_ind = rowSums(trans_5mb_snp * coeff_trans)

gene_meta$cis_1mb_h2 = cis_h2
gene_meta$trans_1mb_h2 = trans_1mb_h2
gene_meta$trans_1mb_h2_ind = trans_1mb_h2_ind
gene_meta$trans_5mb_h2 = trans_5mb_h2
gene_meta$trans_5mb_h2_ind = trans_5mb_h2_ind
fwrite(gene_meta, 'pqtl/05_h2/04_h2_summ/gene_h2_dgn_1mb_5mb.txt', sep = '\t')


