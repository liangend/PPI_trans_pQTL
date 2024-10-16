library(data.table)
setwd('/project/xuanyao/jinghui')

dgn_meta = fread('gtex/06_qtl_z/dgn_gene_meta.txt')
dgn_meta$GeneNameConv = sapply(strsplit(dgn_meta$GeneNameConv, '.', fixed = T), '[', 1)

ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_small_p = fread('pqtl/04_fdr/ukb/ukb_univar_small_p.txt')
ukb_sig = ukb_small_p[ukb_small_p$LOG10P > -log10(5e-8/2923), ]
ukb_sig$gene_id = ukb_prot$gene_id[match(ukb_sig$file, ukb_prot$file)]
ukb_sig = ukb_sig[ukb_sig$gene_id %in% dgn_meta$GeneNameConv, ]

args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

dgn_bim = fread('gtex/00_ref/dgn_geno/all_snp.bim')
dgn_bim = dgn_bim[dgn_bim$V1 == i, ]
ukb_sig_sub = ukb_sig[ukb_sig$CHROM == i, ]
uniq_gene = unique(ukb_sig_sub$gene_id)
for (j in uniq_gene) {
  ukb_sig_j = ukb_sig_sub[ukb_sig_sub$gene_id == j, ]
  gene_set_j = dgn_meta$gene_set[which(dgn_meta$GeneNameConv == j)]
  gene_name_j = dgn_meta$gene_id[which(dgn_meta$GeneNameConv == j)]
  file_dgn_j = paste0('gtex/06_qtl_z/dgn_z/geneSet', gene_set_j, 
                      '.chr', i, '.trans_qtl_pairs.txt.gz')
  sum_stat_j = fread(file_dgn_j)
  sum_stat_j = sum_stat_j[sum_stat_j$phenotype_id == gene_name_j, ]
  
  dgn_bp = as.numeric(sapply(strsplit(sum_stat_j$variant_id, ":", fixed = T), '[', 2))
  ukb_sig_j = cbind(ukb_sig_j, sum_stat_j[match(ukb_sig_j$bp_hg19, dgn_bp), 3:5])
  ukb_sig_j = cbind(ukb_sig_j, dgn_bim[match(ukb_sig_j$bp_hg19, dgn_bim$V4), 5:6])
  colnames(ukb_sig_j)[17:21] = c('dgn_p', 'dgn_beta', 'dgn_se', 'dgn_A0', 'dgn_A1')
  fwrite(ukb_sig_j, paste0('pqtl/12_beta_across_two_data/dgn_beta_ukb_sig/chr', i, 
                           '_', j, '.txt'), sep = '\t')
}






