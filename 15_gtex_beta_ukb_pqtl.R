library(data.table)
library(openxlsx)
setwd('/project/xuanyao/jinghui')
gtex_meta = fread('gtex/06_qtl_z/gtex_lcl_meta.txt')

ukb_pqtl = read.xlsx('pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                     startRow = 5, cols = c(1:3, 7, 9:15, 20), sheet = 10)
colnames(ukb_pqtl)[1] = 'var_id_hg19'
colnames(ukb_pqtl)[12] = 'cis_trans'
ukb_pqtl_sub = ukb_pqtl[ukb_pqtl$Assay.Target %in% gtex_meta$gene_name, ]
ukb_pqtl_sub = ukb_pqtl_sub[ukb_pqtl_sub$CHROM < 23, ]

args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

prot_i = ukb_pqtl_sub$Assay.Target[i]
gene_i = gtex_meta$gene_id[which(gtex_meta$gene_name == prot_i)][1]
chr_i = ukb_pqtl_sub$CHROM[i]

# hg38 bp
bp_i = as.numeric(ukb_pqtl_sub$`GENPOS.(hg38)`[i])
gene_set_i = gtex_meta$gene_set[which(gtex_meta$gene_name == prot_i)]
file_gtex_i = paste0('gtex/06_qtl_z/gtex_lcl_z/geneSet', gene_set_i, 
                    '.chr', chr_i, '.trans_qtl_pairs.txt.gz')
sum_stat_i = fread(file_gtex_i)
sum_stat_i = sum_stat_i[sum_stat_i$phenotype_id == gene_i, ]
gtex_bp = as.numeric(sapply(strsplit(sum_stat_i$variant_id, "_", fixed = T), '[', 2))
gtex_index = which(gtex_bp == bp_i)
if (length(gtex_index) > 0) {
  var_i = sum_stat_i$variant_id[gtex_index[1]]
  beta_i = sum_stat_i$b[gtex_index[1]]
  se_i = sum_stat_i$b_se[gtex_index[1]]
  p_i = sum_stat_i$pval[gtex_index[1]]
} else {
  var_i = 'na'
  beta_i = 0
  se_i = 0
  p_i = 1
}

beta_save = data.frame(ukb_pqtl_sub[i,],
                       corr_var = var_i,
                       corr_beta = beta_i, 
                       corr_se = se_i,
                       corr_p = p_i)
fwrite(beta_save, paste0('pqtl/12_beta_across_two_data/gtex_lcl_beta_ukb_pqtl/', i, '.txt'), 
       col.names = F, sep = '\t')





