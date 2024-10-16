library(data.table)
setwd('/project/xuanyao/jinghui')
gtex_trans = fread('gtex/06_qtl_z/gtex_lcl_trans/all.trans_qtl_pairs.txt.gz')
gtex_trans$phenotype_id = sapply(strsplit(gtex_trans$phenotype_id, '.', fixed = T), 
                                 '[', 1)
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_prot$gene_name = sapply(strsplit(ukb_prot$file, '_', fixed = T), '[', 1)
sum(gtex_trans$phenotype_id %in% ukb_prot$gene_id)

gtex_trans = gtex_trans[gtex_trans$phenotype_id %in% ukb_prot$gene_id, ]
uniq_gene = unique(gtex_trans$phenotype_id)

args = commandArgs(trailingOnly = T)
k = as.numeric(args[1])
i = uniq_gene[k]

gtex_trans_sub = gtex_trans[gtex_trans$phenotype_id == i, ]
gtex_trans_sub$chr = sub('chr', '', sapply(strsplit(gtex_trans_sub$variant_id, '_', fixed = T), '[', 1))
gtex_trans_sub$bp = as.numeric(sapply(strsplit(gtex_trans_sub$variant_id, '_', fixed = T), '[', 2))

prot_file = ukb_prot$file[which(ukb_prot$gene_id == i)][1]
file_chr_all = list.files(paste0('pqtl/UKB_PPP/', prot_file))

ukb_beta_gtex_trans = c()
for (j in unique(gtex_trans_sub$chr)) {
  gtex_trans_sub_j = gtex_trans_sub[gtex_trans_sub$chr == j, ]
  chr_j = paste0('chr', j, '_')
  file_chr_read = file_chr_all[grep(chr_j, file_chr_all)]
  summ_stats_j = fread(paste0('pqtl/UKB_PPP/', prot_file, '/', file_chr_read))
  
  ukb_hg38_j = paste0(summ_stats_j$CHROM, ':', summ_stats_j$GENPOS)
  gtex_hg38_j = paste0(gtex_trans_sub_j$chr, ':', gtex_trans_sub_j$bp)
  
  gtex_trans_sub_j = cbind(gtex_trans_sub_j, 
                          summ_stats_j[match(gtex_hg38_j, ukb_hg38_j), 
                                       c('ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE', 'LOG10P')])
  ukb_beta_gtex_trans = rbind(ukb_beta_gtex_trans, gtex_trans_sub_j)
  #print(j)
}
fwrite(ukb_beta_gtex_trans, paste0('pqtl/12_beta_across_two_data/ukb_beta_gtex_lcl_trans/', i, '.txt'),
       sep = '\t', col.names = F)



