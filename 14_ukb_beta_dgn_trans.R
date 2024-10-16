library(data.table)
setwd('/project/xuanyao/jinghui')
dgn_trans = fread('gtex/06_qtl_z/dgn_trans/fdr01.trans_qtl_pairs.txt.gz')
dgn_trans = dgn_trans[, -(7:11)]
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_prot$gene_name = sapply(strsplit(ukb_prot$file, '_', fixed = T), '[', 1)
sum(dgn_trans$phenotype_id %in% ukb_prot$gene_name)

dgn_trans = dgn_trans[dgn_trans$phenotype_id %in% ukb_prot$gene_name, ]
uniq_gene = unique(dgn_trans$phenotype_id)

args = commandArgs(trailingOnly = T)
k = as.numeric(args[1])
i = uniq_gene[k]

dgn_trans_sub = dgn_trans[dgn_trans$phenotype_id == i, ]
dgn_trans_sub$chr = sapply(strsplit(dgn_trans_sub$variant_id, ':', fixed = T), '[', 1)
dgn_trans_sub$bp = as.numeric(sapply(strsplit(dgn_trans_sub$variant_id, ':', fixed = T), '[', 2))

prot_file = ukb_prot$file[which(ukb_prot$gene_name == i)][1]
file_chr_all = list.files(paste0('pqtl/UKB_PPP/', prot_file))

ukb_beta_dgn_trans = c()
for (j in unique(dgn_trans_sub$chr)) {
  dgn_trans_sub_j = dgn_trans_sub[dgn_trans_sub$chr == j, ]
  chr_j = paste0('chr', j, '_')
  file_chr_read = file_chr_all[grep(chr_j, file_chr_all)]
  summ_stats_j = fread(paste0('pqtl/UKB_PPP/', prot_file, '/', file_chr_read))
  
  hg19_bp_j = sapply(strsplit(summ_stats_j$ID, ':', fixed = T), '[', 2)
  ukb_hg19_j = paste0(summ_stats_j$CHROM, ':', hg19_bp_j)
  
  dgn_trans_sub_j = cbind(dgn_trans_sub_j, 
                         summ_stats_j[match(dgn_trans_sub_j$variant_id, ukb_hg19_j), 
                                      c('ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE', 'LOG10P')])
  ukb_beta_dgn_trans = rbind(ukb_beta_dgn_trans, dgn_trans_sub_j)
  #print(j)
}
fwrite(ukb_beta_dgn_trans, paste0('pqtl/12_beta_across_two_data/ukb_beta_dgn_trans/', i, '.txt'),
       sep = '\t', col.names = F)



