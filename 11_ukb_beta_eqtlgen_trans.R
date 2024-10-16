library(data.table)
setwd('/project/xuanyao/jinghui')
eqtlgen_trans = fread('/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz')
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_prot$gene_name = sapply(strsplit(ukb_prot$file, '_', fixed = T), '[', 1)

eqtlgen_trans = eqtlgen_trans[eqtlgen_trans$Gene %in% ukb_prot$gene_id, ]
uniq_gene = unique(eqtlgen_trans$Gene)

args = commandArgs(trailingOnly = T)
k = as.numeric(args[1])
i = uniq_gene[k]

eqtlgen_trans_sub = eqtlgen_trans[eqtlgen_trans$Gene == i, ]
eqtlgen_trans_sub = eqtlgen_trans_sub[order(eqtlgen_trans_sub$SNPChr), ]

prot_file = ukb_prot$file[which(ukb_prot$gene_id == i)][1]
file_chr_all = list.files(paste0('pqtl/UKB_PPP/', prot_file))

ukb_beta_eqtlgen_trans_sig = c()
for (j in unique(eqtlgen_trans_sub$SNPChr)) {
  eqtlgen_trans_sub_j = eqtlgen_trans_sub[eqtlgen_trans_sub$SNPChr == j, ]
  chr_j = paste0('chr', j, '_')
  file_chr_read = file_chr_all[grep(chr_j, file_chr_all)]
  summ_stats_j = fread(paste0('pqtl/UKB_PPP/', prot_file, '/', file_chr_read))
  
  hg19_bp_j = sapply(strsplit(summ_stats_j$ID, ':', fixed = T), '[', 2)
  hg19_bp_j = paste0(j, ':', hg19_bp_j)
  
  eqtlgen_hg19_j = paste0(eqtlgen_trans_sub_j$SNPChr, ':', eqtlgen_trans_sub_j$SNPPos)
  
  eqtlgen_trans_sub_j = cbind(eqtlgen_trans_sub_j, 
                              summ_stats_j[match(eqtlgen_hg19_j, hg19_bp_j), 
                                           c('ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE', 'LOG10P')])
  ukb_beta_eqtlgen_trans_sig = rbind(ukb_beta_eqtlgen_trans_sig, eqtlgen_trans_sub_j)
  #print(j)
}
fwrite(ukb_beta_eqtlgen_trans_sig, paste0('pqtl/12_beta_across_two_data/ukb_beta_eqtlgen_sig_trans/', i, '.txt'),
       sep = '\t', col.names = F)









