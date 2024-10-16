library(data.table)
setwd('/project/xuanyao/jinghui')
gtex_cis = fread('gtex/00_ref/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz')
gtex_cis = gtex_cis[, c(1,2,3,6,7,8,9)]
colnames(gtex_cis)[c(1,5)] = c('SNP', 'P') 
gtex_cis$gene_id = sapply(strsplit(gtex_cis$gene_id, '.', fixed = T), '[', 1)

ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_prot$gene_name = sapply(strsplit(ukb_prot$file, '_', fixed = T), '[', 1)

gtex_cis = gtex_cis[gtex_cis$gene_id %in% ukb_prot$gene_id, ]
uniq_gene = unique(gtex_cis$gene_id)

args = commandArgs(trailingOnly = T)
k = as.numeric(args[1])
i = uniq_gene[k]

gtex_cis_sub = gtex_cis[gtex_cis$gene_id == i, ]
gtex_cis_sub$chr = sapply(strsplit(gtex_cis_sub$SNP, '_', fixed = T), '[', 1)
gtex_cis_sub$chr = sub('chr', '', gtex_cis_sub$chr)
gtex_cis_sub$bp = as.numeric(sapply(strsplit(gtex_cis_sub$SNP, '_', fixed = T), '[', 2))
gtex_cis_sub$A0 = sapply(strsplit(gtex_cis_sub$SNP, '_', fixed = T), '[', 3)
gtex_cis_sub$A1 = sapply(strsplit(gtex_cis_sub$SNP, '_', fixed = T), '[', 4)

prot_file = ukb_prot$file[which(ukb_prot$gene_id == i)][1]
file_chr_all = list.files(paste0('pqtl/UKB_PPP/', prot_file))

ukb_beta_gtex_cis = c()
for (j in unique(gtex_cis_sub$chr)) {
  gtex_cis_sub_j = gtex_cis_sub[gtex_cis_sub$chr == j, ]
  chr_j = paste0('chr', j, '_')
  file_chr_read = file_chr_all[grep(chr_j, file_chr_all)]
  summ_stats_j = fread(paste0('pqtl/UKB_PPP/', prot_file, '/', file_chr_read))
  
  hg38_bp_j = paste0(summ_stats_j$CHROM, ':', summ_stats_j$GENPOS)
  gtex_hg38_j = paste0(gtex_cis_sub_j$chr, ':', gtex_cis_sub_j$bp)
  
  gtex_cis_sub_j = cbind(gtex_cis_sub_j, 
                         summ_stats_j[match(gtex_hg38_j, hg38_bp_j), 
                                      c('ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE', 'LOG10P')])
  ukb_beta_gtex_cis = rbind(ukb_beta_gtex_cis, gtex_cis_sub_j)
  #print(j)
}
fwrite(ukb_beta_gtex_cis, paste0('pqtl/12_beta_across_two_data/ukb_beta_gtex_liver_cis/', i, '.txt'),
       sep = '\t', col.names = F)



