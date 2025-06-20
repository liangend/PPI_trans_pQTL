setwd('/project/xuanyao/jinghui/pqtl/')
library(data.table)
#library(openxlsx)
args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

## annot each protein target with gene name and corresponding chr for Sun et al., 2018
# prot_list = read.csv('/project/xuanyao/jinghui/pqtl/00_ref/Sun_protein_info.csv')
# prot_to_gene = fread('/project/xuanyao/jinghui/pqtl/00_ref/Sun_prot_to_gene.tsv')
# prot_list$gene = prot_to_gene$To[match(prot_list$UniProt, prot_to_gene$From)]
# gene_meta = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/gene.meta.hg37.gft')
# prot_list$chr = gene_meta$chr[match(prot_list$gene, gene_meta$gene_name)]
# prot_list = prot_list[prot_list$chr %in% paste0('chr', 1:22), ]
# prot_list$tss = gene_meta$start[match(prot_list$gene, gene_meta$gene_name)]

# fwrite(prot_list, '/project/xuanyao/jinghui/pqtl/05_h2/00_ref/Sun_prot_w_coor.txt', sep = '\t')

gtex_meta = fread('/project/xuanyao/jinghui/gtex/06_qtl_z/gtex_blood_meta.txt')
snp_list = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/snplist_w_coor_hg38.txt')

cis_window = 1000000

gene_i = gtex_meta$gene_id[i]
tss_i = gtex_meta$start[i]
chr_i = gtex_meta$chr[i]
sum_stats = fread(paste0('05_h2/01_sumstats/gtex_blood/all/', gene_i, '.sumstats.gz'))
cis_snp = snp_list$SNP[which(paste0('chr', snp_list$CHR) == chr_i & 
                               snp_list$BP > tss_i - cis_window &
                               snp_list$BP < tss_i + cis_window)]
trans_snp1 = snp_list$SNP[which(paste0('chr', snp_list$CHR) != chr_i)]
trans_snp2 = snp_list$SNP[which(paste0('chr', snp_list$CHR) == chr_i & 
                                  (snp_list$BP < tss_i - cis_window |
                                     snp_list$BP > tss_i + cis_window))]
trans_snp = c(trans_snp1, trans_snp2)
cis_sum_stats = sum_stats[sum_stats$SNP %in% cis_snp, ]
cis_sum_stats = na.omit(cis_sum_stats)

trans_sum_stats = sum_stats[sum_stats$SNP %in% trans_snp, ]
trans_sum_stats = na.omit(trans_sum_stats)

fwrite(cis_sum_stats, paste0('05_h2/01_sumstats/gtex_blood/cis/', gene_i, '.sumstats.gz'), sep = '\t')
fwrite(trans_sum_stats, paste0('05_h2/01_sumstats/gtex_blood/trans/', gene_i, '.sumstats.gz'), sep = '\t')


