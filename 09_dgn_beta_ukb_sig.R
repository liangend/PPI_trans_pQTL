library(data.table)
setwd('/project/xuanyao/jinghui')
args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

dgn_meta = fread('gtex/06_qtl_z/dgn_gene_meta.txt')
dgn_meta$GeneNameConv = sapply(strsplit(dgn_meta$GeneNameConv, '.', fixed = T), 
                               '[', 1)

prot_list = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)

prot_list$tss_hg19 = gene_meta$start[match(prot_list$gene_id, gene_meta$gene_id)]
prot_list = na.omit(prot_list)
prot_list$gene_name = sapply(strsplit(prot_list$file, '_', fixed = T), '[', 1)

prot_list = prot_list[which(prot_list$gene_name %in% dgn_meta$gene_id), ]

dgn_snp = fread('gtex/00_ref/dgn_geno/all_snp.bim')

target_i = prot_list$file[i]
chr_i = as.numeric(sub('chr', '', prot_list$chr[i]))
tss_i = prot_list$tss_hg19[i]

cis_window = 1000000
trans_window = 1000000

file_i = list.files(paste0('pqtl/UKB_PPP/', target_i))
file_i = file_i[-grep('chrX', file_i)]
sum_stats_i = c()
for (j in file_i) {
  file_chr_j = fread(paste0('pqtl/UKB_PPP/', target_i, '/', j), 
                     select = c('CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 
                                'ID', 'BETA', 'SE', 'LOG10P'))
  sum_stats_i = rbind(sum_stats_i, file_chr_j)
  #print(j)
}
sum_stats_i$CHROM = as.numeric(sum_stats_i$CHROM)
sum_stats_i$hg19_bp = sapply(strsplit(sum_stats_i$ID, ':', fixed = T), '[', 2)
sum_stats_i$hg19_bp = as.numeric(sum_stats_i$hg19_bp)
sum_stats_i$ID = paste0(sum_stats_i$CHROM, ':', sum_stats_i$hg19_bp)

sum_stats_i = sum_stats_i[sum_stats_i$ID %in% dgn_snp$V2, ]

cis_region = sum_stats_i[sum_stats_i$CHROM == chr_i &
                           sum_stats_i$hg19_bp > tss_i - cis_window &
                           sum_stats_i$hg19_bp < tss_i + cis_window, ]

trans1 = which(sum_stats_i$CHROM != chr_i)
trans2 = which(sum_stats_i$CHROM == chr_i &
                 (sum_stats_i$hg19_bp < tss_i - trans_window |
                    sum_stats_i$hg19_bp > tss_i + trans_window))
trans_region = sum_stats_i[c(trans1, trans2), ]

cis_most_sig = which.max(cis_region$LOG10P)[1]
trans_most_sig = which.max(trans_region$LOG10P)[1]

most_sig = rbind(cis_region[cis_most_sig, ], trans_region[trans_most_sig, ])
most_sig$cis_trans = c('cis', 'trans')
most_sig$target_gene_id = prot_list$gene_id[i]
most_sig$target_gene_name = prot_list$gene_name[i]
most_sig$gene_chr = chr_i
most_sig$gene_tss = tss_i
most_sig$dgn_A1 = dgn_snp$V5[match(most_sig$ID, dgn_snp$V2)]
most_sig$dgn_A2 = dgn_snp$V6[match(most_sig$ID, dgn_snp$V2)]

gene_set_i = dgn_meta$gene_set[which(dgn_meta$gene_id == prot_list$gene_name[i])][1]
file_cis = paste0('gtex/06_qtl_z/dgn_z/geneSet', gene_set_i, 
                  '.chr', most_sig$CHROM[1], '.trans_qtl_pairs.txt.gz')
file_trans = paste0('gtex/06_qtl_z/dgn_z/geneSet', gene_set_i, 
                    '.chr', most_sig$CHROM[2], '.trans_qtl_pairs.txt.gz')

sum_stat_cis = fread(file_cis)
sum_stat_cis = sum_stat_cis[sum_stat_cis$phenotype_id == prot_list$gene_name[i], ]
sun_cis_bp = as.numeric(sapply(strsplit(sum_stat_cis$variant_id, ":", fixed = T), '[', 2))
sun_cis_index = which(sun_cis_bp == most_sig$hg19_bp[1])
if (length(sun_cis_index) > 0) {
  beta_cis = sum_stat_cis$b[sun_cis_index[1]]
  se_cis = sum_stat_cis$b_se[sun_cis_index[1]]
} else {
  beta_cis = 0
  se_cis = 0
}

sum_stat_trans = fread(file_trans)
sum_stat_trans = sum_stat_trans[sum_stat_trans$phenotype_id == prot_list$gene_name[i], ]
sun_trans_bp = as.numeric(sapply(strsplit(sum_stat_trans$variant_id, ":", fixed = T), '[', 2))
sun_trans_index = which(sun_trans_bp == most_sig$hg19_bp[2])
if (length(sun_trans_index) > 0) {
  beta_trans = sum_stat_trans$b[sun_trans_index[1]]
  se_trans = sum_stat_trans$b_se[sun_trans_index[1]]
} else {
  beta_trans = 0
  se_trans = 0
}
most_sig$dgn_beta = c(beta_cis, beta_trans)
most_sig$dgn_se = c(se_cis, se_trans)

fwrite(most_sig, paste0('pqtl/12_beta_across_two_data/dgn_beta_ukb_sig/',
                        i, '.txt'), sep = '\t', col.names = F)



