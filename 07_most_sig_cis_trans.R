library(data.table)
args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

prot_list = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
target_i = prot_list$file[i]
chr_i = as.numeric(sub('chr', '', prot_list$chr[i]))
tss_i = prot_list$gene_start[i]

cis_window = 1000000
trans_window = 5000000
# file_i = list.files(paste0('/project/xuanyao/jinghui/pqtl/AGES_sum_stats/',
#                            target_i, '/harmonised'))
# sum_stats_i = fread(paste0('/project/xuanyao/jinghui/pqtl/AGES_sum_stats/',
#                            target_i, '/harmonised/', file_i))
# sum_stats_i = sum_stats_i[sum_stats_i$chromosome != 'X', ]

file_i = list.files(paste0('/project/xuanyao/jinghui/pqtl/UKB_PPP/', target_i))
file_i = file_i[-grep('chrX', file_i)]
sum_stats_i = c()
for (j in file_i) {
  file_chr_j = fread(paste0('/project/xuanyao/jinghui/pqtl/UKB_PPP/', target_i, '/', j), 
                     select = c('CHROM', 'GENPOS', 'BETA', 'SE', 'LOG10P'))
  sum_stats_i = rbind(sum_stats_i, file_chr_j)
  #print(j)
}
sum_stats_i$CHROM = as.numeric(sum_stats_i$CHROM)

cis_region = sum_stats_i[sum_stats_i$CHROM == chr_i &
                           sum_stats_i$GENPOS > tss_i - cis_window &
                           sum_stats_i$GENPOS < tss_i + cis_window, ]

trans1 = which(sum_stats_i$CHROM != chr_i)
trans2 = which(sum_stats_i$CHROM == chr_i &
                 (sum_stats_i$GENPOS < tss_i - trans_window |
                 sum_stats_i$GENPOS > tss_i + trans_window))
trans_region = sum_stats_i[c(trans1, trans2), ]

cis_most_sig = which.max(cis_region$LOG10P)
trans_most_sig = which.max(trans_region$LOG10P)

most_sig = rbind(cis_region[cis_most_sig, ], trans_region[trans_most_sig, ])
most_sig$cis_trans = c('cis', 'trans')
most_sig$target_gene = prot_list$gene_id[i]
most_sig$gene_chr = chr_i
most_sig$gene_tss = tss_i

fwrite(most_sig, paste0('/project/xuanyao/jinghui/pqtl/10_cis_trans_beta/ukb/',
                        target_i, '.txt'), sep = '\t', col.names = F)

# res_list = list.files('/project/xuanyao/jinghui/pqtl/10_cis_trans_beta/most_sig_results/')
# res_list = sub('.txt', '', res_list)
# setdiff(prot_list$target, res_list)



