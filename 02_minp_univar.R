library(data.table)
setwd('/project/xuanyao/jinghui/')
args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

all_prot = fread('pqtl/05_h2/00_ref/Sun_2018_prot_w_coor.txt')
chr_i = paste0('chrom_', i, '_')
for (j in 1:nrow(all_prot)) {
  prot_j = all_prot$target[j]
  all_file = list.files(paste0('/project2/xuanyao/pQTL/', prot_j))
  file_j = all_file[grep(chr_i, all_file)]
  sum_stats_j = fread(paste0('/project2/xuanyao/pQTL/', prot_j, '/', file_j))
  sum_stats_j = sum_stats_j[sum_stats_j$`log(P)` < log10(5e-8), ]
  sum_stats_j$file = prot_j
  
  fwrite(sum_stats_j, paste0('pqtl/04_fdr/sun_2018/each_prot/', prot_j, '_chr',
                             i, '.txt'), sep = '\t')
}



