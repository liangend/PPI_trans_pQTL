library(data.table)
library(qqman)
n = 42
plot_row = univar_trans_minp[n, ]
target = plot_row$V4[1]
chr = as.numeric(sub('chr', '', plot_row$V1[1]))
bp = plot_row$V2[1]
window_size = 1000000

file_plot = list.files(paste0('AGES_sum_stats/', target, '/harmonised/'))
univar_test = fread(paste0('AGES_sum_stats/', target, '/harmonised/', file_plot))

univar_test = univar_test[, 1:4]
univar_test = univar_test[univar_test$chromosome == chr, ]
univar_test_sub = univar_test[univar_test$base_pair_location > bp - window_size &
                                univar_test$base_pair_location < bp + window_size, ]
colnames(univar_test_sub) = c('SNP', 'P', 'CHR', 'BP')
univar_test_sub$CHR = as.numeric(univar_test_sub$CHR)

mod = plot_row$pco_mod[1]
pco_p_test = readRDS(paste0('03_p/p.mod', mod, '.chr', chr, '.rds'))
pco_test = data.frame(SNP = names(pco_p_test), P = unname(pco_p_test))
pco_test$CHR = as.numeric(sapply(strsplit(pco_test$SNP, ':', fixed = T), '[', 1))
pco_test$BP = as.numeric(sapply(strsplit(pco_test$SNP, ':', fixed = T), '[', 2))
pco_test_sub = pco_test[pco_test$BP > bp - window_size &
                          pco_test$BP < bp + window_size, ]

manhattan_data = rbind(univar_test_sub, pco_test_sub)
manhattan(manhattan_data, xlim = c(min(manhattan_data$BP), max(manhattan_data$BP)),
          highlight = pco_test_sub$SNP)

par(mfrow = c(2,2))

