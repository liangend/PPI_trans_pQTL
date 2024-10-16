library(data.table)
source('/project/xuanyao/jinghui/software/manhattan_fun.R')
setwd('/project/xuanyao/jinghui/pqtl')
pco_univar = fread('04_fdr/ukb/pco_univar_comp_result.txt')
prot_module = fread('11_prot_complex/prot_module_corum_2022.txt')
plot(-log10(pco_univar$pco_p), -log10(pco_univar$univar_p), xlab = '-log(P), PCO', 
     ylab  = '-log(P), univariate', col = 'steelblue')
abline(0, 1, lty = 2)
sum(pco_univar$univar_p > 5e-8/2837)
sum(pco_univar$univar_p > pco_univar$pco_p)
which(pco_univar$univar_p > 1e-6)

par(mfrow = c(1,2))
n_pool = c(80, 2491, 4200, 8520)
n = n_pool[1]
plot_row = pco_univar[n, ]
target = plot_row$univar_gene[1]
mod = as.numeric(sub('mod', '', plot_row$mod[1]))
file = prot_module$file[which(prot_module$mod == mod & prot_module$prot == target)]
chr = as.numeric(sub('chr', '', plot_row$chr[1]))
bp = plot_row$bp[1]
window_size = 1000000

file_plot = list.files(paste0('UKB_PPP/', file))
file_plot = file_plot[grep(paste0('chr', chr, '_'), file_plot)]
univar_test = fread(paste0('UKB_PPP/', file, '/', file_plot))

univar_test = univar_test[, c(1:3, 13)]
univar_test_sub = univar_test[univar_test$GENPOS > bp - window_size &
                                univar_test$GENPOS < bp + window_size, ]
colnames(univar_test_sub) = c('CHR', 'BP', 'SNP', 'P')
univar_test_sub$P = 10 ^ (-univar_test_sub$P)
univar_test_sub$CHR = as.numeric(univar_test_sub$CHR)
univar_test_sub$SNP = 1:nrow(univar_test_sub)

pco_p_test = readRDS(paste0('03_p/ukb_corum_2022/p.mod', mod, '.chr', chr, '.rds'))
pco_test = data.frame(SNP = names(pco_p_test), P = unname(pco_p_test))
pco_test$CHR = as.numeric(sapply(strsplit(pco_test$SNP, ':', fixed = T), '[', 1))
pco_test$BP = as.numeric(sapply(strsplit(pco_test$SNP, ':', fixed = T), '[', 2))
pco_test_sub = pco_test[pco_test$BP > bp - window_size &
                          pco_test$BP < bp + window_size, ]


manhattan_data = rbind(univar_test_sub, pco_test_sub)
manhattan(manhattan_data, xlim = c(min(manhattan_data$BP), max(manhattan_data$BP)),
          highlight = pco_test_sub$SNP, suggestiveline = F, genomewideline = F, 
          col = 'darkgreen')






