setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(qvalue)
library(RColorBrewer)
mcl_pqtl = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
interval_mod = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_in_interval.txt')
interval_mod$mod = paste0('mod', interval_mod$mod)

mcl_pqtl_sub = mcl_pqtl[mcl_pqtl$mod %in% interval_mod$mod, ]
mcl_pqtl_sub = mcl_pqtl_sub[order(mcl_pqtl_sub$mod, mcl_pqtl_sub$chr), ]
# uniq_mod = unique(mcl_pqtl_sub$mod)
# ukb_p = c()
# interval_p = c()
# for (i in uniq_mod) {
#   mcl_pqtl_sub_i = mcl_pqtl_sub[mcl_pqtl_sub$mod == i, ]
#   uniq_chr_i = unique(mcl_pqtl_sub_i$chr)
#   mod_interval_i = interval_mod$interval_mod[which(interval_mod$mod == i)[1]]
#   for (j in uniq_chr_i) {
#     mcl_pqtl_sub_ij = mcl_pqtl_sub_i[mcl_pqtl_sub_i$chr == j, ]
#     interval_p_ij = readRDS(paste0('pqtl/03_p/interval/p.mod', mod_interval_i, 
#                                 '.', j, '.rds'))
#     interval_bp = as.numeric(sapply(strsplit(names(interval_p_ij), '_', fixed = T), 
#                                     '[', 2))
#     interval_p_ij = interval_p_ij[match(mcl_pqtl_sub_ij$bp_hg37, interval_bp)]
#     ukb_p = c(ukb_p, mcl_pqtl_sub_ij$pco_p)
#     interval_p = c(interval_p, interval_p_ij)
#   }
#   print(i)
# }
# 
# comp_p = data.frame(ukb_loci = mcl_pqtl_sub$loci, ukb_mod = mcl_pqtl_sub$mod, 
#                     ukb_p = ukb_p, interval_p = interval_p)
# comp_p$interval_mod = interval_mod$interval_mod[match(comp_p$ukb_mod, interval_mod$mod)]
# fwrite(comp_p, 'pqtl/21_rep_in_interval/comp_ukb_interval.txt', sep = '\t')

comp_p = fread('pqtl/21_rep_in_interval/comp_ukb_interval.txt')
rm_loci = fread('pqtl/25_cell_prop_expr/rm_loci_spearman.txt')
## Remove PPI clusters with immune GWAS convergence do not have complete proteins in INTERVAL
immune_converge_mod = c(224,776,866,901,923,953,954,960,997,1035)
comp_p_sub = comp_p[!comp_p$ukb_mod %in% paste0('mod', immune_converge_mod), ]
## Remove loci that are likely due to cell composition effects
comp_p_sub = comp_p_sub[!comp_p_sub$ukb_loci %in% rm_loci$loci, ]
comp_p_sub = na.omit(comp_p_sub)

## replication rate
comp_p_sub = comp_p_sub[order(comp_p_sub$ukb_p), ]
sum(comp_p_sub$interval_p < 0.05/nrow(comp_p_sub))
## pi1 value
qval_all = qvalue(comp_p_sub$interval_p)
1 - qval_all$pi0

top_sig = seq(50,2050,50)
rep_rate = c()
pi1 = c()
for (i in top_sig) {
  comp_p_sub_i = comp_p_sub[1:i, ]
  rep_rate = c(rep_rate, sum(comp_p_sub_i$interval_p < 0.05/nrow(comp_p_sub)) / 
                 nrow(comp_p_sub_i))
  pi1 = c(pi1, 1-qvalue(comp_p_sub_i$interval_p)$pi0)
}
comp_tab = data.frame(n_sig = c(top_sig, top_sig), prop = c(pi1, rep_rate), 
                      group = rep(c('pi1', 'Rep rate'), each = length(top_sig)))
ggplot(data = comp_tab[comp_tab$group == 'pi1' & comp_tab$n_sig > 100, ], 
       aes(x = n_sig, y = prop, color = group)) +
  geom_point() +
  xlim(c(0,2050)) + ylim(c(0.3,0.6)) + 
  labs(x = 'Top n significant PPI trans-pQTLs from UKB-PPP', 
       y = 'pi1 in INTERVAL', title = '', color = '') +
  scale_color_manual(values = brewer.pal(3, "Set1")) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 14),
        axis.text.y = element_text(colour = 'black', size = 14),
        #axis.ticks.y = element_blank(),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


