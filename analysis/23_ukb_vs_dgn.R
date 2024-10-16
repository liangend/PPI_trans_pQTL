library(data.table)
library(qvalue)
setwd('/project/xuanyao/jinghui')
dgn_beta_ukb_qtl = fread('pqtl/12_beta_across_two_data/dgn_beta_ukb_pqtl_all.txt')
dgn_beta_ukb_qtl = dgn_beta_ukb_qtl[, -c(4:8,10,16:19)]
dgn_beta_ukb_qtl = dgn_beta_ukb_qtl[dgn_beta_ukb_qtl$dgn_beta != 0, ]
adj_dgn_beta = c()
for (i in 1:nrow(dgn_beta_ukb_qtl)) {
  ukb_A0 = dgn_beta_ukb_qtl$ukb_A0[i]
  ukb_A1 = dgn_beta_ukb_qtl$ukb_A1[i]
  dgn_A0 = dgn_beta_ukb_qtl$dgn_A0[i]
  dgn_A1 = dgn_beta_ukb_qtl$dgn_A1[i]
  if (ukb_A0 == dgn_A1 & ukb_A1 == dgn_A0) {
    adj_dgn_beta[i] = dgn_beta_ukb_qtl$dgn_beta[i]
  } else if (ukb_A0 == dgn_A0 & ukb_A1 == dgn_A1) {
    adj_dgn_beta[i] = -dgn_beta_ukb_qtl$dgn_beta[i]
  } else {
    adj_dgn_beta[i] = 0
  }
}
dgn_beta_ukb_qtl$adj_dgn_beta = adj_dgn_beta

### comparing ukb cis-pqtls with DGN cis-eqtls
# UKB pQTL to DGN gene expression effect
dgn_beta_ukb_cis = dgn_beta_ukb_qtl[dgn_beta_ukb_qtl$ukb_cis_trans == 'cis', ]

# DGN eQTL to UKB protein effect
ukb_beta_dgn_cis = fread('pqtl/12_beta_across_two_data/ukb_beta_dgn_cis_all.txt')
colnames(ukb_beta_dgn_cis) = c("dgn_var", "dgn_gene", "dgn_p", "dgn_beta",
                               "dgn_se", "dgn_af", "dgn_chr", "dgn_bp",
                               "dgn_A0", "dgn_A1", "ukb_A0", "ukb_A1", "ukb_A1_frq", 
                               "ukb_beta", "ukb_se", "ukb_logP")
# clump SNPs
uniq_gene = unique(ukb_beta_dgn_cis$dgn_gene)
ukb_beta_dgn_cis_loci = c()
window_size = 250000
for (i in uniq_gene) {
  ukb_beta_sub = ukb_beta_dgn_cis[ukb_beta_dgn_cis$dgn_gene == i, ]
  while(nrow(ukb_beta_sub > 0)) {
    lead_snp = which.min(ukb_beta_sub$dgn_p)[1]
    lead_chr = ukb_beta_sub$dgn_chr[lead_snp]
    lead_pos = ukb_beta_sub$dgn_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(ukb_beta_sub$dgn_chr == lead_chr & 
                       ukb_beta_sub$dgn_bp >= start_pos &
                       ukb_beta_sub$dgn_bp <= end_pos)
    ukb_beta_dgn_cis_loci = rbind(ukb_beta_dgn_cis_loci, ukb_beta_sub[lead_snp, ])
    ukb_beta_sub = ukb_beta_sub[-rm_index, ]
  }
}
ukb_beta_dgn_cis_loci = na.omit(ukb_beta_dgn_cis_loci)
adj_ukb_beta = c()
for (i in 1:nrow(ukb_beta_dgn_cis_loci)) {
  ukb_A0 = ukb_beta_dgn_cis_loci$ukb_A0[i]
  ukb_A1 = ukb_beta_dgn_cis_loci$ukb_A1[i]
  dgn_A0 = ukb_beta_dgn_cis_loci$dgn_A0[i]
  dgn_A1 = ukb_beta_dgn_cis_loci$dgn_A1[i]
  if (ukb_A0 == dgn_A1 & ukb_A1 == dgn_A0) {
    adj_ukb_beta[i] = ukb_beta_dgn_cis_loci$ukb_beta[i]
  } else if (ukb_A0 == dgn_A0 & ukb_A1 == dgn_A1) {
    adj_ukb_beta[i] = -ukb_beta_dgn_cis_loci$ukb_beta[i]
  } else {
    adj_ukb_beta[i] = 0
  }
}
ukb_beta_dgn_cis_loci$adj_ukb_beta = adj_ukb_beta

# cis-pQTL replicated as cis-eQTL
dgn_beta_ukb_cis = na.omit(dgn_beta_ukb_cis)
sum(dgn_beta_ukb_cis$dgn_p < 0.05/nrow(dgn_beta_ukb_cis)) / nrow(dgn_beta_ukb_cis)
qval_dgn_cis = qvalue(dgn_beta_ukb_cis$dgn_p)
1 - qval_dgn_cis$pi0

# cis-eQTL replicated as cis-pQTL
sum(ukb_beta_dgn_cis_loci$ukb_logP > -log10(0.05/nrow(ukb_beta_dgn_cis_loci))) /
  nrow(ukb_beta_dgn_cis_loci)
qval_ukb_cis = qvalue(10^-(ukb_beta_dgn_cis_loci$ukb_logP))
1 - qval_ukb_cis$pi0

#use more strict threshold (threshold in UKB pQTL) for DGN to make fair comparison
ukb_beta_dgn_cis_str = ukb_beta_dgn_cis_loci[ukb_beta_dgn_cis_loci$dgn_p < 5e-8/12585, ]
sum(ukb_beta_dgn_cis_str$ukb_logP > -log10(0.05/nrow(ukb_beta_dgn_cis_str))) /
  nrow(ukb_beta_dgn_cis_str)
qval_ukb_str = qvalue(10^-(ukb_beta_dgn_cis_str$ukb_logP))
1 - qval_ukb_str$pi0

#use top 808 eQTLs (number of cis-pQTLs tested in DGN)
ukb_beta_dgn_cis_top = ukb_beta_dgn_cis_str[order(ukb_beta_dgn_cis_str$dgn_p)[1:808], ]
sum(ukb_beta_dgn_cis_top$ukb_logP > -log10(0.05/nrow(ukb_beta_dgn_cis_top))) /
  nrow(ukb_beta_dgn_cis_top)
qval_ukb_top = qvalue(10^-(ukb_beta_dgn_cis_top$ukb_logP))
1 - qval_ukb_top$pi0

library(ggplot2)
library(gridExtra)
qval_plot = data.frame(qvalue = c(1 - qval_ukb_cis$pi0, 1 - qval_ukb_str$pi0,
                                  1 - qval_ukb_top$pi0),
                       threshold = c('fdr < 0.05', 'p < UKB thre', 'top'))

ggplot(qval_plot, aes(x = threshold, y = qvalue)) + 
  geom_point(size = 2, color = 'steelblue') +
  labs(x = "", y = "pi1", title = 'pi1 of UKB pval') + 
  ylim(c(0.5,1)) + 
  geom_hline(yintercept = 1 - qval_dgn_cis$pi0, lty = 2) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# qqplot of ukb pvalues under different cutoff
p_sort_all = sort(10^-(ukb_beta_dgn_cis_loci$ukb_logP))
p_sort_str = sort(10^-(ukb_beta_dgn_cis_str$ukb_logP))
p_sort_top = sort(10^-(ukb_beta_dgn_cis_top$ukb_logP))
p_sort = data.frame(p_exp = c(-log10(ppoints(p_sort_all)), 
                              -log10(ppoints(p_sort_str)), 
                              -log10(ppoints(p_sort_top))),
                    p_obs = c(-log10(p_sort_all), 
                              -log10(p_sort_str), 
                              -log10(p_sort_top)),
                    threshold = c(rep('fdr < 0.05', length(p_sort_all)),
                                  rep('p < UKB thre', length(p_sort_str)),
                                  rep('top', length(p_sort_top))))

ggplot(p_sort, aes(x = p_exp, y = p_obs, color = threshold)) + 
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(x = "expected", y = "observed", title = 'qqplot of UKB pval') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# cis effect concordance
dgn_cis_ukb_rep = dgn_beta_ukb_cis[qval_dgn_cis$qvalues < 0.05, ]
dgn_cis_ukb_rep = na.omit(dgn_cis_ukb_rep)
sum(sign(dgn_cis_ukb_rep$adj_dgn_beta) == sign(dgn_cis_ukb_rep$ukb_beta)) / 
  nrow(dgn_cis_ukb_rep)

ukb_cis_dgn_rep = ukb_beta_dgn_cis_loci[qval_ukb_cis$qvalues < 0.05, ]
sum(sign(ukb_cis_dgn_rep$adj_ukb_beta) == sign(ukb_cis_dgn_rep$dgn_beta)) / 
  nrow(ukb_cis_dgn_rep)

ukb_cis_dgn_str_rep = ukb_beta_dgn_cis_str[qval_ukb_str$qvalues < 0.05, ]
sum(sign(ukb_cis_dgn_str_rep$adj_ukb_beta) == sign(ukb_cis_dgn_str_rep$dgn_beta)) / 
  nrow(ukb_cis_dgn_str_rep)

ukb_cis_dgn_top_rep = ukb_beta_dgn_cis_top[qval_ukb_top$qvalues < 0.05, ]
sum(sign(ukb_cis_dgn_top_rep$adj_ukb_beta) == sign(ukb_cis_dgn_top_rep$dgn_beta)) / 
  nrow(ukb_cis_dgn_top_rep)

### comparing ukb trans-pqtls with DGN trans-eqtls 
dgn_beta_ukb_trans = dgn_beta_ukb_qtl[dgn_beta_ukb_qtl$ukb_cis_trans == 'trans', ]
ukb_beta_dgn_trans = fread('pqtl/12_beta_across_two_data/ukb_beta_dgn_trans_fdr01.txt')
colnames(ukb_beta_dgn_trans) = c("dgn_var", "dgn_gene", "dgn_p", "dgn_beta",
                                 "dgn_se", "dgn_af", "dgn_A0", "dgn_A1", 
                                 "dgn_chr", "dgn_bp", "ukb_A0", "ukb_A1", "ukb_A1_frq", 
                                 "ukb_beta", "ukb_se", "ukb_logP")
ukb_beta_dgn_trans$ukb_z = ukb_beta_dgn_trans$ukb_beta / ukb_beta_dgn_trans$ukb_se

# clump SNPs
uniq_gene = unique(ukb_beta_dgn_trans$dgn_gene)
ukb_beta_dgn_trans_loci = c()
window_size = 250000
for (i in uniq_gene) {
  ukb_beta_sub = ukb_beta_dgn_trans[ukb_beta_dgn_trans$dgn_gene == i, ]
  while(nrow(ukb_beta_sub) > 0) {
    lead_snp = which.min(ukb_beta_sub$dgn_p)[1]
    lead_chr = ukb_beta_sub$dgn_chr[lead_snp]
    lead_pos = ukb_beta_sub$dgn_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(ukb_beta_sub$dgn_chr == lead_chr & 
                       ukb_beta_sub$dgn_bp >= start_pos &
                       ukb_beta_sub$dgn_bp <= end_pos)
    ukb_beta_dgn_trans_loci = rbind(ukb_beta_dgn_trans_loci, ukb_beta_sub[lead_snp, ])
    ukb_beta_sub = ukb_beta_sub[-rm_index, ]
  }
}
ukb_beta_dgn_trans_loci = na.omit(ukb_beta_dgn_trans_loci)

adj_ukb_z = c()
for (i in 1:nrow(ukb_beta_dgn_trans_loci)) {
  ukb_A0 = ukb_beta_dgn_trans_loci$ukb_A0[i]
  ukb_A1 = ukb_beta_dgn_trans_loci$ukb_A1[i]
  dgn_A0 = ukb_beta_dgn_trans_loci$dgn_A0[i]
  dgn_A1 = ukb_beta_dgn_trans_loci$dgn_A1[i]
  if (ukb_A0 == dgn_A1 & ukb_A1 == dgn_A0) {
    adj_ukb_z[i] = ukb_beta_dgn_trans_loci$ukb_z[i]
  } else if (ukb_A0 == dgn_A0 & ukb_A1 == dgn_A1) {
    adj_ukb_z[i] = -ukb_beta_dgn_trans_loci$ukb_z[i]
  } else {
    adj_ukb_z[i] = 0
  }
}
ukb_beta_dgn_trans_loci$adj_ukb_z = adj_ukb_z
ukb_beta_dgn_trans_loci = ukb_beta_dgn_trans_loci[ukb_beta_dgn_trans_loci$adj_ukb_z != 0, ]
# trans-pQTL replicated as trans-eQTL
sum(dgn_beta_ukb_trans$dgn_p < 0.05/nrow(dgn_beta_ukb_trans)) / nrow(dgn_beta_ukb_trans)
qval_dgn_trans = qvalue(dgn_beta_ukb_trans$dgn_p)
1 - qval_dgn_trans$pi0
# top 124 pQTLs (the number as DGN eQTL investigated in UKB)
dgn_beta_ukb_trans_top = dgn_beta_ukb_trans[order(dgn_beta_ukb_trans$ukb_logP, decreasing = T)[1:124], ]
sum(dgn_beta_ukb_trans_top$dgn_p < 0.05/nrow(dgn_beta_ukb_trans_top)) / nrow(dgn_beta_ukb_trans_top)
qval_dgn_trans_top = qvalue(dgn_beta_ukb_trans_top$dgn_p)
1 - qval_dgn_trans_top$pi0

# trans-eQTL replicated as trans-pQTL
sum(ukb_beta_dgn_trans_loci$ukb_logP > -log10(0.05/nrow(ukb_beta_dgn_trans_loci))) / 
  nrow(ukb_beta_dgn_trans_loci)
qval_ukb_trans = qvalue(10^-(ukb_beta_dgn_trans_loci$ukb_logP))
1 - qval_ukb_trans$pi0

library(ggplot2)
p_sort_ukb_trans = sort(10^-(ukb_beta_dgn_trans_loci$ukb_logP))
p_sort_table = data.frame(p_exp = -log10(ppoints(p_sort_ukb_trans)),
                          p_obs = -log10(p_sort_ukb_trans))
ggplot(p_sort_table, aes(x = p_exp, y = p_obs)) + 
  geom_point(size = 2, color = 'steelblue') +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(x = "expected", y = "observed", title = 'qqplot of DGN trans-eQTLs in UKB') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#use more strict threshold (threshold in UKB pQTL) for DGN to make fair comparison
ukb_beta_dgn_trans_str = ukb_beta_dgn_trans_loci[ukb_beta_dgn_trans_loci$dgn_p < 5e-8/12585, ]
sum(ukb_beta_dgn_trans_str$ukb_logP > -log10(0.05/nrow(ukb_beta_dgn_trans_str))) /
  nrow(ukb_beta_dgn_trans_str)
qval_ukb_trans_str = qvalue(10^-(ukb_beta_dgn_trans_str$ukb_logP))
1 - qval_ukb_trans_str$pi0

# trans effect concordance
dgn_trans_ukb_rep = dgn_beta_ukb_trans[qval_dgn_trans$qvalues < 0.05, ]
dgn_trans_ukb_rep = na.omit(dgn_trans_ukb_rep)
sum(sign(dgn_trans_ukb_rep$adj_dgn_beta) == sign(dgn_trans_ukb_rep$ukb_beta)) / 
  nrow(dgn_trans_ukb_rep)

dgn_trans_ukb_top_rep = dgn_beta_ukb_trans_top[qval_dgn_trans_top$qvalues < 0.05, ]
dgn_trans_ukb_top_rep = na.omit(dgn_trans_ukb_top_rep)
sum(sign(dgn_trans_ukb_top_rep$adj_dgn_beta) == sign(dgn_trans_ukb_top_rep$ukb_beta)) / 
  nrow(dgn_trans_ukb_top_rep)

ukb_beta_dgn_trans_rep = ukb_beta_dgn_trans_loci[qval_ukb_trans$qvalues < 0.05, ]
sum(sign(ukb_beta_dgn_trans_rep$dgn_beta) == sign(ukb_beta_dgn_trans_rep$adj_ukb_z)) / 
  nrow(ukb_beta_dgn_trans_rep)

ukb_beta_dgn_trans_str_rep = ukb_beta_dgn_trans_str[qval_ukb_trans_str$qvalues < 0.05, ]
sum(sign(ukb_beta_dgn_trans_str_rep$dgn_beta) == sign(ukb_beta_dgn_trans_str_rep$adj_ukb_z)) / 
  nrow(ukb_beta_dgn_trans_str_rep)


## cis trans pi1 summary
pi_table = data.frame(pi1 = c(0.84, 0.11, 0.80, 0.77), 
                      rep_group = c('pQTL in DGN',
                                    'pQTL in DGN',
                                    'eQTL in UKB',
                                    'eQTL in UKB'),
                      cis_trans = c('cis', 'trans',
                                    'cis', 'trans'))
ggplot(pi_table, aes(x = rep_group, 
                          y = pi1, fill = cis_trans, label = pi1)) + 
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.4) +
  geom_text(position=position_dodge(0.4), vjust=-0.2, size = 4.5) + 
  labs(x = "", fill = '', title = 'UKB vs DGN') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## plot ukb vs dgn and eQTLGen together
library(ggplot2)
library(RColorBrewer)
pi_table_eqtlgen = data.frame(pi1 = c(0.86, 0.69, 0.24, 0.50), 
                              rep_group = c('UKB-PPP in eQTLGen', 'eQTLGen in UKB-PPP',
                                            'UKB-PPP in eQTLGen', 'eQTLGen in UKB-PPP'),
                              cis_trans = c('cis', 'cis', 'trans', 'trans'))
pi_table_dgn = data.frame(pi1 = c(0.84, 0.11, 0.80, 0.77), 
                          rep_group = c('UKB-PPP in DGN', 'UKB-PPP in DGN',
                                        'DGN in UKB-PPP', 'DGN in UKB-PPP'),
                          cis_trans = c('cis', 'trans', 'cis', 'trans'))

pi_table = rbind(pi_table_eqtlgen, pi_table_dgn)
ggplot(pi_table, aes(x = factor(rep_group, levels = c('DGN in UKB-PPP', 'eQTLGen in UKB-PPP', 
                                                      'UKB-PPP in DGN', 'UKB-PPP in eQTLGen')), 
                     y = pi1, fill = cis_trans, label = pi1)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.4) +
  geom_text(position=position_dodge(0.4), vjust=-0.2, size = 3.5) + 
  labs(x = "", y = ~ paste(pi[1]), fill = '', title = '') +
  ylim(c(0,1))+
  scale_fill_manual(values = brewer.pal(3,"Dark2")[2:3]) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

