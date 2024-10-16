library(data.table)
library(qvalue)
library(openxlsx)
setwd('/project/xuanyao/jinghui')

### comparing ukb cis-pqtls with eqtlgen cis-eqtls
# UKB pQTL to eQTLGen gene expression effect
eqtlgen_z_ukb_cis = fread('pqtl/12_beta_across_two_data/eqtlgen_z_ukb_cis_pqtl.txt')
adj_eqtlgen_z = c()
for (i in 1:nrow(eqtlgen_z_ukb_cis)) {
  ukb_A0 = eqtlgen_z_ukb_cis$ukb_A0[i]
  ukb_A1 = eqtlgen_z_ukb_cis$ukb_A1[i]
  eqtlgen_A0 = eqtlgen_z_ukb_cis$eqtlgen_A0[i]
  eqtlgen_A1 = eqtlgen_z_ukb_cis$eqtlgen_A1[i]
  if (ukb_A0 == eqtlgen_A0 & ukb_A1 == eqtlgen_A1) {
    adj_eqtlgen_z[i] = eqtlgen_z_ukb_cis$eqtlgen_z[i]
  } else if (ukb_A0 == eqtlgen_A1 & ukb_A1 == eqtlgen_A0) {
    adj_eqtlgen_z[i] = -eqtlgen_z_ukb_cis$eqtlgen_z[i]
  } else {
    adj_eqtlgen_z[i] = 0
  }
}
eqtlgen_z_ukb_cis$adj_eqtlgen_z = adj_eqtlgen_z

# eQTLGen eQTL to UKB protein effect
ukb_beta_eqtlgen_cis = fread('pqtl/12_beta_across_two_data/ukb_beta_eqtlgen_cis_all.txt')
colnames(ukb_beta_eqtlgen_cis) = c("eqtlgen_p", "eqtlgen_snp", "eqtlgen_chr", "eqtlgen_pos",
                                   "eqtlgen_A1", "eqtlgen_A0", "eqtlgen_z", "eqtlgen_gene_id",
                                   "eqtlgen_gene_name", "eqtlgen_gene_chr", "eqtlgen_gene_pos", 
                                   "eqtlgen_NrCohorts", "eqtlgen_NrSamples", "eqtlgen_fdr",
                                   "eqtlgen_bon_p", "ukb_A0", "ukb_A1", "ukb_A1_frq",        
                                   "ukb_beta", "ukb_se", "ukb_logP")
# find leading SNP of eQTLGen cis-eQTL
uniq_gene = unique(ukb_beta_eqtlgen_cis$eqtlgen_gene_id)
ukb_beta_eqtlgen_cis_loci = c()
window_size = 250000
for (i in uniq_gene) {
  ukb_beta_sub = ukb_beta_eqtlgen_cis[ukb_beta_eqtlgen_cis$eqtlgen_gene_id == i, ]
  while(nrow(ukb_beta_sub > 0)) {
    lead_snp = which.min(ukb_beta_sub$eqtlgen_p)[1]
    lead_chr = ukb_beta_sub$eqtlgen_chr[lead_snp]
    lead_pos = ukb_beta_sub$eqtlgen_pos[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(ukb_beta_sub$eqtlgen_chr == lead_chr & 
                       ukb_beta_sub$eqtlgen_pos >= start_pos &
                       ukb_beta_sub$eqtlgen_pos <= end_pos)
    ukb_beta_eqtlgen_cis_loci = rbind(ukb_beta_eqtlgen_cis_loci, ukb_beta_sub[lead_snp, ])
    ukb_beta_sub = ukb_beta_sub[-rm_index, ]
  }
}
ukb_beta_eqtlgen_cis_loci = na.omit(ukb_beta_eqtlgen_cis_loci)
adj_eqtlgen_z = c()
for (i in 1:nrow(ukb_beta_eqtlgen_cis_loci)) {
  ukb_A0 = ukb_beta_eqtlgen_cis_loci$ukb_A0[i]
  ukb_A1 = ukb_beta_eqtlgen_cis_loci$ukb_A1[i]
  eqtlgen_A0 = ukb_beta_eqtlgen_cis_loci$eqtlgen_A0[i]
  eqtlgen_A1 = ukb_beta_eqtlgen_cis_loci$eqtlgen_A1[i]
  if (ukb_A0 == eqtlgen_A0 & ukb_A1 == eqtlgen_A1) {
    adj_eqtlgen_z[i] = ukb_beta_eqtlgen_cis_loci$eqtlgen_z[i]
  } else if (ukb_A0 == eqtlgen_A1 & ukb_A1 == eqtlgen_A0) {
    adj_eqtlgen_z[i] = -ukb_beta_eqtlgen_cis_loci$eqtlgen_z[i]
  } else {
    adj_eqtlgen_z[i] = 0
  }
}
ukb_beta_eqtlgen_cis_loci$adj_eqtlgen_z = adj_eqtlgen_z

# cis-pQTL replicated as cis-eQTL
sum(eqtlgen_z_ukb_cis$eqtlgen_p < 0.05/nrow(eqtlgen_z_ukb_cis)) / nrow(eqtlgen_z_ukb_cis)
qval_eqtlgen_cis = qvalue(eqtlgen_z_ukb_cis$eqtlgen_p)
1 - qval_eqtlgen_cis$pi0

# cis-eQTL replicated as cis-pQTL
sum(ukb_beta_eqtlgen_cis_loci$ukb_logP > -log10(0.05/nrow(ukb_beta_eqtlgen_cis_loci))) /
  nrow(ukb_beta_eqtlgen_cis_loci)
qval_ukb_cis = qvalue(10^-(ukb_beta_eqtlgen_cis_loci$ukb_logP))
1 - qval_ukb_cis$pi0

#use more strict threshold (Bonferroni threshold) for eQTLGen to make fair comparison
ukb_beta_eqtlgen_cis_str = ukb_beta_eqtlgen_cis_loci[ukb_beta_eqtlgen_cis_loci$eqtlgen_p < 5e-8/19942, ]
sum(ukb_beta_eqtlgen_cis_str$ukb_logP > -log10(0.05/nrow(ukb_beta_eqtlgen_cis_str))) /
  nrow(ukb_beta_eqtlgen_cis_str)
qval_ukb_cis_str = qvalue(10^-(ukb_beta_eqtlgen_cis_str$ukb_logP))
1 - qval_ukb_cis_str$pi0
# top 1150 (number of cis-pQTLs tested in eQTLGen)
ukb_beta_eqtlgen_cis_top = ukb_beta_eqtlgen_cis_loci[order(ukb_beta_eqtlgen_cis_loci$eqtlgen_p)[1:1150], ]
sum(ukb_beta_eqtlgen_cis_top$ukb_logP > -log10(0.05/nrow(ukb_beta_eqtlgen_cis_top))) /
  nrow(ukb_beta_eqtlgen_cis_top)
qval_ukb_cis_top = qvalue(10^-(ukb_beta_eqtlgen_cis_top$ukb_logP))
1 - qval_ukb_cis_top$pi0

library(ggplot2)
library(gridExtra)
qval_plot = data.frame(qvalue = c(1 - qval_ukb_cis$pi0, 1 - qval_ukb_cis_str$pi0,
                                  1 - qval_ukb_cis_top$pi0),
                       threshold = c('fdr < 0.05', 'p < UKB thre', 'top'))

ggplot(qval_plot, aes(x = threshold, y = qvalue)) + 
  geom_point(size = 2, color = 'steelblue') +
  labs(x = "", y = "pi1", title = 'pi1 of UKB pval') + 
  ylim(c(0.5,1)) + 
  geom_hline(yintercept = 1 - qval_eqtlgen_cis$pi0, lty = 2) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# qqplot of ukb pvalues under different cutoff
p_sort_all = sort(10^-(ukb_beta_eqtlgen_cis_loci$ukb_logP))
p_sort_str = sort(10^-(ukb_beta_eqtlgen_cis_str$ukb_logP))
p_sort_top = sort(10^-(ukb_beta_eqtlgen_cis_top$ukb_logP))
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
eqtlgen_cis_ukb_rep = eqtlgen_z_ukb_cis[qval_eqtlgen_cis$qvalues < 0.05, ]
sum(sign(eqtlgen_cis_ukb_rep$adj_eqtlgen_z) == sign(eqtlgen_cis_ukb_rep$ukb_beta)) / 
  nrow(eqtlgen_cis_ukb_rep)

ukb_cis_eqtlgen_rep = ukb_beta_eqtlgen_cis_loci[qval_ukb_cis$qvalues < 0.05, ]
sum(sign(ukb_beta_eqtlgen_cis_rep$adj_eqtlgen_z) == sign(ukb_beta_eqtlgen_cis_rep$ukb_beta)) / 
  nrow(ukb_beta_eqtlgen_cis_rep)

ukb_cis_eqtlgen_str_rep = ukb_beta_eqtlgen_cis_str[qval_ukb_cis_str$qvalues < 0.05, ]
sum(sign(ukb_cis_eqtlgen_str_rep$adj_eqtlgen_z) == sign(ukb_cis_eqtlgen_str_rep$ukb_beta)) / 
  nrow(ukb_cis_eqtlgen_str_rep)

ukb_cis_eqtlgen_top_rep = ukb_beta_eqtlgen_cis_top[qval_ukb_cis_top$qvalues < 0.05, ]
sum(sign(ukb_cis_eqtlgen_top_rep$adj_eqtlgen_z) == sign(ukb_cis_eqtlgen_top_rep$ukb_beta)) / 
  nrow(ukb_cis_eqtlgen_top_rep)

### comparing ukb trans-pqtls with eqtlgen trans-eqtls 
eqtlgen_z_ukb_trans = fread('pqtl/12_beta_across_two_data/eqtlgen_z_ukb_trans_pqtl.txt')
eqtlgen_z_ukb_trans = eqtlgen_z_ukb_trans[eqtlgen_z_ukb_trans$eqtlgen_A0 != '', ]

ukb_beta_eqtlgen_trans = fread('pqtl/12_beta_across_two_data/ukb_beta_eqtlgen_trans_all.txt')
colnames(ukb_beta_eqtlgen_trans) = c("eqtlgen_p", "eqtlgen_snp", "eqtlgen_chr", "eqtlgen_pos",
                                     "eqtlgen_A1", "eqtlgen_A0", "eqtlgen_z", "eqtlgen_gene_id",
                                     "eqtlgen_gene_name", "eqtlgen_gene_chr", "eqtlgen_gene_pos", 
                                     "eqtlgen_NrCohorts", "eqtlgen_NrSamples", "eqtlgen_fdr",
                                     "eqtlgen_bon_p", "ukb_A0", "ukb_A1", "ukb_A1_frq",        
                                     "ukb_beta", "ukb_se", "ukb_logP")
ukb_beta_eqtlgen_trans$ukb_z = ukb_beta_eqtlgen_trans$ukb_beta / 
  ukb_beta_eqtlgen_trans$ukb_se

# clump eQTLGen SNPs
uniq_gene = unique(ukb_beta_eqtlgen_trans$eqtlgen_gene_id)
ukb_beta_eqtlgen_trans_loci = c()
window_size = 250000
for (i in uniq_gene) {
  ukb_beta_sub = ukb_beta_eqtlgen_trans[ukb_beta_eqtlgen_trans$eqtlgen_gene_id == i, ]
  while(nrow(ukb_beta_sub > 0)) {
    lead_snp = which.min(ukb_beta_sub$eqtlgen_p)[1]
    lead_chr = ukb_beta_sub$eqtlgen_chr[lead_snp]
    lead_pos = ukb_beta_sub$eqtlgen_pos[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(ukb_beta_sub$eqtlgen_chr == lead_chr & 
                       ukb_beta_sub$eqtlgen_pos >= start_pos &
                       ukb_beta_sub$eqtlgen_pos <= end_pos)
    ukb_beta_eqtlgen_trans_loci = rbind(ukb_beta_eqtlgen_trans_loci, ukb_beta_sub[lead_snp, ])
    ukb_beta_sub = ukb_beta_sub[-rm_index, ]
  }
}
ukb_beta_eqtlgen_trans_loci = na.omit(ukb_beta_eqtlgen_trans_loci)

adj_ukb_z = c()
for (i in 1:nrow(ukb_beta_eqtlgen_trans_loci)) {
  ukb_A0 = ukb_beta_eqtlgen_trans_loci$ukb_A0[i]
  ukb_A1 = ukb_beta_eqtlgen_trans_loci$ukb_A1[i]
  eqtlgen_A0 = ukb_beta_eqtlgen_trans_loci$eqtlgen_A0[i]
  eqtlgen_A1 = ukb_beta_eqtlgen_trans_loci$eqtlgen_A1[i]
  if (ukb_A0 == eqtlgen_A0 & ukb_A1 == eqtlgen_A1) {
    adj_ukb_z[i] = ukb_beta_eqtlgen_trans_loci$ukb_z[i]
  } else if (ukb_A0 == eqtlgen_A1 & ukb_A1 == eqtlgen_A0) {
    adj_ukb_z[i] = -ukb_beta_eqtlgen_trans_loci$ukb_z[i]
  } else {
    adj_ukb_z[i] = 0
  }
}
ukb_beta_eqtlgen_trans_loci$adj_ukb_z = adj_ukb_z

# trans-pQTL replicated as trans-eQTL
sum(eqtlgen_z_ukb_trans$eqtlgen_p < 0.05/nrow(eqtlgen_z_ukb_trans)) / nrow(eqtlgen_z_ukb_trans)
qval_eqtlgen_trans = qvalue(eqtlgen_z_ukb_trans$eqtlgen_p)
1 - qval_eqtlgen_trans$pi0

# top 500 trans-pQTL replicated as trans-eQTL
ukb_qtl = read.xlsx('pqtl/00_ref/2022_ukbiobank_prot.xlsx', sheet = 10, startRow = 5)
ukb_qtl$snp_gene_pair = paste0(ukb_qtl$`Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)`, ':', ukb_qtl$Assay.Target)
eqtlgen_z_ukb_trans$snp_gene_pair = paste0(eqtlgen_z_ukb_trans$var_id_hg19, ':', eqtlgen_z_ukb_trans$Assay.Target)

eqtlgen_z_ukb_trans$ukb_log10p = ukb_qtl$`log10(p).(discovery)`[match(eqtlgen_z_ukb_trans$snp_gene_pair, 
                                                                      ukb_qtl$snp_gene_pair)]
eqtlgen_z_ukb_trans = eqtlgen_z_ukb_trans[order(eqtlgen_z_ukb_trans$ukb_log10p, decreasing = T), ]
qval_eqtlgen_trans_top500 = qvalue(eqtlgen_z_ukb_trans$eqtlgen_p[1:500])
1 - qval_eqtlgen_trans_top500$pi0

# trans-eQTL replicated as trans-pQTL
sum(ukb_beta_eqtlgen_trans_loci$ukb_logP > -log10(0.05/nrow(ukb_beta_eqtlgen_trans_loci))) / 
  nrow(ukb_beta_eqtlgen_trans_loci)
qval_ukb_trans = qvalue(10^-(ukb_beta_eqtlgen_trans_loci$ukb_logP))
1-qval_ukb_trans$pi0

ukb_beta_eqtlgen_trans_str = ukb_beta_eqtlgen_trans_loci[ukb_beta_eqtlgen_trans_loci$eqtlgen_p < 5e-8/19942, ]
ukb_beta_eqtlgen_trans_str = ukb_beta_eqtlgen_trans_str[order(ukb_beta_eqtlgen_trans_str$eqtlgen_p), ]
sum(ukb_beta_eqtlgen_trans_str$ukb_logP > -log10(0.05/nrow(ukb_beta_eqtlgen_trans_str))) / 
  nrow(ukb_beta_eqtlgen_trans_str)
qval_ukb_trans_str = qvalue(10^-(ukb_beta_eqtlgen_trans_str$ukb_logP))
1-qval_ukb_trans_str$pi0
# top 500 trans-eQTL replicated as trans-pQTL
qval_ukb_trans_top500 = qvalue(10^-(ukb_beta_eqtlgen_trans_str$ukb_logP[1:500]))
1-qval_ukb_trans_top500$pi0

# trans effect concordance
ukb_beta_eqtlgen_trans_rep = ukb_beta_eqtlgen_trans_loci[qval_ukb_trans$qvalues < 0.05, ]
sum(sign(ukb_beta_eqtlgen_trans_rep$eqtlgen_z) == sign(ukb_beta_eqtlgen_trans_rep$adj_ukb_z)) / 
  nrow(ukb_beta_eqtlgen_trans_rep)

ukb_beta_eqtlgen_trans_str_rep = ukb_beta_eqtlgen_trans_str[qval_ukb_trans_str$qvalues < 0.05, ]
sum(sign(ukb_beta_eqtlgen_trans_str_rep$eqtlgen_z) == sign(ukb_beta_eqtlgen_trans_str_rep$adj_ukb_z)) / 
  nrow(ukb_beta_eqtlgen_trans_str_rep)

## cis trans pi1 summary
# pi_table = data.frame(pi1 = c(1 - qval_eqtlgen_cis$pi0, 
#                               1 - qval_ukb_cis_str$pi0, 
#                               1 - qval_eqtlgen_trans$pi0, 
#                               1-qval_ukb_trans_str$pi0), 
#                       rep_group = c('pQTL in eQTLGen',
#                                     'eQTL in UKB',
#                                     'pQTL in eQTLGen',
#                                     'eQTL in UKB'),
#                       cis_trans = c('cis', 'cis', 'trans', 'trans'))

pi_table = data.frame(pi1 = c(0.86, 0.69, 0.24, 0.50), 
                      rep_group = c('pQTL in eQTLGen',
                                    'eQTL in UKB',
                                    'pQTL in eQTLGen',
                                    'eQTL in UKB'),
                      cis_trans = c('cis', 'cis', 'trans', 'trans'))

ggplot(pi_table, aes(x = rep_group, y = pi1, fill = cis_trans, label = round(pi1, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.4) +
  geom_text(position=position_dodge(0.4), vjust=-0.2, size = 4.5) + 
  labs(x = "", fill = '', title = 'UKB vs eQTLGen') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')





