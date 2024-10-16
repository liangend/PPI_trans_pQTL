library(data.table)
library(qvalue)
setwd('/project/xuanyao/jinghui')
gtex_beta_ukb_qtl = fread('pqtl/12_beta_across_two_data/gtex_lcl_beta_ukb_pqtl_all.txt')
gtex_beta_ukb_qtl = gtex_beta_ukb_qtl[, -c(1,4,6,7)]
gtex_beta_ukb_qtl = gtex_beta_ukb_qtl[gtex_beta_ukb_qtl$gtex_beta != 0, ]
adj_gtex_beta = c()
for (i in 1:nrow(gtex_beta_ukb_qtl)) {
  ukb_A0 = gtex_beta_ukb_qtl$ukb_A0[i]
  ukb_A1 = gtex_beta_ukb_qtl$ukb_A1[i]
  gtex_A0 = gtex_beta_ukb_qtl$gtex_A0[i]
  gtex_A1 = gtex_beta_ukb_qtl$gtex_A1[i]
  if (ukb_A0 == gtex_A1 & ukb_A1 == gtex_A0) {
    adj_gtex_beta[i] = gtex_beta_ukb_qtl$gtex_beta[i]
  } else if (ukb_A0 == gtex_A0 & ukb_A1 == gtex_A1) {
    adj_gtex_beta[i] = -gtex_beta_ukb_qtl$gtex_beta[i]
  } else {
    adj_gtex_beta[i] = 0
  }
}

gtex_beta_ukb_qtl$adj_gtex_beta = adj_gtex_beta

### comparing ukb cis-pqtls with DGN cis-eqtls
# UKB pQTL to DGN gene expression effect
gtex_beta_ukb_cis = gtex_beta_ukb_qtl[gtex_beta_ukb_qtl$ukb_cis_trans == 'cis', ]

# gtex eQTL to UKB protein effect
ukb_beta_gtex_cis = fread('pqtl/12_beta_across_two_data/ukb_beta_gtex_lcl_cis_all.txt')
colnames(ukb_beta_gtex_cis) = c("gtex_var", "gtex_gene_id", "gtex_tss_dist", "gtex_maf",
                                "gtex_p", "gtex_beta", "gtex_se", "gtex_chr",
                                "gtex_bp", "gtex_A0", "gtex_A1", "ukb_A0", "ukb_A1",
                                "ukb_A1_frq", "ukb_beta", "ukb_se", "ukb_logP")
# clump SNPs
uniq_gene = unique(ukb_beta_gtex_cis$gtex_gene_id)
ukb_beta_gtex_cis_loci = c()
window_size = 250000
for (i in uniq_gene) {
  ukb_beta_sub = ukb_beta_gtex_cis[ukb_beta_gtex_cis$gtex_gene_id == i, ]
  while(nrow(ukb_beta_sub > 0)) {
    lead_snp = which.min(ukb_beta_sub$gtex_p)[1]
    lead_chr = ukb_beta_sub$gtex_chr[lead_snp]
    lead_pos = ukb_beta_sub$gtex_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(ukb_beta_sub$gtex_chr == lead_chr & 
                       ukb_beta_sub$gtex_bp >= start_pos &
                       ukb_beta_sub$gtex_bp <= end_pos)
    ukb_beta_gtex_cis_loci = rbind(ukb_beta_gtex_cis_loci, ukb_beta_sub[lead_snp, ])
    ukb_beta_sub = ukb_beta_sub[-rm_index, ]
  }
}
ukb_beta_gtex_cis_loci = na.omit(ukb_beta_gtex_cis_loci)
adj_ukb_beta = c()
for (i in 1:nrow(ukb_beta_gtex_cis_loci)) {
  ukb_A0 = ukb_beta_gtex_cis_loci$ukb_A0[i]
  ukb_A1 = ukb_beta_gtex_cis_loci$ukb_A1[i]
  gtex_A0 = ukb_beta_gtex_cis_loci$gtex_A0[i]
  gtex_A1 = ukb_beta_gtex_cis_loci$gtex_A1[i]
  if (ukb_A0 == gtex_A1 & ukb_A1 == gtex_A0) {
    adj_ukb_beta[i] = -ukb_beta_gtex_cis_loci$ukb_beta[i]
  } else if (ukb_A0 == gtex_A0 & ukb_A1 == gtex_A1) {
    adj_ukb_beta[i] = ukb_beta_gtex_cis_loci$ukb_beta[i]
  } else {
    adj_ukb_beta[i] = 0
  }
}
ukb_beta_gtex_cis_loci$adj_ukb_beta = adj_ukb_beta

## cis-pQTL replicated as cis-eQTL
gtex_beta_ukb_cis = na.omit(gtex_beta_ukb_cis)
sum(gtex_beta_ukb_cis$gtex_p < 0.05/nrow(gtex_beta_ukb_cis)) / nrow(gtex_beta_ukb_cis)
length(unique(gtex_beta_ukb_cis$ukb_prot))
qval_gtex_cis = qvalue(gtex_beta_ukb_cis$gtex_p)
1 - qval_gtex_cis$pi0
# top 100 cis-pQTLs
gtex_beta_ukb_cis_top = gtex_beta_ukb_cis[order(gtex_beta_ukb_cis$ukb_logP, 
                                                decreasing = T)[1:100], ]
sum(gtex_beta_ukb_cis_top$gtex_p < 0.05/nrow(gtex_beta_ukb_cis_top)) / 
  nrow(gtex_beta_ukb_cis_top)

## cis-eQTL replicated as cis-pQTL
sum(ukb_beta_gtex_cis_loci$ukb_logP > -log10(0.05/nrow(ukb_beta_gtex_cis_loci))) /
  nrow(ukb_beta_gtex_cis_loci)
qval_ukb_cis = qvalue(10^-(ukb_beta_gtex_cis_loci$ukb_logP))
1 - qval_ukb_cis$pi0

#use more strict threshold (threshold in UKB pQTL) for gtex to make fair comparison
ukb_beta_gtex_cis_str = ukb_beta_gtex_cis_loci[ukb_beta_gtex_cis_loci$gtex_p < 5e-8/2922, ]
length(unique(ukb_beta_gtex_cis_str$gtex_gene_id))
sum(ukb_beta_gtex_cis_str$ukb_logP > -log10(0.05/nrow(ukb_beta_gtex_cis_str))) /
  nrow(ukb_beta_gtex_cis_str)
qval_ukb_str = qvalue(10^-(ukb_beta_gtex_cis_str$ukb_logP))
1 - qval_ukb_str$pi0

# top 100 cis-eQTLs
ukb_beta_gtex_cis_top = ukb_beta_gtex_cis_loci[order(ukb_beta_gtex_cis_loci$gtex_p)[1:100], ]
sum(ukb_beta_gtex_cis_top$ukb_logP > -log10(0.05/nrow(ukb_beta_gtex_cis_top))) /
  nrow(ukb_beta_gtex_cis_top)


library(ggplot2)
library(gridExtra)
# number of signals
n_signal_plot = data.frame(n_loci_gtex = c(755, 167, 91),
                           n_gene_gtex = c(647, 153, 87),
                           n_loci_ukb = c(925, 1035, 850),
                           n_gene_ukb = c(925, 1035, 850),
                           tissue = c('blood', 'liver', 'LCL'))
n_signal_plot$n_signal_gtex = paste0(n_signal_plot$n_loci_gtex, ' (', 
                                     n_signal_plot$n_gene_gtex, ')')
n_signal_plot$n_signal_ukb = paste0(n_signal_plot$n_loci_ukb, ' (', 
                                    n_signal_plot$n_gene_ukb, ')')

p1 = ggplot(n_signal_plot, aes(x=reorder(tissue, -n_loci_gtex), y=n_loci_gtex, 
                               label = n_signal_gtex)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue', width = 0.5) +
  labs(title = "# cis-eQTLs (# genes) analyzed", x = '', y = 'N') +
  geom_text(vjust=-0.2, size = 4.5) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2 = ggplot(n_signal_plot, aes(x=reorder(tissue, -n_loci_ukb), y=n_loci_ukb, 
                               label = n_signal_ukb)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue', width = 0.5) +
  labs(title = "# cis-pQTLs (# proteins) analyzed", x = '', y = 'N') +
  geom_text(vjust=-0.2, size = 4.5) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)

## replication and pi1
# GTEx cis-eQTL p < 5e-8/2922
gtex_rep_in_ukb = data.frame(rep = c(0.83, 0.66, NA, 0.77, 0.88, 0.71), 
                             tissue = rep(c('blood', 'liver', 'LCL'), each = 2),
                             group = rep(c('pi1', 'rep rate'), 3))
# All GTEx cis-eQTL
gtex_all_rep_in_ukb = data.frame(rep = c(0.59, 0.45, 0.71, 0.57, 0.59, 0.49), 
                                 tissue = rep(c('blood', 'liver', 'LCL'), each = 2),
                                 group = rep(c('pi1', 'rep rate'), 3))

# Top 100 GTEx cis-eQTLs
gtex_top_rep_in_ukb = data.frame(rep = c(0.79, 0.81, 0.73), 
                                 tissue = rep(c('blood', 'liver', 'LCL')))

# All UKB cis-pQTLs
ukb_rep_in_gtex = data.frame(rep = c(0.63, 0.29, 0.4, 0.12, 0.36, 0.08), 
                             tissue = rep(c('blood', 'liver', 'LCL'), each = 2),
                             group = rep(c('pi1', 'rep rate'), 3))
# Top 100 UKB cis-pQTLs
ukb_top_rep_in_gtex = data.frame(rep = c(0.4, 0.27, 0.18), 
                                 tissue = rep(c('blood', 'liver', 'LCL')))

p1 = ggplot(gtex_rep_in_ukb, aes(x=reorder(tissue, -rep), y=rep, label = round(rep, 2), fill = group)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  labs(title = "Replication of GTEx cis-eQTL in UKB", x = '', y = 'rate', fill = '') +
  geom_text(vjust=-0.2, size = 4.5, position = position_dodge(0.5)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2 = ggplot(ukb_rep_in_gtex, aes(x=reorder(tissue, -rep), y=rep, label = round(rep, 2), fill = group)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  labs(title = "Replication of UKB cis-pQTL in GTEx", x = '', y = 'rate', fill = '') +
  geom_text(vjust=-0.2, size = 4.5, position = position_dodge(0.5)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, nrow = 1)

p3 = ggplot(gtex_top_rep_in_ukb, aes(x=reorder(tissue, -rep), y=rep, label = round(rep, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5, fill = 'steelblue') +
  labs(title = "Replication of top 100 GTEx cis-eQTL in UKB", x = '', y = 'rate', fill = '') +
  geom_text(vjust=-0.2, size = 4.5, position = position_dodge(0.5)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p4 = ggplot(ukb_top_rep_in_gtex, aes(x=reorder(tissue, -rep), y=rep, label = round(rep, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5, fill = 'steelblue') +
  labs(title = "Replication of top 100 UKB cis-pQTL in GTEx", x = '', y = 'rate', fill = '') +
  geom_text(vjust=-0.2, size = 4.5, position = position_dodge(0.5)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p3, p4, nrow = 1)

# cis effect concordance
gtex_cis_ukb_rep = gtex_beta_ukb_cis[gtex_beta_ukb_cis$gtex_p < 
                                       0.05/nrow(gtex_beta_ukb_cis), ]
gtex_cis_ukb_rep = na.omit(gtex_cis_ukb_rep)
sum(sign(gtex_cis_ukb_rep$adj_gtex_beta) == sign(gtex_cis_ukb_rep$ukb_beta)) / 
  nrow(gtex_cis_ukb_rep)

ukb_cis_gtex_rep = ukb_beta_gtex_cis_loci[ukb_beta_gtex_cis_loci$ukb_logP > 
                                            -log10(0.05/nrow(ukb_beta_gtex_cis_loci)), ]
sum(sign(ukb_cis_gtex_rep$adj_ukb_beta) == sign(ukb_cis_gtex_rep$gtex_beta)) / 
  nrow(ukb_cis_gtex_rep)

ukb_cis_gtex_str_rep = ukb_beta_gtex_cis_str[ukb_beta_gtex_cis_str$ukb_logP > 
                                               -log10(0.05/nrow(ukb_beta_gtex_cis_loci)), ]
sum(sign(ukb_cis_gtex_str_rep$adj_ukb_beta) == sign(ukb_cis_gtex_str_rep$gtex_beta)) / 
  nrow(ukb_cis_gtex_str_rep)

### comparing ukb trans-pqtls with gtex trans-eqtls 
gtex_beta_ukb_trans = gtex_beta_ukb_qtl[gtex_beta_ukb_qtl$ukb_cis_trans == 'trans', ]
ukb_beta_gtex_trans = fread('pqtl/12_beta_across_two_data/ukb_beta_gtex_blood_trans_all.txt')
colnames(ukb_beta_gtex_trans) = c("gtex_var", "gtex_gene", "gtex_p", "gtex_beta",
                                  "gtex_se", "gtex_af", "gtex_chr", "gtex_bp", 
                                  "ukb_A0", "ukb_A1", "ukb_A1_frq", 
                                  "ukb_beta", "ukb_se", "ukb_logP")
gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.gtf')
gene_meta = gene_meta[, c("chr", "start", "end", "gene_id")]
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
ukb_beta_gtex_trans$gene_chr = gene_meta$chr[match(ukb_beta_gtex_trans$gtex_gene, gene_meta$gene_id)]
ukb_beta_gtex_trans$gene_chr = as.numeric(sub('chr', '', ukb_beta_gtex_trans$gene_chr))
ukb_beta_gtex_trans$gene_tss = gene_meta$start[match(ukb_beta_gtex_trans$gtex_gene, gene_meta$gene_id)]
ukb_keep = (ukb_beta_gtex_trans$gene_chr != ukb_beta_gtex_trans$gtex_chr |
              (ukb_beta_gtex_trans$gene_chr == ukb_beta_gtex_trans$gtex_chr &
                 abs(ukb_beta_gtex_trans$gtex_bp - ukb_beta_gtex_trans$gene_tss) > 5000000))
ukb_beta_gtex_trans = ukb_beta_gtex_trans[ukb_keep, ]
ukb_beta_gtex_trans$gtex_A0 = sapply(strsplit(ukb_beta_gtex_trans$gtex_var, '_', fixed = T), 
                                     '[', 3)
ukb_beta_gtex_trans$gtex_A1 = sapply(strsplit(ukb_beta_gtex_trans$gtex_var, '_', fixed = T), 
                                     '[', 4)
ukb_beta_gtex_trans$ukb_z = ukb_beta_gtex_trans$ukb_beta / ukb_beta_gtex_trans$ukb_se

# clump SNPs
uniq_gene = unique(ukb_beta_gtex_trans$gtex_gene)
ukb_beta_gtex_trans_loci = c()
window_size = 250000
for (i in uniq_gene) {
  ukb_beta_sub = ukb_beta_gtex_trans[ukb_beta_gtex_trans$gtex_gene == i, ]
  while(nrow(ukb_beta_sub > 0)) {
    lead_snp = which.min(ukb_beta_sub$gtex_p)[1]
    lead_chr = ukb_beta_sub$gtex_chr[lead_snp]
    lead_pos = ukb_beta_sub$gtex_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(ukb_beta_sub$gtex_chr == lead_chr & 
                       ukb_beta_sub$gtex_bp >= start_pos &
                       ukb_beta_sub$gtex_bp <= end_pos)
    ukb_beta_gtex_trans_loci = rbind(ukb_beta_gtex_trans_loci, ukb_beta_sub[lead_snp, ])
    ukb_beta_sub = ukb_beta_sub[-rm_index, ]
  }
}
ukb_beta_gtex_trans_loci = na.omit(ukb_beta_gtex_trans_loci)

adj_ukb_z = c()
for (i in 1:nrow(ukb_beta_gtex_trans_loci)) {
  ukb_A0 = ukb_beta_gtex_trans_loci$ukb_A0[i]
  ukb_A1 = ukb_beta_gtex_trans_loci$ukb_A1[i]
  gtex_A0 = ukb_beta_gtex_trans_loci$gtex_A0[i]
  gtex_A1 = ukb_beta_gtex_trans_loci$gtex_A1[i]
  if (ukb_A0 == gtex_A1 & ukb_A1 == gtex_A0) {
    adj_ukb_z[i] = -ukb_beta_gtex_trans_loci$ukb_z[i]
  } else if (ukb_A0 == gtex_A0 & ukb_A1 == gtex_A1) {
    adj_ukb_z[i] = ukb_beta_gtex_trans_loci$ukb_z[i]
  } else {
    adj_ukb_z[i] = 0
  }
}
ukb_beta_gtex_trans_loci$adj_ukb_z = adj_ukb_z

# trans-pQTL replicated as trans-eQTL
sum(gtex_beta_ukb_trans$gtex_p < 0.05/nrow(gtex_beta_ukb_trans)) / nrow(gtex_beta_ukb_trans)
qval_gtex_trans = qvalue(gtex_beta_ukb_trans$gtex_p)
1 - qval_gtex_trans$pi0
# top pQTLs (the number as gtex eQTL investigated in UKB)
gtex_beta_ukb_trans_top = gtex_beta_ukb_trans[order(gtex_beta_ukb_trans$ukb_logP, decreasing = T)[1:100], ]
sum(gtex_beta_ukb_trans_top$gtex_p < 0.05/nrow(gtex_beta_ukb_trans_top)) / nrow(gtex_beta_ukb_trans_top)
qval_gtex_trans_top = qvalue(gtex_beta_ukb_trans_top$gtex_p)
1 - qval_gtex_trans_top$pi0

# trans-eQTL replicated as trans-pQTL
sum(ukb_beta_gtex_trans_loci$ukb_logP > -log10(0.05/nrow(ukb_beta_gtex_trans_loci))) / 
  nrow(ukb_beta_gtex_trans_loci)
qval_ukb_trans = qvalue(10^-(ukb_beta_gtex_trans_loci$ukb_logP))
1 - qval_ukb_trans$pi0

#use more strict threshold (threshold in UKB pQTL) for gtex to make fair comparison
ukb_beta_gtex_trans_str = ukb_beta_gtex_trans_loci[ukb_beta_gtex_trans_loci$gtex_p < 5e-8/2922, ]
sum(ukb_beta_gtex_trans_str$ukb_logP > -log10(0.05/nrow(ukb_beta_gtex_trans_str))) /
  nrow(ukb_beta_gtex_trans_str)
qval_ukb_trans_str = qvalue(10^-(ukb_beta_gtex_trans_str$ukb_logP))
1 - qval_ukb_trans_str$pi0

## replication and pi1
# All UKB trans-pQTLs
ukb_rep_in_gtex = data.frame(rep = c(0.039, 0.00017, 0, 0, 0, 0), 
                             tissue = rep(c('blood', 'liver', 'LCL'), each = 2),
                             group = rep(c('pi1', 'rep rate'), 3))
# Top 100 UKB trans-pQTLs
ukb_top_rep_in_gtex = data.frame(rep = c(0, 0.01, 0.32, 0, 0.29, 0), 
                                 tissue = rep(c('blood', 'liver', 'LCL'), each = 2),
                                 group = rep(c('pi1', 'rep rate'), 3))

ggplot(ukb_top_rep_in_gtex[ukb_top_rep_in_gtex$group == 'pi1', ], 
       aes(x=reorder(tissue, -rep), y=rep, label = round(rep, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5, fill = 'steelblue') +
  labs(title = "pi1 of top 100 UKB trans-pQTLs in GTEx", x = '', y = 'rate') +
  geom_text(vjust=-0.2, size = 4.5, position = position_dodge(0.5)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

all_ukb_gtex_blood_p = gtex_beta_ukb_trans$gtex_p
all_ukb_gtex_liver_p = gtex_beta_ukb_trans$gtex_p
all_ukb_gtex_lcl_p = gtex_beta_ukb_trans$gtex_p

p_sort_blood = sort(all_ukb_gtex_blood_p)
p_sort_liver = sort(all_ukb_gtex_liver_p)
p_sort_lcl = sort(all_ukb_gtex_lcl_p)
p_sort = data.frame(p_exp = c(-log10(ppoints(p_sort_blood)), 
                              -log10(ppoints(p_sort_liver)), 
                              -log10(ppoints(p_sort_lcl))),
                    p_obs = c(-log10(p_sort_blood), 
                              -log10(p_sort_liver), 
                              -log10(p_sort_lcl)),
                    tissue = c(rep('blood', length(p_sort_blood)),
                               rep('liver', length(p_sort_liver)),
                               rep('lcl', length(p_sort_lcl))))

ggplot(p_sort, aes(x = p_exp, y = p_obs, color = tissue)) + 
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(x = "expected", y = "observed", title = 'qqplot of all UKB pQTL in GTEx') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())





