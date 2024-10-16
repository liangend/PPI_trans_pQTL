library(data.table)
library(qvalue)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')
### Reading data
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
dgn_beta_ukb_qtl = na.omit(dgn_beta_ukb_qtl)
dgn_beta_ukb_qtl = dgn_beta_ukb_qtl[dgn_beta_ukb_qtl$adj_dgn_beta != 0, ]

gtex_blood_beta_ukb_qtl = fread('pqtl/12_beta_across_two_data/gtex_blood_beta_ukb_pqtl_all.txt')
gtex_blood_beta_ukb_qtl = gtex_blood_beta_ukb_qtl[, -c(1,4,6,7)]
adj_gtex_beta = c()
for (i in 1:nrow(gtex_blood_beta_ukb_qtl)) {
  ukb_A0 = gtex_blood_beta_ukb_qtl$ukb_A0[i]
  ukb_A1 = gtex_blood_beta_ukb_qtl$ukb_A1[i]
  gtex_A0 = gtex_blood_beta_ukb_qtl$gtex_A0[i]
  gtex_A1 = gtex_blood_beta_ukb_qtl$gtex_A1[i]
  if (ukb_A0 == gtex_A1 & ukb_A1 == gtex_A0) {
    adj_gtex_beta[i] = gtex_blood_beta_ukb_qtl$gtex_beta[i]
  } else if (ukb_A0 == gtex_A0 & ukb_A1 == gtex_A1) {
    adj_gtex_beta[i] = -gtex_blood_beta_ukb_qtl$gtex_beta[i]
  } else {
    adj_gtex_beta[i] = 0
  }
}

gtex_blood_beta_ukb_qtl$adj_gtex_beta = adj_gtex_beta
gtex_blood_beta_ukb_qtl = na.omit(gtex_blood_beta_ukb_qtl)
gtex_blood_beta_ukb_qtl = gtex_blood_beta_ukb_qtl[gtex_blood_beta_ukb_qtl$adj_gtex_beta != 0, ]


gtex_liver_beta_ukb_qtl = fread('pqtl/12_beta_across_two_data/gtex_liver_beta_ukb_pqtl_all.txt')
gtex_liver_beta_ukb_qtl = gtex_liver_beta_ukb_qtl[, -c(1,4,6,7)]
adj_gtex_beta = c()
for (i in 1:nrow(gtex_liver_beta_ukb_qtl)) {
  ukb_A0 = gtex_liver_beta_ukb_qtl$ukb_A0[i]
  ukb_A1 = gtex_liver_beta_ukb_qtl$ukb_A1[i]
  gtex_A0 = gtex_liver_beta_ukb_qtl$gtex_A0[i]
  gtex_A1 = gtex_liver_beta_ukb_qtl$gtex_A1[i]
  if (ukb_A0 == gtex_A1 & ukb_A1 == gtex_A0) {
    adj_gtex_beta[i] = gtex_liver_beta_ukb_qtl$gtex_beta[i]
  } else if (ukb_A0 == gtex_A0 & ukb_A1 == gtex_A1) {
    adj_gtex_beta[i] = -gtex_liver_beta_ukb_qtl$gtex_beta[i]
  } else {
    adj_gtex_beta[i] = 0
  }
}

gtex_liver_beta_ukb_qtl$adj_gtex_beta = adj_gtex_beta
gtex_liver_beta_ukb_qtl = na.omit(gtex_liver_beta_ukb_qtl)
gtex_liver_beta_ukb_qtl = gtex_liver_beta_ukb_qtl[gtex_liver_beta_ukb_qtl$adj_gtex_beta != 0, ]

gtex_lcl_beta_ukb_qtl = fread('pqtl/12_beta_across_two_data/gtex_lcl_beta_ukb_pqtl_all.txt')
gtex_lcl_beta_ukb_qtl = gtex_lcl_beta_ukb_qtl[, -c(1,4,6,7)]
adj_gtex_beta = c()
for (i in 1:nrow(gtex_lcl_beta_ukb_qtl)) {
  ukb_A0 = gtex_lcl_beta_ukb_qtl$ukb_A0[i]
  ukb_A1 = gtex_lcl_beta_ukb_qtl$ukb_A1[i]
  gtex_A0 = gtex_lcl_beta_ukb_qtl$gtex_A0[i]
  gtex_A1 = gtex_lcl_beta_ukb_qtl$gtex_A1[i]
  if (ukb_A0 == gtex_A1 & ukb_A1 == gtex_A0) {
    adj_gtex_beta[i] = gtex_lcl_beta_ukb_qtl$gtex_beta[i]
  } else if (ukb_A0 == gtex_A0 & ukb_A1 == gtex_A1) {
    adj_gtex_beta[i] = -gtex_lcl_beta_ukb_qtl$gtex_beta[i]
  } else {
    adj_gtex_beta[i] = 0
  }
}

gtex_lcl_beta_ukb_qtl$adj_gtex_beta = adj_gtex_beta
gtex_lcl_beta_ukb_qtl = na.omit(gtex_lcl_beta_ukb_qtl)
gtex_lcl_beta_ukb_qtl = gtex_lcl_beta_ukb_qtl[gtex_lcl_beta_ukb_qtl$adj_gtex_beta != 0, ]


## cis- and trans-pQTL replication
dgn_beta_ukb_cis = dgn_beta_ukb_qtl[dgn_beta_ukb_qtl$ukb_cis_trans == 'cis', ]
gtex_blood_beta_ukb_cis = gtex_blood_beta_ukb_qtl[gtex_blood_beta_ukb_qtl$ukb_cis_trans == 'cis', ]
gtex_liver_beta_ukb_cis = gtex_liver_beta_ukb_qtl[gtex_liver_beta_ukb_qtl$ukb_cis_trans == 'cis', ]
gtex_lcl_beta_ukb_cis = gtex_lcl_beta_ukb_qtl[gtex_lcl_beta_ukb_qtl$ukb_cis_trans == 'cis', ]

dgn_beta_ukb_trans = dgn_beta_ukb_qtl[dgn_beta_ukb_qtl$ukb_cis_trans == 'trans', ]
gtex_blood_beta_ukb_trans = gtex_blood_beta_ukb_qtl[gtex_blood_beta_ukb_qtl$ukb_cis_trans == 'trans', ]
gtex_liver_beta_ukb_trans = gtex_liver_beta_ukb_qtl[gtex_liver_beta_ukb_qtl$ukb_cis_trans == 'trans', ]
gtex_lcl_beta_ukb_trans = gtex_lcl_beta_ukb_qtl[gtex_lcl_beta_ukb_qtl$ukb_cis_trans == 'trans', ]

p_sort_dgn_cis = sort(dgn_beta_ukb_cis$dgn_p)
p_sort_blood_cis = sort(gtex_blood_beta_ukb_cis$gtex_p)
p_sort_liver_cis = sort(gtex_liver_beta_ukb_cis$gtex_p)
p_sort_lcl_cis = sort(gtex_lcl_beta_ukb_cis$gtex_p)

p_sort_dgn_trans = sort(dgn_beta_ukb_trans$dgn_p)
p_sort_blood_trans = sort(gtex_blood_beta_ukb_trans$gtex_p)
p_sort_liver_trans = sort(gtex_liver_beta_ukb_trans$gtex_p)
p_sort_lcl_trans = sort(gtex_lcl_beta_ukb_trans$gtex_p)

## pi1 values
pi1_table = data.frame(pi1 = c(1 - qvalue(p_sort_dgn_cis)$pi0,
                               1 - qvalue(p_sort_blood_cis)$pi0,
                               1 - qvalue(p_sort_liver_cis)$pi0,
                               1 - qvalue(p_sort_lcl_cis)$pi0,
                               1 - qvalue(p_sort_dgn_trans)$pi0,
                               1 - qvalue(p_sort_blood_trans)$pi0,
                               1 - qvalue(p_sort_liver_trans)$pi0,
                               1 - qvalue(p_sort_lcl_trans)$pi0),
                       cis_trans = rep(c('cis', 'trans'), each = 4),
                       tissue = rep(c('dgn', 'blood', 'liver', 'lcl'), 2))
ggplot(pi1_table[pi1_table$tissue != 'dgn', ], 
       aes(x=reorder(cis_trans, -pi1), y=pi1, label = round(pi1, 2), fill = reorder(tissue, -pi1))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  labs(title = "", x = '', y = 'pi1', fill = '') +
  scale_fill_brewer(palette="Set1") +
  geom_text(vjust=-0.2, size = 4.5, position=position_dodge(0.5)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## qq plot
p_sort_cis = data.frame(p_exp = c(-log10(ppoints(p_sort_dgn_cis)),
                                  -log10(ppoints(p_sort_blood_cis)), 
                                    -log10(ppoints(p_sort_liver_cis)), 
                                    -log10(ppoints(p_sort_lcl_cis))),
                          p_obs = c(-log10(p_sort_dgn_cis),
                                    -log10(p_sort_blood_cis), 
                                    -log10(p_sort_liver_cis), 
                                    -log10(p_sort_lcl_cis)),
                          tissue = c(rep('dgn', length(p_sort_dgn_cis)),
                                     rep('blood', length(p_sort_blood_cis)),
                                     rep('liver', length(p_sort_liver_cis)),
                                     rep('lcl', length(p_sort_lcl_cis))))

p1 = ggplot(p_sort_cis, aes(x = p_exp, y = p_obs, color = tissue)) + 
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(x = "expected", y = "observed", title = 'qqplot of all UKB cis-pQTLs') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p_sort_trans = data.frame(p_exp = c(-log10(ppoints(p_sort_dgn_trans)), 
                                    -log10(ppoints(p_sort_blood_trans)), 
                              -log10(ppoints(p_sort_liver_trans)), 
                              -log10(ppoints(p_sort_lcl_trans))),
                    p_obs = c(-log10(p_sort_dgn_trans), 
                              -log10(p_sort_blood_trans), 
                              -log10(p_sort_liver_trans), 
                              -log10(p_sort_lcl_trans)),
                    tissue = c(rep('dgn', length(p_sort_dgn_trans)),
                               rep('blood', length(p_sort_blood_trans)),
                               rep('liver', length(p_sort_liver_trans)),
                               rep('lcl', length(p_sort_lcl_trans))))

p2 = ggplot(p_sort_trans, aes(x = p_exp, y = p_obs, color = tissue)) + 
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(x = "expected", y = "observed", title = 'qqplot of all UKB trans-pQTLs') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, nrow = 1)


## number of genes and qtls in each dataset
n_signal_table = data.frame(n_loci = c(length(p_sort_dgn_cis),
                                       length(p_sort_blood_cis),
                                       length(p_sort_liver_cis),
                                       length(p_sort_lcl_cis),
                                       length(p_sort_dgn_trans),
                                       length(p_sort_blood_trans),
                                       length(p_sort_liver_trans),
                                       length(p_sort_lcl_trans)),
                            n_gene = c(length(unique(dgn_beta_ukb_cis$ukb_gene)),
                                       length(unique(gtex_blood_beta_ukb_cis$ukb_prot)),
                                       length(unique(gtex_liver_beta_ukb_cis$ukb_prot)),
                                       length(unique(gtex_lcl_beta_ukb_cis$ukb_prot)),
                                       length(unique(dgn_beta_ukb_trans$ukb_gene)),
                                       length(unique(gtex_blood_beta_ukb_trans$ukb_prot)),
                                       length(unique(gtex_liver_beta_ukb_trans$ukb_prot)),
                                       length(unique(gtex_lcl_beta_ukb_trans$ukb_prot))),
                       cis_trans = rep(c('cis', 'trans'), each = 4),
                       tissue = rep(c('dgn', 'blood', 'liver', 'lcl'), 2))
n_signal_table$n_signal = paste0(n_signal_table$n_loci, ' (', n_signal_table$n_gene, ')')
ggplot(n_signal_table, aes(x=cis_trans, y=n_loci, 
                           label = n_signal, fill = reorder(tissue, -n_loci))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  labs(title = "# of UKB pQTLs (# genes) analyzed in each dataset", x = '', y = 'N', fill = '') +
  geom_text(vjust=-0.2, size = 4, position=position_dodge(0.5)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



