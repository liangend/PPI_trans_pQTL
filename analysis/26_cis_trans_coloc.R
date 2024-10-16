library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(boot)
setwd('/project/xuanyao/jinghui')
cis_trans_all = fread('pqtl/14_cis_trans_coloc/ukb_trans_coloc_cis.txt')
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
cis_trans_all$cis_gene = ukb_prot$gene_name[match(cis_trans_all$cis_prot, ukb_prot$file)]
cis_trans_all$cis_trans_pair = paste0(cis_trans_all$cis_gene, ',', cis_trans_all$target_gene)
cis_trans_all$cis_pval = pnorm(abs(cis_trans_all$cis_beta/cis_trans_all$cis_se), lower.tail = F) * 2
### protein-protein interaction and complex
ppi_list = readRDS('pqtl/11_prot_complex/ppi_list.rds')

cis_trans_all$is_ppi = (cis_trans_all$cis_trans_pair %in% ppi_list)


cis_trans_coloc = cis_trans_all[cis_trans_all$pp4 > 0.75, ]
cis_trans_most_coloc = c()
uniq_trans_loci = unique(cis_trans_coloc$loci)
for (i in uniq_trans_loci) {
  cis_trans_coloc_i = cis_trans_coloc[cis_trans_coloc$loci == i, ]
  cis_trans_most_coloc = rbind(cis_trans_most_coloc, cis_trans_coloc_i[which.max(cis_trans_coloc_i$pp4), ])
}

cis_trans_ppi = cis_trans_all[cis_trans_all$is_ppi, ]
# ppi cor
# mcl_prot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')
# cis_trans_cor = c()
# for (i in 1:nrow(cis_trans_ppi)) {
#   prot1 = cis_trans_ppi$target_gene[i]
#   prot2 = cis_trans_ppi$cis_gene[i]
#   mod1 = mcl_prot$mod[mcl_prot$prot == prot1]
#   mod2 = mcl_prot$mod[mcl_prot$prot == prot2]
#   mod_i = intersect(mod1, mod2)[1]
#   if (is.na(mod_i)) {
#     next
#   }
#   mod_cor_i = readRDS(paste0('pqtl/02_sigma/ukb/mod', mod_i, '.rds'))
#   cis_trans_cor[i] = mod_cor_i[prot1, prot2]
# }
# cis_trans_ppi$cis_trans_cor = cis_trans_cor
# kkk = na.omit(cis_trans_ppi)
# mean(sign(kkk$cis_beta) * sign(kkk$trans_beta) == sign(kkk$cis_trans_cor))

cis_trans_coloc_ppi = cis_trans_coloc[cis_trans_coloc$is_ppi, ]
cis_trans_most_coloc_ppi = cis_trans_most_coloc[cis_trans_most_coloc$is_ppi, ]

dir_table = data.frame(all_prop = c(mean(sign(cis_trans_all$trans_beta) == sign(cis_trans_all$cis_beta)),
                                    mean(sign(cis_trans_ppi$trans_beta) == sign(cis_trans_ppi$cis_beta)),
                                    mean(sign(cis_trans_coloc$trans_beta) == sign(cis_trans_coloc$cis_beta)),
                                    mean(sign(cis_trans_most_coloc$trans_beta) == sign(cis_trans_most_coloc$cis_beta)),
                                    mean(sign(cis_trans_coloc_ppi$trans_beta) == sign(cis_trans_coloc_ppi$cis_beta)),
                                    mean(sign(cis_trans_most_coloc_ppi$trans_beta) == 
                                           sign(cis_trans_most_coloc_ppi$cis_beta))),
                       pos_prop = c(mean(cis_trans_all[cis_trans_all$cis_beta > 0, 'trans_beta'] > 0),
                                    mean(cis_trans_ppi[cis_trans_ppi$cis_beta > 0, 'trans_beta'] > 0),
                                    mean(cis_trans_coloc[cis_trans_coloc$cis_beta > 0, 'trans_beta'] > 0),
                                    mean(cis_trans_most_coloc[cis_trans_most_coloc$cis_beta > 0, 'trans_beta'] > 0),
                                    mean(cis_trans_coloc_ppi[cis_trans_coloc_ppi$cis_beta > 0, 'trans_beta'] > 0),
                                    mean(cis_trans_most_coloc_ppi[cis_trans_most_coloc_ppi$cis_beta > 0, 
                                                                  'trans_beta'] > 0)),
                       neg_prop = c(mean(cis_trans_all[cis_trans_all$cis_beta < 0, 'trans_beta'] < 0),
                                    mean(cis_trans_ppi[cis_trans_ppi$cis_beta < 0, 'trans_beta'] < 0),
                                    mean(cis_trans_coloc[cis_trans_coloc$cis_beta < 0, 'trans_beta'] < 0),
                                    mean(cis_trans_most_coloc[cis_trans_most_coloc$cis_beta < 0, 'trans_beta'] < 0),
                                    mean(cis_trans_coloc_ppi[cis_trans_coloc_ppi$cis_beta < 0, 'trans_beta'] < 0),
                                    mean(cis_trans_most_coloc_ppi[cis_trans_most_coloc_ppi$cis_beta < 0, 
                                                                  'trans_beta'] < 0)),
                       all_n = c(nrow(cis_trans_all), nrow(cis_trans_ppi), nrow(cis_trans_coloc),
                                 nrow(cis_trans_most_coloc), nrow(cis_trans_coloc_ppi), nrow(cis_trans_most_coloc_ppi)),
                       pos_n = c(sum(cis_trans_all$cis_beta > 0), sum(cis_trans_ppi$cis_beta > 0), 
                                 sum(cis_trans_coloc$cis_beta > 0), sum(cis_trans_most_coloc$cis_beta > 0),
                                 sum(cis_trans_coloc_ppi$cis_beta > 0), sum(cis_trans_most_coloc_ppi$cis_beta > 0)),
                       neg_n = c(sum(cis_trans_all$cis_beta < 0), sum(cis_trans_ppi$cis_beta < 0), 
                                 sum(cis_trans_coloc$cis_beta < 0), sum(cis_trans_most_coloc$cis_beta < 0),
                                 sum(cis_trans_coloc_ppi$cis_beta < 0), sum(cis_trans_most_coloc_ppi$cis_beta < 0)),
                       group = c('all', 'PPI', 'coloc', 'coloc most sig', 'coloc & PPI', 'coloc most sig & PPI'))

qtl_prop_plot = reshape(dir_table, direction = "long", idvar = c('group'), 
                        varying = list(1:3, 4:6), v.names = c('prop', 'N'),
                        timevar = 'cat', times = c('all', 'pos', 'neg'))
qtl_prop_plot$group = factor(qtl_prop_plot$group, levels = c('all', 'PPI', 'coloc', 'coloc & PPI',
                                                             'coloc most sig', 'coloc most sig & PPI'))
qtl_prop_plot = rbind(qtl_prop_plot, 
                      data.frame(group = rep('coloc & CORUM', 3),
                                 cat = c('all', 'pos', 'neg'),
                                 prop = c(0.8095238, 0.9285714, 0.5714286),
                                 N = c(21, 14, 7)))
ggplot(qtl_prop_plot, aes(x = group, y = prop, fill = cat, label = N)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = 'Proportion', fill = 'cis direction',
       title = 'Proportion of consistent direction of trans and cis pQTL') +
  geom_text(position=position_dodge(0.9), vjust=-0.2) + 
  scale_fill_manual(values = brewer.pal(3, "Set1")) +
  coord_cartesian(ylim = c(0.25, 1)) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

qtl_prop_plot2 = qtl_prop_plot[qtl_prop_plot$cat == 'all', ]
qtl_prop_plot2 = qtl_prop_plot2[qtl_prop_plot2$group %in% c('coloc', 'coloc & PPI'), ]
ggplot(qtl_prop_plot2, aes(x = group, y = prop)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue', width = 0.5) +
  labs(x = "", y = 'Proportion',
       title = 'Consistent direction of cis- and trans-pQTLs') +
  scale_fill_manual(values = brewer.pal(3, "Set1")) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#### pco pleiotropic loci effect direction
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
pco_pleio = pco_all[pco_all$pleio, ]
univar_z = sapply(strsplit(pco_pleio$univar_z, ',', fixed = T), as.numeric)
n_prot = sapply(strsplit(pco_pleio$univar_trans_prot, ',', fixed = T), length)
# some univar_z have more z values than number of proteins due to different SNPs at the same genomic position
for (i in 1:length(n_prot)) {
  univar_z[[i]] = univar_z[[i]][1:n_prot[i]]
}
univar_p = sapply(univar_z, function(x){pnorm(abs(x), lower.tail = F)*2})
univar_sig_z = sapply(1:length(univar_z), 
                      function(x){univar_z[[x]][which(univar_p[[x]] < 0.05/length(univar_p[[x]]))]})
sig_n = sapply(univar_sig_z, length)
sign_z = sapply(univar_sig_z, function(x){unique(sign(x))})
sign_length = sapply(sign_z, length)
sig_plot = as.data.frame(table(sig_n))
colnames(sig_plot) = c('n_prot', 'n_qtl')
concor_dat = as.data.frame(table(sig_n[which(sign_length == 1)]))
sig_plot$n_concor_qtl = concor_dat$Freq[match(sig_plot$n_prot, concor_dat$Var1)]
sig_plot[is.na(sig_plot)] = 0 
sig_plot$concor_prop = sig_plot$n_concor_qtl/sig_plot$n_qtl
sig_plot$expect_same_dir = 2/2^(c(2:22,26,27))
# proportion of all sig univar z having the same effect direction
ggplot(sig_plot, aes(x = n_prot, y = concor_prop, label = n_qtl)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue', width = 0.5) +
  geom_text(position=position_dodge(0.9), vjust=-0.4, size = 3) + 
  geom_point(aes(x = n_prot, y = expect_same_dir), color = 'red') + 
  labs(x = "PPI size", y = 'Proportion of concordance',
       title = '') +
  scale_fill_manual(values = brewer.pal(3, "Set1")) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# random effect direction
rand_dat = data.frame(n_prot = 2:10)
for (i in 1:100) {
  rand_dir = sapply(sig_n, function(x){length(unique(sample(c(-1,1), size = x, replace = T)))})
  rand_ppi_concor = as.data.frame(table(sig_n[which(rand_dir == 1)]))
  rand_dat = cbind(rand_dat, n = rand_ppi_concor$Freq[match(rand_dat$n_prot, rand_ppi_concor$Var1)])
  rand_concor[i] = table(rand_dir)[1] / length(rand_dir)
}
table(sign_length)[1] / length(sign_length) / mean(rand_concor)

sig_enrich = apply(rand_dat[,-1], 2, function(x){sig_vs_rand$sig_concor/x})
sig_enrich[sig_enrich > 100] = NA
sig_vs_rand = data.frame(n_prot = 2:10)
sig_vs_rand$enrich = apply(sig_enrich, 1, function(x){mean(x, na.rm = T)})
sig_vs_rand$enrich_lower = apply(sig_enrich, 1, function(x){quantile(x, 0.025, na.rm = T)})
sig_vs_rand$enrich_upper = apply(sig_enrich, 1, function(x){quantile(x, 0.975, na.rm = T)})

ggplot(sig_vs_rand[1:6, ], aes(x = n_prot, y = enrich)) + 
  geom_point(stat="identity", size = 3, color = 'black') +
  geom_errorbar(aes(ymin=enrich_lower, ymax=enrich_upper), width=.08) +
  geom_hline(yintercept = 1, lty = 2) + 
  labs(title = '', x = "PPI size", y = "Concordance enrichment", fill = '') +
  theme(text = element_text(size=14, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# example of concordant effect
pco_ex = pco_pleio[2307, ]
plot_ex = data.frame(prot = unlist(strsplit(pco_ex$univar_trans_prot, ',', fixed = T)),
                     univar_z = as.numeric(unlist(strsplit(pco_ex$univar_z, ',', fixed = T))))
plot_ex$univar_p = pnorm(abs(plot_ex$univar_z), lower.tail = F)*2
plot_ex = plot_ex[plot_ex$univar_p < 0.05/nrow(plot_ex), ]
ggplot(plot_ex, aes(x = prot, y = univar_z)) + 
  geom_point(color = 'steelblue', size = 4) +
  labs(x = "Proteins affected in PPI", y = 'Univariate z score',
       title = '') +
  ylim(c(-22,4)) + geom_abline(slope = 0, intercept = 0, linetype="dashed") +
  scale_fill_manual(values = brewer.pal(3, "Set1")) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

concor_prop = sapply(univar_sig_z, function(x){max(table(sign(x))/length(x))})
concor_qtl = data.frame(n_prot = sig_n, concor_prop = concor_prop)
concor_qtl = concor_qtl[concor_qtl$n_prot > 3, ]
concor_qtl$n_prot = as.factor(concor_qtl$n_prot)
# proportion of sig univar z having the same effect direction in each PPI trans-pQTL
ggplot(concor_qtl, aes(x = n_prot, y = concor_prop)) + 
  geom_boxplot(fill = 'steelblue', width = 0.5) +
  labs(x = "n sig univar z", y = 'Proportion',
       title = 'Consistent direction of univar prot') +
  scale_fill_manual(values = brewer.pal(3, "Set1")) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


## trans-pQTLs that are not trans-eQTLs but coloc with cis-pQTLs
dgn_beta_ukb_sig = fread('pqtl/12_beta_across_two_data/dgn_beta_ukb_combined_sig.txt')
dgn_beta_ukb_sig = na.omit(dgn_beta_ukb_sig)
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
dgn_beta_ukb_sig$gene_chr = ukb_prot$chr[match(dgn_beta_ukb_sig$file, 
                                               ukb_prot$file)]
dgn_beta_ukb_sig$gene_tss = ukb_prot$gene_start[match(dgn_beta_ukb_sig$file, 
                                                      ukb_prot$file)]
dgn_beta_ukb_sig$gene_name = ukb_prot$gene_name[match(dgn_beta_ukb_sig$file, 
                                                      ukb_prot$file)]
dgn_beta_ukb_sig$cis_trans = 'trans'
dgn_beta_ukb_sig$cis_trans[which(dgn_beta_ukb_sig$CHROM == dgn_beta_ukb_sig$gene_chr &
                                   abs(dgn_beta_ukb_sig$bp_hg38 - dgn_beta_ukb_sig$gene_tss) < 1000000)] = 'cis'
dgn_beta_ukb_trans = dgn_beta_ukb_sig[dgn_beta_ukb_sig$cis_trans == 'trans', ]
cis_trans_coloc$snp_prot = paste0(cis_trans_coloc$ID, ':', cis_trans_coloc$target_gene)
dgn_beta_ukb_trans$snp_prot = paste0(dgn_beta_ukb_trans$ID, ':', dgn_beta_ukb_trans$gene_name)
cis_trans_coloc$dgn_p = dgn_beta_ukb_trans$dgn_p[match(cis_trans_coloc$snp_prot, dgn_beta_ukb_trans$snp_prot)]
cis_trans_coloc_dgn = na.omit(cis_trans_coloc)

cis_trans_coloc_dgn$dgn_p_interval = cut(cis_trans_coloc_dgn$dgn_p, c(0, 0.01, 0.05, 0.1, 1), 
                                         include.lowest = T)
cis_trans_coloc_dgn$sign_cons = sign(cis_trans_coloc_dgn$cis_beta) * sign(cis_trans_coloc_dgn$trans_beta)
cis_trans_coloc_dgn$sign_cons[cis_trans_coloc_dgn$sign_cons == -1] = 0

boot_prop = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, mean)
  return(aggre_mean[,2])
}

cons_prop = boot(cis_trans_coloc_dgn, boot_prop, 1000, formula = sign_cons~dgn_p_interval)
cons_prop_plot = data.frame(prop = apply(cons_prop$t, 2, mean),
                      quant_lower = apply(cons_prop$t, 2, function(x){quantile(x, 0.025)}), 
                      quant_higher = apply(cons_prop$t, 2, function(x){quantile(x, 0.975)}),
                      DGN_p = c('0 ~ 0.01', '0.01 ~ 0.05', '0.05 ~ 0.1', '> 0.1'))
ggplot(cons_prop_plot, aes(x = factor(DGN_p, levels = c('0 ~ 0.01', '0.01 ~ 0.05', '0.05 ~ 0.1', '> 0.1')), 
                     y = prop)) + 
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  labs(x = "DGN p", y = "Proportion of consistent cis-trans", color = '') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# for each trans-pqtl, select one cis-pqtl with the largest pp4
cis_trans_most_coloc_dgn = c()
uniq_trans_loci = unique(cis_trans_coloc_dgn$loci)
for (i in uniq_trans_loci) {
  cis_trans_coloc_i = cis_trans_coloc_dgn[cis_trans_coloc_dgn$loci == i, ]
  cis_trans_most_coloc_dgn = rbind(cis_trans_most_coloc_dgn, cis_trans_coloc_i[which.max(cis_trans_coloc_i$pp4), ])
}

cons_prop_most_coloc = boot(cis_trans_most_coloc_dgn, boot_prop, 1000, formula = sign_cons~dgn_p_interval)
cons_prop_most_coloc_plot = data.frame(prop = apply(cons_prop_most_coloc$t, 2, mean),
                            quant_lower = apply(cons_prop_most_coloc$t, 2, function(x){quantile(x, 0.025)}), 
                            quant_higher = apply(cons_prop_most_coloc$t, 2, function(x){quantile(x, 0.975)}),
                            DGN_p = c('0 ~ 0.01', '0.01 ~ 0.05', '0.05 ~ 0.1', '> 0.1'))
ggplot(cons_prop_most_coloc_plot, aes(x = factor(DGN_p, levels = c('0 ~ 0.01', '0.01 ~ 0.05', '0.05 ~ 0.1', '> 0.1')), 
                           y = prop)) + 
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  labs(x = "DGN p", y = "Proportion of consistent cis-trans", color = '') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# look at unique loci of cis-trans pairs
cis_trans_uniq = c()
temp = cis_trans_most_coloc_dgn
window_size = 250000
while (nrow(temp) > 0) {
  lead_snp = which.max(temp$LOG10P)
  add_loci = temp[lead_snp, ]
  loci_bp = add_loci$bp_hg38
  loci_chr = add_loci$CHROM
  rm_index = which(temp$CHROM == loci_chr &
                     abs(temp$bp_hg38 - loci_bp) <= window_size)

  cis_trans_uniq = rbind(cis_trans_uniq, add_loci)
  temp = temp[-rm_index, ]
}

#cons_prop_uniq = boot(cis_trans_uniq, boot_prop, 1000, formula = sign_cons~dgn_p_interval)
cons_prop_uniq_plot = aggregate(sign_cons~dgn_p_interval, data = cis_trans_uniq, mean)
colnames(cons_prop_uniq_plot) = c('DGN_p', 'prop')
cons_prop_uniq_plot$DGN_p = c('0 ~ 0.01', '0.01 ~ 0.05', '0.05 ~ 0.1', '> 0.1')
ggplot(cons_prop_uniq_plot, aes(x = factor(DGN_p, levels = c('0 ~ 0.01', '0.01 ~ 0.05', '0.05 ~ 0.1', '> 0.1')), 
                                      y = prop)) + 
  #geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  labs(x = "DGN p", y = "Proportion of consistent cis-trans", color = '') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




