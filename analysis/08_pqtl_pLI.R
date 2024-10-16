library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(boot)
setwd('/project/xuanyao/jinghui')
lof = fread('pqtl/06_LoF/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz')
gene_bayes = fread('pqtl/06_LoF/gene_bayes.tsv')
lof$gene_bayes = gene_bayes$post_mean[match(lof$gene_id, gene_bayes$ensg)]

## GWAS gene (Mostafavi et al., 2023)
gwas_hit = read.xlsx('pqtl/06_LoF/Mostafavi_gwas_eqtl.xlsx', sheet = 3)
gene_list = read.xlsx('pqtl/06_LoF/Mostafavi_gwas_eqtl.xlsx', sheet = 2)

gwas_hit$chr = sapply(strsplit(gwas_hit$Variant, ':', fixed = T), '[', 1)
gwas_hit$chr = paste0('chr', gwas_hit$chr)
gwas_hit$bp = as.numeric(sapply(strsplit(gwas_hit$Variant, ':', fixed = T), '[', 2))
gwas_gene_id = c()
gwas_gene_name = c()
for (i in 1:nrow(gwas_hit)) {
  chr_i = gwas_hit$chr[i]
  bp_i = gwas_hit$bp[i]
  gene_i = gene_list[gene_list$Chr == chr_i, ]
  near_index = which.min(abs(gene_i$TSS - bp_i))
  gwas_gene_id[i] = gene_i$Ensembl_ID[near_index]
  gwas_gene_name[i] = gene_i$Name[near_index]
}
gwas_hit$gene_id = gwas_gene_id
gwas_hit$gene_name = gwas_gene_name
gwas_hit$pLI = lof$pLI[match(gwas_hit$gene_name, lof$gene)]
gwas_hit$loeuf = lof$oe_lof_upper[match(gwas_hit$gene_name, lof$gene)]
gwas_hit_uniq = na.omit(unique(gwas_hit[,6:9]))

boot_gwas = function(data, i){
  df = data[i, ]
  c(sum(df$pLI > 0.9)/nrow(df), median(df$pLI), 
    mean(df$loeuf), median(df$loeuf))
}

gwas_res = boot(gwas_hit_uniq, boot_gwas, R = 1000)

## 2018 Sun
sun_pqtl = read.xlsx('pqtl/00_ref/2018_Sun_supplement.xlsx', 
                     sheet = 4, startRow = 5, cols = 1:16)
sun_prot_list = fread('pqtl/05_h2/00_ref/Sun_2018_prot_w_coor.txt')

sun_prot_list$pLI = lof$pLI[match(sun_prot_list$gene, lof$gene)]
sun_prot_list$loeuf = lof$oe_lof_upper[match(sun_prot_list$gene, lof$gene)]
sun_prot_list$gene_bayes = lof$gene_bayes[match(sun_prot_list$gene, lof$gene)]

uniq_targe = unique(sun_pqtl$SOMAmer.ID)
signal = c()
for (i in 1:length(uniq_targe)) {
  signal_i = unique(sun_pqtl$`cis/.trans`[which(sun_pqtl$SOMAmer.ID == uniq_targe[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
sun_prot_list$signal = signal[match(sun_prot_list$target, uniq_targe)]
sun_prot_list$signal[which(is.na(sun_prot_list$signal))] = 'no'
sun_prot_list$signal2 = sun_prot_list$signal
sun_prot_list$signal2[sun_prot_list$signal2 == 'both'] = 'cis'
sun_prot_list = na.omit(sun_prot_list)

## 2022 Gudjonsson
gud_pqtl = read.xlsx('pqtl/00_ref/2022_Gudjonsson_pqtl.xlsx', 
                     startRow = 4, cols = 1:12)
gud_prot_info = read.xlsx('pqtl/00_ref/2022_Gudjonsson_prot_info.xlsx', 
                          startRow = 3)
gud_pqtl$target = gud_prot_info$Study.Accession[match(gud_pqtl$SOMAmer, gud_prot_info$Study.tag)]
gud_prot_list = fread('pqtl/05_h2/04_h2_summ/prot_h2_gud_1mb_5mb.txt', select = 1:4)
gud_prot_list$pLI = lof$pLI[match(gud_prot_list$prot, lof$gene)]
gud_prot_list$loeuf = lof$oe_lof_upper[match(gud_prot_list$prot, lof$gene)]
gud_prot_list$gene_bayes = lof$gene_bayes[match(gud_prot_list$prot, lof$gene)]

uniq_targe = unique(gud_pqtl$`Protein.(Entrez.symbol)`)
signal = c()
for (i in 1:length(uniq_targe)) {
  signal_i = unique(gud_pqtl$`cis/trans`[which(gud_pqtl$`Protein.(Entrez.symbol)` == uniq_targe[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
gud_prot_list$signal = signal[match(gud_prot_list$prot, uniq_targe)]
gud_prot_list$signal[which(is.na(gud_prot_list$signal))] = 'no'
gud_prot_list$signal2 = gud_prot_list$signal
gud_prot_list$signal2[gud_prot_list$signal2 == 'both'] = 'cis'
gud_prot_list = na.omit(gud_prot_list)

## 2022 ukbiobank
ukb_pqtl = read.xlsx('pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                     startRow = 5, sheet = 11)
ukb_prot_list = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_prot_list$pLI = lof$pLI[match(ukb_prot_list$gene_name, lof$gene)]
ukb_prot_list$loeuf = lof$oe_lof_upper[match(ukb_prot_list$gene_name, lof$gene)]
ukb_prot_list$gene_bayes = lof$gene_bayes[match(ukb_prot_list$gene_name, lof$gene)]

uniq_target = ukb_prot_list$gene_name
signal = c()
for (i in 1:length(uniq_target)) {
  signal_i = unique(ukb_pqtl$`cis/trans`[which(ukb_pqtl$Assay.Target == uniq_target[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else if (length(signal_i) == 1) {
    signal[i] = signal_i
  } else {
    signal[i] = 'no'
  }
}
ukb_prot_list$signal = signal
ukb_prot_list$signal2 = signal
ukb_prot_list$signal2[which(ukb_prot_list$signal2 == 'both')] = 'cis'
ukb_prot_list = na.omit(ukb_prot_list)

# fwrite(ukb_prot_list, 'pqtl/06_LoF/ukb_gene_list.txt', sep = '\t')

## signal proportion
all_data = rbind(as.data.frame(table(sun_prot_list$signal) / nrow(sun_prot_list)), 
                 as.data.frame(table(gud_prot_list$signal) / nrow(gud_prot_list)),
                 as.data.frame(table(ukb_prot_list$signal) / nrow(ukb_prot_list)))
colnames(all_data) = c('signal', 'proportion')
all_data$data = rep(c('Sun', 'Gudjonsson', 'UKB'), each = 4)

ggplot(all_data, aes(x = factor(data, levels = c('Sun', 'Gudjonsson', 'UKB')), 
                          y = proportion, 
                          fill = factor(signal, levels = c('cis', 'trans', 'both', 'no')))) + 
  geom_bar(stat="identity", color="black", width = 0.4) +
  labs(x = "", fill = 'signal') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## pLI, LOEUF difference in different pQTL genes
boot_median = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, function(x){sum(x > 0.9) / length(x)})
  return(aggre_mean[,2])
}

sun_boot = boot(sun_prot_list, boot_median, 1000, formula = pLI~signal)
gud_boot = boot(gud_prot_list, boot_median, 1000, formula = pLI~signal)
ukb_boot = boot(ukb_prot_list, boot_median, 1000, formula = pLI~signal)

pLI_plot = data.frame(pLI = c(apply(sun_boot$t, 2, mean), apply(gud_boot$t, 2, mean), 
                              apply(ukb_boot$t, 2, mean), mean(gwas_res$t[, 1])),
                      quant_lower = c(apply(sun_boot$t, 2, function(x){quantile(x, 0.025)}), 
                                      apply(gud_boot$t, 2, function(x){quantile(x, 0.025)}),
                                      apply(ukb_boot$t, 2, function(x){quantile(x, 0.025)}),
                                      quantile(gwas_res$t[, 1], 0.025)), 
                      quant_higher = c(apply(sun_boot$t, 2, function(x){quantile(x, 0.975)}), 
                                       apply(gud_boot$t, 2, function(x){quantile(x, 0.975)}),
                                       apply(ukb_boot$t, 2, function(x){quantile(x, 0.975)}),
                                       quantile(gwas_res$t[, 1], 0.975)),
                      signal = c(rep(c('cis & trans', 'cis only', 'no signal', 'trans only'), 3), 'GWAS'), 
                      data = c(rep(c('Sun et al., 2018', 'Gudjonsson et al., 2022', 'UKB-PPP, 2023'), each = 4), 
                               'Mostafavi et al., 2023'))


# pval_comp = data.frame(pLI = c(unlist(ukb_boot$t), gwas_res$t[, 1]),
#                        signal = rep(c('cis & trans', 'cis only', 'no signal', 
#                                       'trans only', 'GWAS'), each = 1000))
library(RColorBrewer)
library(ggpubr)
library(grid)
## p value of cis & trans vs. whole genome
pLI_dif_both = ukb_boot$t[, 1] - 0.1554507
empirical_p0 = mean(pLI_dif_both > 0) 
2*min(c(empirical_p1, 1-empirical_p1))

## p value of cis only vs. cis & trans
pLI_dif_cis_both = ukb_boot$t[,2] - ukb_boot$t[,1]
empirical_p1 = mean(pLI_dif_cis_both > 0) 
2*min(c(empirical_p1, 1-empirical_p1))

## p value of trans only vs. no signal
pLI_dif_trans_no = ukb_boot$t[,4] - ukb_boot$t[,3]
empirical_p2 = mean(pLI_dif_trans_no > 0) 
2*min(c(empirical_p2, 1-empirical_p2))

## p value of trans only vs. GWAS
pLI_dif_trans_gwas = ukb_boot$t[,4] - gwas_res$t[, 1]
empirical_p3 = mean(pLI_dif_trans_gwas > 0) 
2*min(c(empirical_p3, 1-empirical_p3))

## p value of trans only vs. cis (only + both)
ukb_boot2 = boot(ukb_prot_list, boot_median, 1000, formula = pLI~signal2)
pLI_dif_trans_cis = ukb_boot2$t[,3] - ukb_boot2$t[,1]
empirical_p4 = mean(pLI_dif_trans_cis > 0) 
2*min(c(empirical_p4, 1-empirical_p4))

ggplot(pLI_plot[pLI_plot$data %in% c('UKB-PPP, 2023', 'Mostafavi et al., 2023'), ], 
       aes(x = factor(signal, levels = c('cis only', 'cis & trans', 'trans only', 'no signal', 'GWAS')), y = pLI,
           color = factor(signal, levels = c('cis only', 'cis & trans', 'trans only', 'no signal', 'GWAS')))) + 
  geom_point(position=position_dodge(0.3), size = 5) +
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, 
                position=position_dodge(0.3)) +
  geom_hline(yintercept = sum(lof$pLI > 0.9, na.rm = T)/nrow(lof), lty = 2) + 
  geom_segment(aes(x = 1, y = 0.33, xend = 2, yend = 0.33), color = "black", linewidth=0.18) +
  geom_segment(aes(x = 1.5, y = 0.34, xend = 1.5, yend = 0.35), color = "black", linewidth=0.18) +
  geom_segment(aes(x = 3, y = 0.34, xend = 3, yend = 0.35), color = "black", linewidth=0.18) +
  geom_segment(aes(x = 1.5, y = 0.35, xend = 3, yend = 0.35), color = "black", linewidth=0.18) +
  annotate(geom="text", x=2.25, y=0.36, label="p < 0.001", size = 4) + 
  geom_segment(aes(x = 3, y = 0.36, xend = 3, yend = 0.37), color = "black", linewidth=0.18) +
  geom_segment(aes(x = 4, y = 0.36, xend = 4, yend = 0.37), color = "black", linewidth=0.18) +
  geom_segment(aes(x = 3, y = 0.37, xend = 4, yend = 0.37), color = "black", linewidth=0.18) +
  annotate(geom="text", x=3.5, y=0.38, label="p = 0.46", size = 4) +
  geom_segment(aes(x = 3, y = 0.39, xend = 3, yend = 0.4), color = "black", linewidth=0.18) +
  geom_segment(aes(x = 5, y = 0.39, xend = 5, yend = 0.4), color = "black", linewidth=0.18) +
  geom_segment(aes(x = 3, y = 0.4, xend = 5, yend = 0.4), color = "black", linewidth=0.18) +
  annotate(geom="text", x=4, y=0.41, label="p = 0.052", size = 4) +
  labs(x = "", y = "Proportion of genes with pLI > 0.9", color = 'signal', 
       title = 'pLI of pQTL target genes') +
  scale_color_manual(values=c(brewer.pal(4,"Dark2"), 'steelblue')) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


ggplot(pLI_plot, aes(x = factor(data, levels = c('Sun et al., 2018', 'Gudjonsson et al., 2022', 
                                                 'UKB-PPP, 2023', 'Mostafavi et al., 2023')), 
                          y = pLI, 
                          color = factor(signal, levels = c('cis only', 'cis & trans', 'trans only', 
                                                           'no signal', 'GWAS')))) + 
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  geom_hline(yintercept = sum(lof$pLI > 0.9, na.rm = T)/nrow(lof), lty = 2) + 
  labs(x = "", y = "Proportion of genes with pLI > 0.9", color = '') +
  scale_color_manual(values=c(brewer.pal(4,"Dark2"), 'steelblue')) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## pLI of pQTL nearby genes
ukb_pqtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
ukb_pqtl$pLI = lof$pLI[match(ukb_pqtl$nearest_gene, lof$gene)]

## eqtl coloc
eqtl_pp4_gen = fread('pqtl/14_cis_trans_coloc/ukb_trans_eqtlgen_eqtl_coloc.txt')
colnames(eqtl_pp4_gen) = c('gene', 'ukb_loci', 'pp4', 'eqtlgen_minp')
eqtl_coloc_gen = eqtl_pp4_gen[eqtl_pp4_gen$pp4 > 0.75, ]
is_eqtl_gen = rep(F, nrow(ukb_pqtl))
is_eqtl_gen[which(ukb_pqtl$loci %in% unique(eqtl_coloc_gen$ukb_loci))] = T

eqtl_pp4_interval = fread('pqtl/14_cis_trans_coloc/ukb_trans_interval_eqtl_coloc.txt')
colnames(eqtl_pp4_interval) = c('gene', 'ukb_loci', 'pp4', 'interval_minp')
eqtl_coloc_interval = eqtl_pp4_interval[eqtl_pp4_interval$pp4 > 0.75, ]
is_eqtl_interval = rep(F, nrow(ukb_pqtl))
is_eqtl_interval[which(ukb_pqtl$loci %in% unique(eqtl_coloc_interval$ukb_loci))] = T

ukb_pqtl$trans_is_eqtl = (is_eqtl_gen | is_eqtl_interval)
ukb_pqtl$nearest_gene_group = ukb_pqtl$cis_trans
ukb_pqtl$nearest_gene_group[which(ukb_pqtl$cis_trans == 'trans' & ukb_pqtl$trans_is_eqtl)] = 'trans coloc eqtl'
ukb_pqtl$nearest_gene_group[which(ukb_pqtl$cis_trans == 'trans' & !ukb_pqtl$trans_is_eqtl)] = 'trans not coloc eqtl'

table(ukb_pqtl$nearest_gene_group)

ukb_pqtl_boot = boot(ukb_pqtl, boot_median, 1000, formula = pLI ~ nearest_gene_group)


pLI_dif_target = ukb_pqtl_boot$t[,2] - ukb_pqtl_boot$t[,1]
empirical_p1 = mean(pLI_dif_target > 0)
2*min(c(empirical_p1, 1-empirical_p1 ))

pLI_dif_trans = ukb_pqtl_boot$t[,3] - ukb_pqtl_boot$t[,2]
empirical_p2 = mean(pLI_dif_trans > 0)
2*min(c(empirical_p2, 1-empirical_p2 ))

pLI_plot2 = data.frame(pLI = c(apply(ukb_pqtl_boot$t, 2, mean), 
                              mean(gwas_res$t[, 1])),
                      quant_lower = c(apply(ukb_pqtl_boot$t, 2, function(x){quantile(x, 0.025)}), 
                                      quantile(gwas_res$t[, 1], 0.025)), 
                      quant_higher = c(apply(ukb_pqtl_boot$t, 2, function(x){quantile(x, 0.975)}), 
                                       quantile(gwas_res$t[, 1], 0.975)),
                      group = c('cis-pQTL', 'trans-pQTL coloc \n with cis-eQTL',  
                                'trans-pQTL not coloc \n with cis-eQTL', 'GWAS'))

ggplot(pLI_plot2, aes(x = factor(group, levels = c('cis-pQTL', 'trans-pQTL coloc \n with cis-eQTL',  
                                                       'trans-pQTL not coloc \n with cis-eQTL', 'GWAS')), 
                     y = pLI, 
                     color = factor(group, levels = c('cis-pQTL', 'trans-pQTL coloc \n with cis-eQTL',  
                                                          'trans-pQTL not coloc \n with cis-eQTL', 'GWAS')))) + 
  geom_point(position=position_dodge(0.5), size = 5) +
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, position=position_dodge(0.5)) +
  geom_hline(yintercept = sum(lof$pLI > 0.9, na.rm = T)/nrow(lof), lty = 2) + 
  labs(x = "", y = "Proportion of genes with pLI > 0.9", color = '', 
       title = 'pLI of pQTL nearby genes') + ylim(c(0.08, 0.28)) +
  scale_color_manual(values=c(brewer.pal(3,"Dark2"), 'steelblue')) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


## constraint quantiles vs. number of different genes
ukb_prot_sub = na.omit(ukb_prot_list)
ukb_prot_sub$pLI_bin = cut(ukb_prot_sub$pLI, breaks=quantile(ukb_prot_sub$pLI, seq(0,1,0.1)), right = T)
ukb_prot_sub$loeuf_bin = cut(ukb_prot_sub$loeuf, breaks=quantile(ukb_prot_sub$loeuf, seq(0,1,0.1)), right = T)
ukb_prot_sub$bayes_bin = cut(ukb_prot_sub$gene_bayes, breaks=quantile(ukb_prot_sub$gene_bayes, seq(0,1,0.1)), 
                             right = T)

ukb_pLI_gene = table(ukb_prot_sub$pLI_bin, ukb_prot_sub$signal)
ukb_pLI_gene = as.data.frame(ukb_pLI_gene)
colnames(ukb_pLI_gene) = c('pLI', 'signal', 'n_gene')
ukb_pLI_gene$signal = as.character(ukb_pLI_gene$signal)
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'trans'] = 'trans only'
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'cis'] = 'cis only'
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'both'] = 'cis & trans'
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'no'] = 'no signal'

ggplot(ukb_pLI_gene, aes(x = pLI, y = n_gene, 
                         fill = factor(signal, levels = c('cis & trans', 'cis only', 'trans only', 'no signal')))) + 
  geom_bar(stat="identity", width = 0.8) +
  labs(x = "pLI bin", y = "Number of genes", fill = '', title = '') +
  scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels = c('low', rep('', 8), 'high')) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())





