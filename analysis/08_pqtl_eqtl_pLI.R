library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
lof = fread('pqtl/06_LoF/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz')
gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)

# gene_bayes = fread('pqtl/06_LoF/gene_bayes.tsv')
# lof$gene_bayes = gene_bayes$post_mean[match(lof$gene_id, gene_bayes$ensg)]

## eQTLGen
eqtlgen_cis = fread('/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
cis_egene = unique(eqtlgen_cis$Gene)
eqtlgen_trans = fread('/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz')
trans_egene = unique(eqtlgen_trans$Gene)
eqtlgen_gene_list = fread('/project2/xuanyao/data/eQTLGen/eQTLGen.gene.txt', header = F)
colnames(eqtlgen_gene_list) = 'gene_id'

eqtlgen_gene_list$pLI = lof$pLI[match(eqtlgen_gene_list$gene_id, lof$gene_id)]
eqtlgen_gene_list$loeuf = lof$oe_lof_upper[match(eqtlgen_gene_list$gene_id, lof$gene_id)]
#eqtlgen_gene_list$gene_bayes = lof$gene_bayes[match(eqtlgen_gene_list$gene_id, lof$gene_id)]
eqtlgen_gene_list$is_cis = (eqtlgen_gene_list$gene_id %in% cis_egene)
eqtlgen_gene_list$is_trans = (eqtlgen_gene_list$gene_id %in% trans_egene)

signal = c()
for (i in 1:nrow(eqtlgen_gene_list)) {
  if (eqtlgen_gene_list$is_cis[i] & eqtlgen_gene_list$is_trans[i]) {
    signal[i] = 'both'
  } else if (eqtlgen_gene_list$is_cis[i] & !eqtlgen_gene_list$is_trans[i]){
    signal[i] = 'cis'
  } else if (!eqtlgen_gene_list$is_cis[i] & eqtlgen_gene_list$is_trans[i]) {
    signal[i] = 'trans'
  } else if (!eqtlgen_gene_list$is_cis[i] & !eqtlgen_gene_list$is_trans[i]) {
    signal[i] = 'no'
  }
}

eqtlgen_gene_list$signal = signal
eqtlgen_gene_list = na.omit(eqtlgen_gene_list)

## 2022 ukbiobank
ukb_pqtl = read.xlsx('pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                     startRow = 5, cols = c(1:11, 20), sheet = 10)
ukb_prot_list = read.xlsx('pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                          startRow = 3, cols = 1:9, sheet = 4)
ukb_prot_list$pLI = lof$pLI[match(ukb_prot_list$Gene.symbol, lof$gene)]
ukb_prot_list$loeuf = lof$oe_lof_upper[match(ukb_prot_list$Gene.symbol, lof$gene)]
#ukb_prot_list$gene_bayes = lof$gene_bayes[match(ukb_prot_list$Gene.symbol, lof$gene)]

uniq_target = unique(ukb_pqtl$Assay.Target)
signal = c()
for (i in 1:length(uniq_target)) {
  signal_i = unique(ukb_pqtl$`cis/trans`[which(ukb_pqtl$Assay.Target == uniq_target[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
ukb_prot_list$signal = signal[match(ukb_prot_list$Assay.Target, uniq_target)]
ukb_prot_list$signal[which(is.na(ukb_prot_list$signal))] = 'no'
ukb_prot_list = na.omit(ukb_prot_list)
ukb_prot_list$gene_id = gene_meta$gene_id[match(ukb_prot_list$Assay.Target, gene_meta$gene_name)]

## signal proportion
all_data = rbind(as.data.frame(table(eqtlgen_gene_list$signal) / nrow(eqtlgen_gene_list)), 
                 as.data.frame(table(ukb_prot_list$signal) / nrow(ukb_prot_list)))
colnames(all_data) = c('signal', 'proportion')
all_data$data = rep(c('eQTLGen', 'UKB prot'), each = 4)

p1 = ggplot(all_data, aes(x = data, y = proportion, 
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

## pLI, LOEUF, gene_bayes difference in different pQTL genes
library(boot)
boot_mean = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, mean)
  return(aggre_mean[,2])
}
eqtlgen_boot = boot(eqtlgen_gene_list, boot_mean, 1000, formula = pLI~signal)
ukb_boot = boot(ukb_prot_list, boot_mean, 1000, formula = pLI~signal)

eqtlgen_boot2 = boot(eqtlgen_gene_list, boot_mean, 1000, formula = loeuf~signal)
ukb_boot2 = boot(ukb_prot_list, boot_mean, 1000, formula = loeuf~signal)

# eqtlgen_boot3 = boot(eqtlgen_gene_list, boot_mean, 1000, formula = gene_bayes~signal)
# ukb_boot3 = boot(ukb_prot_list, boot_mean, 1000, formula = gene_bayes~signal)

pLI_plot = data.frame(pLI = c(apply(eqtlgen_boot$t, 2, mean), apply(ukb_boot$t, 2, mean)),
                      loeuf = c(apply(eqtlgen_boot2$t, 2, mean), apply(ukb_boot2$t, 2, mean)),
                      quant_lower = c(apply(eqtlgen_boot$t, 2, function(x){quantile(x, 0.025)}), 
                                      apply(ukb_boot$t, 2, function(x){quantile(x, 0.025)})), 
                      quant_higher = c(apply(eqtlgen_boot$t, 2, function(x){quantile(x, 0.975)}), 
                                       apply(ukb_boot$t, 2, function(x){quantile(x, 0.975)})),
                      quant_lower2 = c(apply(eqtlgen_boot2$t, 2, function(x){quantile(x, 0.025)}), 
                                       apply(ukb_boot2$t, 2, function(x){quantile(x, 0.025)})), 
                      quant_higher2 = c(apply(eqtlgen_boot2$t, 2, function(x){quantile(x, 0.975)}), 
                                        apply(ukb_boot2$t, 2, function(x){quantile(x, 0.975)})),
                      group = rep(c('cis & trans', 'cis only', 'no signal', 'trans only'), 2), 
                      data = rep(c('eQTLGen', 'UKB prot'), each = 4))

p2 = ggplot(pLI_plot, aes(x = data, y = pLI,
                          color = factor(group, levels = c('cis & trans', 'cis only', 'trans only', 'no signal')))) + 
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  labs(x = "", y = "pLI", color = 'signal') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p3 = ggplot(pLI_plot, aes(x = data, y = loeuf, 
                          color = factor(group, levels = c('cis & trans', 'cis only', 'trans only', 'no signal')))) + 
  geom_errorbar(aes(ymin=quant_lower2, ymax=quant_higher2), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  labs(x = "", y = "LOEUF", color = 'signal') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, p3, nrow = 2)


## contraint quantiles vs. number of different genes
ukb_prot_sub = na.omit(ukb_prot_list)
ukb_prot_sub$pLI_bin = cut(ukb_prot_sub$pLI, breaks=quantile(ukb_prot_sub$pLI, seq(0,1,0.1)), right = T)
ukb_prot_sub$loeuf_bin = cut(ukb_prot_sub$loeuf, breaks=quantile(ukb_prot_sub$loeuf, seq(0,1,0.1)), right = T)

ukb_pLI_gene = table(ukb_prot_sub$pLI_bin, ukb_prot_sub$signal)
ukb_pLI_gene = as.data.frame(ukb_pLI_gene)
colnames(ukb_pLI_gene) = c('pLI', 'signal', 'n_gene')
ukb_pLI_gene$signal = as.character(ukb_pLI_gene$signal)
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'trans'] = 'trans only'
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'cis'] = 'cis only'
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'both'] = 'cis & trans'
ukb_pLI_gene$signal[ukb_pLI_gene$signal == 'no'] = 'no signal'

p5=ggplot(ukb_pLI_gene, aes(x = pLI, y = n_gene, 
                            color = factor(signal, levels = c('cis & trans', 'cis only', 'trans only', 'no signal')))) + 
  geom_point(size = 2) +
  labs(x = "pLI", y = "# genes", color = 'signal', title = 'UKB, 2023') +
  scale_x_discrete(labels = c('low', rep('', 8), 'high')) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

ukb_loeuf_gene = table(ukb_prot_sub$loeuf_bin, ukb_prot_sub$signal)
ukb_loeuf_gene = as.data.frame(ukb_loeuf_gene)
colnames(ukb_loeuf_gene) = c('LOEUF', 'signal', 'n_gene')
ukb_loeuf_gene$signal = as.character(ukb_loeuf_gene$signal)
ukb_loeuf_gene$signal[ukb_loeuf_gene$signal == 'trans'] = 'trans only'
ukb_loeuf_gene$signal[ukb_loeuf_gene$signal == 'cis'] = 'cis only'
ukb_loeuf_gene$signal[ukb_loeuf_gene$signal == 'both'] = 'cis & trans'
ukb_loeuf_gene$signal[ukb_loeuf_gene$signal == 'no'] = 'no signal'
p6=ggplot(ukb_loeuf_gene, aes(x = LOEUF, y = n_gene, 
                              color = factor(signal, levels = c('cis & trans', 'cis only', 'trans only', 'no signal')))) + 
  geom_point(size = 2) +
  labs(x = "LOEUF", y = "# genes", color = 'signal', title = 'UKB, 2023') +
  scale_x_discrete(labels = c('low', rep('', 8), 'high')) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

eqtlgen_sub = na.omit(eqtlgen_gene_list)
eqtlgen_sub$pLI_bin = cut(eqtlgen_sub$pLI, breaks=quantile(eqtlgen_sub$pLI, seq(0,1,0.1)), right = T)
eqtlgen_sub$loeuf_bin = cut(eqtlgen_sub$loeuf, breaks=quantile(eqtlgen_sub$loeuf, seq(0,1,0.1)), right = T)

eqtlgen_pLI_gene = table(eqtlgen_sub$pLI_bin, eqtlgen_sub$signal)
eqtlgen_pLI_gene = as.data.frame(eqtlgen_pLI_gene)
colnames(eqtlgen_pLI_gene) = c('pLI', 'signal', 'n_gene')
eqtlgen_pLI_gene$signal = as.character(eqtlgen_pLI_gene$signal)
eqtlgen_pLI_gene$signal[eqtlgen_pLI_gene$signal == 'trans'] = 'trans only'
eqtlgen_pLI_gene$signal[eqtlgen_pLI_gene$signal == 'cis'] = 'cis only'
eqtlgen_pLI_gene$signal[eqtlgen_pLI_gene$signal == 'both'] = 'cis & trans'
eqtlgen_pLI_gene$signal[eqtlgen_pLI_gene$signal == 'no'] = 'no signal'

p7=ggplot(eqtlgen_pLI_gene, aes(x = pLI, y = n_gene, 
                                color = factor(signal, levels = c('cis & trans', 'cis only', 'trans only', 'no signal')))) + 
  geom_point(size = 2) +
  labs(x = "pLI", y = "# genes", color = 'signal', title = 'eQTLGene') +
  scale_x_discrete(labels = c('low', rep('', 8), 'high')) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

eqtlgen_loeuf_gene = table(eqtlgen_sub$loeuf_bin, eqtlgen_sub$signal)
eqtlgen_loeuf_gene = as.data.frame(eqtlgen_loeuf_gene)
colnames(eqtlgen_loeuf_gene) = c('LOEUF', 'signal', 'n_gene')
eqtlgen_loeuf_gene$signal = as.character(eqtlgen_loeuf_gene$signal)
eqtlgen_loeuf_gene$signal[eqtlgen_loeuf_gene$signal == 'trans'] = 'trans only'
eqtlgen_loeuf_gene$signal[eqtlgen_loeuf_gene$signal == 'cis'] = 'cis only'
eqtlgen_loeuf_gene$signal[eqtlgen_loeuf_gene$signal == 'both'] = 'cis & trans'
eqtlgen_loeuf_gene$signal[eqtlgen_loeuf_gene$signal == 'no'] = 'no signal'

p8=ggplot(eqtlgen_loeuf_gene, aes(x = LOEUF, y = n_gene, 
                                  color = factor(signal, levels = c('cis & trans', 'cis only', 'trans only', 'no signal')))) + 
  geom_point(size = 2) +
  labs(x = "LOEUF", y = "# genes", color = 'signal', title = 'eQTLGene') +
  scale_x_discrete(labels = c('low', rep('', 8), 'high')) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
grid.arrange(p2, p3, p5, p6, p7, p8, nrow = 3)





