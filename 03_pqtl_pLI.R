library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(boot)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')
lof = fread('pqtl/06_LoF/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz')

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
gwas_hit_uniq = na.omit(unique(gwas_hit[,6:9]))

boot_gwas = function(data, i){
  df = data[i, ]
  c(sum(df$pLI > 0.9)/nrow(df), median(df$pLI), 
    mean(df$loeuf), median(df$loeuf))
}

gwas_res = boot(gwas_hit_uniq, boot_gwas, R = 1000)

## UKB-PPP
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
ukb_prot_list = na.omit(ukb_prot_list)

## pLI difference in different pQTL genes
boot_median = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, function(x){sum(x > 0.9) / length(x)})
  return(aggre_mean[,2])
}

ukb_boot = boot(ukb_prot_list, boot_median, 1000, formula = pLI~signal)

pLI_plot = data.frame(pLI = c(apply(ukb_boot$t, 2, mean), mean(gwas_res$t[, 1])),
                      quant_lower = c(apply(ukb_boot$t, 2, function(x){quantile(x, 0.025)}),
                                      quantile(gwas_res$t[, 1], 0.025)), 
                      quant_higher = c(apply(ukb_boot$t, 2, function(x){quantile(x, 0.975)}),
                                       quantile(gwas_res$t[, 1], 0.975)),
                      signal = c(c('cis & trans', 'cis only', 'no signal', 'trans only'), 'GWAS'), 
                      data = c(rep('UKB-PPP, 2023', 4), 'Mostafavi et al., 2023'))

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

ggplot(pLI_plot, 
       aes(x = factor(signal, levels = c('cis only', 'cis & trans', 'trans only', 'no signal', 'GWAS')), y = pLI,
           color = factor(signal, levels = c('cis only', 'cis & trans', 'trans only', 'no signal', 'GWAS')))) + 
  geom_point(position=position_dodge(0.3), size = 5) +
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, 
                position=position_dodge(0.3)) +
  geom_hline(yintercept = sum(lof$pLI > 0.9, na.rm = T)/nrow(lof), lty = 2) + 
  labs(x = "", y = "Proportion of genes with pLI > 0.9", color = 'signal', 
       title = 'pLI of pQTL target genes') +
  scale_color_manual(values=c(brewer.pal(4,"Dark2"), 'steelblue')) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

