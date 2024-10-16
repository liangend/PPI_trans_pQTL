library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui/pqtl')
lof = fread('06_LoF/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz')
gene_meta = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
ukb_prot = fread('05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_prot$gene_name = gene_meta$gene_name[match(ukb_prot$gene_id, gene_meta$gene_id)]
ukb_prot$pLI = lof$pLI[match(ukb_prot$gene_name, lof$gene)]
ukb_prot = na.omit(ukb_prot)

## ukb trans h2
trans_h2 = c()
for (i in 1:nrow(ukb_prot)) {
  file_i = ukb_prot$file[i]
  h2_i = fread(paste0('05_h2/02_h2_results/ukb_h2/trans_5mb/', file_i, '_h2.log'), 
              skip = 29, nrows = 1, header = F)
  trans_h2[i] = h2_i$V5
}

## ukb total h2
total_h2 = c()
for (i in 1:nrow(ukb_prot)) {
  file_i = ukb_prot$file[i]
  h2_i = fread(paste0('05_h2/02_h2_results/ukb_h2/all/', file_i, '_h2.log'), 
               skip = 29, nrows = 1, header = F)
  total_h2[i] = h2_i$V5
}

ukb_prot$trans_h2 = trans_h2
ukb_prot$total_h2 = total_h2
ukb_prot$trans_prop = ukb_prot$trans_h2 / ukb_prot$total_h2
ukb_prot_sub = ukb_prot[ukb_prot$trans_prop > 0 & ukb_prot$trans_prop < 1, ]
ukb_prot_sub$trans_prop_bin = cut(ukb_prot_sub$trans_prop, 
                              breaks=quantile(ukb_prot_sub$trans_prop, seq(0,1,0.1)), 
                              right = T)

ukb_trans_summ = as.data.frame(aggregate(pLI ~ trans_prop_bin, data = ukb_prot_sub, mean))
ukb_trans_summ$se = aggregate(pLI ~ trans_prop_bin, data = ukb_prot_sub, sd)[,2] / 
  table(ukb_prot_sub$trans_prop_bin)^0.5

ggplot(ukb_trans_summ, aes(x = trans_prop_bin, y= pLI)) + 
  geom_point(color = 'steelblue') +
  geom_errorbar(aes(ymin=pLI-se, ymax=pLI+se), width=.05, color = 'steelblue') + 
  labs(x = "trans-h2 prop", y = 'pLI', title = 'UKB, 2022') +
  theme(text = element_text(size=16, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

