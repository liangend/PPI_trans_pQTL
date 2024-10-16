library(data.table)
library(openxlsx)
library(boot)
library(ggplot2)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')

pco_tf = fread('pqtl/04_fdr/ukb/pco_tf_univar_comp.txt')
sum(pco_tf$univar_max_logP < -log10(5e-8/2922)) / nrow(pco_tf)

#### pco pqtl enrichment in TF
tf_gene_id = fread('pqtl/15_TF/humantf_ccbr/TFs_Ensembl_v_1.01.txt', header = F)
hg38_gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.unadj.gtf')
hg38_gene_meta$gene_id = sapply(strsplit(hg38_gene_meta$gene_id, '.', fixed = T), '[', 1)
hg38_gene_meta = hg38_gene_meta[hg38_gene_meta$gene_type == 'protein_coding']
tf_gene_id$gene_name = hg38_gene_meta$gene_name[match(tf_gene_id$V1, hg38_gene_meta$gene_id)]

gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta_sub = gene_meta_hg37[gene_meta_hg37$gene_type == 'protein_coding', ]
gene_meta_sub = gene_meta_sub[!grepl('gene_status', gene_meta_sub$gene_name), ]
gene_meta_sub = gene_meta_sub[!grepl('ENSG', gene_meta_sub$gene_name), ]

near_gene = c()
for (i in 1:nrow(pco_tf)) {
  chr_i = paste0('chr', pco_tf$chr[i])
  gene_meta_i = gene_meta_sub[gene_meta_sub$chr == chr_i, ]
  near_gene[i] = gene_meta_i$gene_name[which.min(abs(gene_meta_i$start - pco_tf$bp_hg37[i]))]
}
pco_tf$near_gene = near_gene

sum(pco_tf$chr == pco_tf$mod_tf_chr & 
      abs(pco_tf$bp_hg37 - pco_tf$mod_tf_tss) < 1000000, na.rm = T)
pco_tf$is_near_gene_tf = (pco_tf$near_gene %in% tf_gene_id$gene_name)


boot_tf = function(trans_all, ukb_p_thre1, ukb_p_thre2, n_sample){
  ukb_sub = trans_all[trans_all$univar_max_logP >= ukb_p_thre1 &
                        trans_all$univar_max_logP < ukb_p_thre2, ]
  print(paste0('Number of cases: ', nrow(ukb_sub)))
  print(paste0('Number of TF: ', sum(ukb_sub$is_tf)))
  #baseline = trans_all[trans_all$ukb_logP > ukb_p_thre, ]
  n_qtl = nrow(ukb_sub)
  enrich = c()
  for (i in 1:n_sample) {
    sample_i = trans_all[sample(1:nrow(trans_all), n_qtl), ]
    enrich[i] = sum(ukb_sub$is_tf) / sum(sample_i$is_tf)
  }
  return(enrich)
}

# TF enrichment
tf_enrich1 = boot_tf(pco_tf, 0, -log10(1e-5), 1000)
tf_enrich2 = boot_tf(pco_tf, -log10(1e-5), -log10(5e-8), 1000)
tf_enrich3 = boot_tf(pco_tf, -log10(5e-8), -log10(5e-8/2922), 1000)
tf_enrich4 = boot_tf(pco_tf, -log10(5e-8/2922), 10000, 100)
tf_enrich5 = boot_tf(pco_tf, 0, -log10(5e-8/2922), 100)

boot_mean = function(data, i) {
  data_i = data[i]
  mean(data_i)
}

tf_boot1 = boot(tf_enrich1, boot_mean, 1000)
tf_boot2 = boot(tf_enrich2, boot_mean, 1000)
tf_boot3 = boot(tf_enrich3, boot_mean, 1000)
tf_boot4 = boot(tf_enrich4, boot_mean, 1000)
tf_boot5 = boot(tf_enrich5, boot_mean, 1000)

tf_plot = data.frame(TF = c(mean(tf_enrich1), mean(tf_enrich2), 
                            mean(tf_enrich3), mean(tf_enrich4)),
                     lower = c(quantile(tf_boot1$t, 0.025), 
                               quantile(tf_boot2$t, 0.025),
                               quantile(tf_boot3$t, 0.025),
                               quantile(tf_boot4$t, 0.025)), 
                     upper = c(quantile(tf_boot1$t, 0.975),
                               quantile(tf_boot2$t, 0.975),
                               quantile(tf_boot3$t, 0.975),
                               quantile(tf_boot4$t, 0.975)),
                     pqtl = c('univar_p1', 'univar_p2', 
                              'univar_p3', 'univar_p4'))

ggplot(tf_plot, aes(x = pqtl, y = TF)) + 
  geom_point(size = 2, color = 'steelblue') +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05, color = 'steelblue') +
  labs(x = "univar minp", y = 'Enrichment', title = 'TF enrichment') + 
  scale_x_discrete(labels = c('> 1e-5', '5e-8 ~ 1e-5', 
                              '5e-8/2922 ~ 5e-8', '< 5e-8/2922')) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# plot ppi and tf enrichment together
tab_tf = data.frame(enrich = c(mean(tf_enrich4), mean(tf_enrich5)),
                    ci_low = c(quantile(tf_boot4$t, 0.025),
                               quantile(tf_boot5$t, 0.025)), 
                    ci_high = c(quantile(tf_boot4$t, 0.975), 
                                quantile(tf_boot5$t, 0.975)),
                    is_nov = c('Shared', 'Novel'))

tf_dif = tf_boot5$t - tf_boot4$t 
pnorm(mean(tf_dif)/sd(tf_dif), lower.tail = F) * 2

ggplot(tab_tf, aes(x=is_nov, y=enrich)) +
  geom_bar(stat="identity", position=position_dodge(0.6), width = 0.3, fill = 'steelblue') +
  geom_errorbar(position=position_dodge(0.6), aes(ymin=ci_low, ymax=ci_high), width=.15) +
  labs(x = "", y = 'Enrichment', title = '', fill = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
