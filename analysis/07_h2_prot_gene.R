library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(boot)
setwd('/project/xuanyao/jinghui')
h2_summ_sun = fread('pqtl/05_h2/04_h2_summ/prot_h2_sun_1mb_5mb.txt')
h2_summ_gud = fread('pqtl/05_h2/04_h2_summ/prot_h2_gud_1mb_5mb.txt')
h2_summ_ukb = fread('pqtl/05_h2/04_h2_summ/prot_h2_ukb_1mb_5mb.txt')
h2_summ_dgn = fread('pqtl/05_h2/04_h2_summ/gene_h2_dgn_1mb_5mb.txt')
h2_summ_dgn = h2_summ_dgn[h2_summ_dgn$cis_1mb_h2 !=0, ]
h2_summ_gtex = fread('pqtl/05_h2/04_h2_summ/gene_h2_gtex_1mb_5mb.txt')

trans_prop_sun = mean(h2_summ_sun$trans_5mb_h2_ind) / 
  mean(h2_summ_sun$trans_5mb_h2_ind + h2_summ_sun$cis_1mb_h2)
trans_prop_gud = mean(h2_summ_gud$trans_5mb_h2_ind) / 
  mean(h2_summ_gud$trans_5mb_h2_ind + h2_summ_gud$cis_1mb_h2)
trans_prop_ukb = mean(h2_summ_ukb$trans_5mb_h2_ind) / 
  mean(h2_summ_ukb$trans_5mb_h2_ind + h2_summ_ukb$cis_1mb_h2)
trans_prop_dgn = mean(h2_summ_dgn$trans_5mb_h2_ind) / 
  mean(h2_summ_dgn$trans_5mb_h2_ind + h2_summ_dgn$cis_1mb_h2)
trans_prop_gtex = mean(h2_summ_gtex$trans_5mb_h2_ind) / 
  mean(h2_summ_gtex$trans_5mb_h2_ind + h2_summ_gtex$cis_1mb_h2)

boot_mean = function(data, i) {
  df = data[i, ]
  mean(df$trans_5mb_h2_ind)/
    mean(df$trans_5mb_h2_ind + df$cis_1mb_h2)
}
sun_boot = boot(h2_summ_sun, boot_mean, 1000)
gud_boot = boot(h2_summ_gud, boot_mean, 1000)
ukb_boot = boot(h2_summ_ukb, boot_mean, 1000)
dgn_boot = boot(h2_summ_dgn, boot_mean, 1000)
gtex_boot = boot(h2_summ_gtex, boot_mean, 1000)

trans_prop_all = data.frame(trans_prop = c(trans_prop_sun, trans_prop_gud, trans_prop_ukb, 
                                           trans_prop_dgn, trans_prop_gtex),
                            lower = c(quantile(sun_boot$t, 0.025), 
                                      quantile(gud_boot$t, 0.025),
                                      quantile(ukb_boot$t, 0.025),
                                      quantile(dgn_boot$t, 0.025),
                                      quantile(gtex_boot$t, 0.025)), 
                            upper = c(quantile(sun_boot$t, 0.975), 
                                      quantile(gud_boot$t, 0.975),
                                      quantile(ukb_boot$t, 0.975),
                                      quantile(dgn_boot$t, 0.975),
                                      quantile(gtex_boot$t, 0.975)),
                            data = c('Sun et al., 2018', 'Gudjonsson et al., 2022',
                                     'UKB-PPP, 2023', 'DGN, 2014', 'GTEx blood, 2020'),
                            type = c(rep('Protein', 3), rep('Gene', 2)))
ggplot(trans_prop_all, aes(x=data,  y=trans_prop, label = round(trans_prop, 2), color = type)) +
  geom_point(position=position_dodge(0.3), size = 1.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.07) +
  geom_text(hjust=0.4, vjust=-2, size = 4.5, color = 'black') + 
  labs(x = "", y = bquote('trans/(cis + trans) h'^2), title = '', color = '') + 
  ylim(c(0.5,1.1)) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## using common genes to calculate h2
gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
h2_summ_gud$gene_id = gene_meta$gene_id[match(h2_summ_gud$prot, gene_meta$gene_name)]
h2_summ_sun$gene_id = gene_meta$gene_id[match(h2_summ_sun$gene, gene_meta$gene_name)]
h2_summ_dgn$GeneNameConv = sapply(strsplit(h2_summ_dgn$GeneNameConv, '.', fixed = T), '[', 1)
h2_summ_gtex$gene_id = sapply(strsplit(h2_summ_gtex$gene_id, '.', fixed = T), '[', 1)

comm_gene = Reduce(intersect, list(h2_summ_gud$gene_id, h2_summ_sun$gene_id, h2_summ_ukb$gene_id,
                                   h2_summ_dgn$GeneNameConv, h2_summ_gtex$gene_id))

h2_summ_gud = h2_summ_gud[h2_summ_gud$gene_id %in% comm_gene, ]
h2_summ_sun = h2_summ_sun[h2_summ_sun$gene_id %in% comm_gene, ]
h2_summ_ukb = h2_summ_ukb[h2_summ_ukb$gene_id %in% comm_gene, ]

h2_summ_dgn = h2_summ_dgn[h2_summ_dgn$GeneNameConv %in% comm_gene, ]
h2_summ_gtex = h2_summ_gtex[h2_summ_gtex$gene_id %in% comm_gene, ]

trans_prop_sun = mean(h2_summ_sun$trans_5mb_h2_ind) / 
  mean(h2_summ_sun$trans_5mb_h2_ind + h2_summ_sun$cis_1mb_h2)
trans_prop_gud = mean(h2_summ_gud$trans_5mb_h2_ind) / 
  mean(h2_summ_gud$trans_5mb_h2_ind + h2_summ_gud$cis_1mb_h2)
trans_prop_ukb = mean(h2_summ_ukb$trans_5mb_h2_ind) / 
  mean(h2_summ_ukb$trans_5mb_h2_ind + h2_summ_ukb$cis_1mb_h2)
trans_prop_dgn = mean(h2_summ_dgn$trans_5mb_h2_ind) / 
  mean(h2_summ_dgn$trans_5mb_h2_ind + h2_summ_dgn$cis_1mb_h2)
trans_prop_gtex = mean(h2_summ_gtex$trans_5mb_h2_ind) / 
  mean(h2_summ_gtex$trans_5mb_h2_ind + h2_summ_gtex$cis_1mb_h2)

boot_mean = function(data, i) {
  df = data[i, ]
  mean(df$trans_5mb_h2_ind)/
    mean(df$trans_5mb_h2_ind + df$cis_1mb_h2)
}
sun_boot = boot(h2_summ_sun, boot_mean, 1000)
gud_boot = boot(h2_summ_gud, boot_mean, 1000)
ukb_boot = boot(h2_summ_ukb, boot_mean, 1000)
dgn_boot = boot(h2_summ_dgn, boot_mean, 1000)
gtex_boot = boot(h2_summ_gtex, boot_mean, 1000)

trans_prop_comm = data.frame(trans_prop = c(trans_prop_sun, trans_prop_gud, trans_prop_ukb, 
                                           trans_prop_dgn, trans_prop_gtex),
                            lower = c(quantile(sun_boot$t, 0.025), 
                                      quantile(gud_boot$t, 0.025),
                                      quantile(ukb_boot$t, 0.025),
                                      quantile(dgn_boot$t, 0.025),
                                      quantile(gtex_boot$t, 0.025)), 
                            upper = c(quantile(sun_boot$t, 0.975), 
                                      quantile(gud_boot$t, 0.975),
                                      quantile(ukb_boot$t, 0.975),
                                      quantile(dgn_boot$t, 0.975),
                                      quantile(gtex_boot$t, 0.975)),
                            data = c('Sun et al., 2018', 'Gudjonsson et al., 2022',
                                     'UKB-PPP, 2023', 'DGN, 2014', 'GTEx blood, 2020'),
                            type = c(rep('protein', 3), rep('mRNA', 2)))
ggplot(trans_prop_comm, aes(x=data,  y=trans_prop, label = round(trans_prop, 2), color = type)) +
  geom_point(position=position_dodge(0.3), size = 2.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.12) +
  geom_text(hjust=0.4, vjust=-2.4, size = 4, color = 'black') + 
  labs(x = "", y = bquote('trans/(cis + trans) h'^2), title = 'h2 explained by trans', color = '') + 
  ylim(c(0.5,1.1)) + 
  scale_color_manual(values = brewer.pal(5,"Set1")[c(2,5)]) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

trans_prop_comm_sub = trans_prop_comm[trans_prop_comm$data %in% 
                                        c('GTEx blood, 2020', 'UKB-PPP, 2023'), ]
ggplot(trans_prop_comm_sub, aes(x=data,  y=trans_prop, label = round(trans_prop, 2), fill = type)) +
  geom_bar(stat="identity", width = 0.5) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  scale_fill_manual(values = brewer.pal(3,"Dark2")[1:2]) +
  labs(x = "", y = bquote('trans/(cis + trans) h'^2), title = '', color = '') + 
  # ylim(c(0.5,1.1)) + 
  theme(text = element_text(size=12, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

