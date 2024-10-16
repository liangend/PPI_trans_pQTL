library(data.table)
library(ggplot2)
setwd('/project/xuanyao/jinghui/')
dgn_beta_ukb_sig = fread('pqtl/12_beta_across_two_data/dgn_beta_ukb_combined_sig.txt')
dgn_beta_ukb_sig = na.omit(dgn_beta_ukb_sig)
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')

dgn_beta_ukb_sig$gene_chr = ukb_prot$chr[match(dgn_beta_ukb_sig$file, 
                                               ukb_prot$file)]
dgn_beta_ukb_sig$gene_tss = ukb_prot$gene_start[match(dgn_beta_ukb_sig$file, 
                                                      ukb_prot$file)]
dgn_beta_ukb_sig$cis_trans = 'trans'
dgn_beta_ukb_sig$cis_trans[which(dgn_beta_ukb_sig$CHROM == dgn_beta_ukb_sig$gene_chr &
                                   abs(dgn_beta_ukb_sig$bp_hg38 - dgn_beta_ukb_sig$gene_tss) < 1000000)] = 'cis'


dgn_beta_ukb_trans = dgn_beta_ukb_sig[dgn_beta_ukb_sig$cis_trans == 'trans', ]
length(unique(dgn_beta_ukb_trans$ID)) / nrow(dgn_beta_ukb_trans)

ukb_trans_all = dgn_beta_ukb_trans
window_size = 250000
uniq_trans_loci = c()
while (nrow(ukb_trans_all) > 0) {
  lead_snp = which.max(ukb_trans_all$LOG10P)
  add_loci = ukb_trans_all[lead_snp, 1:17]
  loci_bp = add_loci$bp_hg38
  loci_chr = add_loci$CHROM
  rm_index = which(ukb_trans_all$CHROM == loci_chr &
                     abs(ukb_trans_all$bp_hg38 - loci_bp) <= window_size)
  ukb_trans_sub = ukb_trans_all[rm_index, ]
  add_loci$dgn_p = min(ukb_trans_sub$dgn_p)
  add_loci$n_prot_affect = length(unique(ukb_trans_sub$gene_id))
  uniq_trans_loci = rbind(uniq_trans_loci, add_loci)
  ukb_trans_all = ukb_trans_all[-rm_index, ]
}

uniq_trans_loci$dgn_p_interval = cut(uniq_trans_loci$dgn_p, c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), 
                                     include.lowest = T)

p1 = ggplot(uniq_trans_loci, aes(x = dgn_p, y = n_prot_affect)) + 
  geom_point(stat="identity", size = 1.5) +
  labs(x = "DGN p value", y = '# prot affected', title = '2065 uniq trans-pQTLs') + 
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

p2 = ggplot(uniq_trans_loci, aes(x = dgn_p_interval, y = n_prot_affect)) + 
  geom_boxplot() +
  labs(x = "DGN p value", y = '# prot affected', title = '2065 uniq trans-pQTLs') + 
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
grid.arrange(p1, p2, nrow = 1)

uniq_trans_loci_one_prot = uniq_trans_loci[uniq_trans_loci$n_prot_affect == 1, ]
uniq_trans_loci_multi_prot = uniq_trans_loci[uniq_trans_loci$n_prot_affect > 1, ]

eqtl_prop_plot = data.frame(prop = c(mean(uniq_trans_loci_one_prot$dgn_p < 0.01),
                                     mean(uniq_trans_loci_multi_prot$dgn_p < 0.01)),
                            N = c(nrow(uniq_trans_loci_one_prot), 
                                  nrow(uniq_trans_loci_multi_prot)),
                            group = c('one prot', 'multiple prot'))

ggplot(eqtl_prop_plot, aes(x = group, y = prop, label = N)) + 
  geom_bar(stat="identity", width = 0.5) +
  labs(x = "", y = 'Proportion', fill = 'cis direction',
       title = 'Proportion of trans-pQTL \n with DGN p < 0.01') +
  geom_text(position=position_dodge(0.9), vjust=-0.2) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




