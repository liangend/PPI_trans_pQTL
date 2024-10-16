library(data.table)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui/pqtl')
sig_ukb_overlap = fread('11_prot_complex/shared_cis_trans.txt')
prot_module = fread('11_prot_complex/prot_module_corum_2022.txt')
prot_module$mod = paste0('mod', prot_module$mod)
colnames(prot_module)[5] = 'gene_chr'

sig_ukb_merge = merge(sig_ukb_overlap, prot_module, by = 'mod')

beta_all = c()
se_all = c()
p_all = c()
for (i in 1:nrow(sig_ukb_merge)) {
  file_i = sig_ukb_merge$file[i]
  bp_i = sig_ukb_merge$bp[i]
  chr_i = sig_ukb_merge$chr[i]
  chr_all = list.files(paste0('UKB_PPP/', file_i))
  chr_i = chr_all[grep(paste0('chr', chr_i, '_'), chr_all)]
  summ_i = fread(paste0('UKB_PPP/', file_i, '/', chr_i))
  pos_i = grep(paste0(':', bp_i, ':'), summ_i$ID)[1]
  if (length(pos_i) > 0) {
    beta_all[i] = summ_i$BETA[pos_i]
    se_all[i] = summ_i$SE[pos_i]
    p_all[i] = summ_i$LOG10P[pos_i]
  } else {
    beta_all[i] = NA
    se_all[i] = NA
    p_all[i] = NA
  }
  print(i)
}
sig_ukb_merge$beta = beta_all
sig_ukb_merge$se = se_all
sig_ukb_merge$logP = p_all

rows = c()
cols = c()
uniq_loci = sig_ukb_overlap$loci
n_loci = 9
n_prot = table(sig_ukb_merge$loci)
for (i in 1:9) {
  loci_group = uniq_loci[((i-1)*9+1):min(i*9, length(uniq_loci))]
  n_i = sum(sig_ukb_merge$loci %in% loci_group)
  rows = c(rows, rep(i, n_i))
  cols = c(cols, rep(1:length(loci_group), n_prot[which(names(n_prot) %in% loci_group)]))
}
sig_ukb_merge$cis_trans = 'cis'
sig_ukb_merge$cis_trans[which(sig_ukb_merge$chr != sig_ukb_merge$gene_chr)] = 'trans'
sig_ukb_merge$rows = rows
sig_ukb_merge$cols = cols
sig_ukb_merge$sig_or_not = 'sig'
sig_ukb_merge$sig_or_not[which(sig_ukb_merge$logP < -log10(0.05/245))] = 'not'
#fwrite(sig_ukb_merge, '11_prot_complex/shared_cis_trans_w_beta.txt', sep = '\t')


sig_ukb_merge = fread('11_prot_complex/shared_cis_trans_w_beta.txt')
sig_ukb_merge$prot_plot = paste0(sig_ukb_merge$cis_trans, sig_ukb_merge$prot)
sig_ukb_merge$sign = as.factor(sign(sig_ukb_merge$beta))
loci_keep = c(25,1919,1956,2121,2188,2292,2332,2380,2940,3349,3389,4243,5000,
              5312,7056,7257,7265,7267,7270, 7272,7276,7287,7310,7456,7463,7465,7467,
              7485,7513,7543,8048,8052,8106,8111,9373,9395,9779,10614,10630,10683,10899)
sig_ukb_merge_plot = sig_ukb_merge[sig_ukb_merge$loci %in% loci_keep, ]
rows = c()
n_loci = 6
n_prot = table(sig_ukb_merge_plot$loci)
for (i in 1:7) {
  loci_group = loci_keep[((i-1)*n_loci+1):min(i*n_loci, length(loci_keep))]
  n_i = sum(sig_ukb_merge_plot$loci %in% loci_group)
  rows = c(rows, rep(i, n_i))
}
i=2
sig_ukb_merge_plot$rows = rows
p1 = ggplot(sig_ukb_merge_plot[sig_ukb_merge_plot$rows == 1, ], 
       aes(x = prot_plot, y = beta, fill = sign, color = sig_or_not)) + 
  geom_bar(stat="identity", width = 0.5) + 
  facet_grid(cols = vars(loci), scales = 'free', drop = F) +
  scale_fill_manual(values=c('lightblue','pink')) +
  scale_color_manual(values=c('lightgray','black')) + 
  scale_x_discrete(name = '', breaks=sig_ukb_merge_plot$prot_plot, labels=sig_ukb_merge_plot$prot) + 
  theme(text = element_text(size=10, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        legend.position = 'none')

p2 = ggplot(sig_ukb_merge_plot[sig_ukb_merge_plot$rows == 2, ], 
            aes(x = prot_plot, y = beta, fill = sign, color = sig_or_not)) + 
  geom_bar(stat="identity", width = 0.5) + 
  facet_grid(cols = vars(loci), scales = 'free', drop = F) +
  scale_fill_manual(values=c('lightblue','pink')) +
  scale_color_manual(values=c('lightgray','black')) + 
  scale_x_discrete(name = '', breaks=sig_ukb_merge_plot$prot_plot, labels=sig_ukb_merge_plot$prot) + 
  theme(text = element_text(size=10, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        legend.position = 'none')

p3 = ggplot(sig_ukb_merge_plot[sig_ukb_merge_plot$rows == 3, ], 
            aes(x = prot_plot, y = beta, fill = sign, color = sig_or_not)) + 
  geom_bar(stat="identity", width = 0.5) + 
  facet_grid(cols = vars(loci), scales = 'free', drop = F) +
  scale_fill_manual(values=c('lightblue','pink')) +
  scale_color_manual(values=c('lightgray','black')) + 
  scale_x_discrete(name = '', breaks=sig_ukb_merge_plot$prot_plot, labels=sig_ukb_merge_plot$prot) + 
  theme(text = element_text(size=10, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        legend.position = 'none')

p4 = ggplot(sig_ukb_merge_plot[sig_ukb_merge_plot$rows == 4, ], 
            aes(x = prot_plot, y = beta, fill = sign, color = sig_or_not)) + 
  geom_bar(stat="identity", width = 0.5) + 
  facet_grid(cols = vars(loci), scales = 'free', drop = F) +
  scale_fill_manual(values=c('lightblue','pink')) +
  scale_color_manual(values=c('lightgray','black')) + 
  scale_x_discrete(name = '', breaks=sig_ukb_merge_plot$prot_plot, labels=sig_ukb_merge_plot$prot) + 
  theme(text = element_text(size=10, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        legend.position = 'none')

p5 = ggplot(sig_ukb_merge_plot[sig_ukb_merge_plot$rows == 5, ], 
            aes(x = prot_plot, y = beta, fill = sign, color = sig_or_not)) + 
  geom_bar(stat="identity", width = 0.5) + 
  facet_grid(cols = vars(loci), scales = 'free', drop = F) +
  scale_fill_manual(values=c('lightblue','pink')) +
  scale_color_manual(values=c('lightgray','black')) + 
  scale_x_discrete(name = '', breaks=sig_ukb_merge_plot$prot_plot, labels=sig_ukb_merge_plot$prot) + 
  theme(text = element_text(size=10, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        legend.position = 'none')

p6 = ggplot(sig_ukb_merge_plot[sig_ukb_merge_plot$rows == 6, ], 
            aes(x = prot_plot, y = beta, fill = sign, color = sig_or_not)) + 
  geom_bar(stat="identity", width = 0.5) + 
  facet_grid(cols = vars(loci), scales = 'free', drop = F) +
  scale_fill_manual(values=c('lightblue','pink')) +
  scale_color_manual(values=c('lightgray','black')) + 
  scale_x_discrete(name = '', breaks=sig_ukb_merge_plot$prot_plot, labels=sig_ukb_merge_plot$prot) + 
  theme(text = element_text(size=10, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        legend.position = 'none')

p7 = ggplot(sig_ukb_merge_plot[sig_ukb_merge_plot$rows == 7, ], 
            aes(x = prot_plot, y = beta, fill = sign, color = sig_or_not)) + 
  geom_bar(stat="identity", width = 0.5) + 
  facet_grid(cols = vars(loci), scales = 'free', drop = F) +
  scale_fill_manual(values=c('lightblue','pink')) +
  scale_color_manual(values=c('lightgray','black')) + 
  scale_x_discrete(name = '', breaks=sig_ukb_merge_plot$prot_plot, labels=sig_ukb_merge_plot$prot) + 
  theme(text = element_text(size=10, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        legend.position = 'none')

grid.arrange(p1,p2,p3, nrow = 3)
grid.arrange(p4,p5,p6,p7, nrow = 4)


