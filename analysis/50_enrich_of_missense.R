library(data.table)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(boot)
setwd('/project/xuanyao/jinghui')
misloc = read.xlsx('pqtl/19_prot_misloc/misloc_supp.xlsx', sheet = 2)
misloc = misloc[!is.na(misloc$dbSNP_ID), ]
snp_list = fread('gtex/10_GWAS_coloc/snp_rsid_hg38/snp.bim')
misloc$chr = snp_list$V1[match(misloc$dbSNP_ID, snp_list$V2)]
misloc$bp_hg38 = snp_list$V4[match(misloc$dbSNP_ID, snp_list$V2)]
misloc_sub = misloc[!is.na(misloc$bp_hg38), ]
misloc_sub = misloc_sub[misloc_sub$chr != 'X', ]

all_ukb_snp = fread('pqtl/00_ref/ukb_geno/all.txt')
all_snp = paste0(all_ukb_snp$CHROM, ':', all_ukb_snp$GENPOS)
all_misloc_snp = paste0(misloc_sub$chr, ':', misloc_sub$bp_hg38)
misloc_sub$ukb_snpID = all_ukb_snp$ID[match(all_misloc_snp, all_snp)]
misloc_sub = misloc_sub[!is.na(misloc_sub$ukb_snpID), ]
misloc_sub$bp_hg37 = sapply(strsplit(misloc_sub$ukb_snpID, ':', fixed = T), '[', 2)
misloc_sub$pco_snp_id = paste0(misloc_sub$chr, ':', misloc_sub$bp_hg37)

common_var = misloc_sub[misloc_sub$Collection == 'hmORFeome 1.1 common variant', ]


ukb_small_p = fread('pqtl/04_fdr/ukb/ukb_univar_small_p.txt')
ukb_sig = ukb_small_p[ukb_small_p$LOG10P > -log10(5e-8/2923), ]
gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')

ukb_sig$prot = sapply(strsplit(ukb_sig$file, '_', fixed = T), '[', 1)
ukb_sig = ukb_sig[ukb_sig$prot %in% gene_meta$gene_name, ]
ukb_sig$CHROM = paste0('chr', ukb_sig$CHROM)
ukb_sig$gene_chr = gene_meta$chr[match(ukb_sig$prot, gene_meta$gene_name)]
ukb_sig$gene_tss = gene_meta$start[match(ukb_sig$prot, gene_meta$gene_name)]
ukb_sig$cis_trans = 'trans'
ukb_sig$cis_trans[which(ukb_sig$gene_chr == ukb_sig$CHROM & 
                          abs(ukb_sig$bp_hg19 - ukb_sig$gene_tss) < 1000000)] = 'cis'

pco_small_p = fread('pqtl/04_fdr/ukb/pco_small_p_mcl.txt')
pco_sig = pco_small_p[pco_small_p$V2 < 5e-8/1088, ]
pco_sig$chr = sapply(strsplit(pco_sig$V1, ':', fixed = T), '[', 2)
pco_sig$bp = sapply(strsplit(pco_sig$V1, ':', fixed = T), '[', 3)
pco_sig$snp = paste0(pco_sig$chr, ':', pco_sig$bp)


n_all_snp = nrow(all_ukb_snp) - 21
n_ukb_sig_snp = length(unique(ukb_sig$ID))
n_pco_sig_snp = length(unique(pco_sig$snp))

n_ukb_sig_snp/n_all_snp
n_pco_sig_snp/n_all_snp

ukb_sig_trans = ukb_sig[ukb_sig$cis_trans == 'trans', ]
ukb_sig_cis = ukb_sig[ukb_sig$cis_trans == 'cis', ]

misloc_sub$is_ukb_sig = misloc_sub$ukb_snpID %in% ukb_sig$ID
misloc_sub$is_ukb_trans = misloc_sub$ukb_snpID %in% ukb_sig_trans$ID
is_ukb_cis = c()
for (i in 1:nrow(misloc_sub)) {
  target_i = misloc_sub$Gene[i]
  snp_i = misloc_sub$ukb_snpID[i]
  cis_i = which(ukb_sig_cis$prot == target_i & ukb_sig_cis$ID == snp_i)
  is_ukb_cis[i] = ifelse(length(cis_i) > 0, T, F)
  print(i)
}
misloc_sub$is_ukb_cis = is_ukb_cis

misloc_sub$ukb_cat = 'none'
misloc_sub$ukb_cat[which(misloc_sub$is_ukb_cis & !misloc_sub$is_ukb_trans)] = 'cis only'
misloc_sub$ukb_cat[which(!misloc_sub$is_ukb_cis & misloc_sub$is_ukb_trans)] = 'trans only'
misloc_sub$ukb_cat[which(misloc_sub$is_ukb_cis & misloc_sub$is_ukb_trans)] = 'both'
misloc_w_score = misloc_sub[!is.na(misloc_sub$Impact.score), ]
n_signal = as.data.frame(table(misloc_w_score$ukb_cat))
n_signal$Var1 = factor(n_signal$Var1, levels = c('none', 'cis only', 'trans only', 'both'))
ggplot(data = n_signal, 
       aes(x = Var1, y = Freq, fill = Var1, label = Freq)) +
  geom_bar(stat="identity", width = 0.8) +
  geom_text(hjust=0.4, vjust=-0.5, size = 4.5, color = 'black') +
  labs(x = '', y = '# SNPs', title = '') +
  scale_fill_manual(values = brewer.pal(4, "Set1")) + 
  ylim(c(0,170)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 14),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

boot_mean = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, mean)
  return(aggre_mean[,2])
}

score_boot = boot(misloc_w_score, boot_mean, 1000, formula = Impact.score~ukb_cat)
score_plot = data.frame(impact = apply(score_boot$t, 2, mean),
                        quant_lower = apply(score_boot$t, 2, function(x){quantile(x, 0.025)}),
                        quant_higher = apply(score_boot$t, 2, function(x){quantile(x, 0.975)}),
                        cat = c('both', 'cis only', 'none', 'trans only'))
score_plot$cat = factor(score_plot$cat, levels = c('none', 'cis only', 'trans only', 'both'))
ggplot(score_plot, aes(x = cat, y = impact, color = cat)) + 
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  labs(x = "", y = "Impact score", color = '') +
  scale_color_manual(values=c(brewer.pal(4,"Set1"))) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

misloc_sub$is_pco_sig = misloc_sub$pco_snp_id %in% pco_sig$snp
sum(misloc_sub$is_ukb_sig) / nrow(misloc_sub) / 
  (n_ukb_sig_snp/n_all_snp)
sum(misloc_sub$is_pco_sig) / nrow(misloc_sub) / 
  (n_pco_sig_snp/n_all_snp)

enrich_plot = data.frame(data = c('Univar', 'PPI'),
                         enrich = c(5.6, 8.3))
ggplot(data = enrich_plot, 
       aes(x = data, y = enrich, fill = data)) +
  geom_bar(stat="identity", width = 0.8) +
  labs(x = '', y = 'Enrichment of snp in significant trans-pQTLs', title = '') +
  scale_fill_manual(values = brewer.pal(3, "Set1")) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 14),
        axis.text.y = element_text(colour = 'black', size = 14),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

misloc_w_score = misloc_sub[!is.na(misloc_sub$Impact.score), ]
aggregate(Impact.score ~ is_ukb_sig, data = misloc_w_score, mean)
aggregate(Impact.score ~ is_pco_sig, data = misloc_w_score, mean)



