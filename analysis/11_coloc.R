library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')
# trait_list = list.files('pqtl/08_gwas_coloc/pco_mcl_coloc/')
# trait_list = trait_list[!grepl('.txt', trait_list)]
# pp4_file = list.files('pqtl/08_gwas_coloc/pco_mcl_coloc/AD/', pattern = 'pp4')
# pp4_all = c()
# for (i in 1:length(trait_list)) {
#   pp4_i = fread(paste0('pqtl/08_gwas_coloc/pco_mcl_coloc/', trait_list[i], '/', pp4_file))
#   pp4_all = cbind(pp4_all, pp4_i$pp4)
#   print(i)
# }
# colnames(pp4_all) = trait_list
# pp4_all = as.data.frame(pp4_all)
# file_info = strsplit(pp4_i$file, '_', fixed = T)
# mod_target = sapply(file_info, function(x){return(x[1])})
# cis_trans = sapply(file_info, function(x){return(x[length(x) - 1])})
# chr_all = sapply(file_info, function(x){return(x[length(x) - 2])})
# loci_all = sapply(file_info, function(x){return(x[length(x)])})
# loci_all = as.numeric(sub('.txt.gz', '', loci_all))
# 
# pp4_all = cbind(data.frame(loci = loci_all),
#                 pp4_all)
# fwrite(pp4_all, 'pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4.txt',
#        sep = '\t')

# pp4_ukb = fread('pqtl/08_gwas_coloc/ukb_coloc/ukb_all_pp4.txt')
# pp4_ukb = as.data.frame(pp4_ukb)
ukb_loci = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')

## keep coloc with gwas loci p < 5e-8
# all_trait = colnames(pp4_ukb)[-(1:3)]
# pp4_sig = c()
# for (i in all_trait) {
#   gwas_loci_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', 
#                              i, '.txt'))
#   gwas_loci_i = gwas_loci_i[gwas_loci_i$P < 5e-8, ]
#   pp4_i = pp4_ukb[, c('file', 'loci', 'cis_trans', i)]
#   colnames(pp4_i)[4] = 'pp4'
#   pp4_i = pp4_i[pp4_i$pp4 > 0.75, ]
#   pp4_i = cbind(pp4_i, ukb_loci[match(pp4_i$loci, ukb_loci$loci), c('CHROM', 'bp_hg38', 'bp_hg19')])
#   pp4_keep = c()
#   for (j in 1:nrow(pp4_i)) {
#     chr_j = pp4_i$CHROM[j]
#     if (i %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {  ## gwas with hg38 bp
#       bp_j = pp4_i$bp_hg38[j]
#     } else {
#       bp_j = pp4_i$bp_hg19[j]
#     }
#     gwas_sub = gwas_loci_i[gwas_loci_i$CHR == chr_j, ]
#     if_in_gwas = which(gwas_sub$loci_start <= bp_j & gwas_sub$loci_end >= bp_j)
#     if (length(if_in_gwas) > 0) {
#       pp4_keep = c(pp4_keep, j)
#     }
#   }
#   pp4_i = pp4_i[pp4_keep, ]
#   pp4_i$trait = i
#   pp4_sig = rbind(pp4_sig, pp4_i)
# }
# fwrite(pp4_sig, 'pqtl/08_gwas_coloc/ukb_coloc/ukb_pp4_gwas_5e8.txt', sep = '\t')

# ukb_uniq_loci = fread('pqtl/04_fdr/ukb/ukb_uniq_500kb_loci.txt')
# ukb_uniq_loci = ukb_uniq_loci[ukb_uniq_loci$CHROM != 23, ]
# pp4_ukb_uniq = pp4_ukb[pp4_ukb$loci %in% ukb_uniq_loci$loci, ]

pp4_sig = fread('pqtl/08_gwas_coloc/ukb_coloc/ukb_pp4_gwas_5e8.txt')
pp4_trans_sig = pp4_sig[pp4_sig$cis_trans == 'trans', ]
pp4_cis_sig = pp4_sig[pp4_sig$cis_trans == 'cis', ]
pp4_ukb_plot = data.frame(n_loci = c(sum(ukb_loci$cis_trans == 'cis'),
                                     sum(ukb_loci$cis_trans == 'trans')),
                          n_coloc_loci = c(length(unique(pp4_cis_sig$loci)),
                                           length(unique(pp4_trans_sig$loci))),
                          total_coloc = c(nrow(pp4_cis_sig), nrow(pp4_trans_sig)),
                          cis_trans = c('cis', 'trans'))


loci_prop = data.frame(prop = c(pp4_ukb_plot$n_coloc_loci/pp4_ukb_plot$n_loci,
                                1-pp4_ukb_plot$n_coloc_loci/pp4_ukb_plot$n_loci),
                       coloc = rep(c('True', 'False'), each = 2),
                       cis_trans = rep(c('cis', 'trans'), 2))

# n_loci_plot$cum_prop = n_loci_plot$loci_prop
# n_loci_plot$cum_prop[4:6] = 1
ggplot(loci_prop, aes(x = cis_trans, y = prop, fill = coloc)) + 
  geom_bar(stat="identity", width = 0.5) +
  #geom_text(aes(y = cum_prop, label = n_loci), vjust=1.6, 
  #          color="black", size=5) +
  labs(x = "", y = "Proportion", title = '') +
  scale_fill_brewer(palette="Paired") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(pp4_ukb_plot, 
       aes(x = cis_trans, y = total_coloc/n_loci, fill = cis_trans)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  #geom_text(aes(label = round(n_causal/n_signal,2)), vjust=-0.3, size=4, position = position_dodge(0.5)) + 
  labs(x = "", y = "Number of coloc traits per locus", title = '') +
  scale_fill_manual(values = brewer.pal(3,"Dark2")[2:3]) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


### gwas explained by cis and trans
immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA', 'asthma')
pp4_pqtl = c()
for (i in immune_trait) {
  pp4_file = list.files(paste0('pqtl/08_gwas_coloc/gwas_coloc_pqtl/', i))
  pp4_all_i = c()
  for (j in pp4_file) {
    pp4_j = fread(paste0('pqtl/08_gwas_coloc/gwas_coloc_pqtl/', i, '/', j))
    pp4_all_i = rbind(pp4_all_i, pp4_j)
  }
  pp4_all_i$trait = i
  gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/gwas_loci_500kb/', 
                        i,'.txt'))
  pp4_all_i$gwas_chr = gwas_i$CHR[match(pp4_all_i$gwas_loci, gwas_i$loci)]
  pp4_pqtl = rbind(pp4_pqtl, pp4_all_i)
  print(i)
}
pp4_pqtl = pp4_pqtl[pp4_pqtl$pp4 > 0.75, ]
pp4_pqtl$gene_name = sapply(strsplit(pp4_pqtl$prot, '_', fixed = T), '[', 1)
pp4_pqtl$signal = ukb_prot_list$signal[match(pp4_pqtl$gene_name, ukb_prot_list$Gene.symbol)]

pp4_cis = pp4_pqtl[pp4_pqtl$pqtl_cis_trans == 'cis', ]
pp4_trans = pp4_pqtl[pp4_pqtl$pqtl_cis_trans == 'trans', ]


n_loci = c()
n_coloc_cis = c()
n_coloc_trans = c()

for (i in immune_trait) {
  gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/gwas_loci_500kb/', 
                        i,'.txt'))
  n_loci[i] = nrow(gwas_i)
  
  pp4_cis_i = pp4_cis[pp4_cis$trait == i, ]
  pp4_trans_i = pp4_trans[pp4_trans$trait == i, ]
  
  cis_coloc_i = unique(pp4_cis_i$gwas_loci)
  trans_coloc_i = unique(pp4_trans_i$gwas_loci)
  
  n_coloc_cis[i] = length(cis_coloc_i)
  n_coloc_trans[i] = length(trans_coloc_i)
}

immune_coloc = data.frame(trait = rep(immune_trait, 2), 
                          coloc_prop = c(n_coloc_trans/n_loci,n_coloc_cis/n_loci),
                          coloc = rep(c('trans', 'cis'), each = length(immune_trait)))
ggplot(data = immune_coloc, aes(x = coloc, y = coloc_prop, 
                                fill = coloc)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(color="black", alpha=0.9, width = 0.2) +
  labs(x = '', y = 'Proportion', fill = 'trans coloc', 
       title = '') +
  scale_fill_manual(values = brewer.pal(4, "Blues")) + 
  # scale_fill_brewer(palette = "Paired") +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


### specific examples
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
pp4_ukb_trait = pp4_ukb_uniq[, c('file', 'loci', 'chr', 'prot', 'cis_trans', 'loci_cat', 'RA')]
pp4_ukb_trait = pp4_ukb_trait[pp4_ukb_trait$RA > 0.75, ]
gwas_locus = fread('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/RA_Ishigaki.tsv.gz')
gwas_locus = gwas_locus[, c('CHR', 'POS', 'P', 'SNP')]


plot_coloc = pp4_ukb_trait[5, ]

# ukb file
ukb_summ_stats = fread(paste0('pqtl/08_gwas_coloc/ukb_qtl_loci/', plot_coloc$file))
ukb_summ_stats$bp_hg37 = as.numeric(sapply(strsplit(ukb_summ_stats$ID, ':', fixed = T), '[', 2))
ukb_summ_stats = ukb_summ_stats[, c(1,10,5,9)]
lead_snp = ukb_summ_stats$snp[which.max(ukb_summ_stats$LOG10P)]
ukb_summ_stats$LOG10P[ukb_summ_stats$LOG10P > 310] = 310
ukb_summ_stats$LOG10P = 10^(-ukb_summ_stats$LOG10P)
colnames(ukb_summ_stats) = c('chrom', 'pos', 'p', 'rsid')
ukb_rsid = ukb_summ_stats$rsid[ukb_summ_stats$rsid != '']

# gwas file
gwas_sub = gwas_locus[gwas_locus$SNP %in% ukb_rsid, ]
colnames(gwas_sub) = c('chrom', 'pos', 'p', 'rsid')
gwas_sub$p = as.numeric(gwas_sub$p)
# gwas_sub$pos = ukb_summ_stats$pos[match(gwas_sub$rsid, ukb_summ_stats$rsid)]

loc_ukb <- locus(data = ukb_summ_stats, index_snp = lead_snp, 
                 flank = 2.5e5,
                 ens_db = "EnsDb.Hsapiens.v75")
#locus_plot(loc_ukb, labels = c('index'))

loc_gwas <- locus(data = gwas_sub, index_snp = lead_snp, 
                  flank = 2.5e5,
                  ens_db = "EnsDb.Hsapiens.v75")
#locus_plot(loc_gwas, labels = c('index'))

multi_layout(nrow = 2,
             plots = {
               locus_plot(loc_ukb, use_layout = FALSE, labels = c('index'), main = 'UKB trans 2')
               locus_plot(loc_gwas, use_layout = FALSE, labels = c('index'), main = 'RA GWAS')
             })








