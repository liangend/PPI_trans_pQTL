library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')

pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
pp4_all = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_all_gwas_5e8.txt')

mcl_mod = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_annot.txt')
mcl_mod$mod = paste0('mod', mcl_mod$mod)

pco_nov = pco_all[pco_all$is_nov & pco_all$nearest_gene_is_ppi, ]
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
pco_nov$nearest_gene_is_ukb = (pco_nov$nearest_gene %in% ukb_prot$gene_name)

mod_n_coloc = aggregate(trait ~ loci, data = pp4_sig, function(x){length(unique(x))})
mod_n_coloc = cbind(mod_n_coloc, mcl_mod[match(mod_n_coloc$mod, mcl_mod$mod), -1])

### pco trans-pQTL pvalue vs univariate p value
pco_ex = pco_all[pco_all$loci == 11518, ]
univar_z = as.numeric(unlist(strsplit(pco_ex$univar_z, ',', fixed = T)))
univar_logp = -((pnorm(abs(univar_z), lower.tail = F, log.p = T) + log(2)) * log10(exp(1)))
prot_name = unlist(strsplit(pco_ex$univar_trans_prot, ',', fixed = T))
logp_tab = data.frame(logp = c(univar_logp, -log10(pco_ex$pco_p)), 
                      target = c(prot_name, paste0('PPI ', sub('mod', 'clu', pco_ex$mod))), 
                      group = c(rep('Univariate', length(univar_z)), 'Trans-PCO'))
logp_tab$target = factor(logp_tab$target, levels = c(paste0('PPI ', sub('mod', 'clu', pco_ex$mod)),
                                                     prot_name[order(univar_logp, decreasing = T)]))
ggplot(logp_tab, aes(x = target, y = logp, fill = group)) + 
  geom_bar(stat="identity", width = 0.7) +
  # coord_flip() +
  labs(x = "", y = "-log10(p)", color = '') +
  scale_fill_manual(values=c('darkred', 'lightblue')) +
  geom_hline(yintercept = -log10(0.05/length(prot_name)), lty = 2) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11, angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

### gwas coloc heatmap
pp4_sig_sub = pp4_all[pp4_all$loci %in% c(35306, 35326), ]
#pp4_sig_sub = pp4_all[pp4_all$loci %in% c(11151, 11156), ]
pp4_sig_sub$loci_trait = paste0(pp4_sig_sub$loci, ':', pp4_sig_sub$trait)
coloc_heatmap = data.frame(trait = rep(unique(pp4_sig_sub$trait), 
                                       length(unique(pp4_sig_sub$loci))),
                           loci = rep(unique(pp4_sig_sub$loci), 
                                      each = length(unique(pp4_sig_sub$trait))))
coloc_heatmap$loci_trait = paste0(coloc_heatmap$loci, ':', coloc_heatmap$trait)
# coloc_heatmap$coloc = as.numeric(coloc_heatmap$loci_trait %in% pp4_sig_sub$loci_trait)
coloc_heatmap$pp4 = pp4_sig_sub$pp4[match(coloc_heatmap$loci_trait, pp4_sig_sub$loci_trait)]
coloc_heatmap[is.na(coloc_heatmap)] = 0
coloc_heatmap$loci = factor(coloc_heatmap$loci, levels = c(35326, 35306))
#coloc_heatmap$loci = factor(coloc_heatmap$loci, levels = c(11151, 11156))
pp4_sum = aggregate(pp4 ~ trait, data = coloc_heatmap, sum)
rm_trait = pp4_sum$trait[which(pp4_sum$pp4 < 0.1)]

half_red = colorRampPalette(c("white","red"))(100)[10]
ggplot(coloc_heatmap[!(coloc_heatmap$trait %in% rm_trait), ], aes(x = trait, y = loci, fill = pp4)) +
  geom_tile(color = "black") + 
  scale_fill_gradientn(colours = c(colorRampPalette(c("white",half_red))(50), 
                                   colorRampPalette(c(half_red, 'red'))(50))) +
  labs(x = '', y = '', fill = '',  title = 'Colocalization with GWAS') + 
  theme(text = element_text(size=13, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 10, angle = 90, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 0),
        axis.ticks.y=element_blank(), 
        axis.line = element_blank(),
        strip.text.y = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

### ppi pqtl coloc near genes that are also in the target ppi
# pp4_sig = pp4_all[pp4_all$pp4 > 0.75, ]
# pp4_sig = cbind(pp4_sig, pco_all[match(pp4_sig$loci, pco_all$loci),
#                                  c(12,13,15,16,17)])
# pp4_sig = pp4_sig[pp4_sig$near_ppi_dist < 1000000, ]
# pp4_sig = cbind(pp4_sig, mcl_mod[match(pp4_sig$mod, mcl_mod$mod), 4:15])
# fwrite(pp4_sig, 'pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_gwas_5e8_coloc_near_ppi.txt', sep = '\t')
pp4_sig = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_gwas_5e8_coloc_near_ppi.txt')
pp4_sig$pleio = pco_all$pleio[match(pp4_sig$loci, pco_all$loci)]
pp4_sig$n_prot = pco_all$n_prot[match(pp4_sig$loci, pco_all$loci)]
pp4_sig$ukb_prot = pco_all$univar_trans_prot[match(pp4_sig$loci, pco_all$loci)]
pp4_sig_pleio = pp4_sig[pp4_sig$pleio, ]
# fwrite(pp4_sig_pleio, 'pqtl/pp4_mcl.txt', sep = '\t')

##### specific example
gwas_locus = fread('gtex/10_GWAS_coloc/GWAS_sum_stats/quant_trait/PLT_30080.tsv.gz')
gwas_locus$P = as.numeric(gwas_locus$P)
# pco p file
pp4_all = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4.txt')
pp4_sig_pleio$file = pp4_all$file[match(pp4_sig_pleio$loci, pp4_all$loci)]
plot_coloc = pp4_sig_pleio[1071, ]
pco_locus = fread(paste0('pqtl/08_gwas_coloc/pco_mcl_loci/', plot_coloc$file))
pco_locus$bp_hg37 = as.numeric(sapply(strsplit(pco_locus$ID, ':', fixed = T), '[', 2))
pco_locus = pco_locus[, c(1,8,6,7)]
colnames(pco_locus) = c('chrom', 'pos', 'p', 'rsid')

common_snp = intersect(pco_locus$rsid, gwas_locus$SNP)
pco_locus = pco_locus[pco_locus$rsid %in% common_snp, ]

# gwas file
gwas_sub = gwas_locus[gwas_locus$SNP %in% common_snp, ]
gwas_sub = gwas_sub[, c('CHR', 'POS', 'P', 'SNP')]
colnames(gwas_sub) = c('chrom', 'pos', 'p', 'rsid')
gwas_sub$pos = pco_locus$pos[match(gwas_sub$rsid, pco_locus$rsid)]

library(locuszoomr)
library(EnsDb.Hsapiens.v75)
lead_snp = pco_locus$rsid[which.min(pco_locus$p)]
lead_pos = pco_locus$pos[which.min(pco_locus$p)]
loc_pco <- locus(data = pco_locus, index_snp = lead_snp, 
                 flank = 2e5,
                 ens_db = "EnsDb.Hsapiens.v75")
locus_plot(loc_pco, 
           # labels = c('index'),
           # xticks = 'none', 
           # filter_gene_name = c('GTF3C5', 'CEL', 'RALGDS','GBGT1','OBP2B','ABO',
           #                      'SURF6','MED22','RPL7A','SURF1','SURF2','SURF4','REXO4',
           #                      'ADAMTS13','CACFD1','SLC2A6','TMEM8C','C9orf96'), showExons = F,
           filter_gene_biotype = 'protein_coding', 
           #cex.text = 0.5,
           cex.axis = 1.2, cex.lab = 1.3, 
           main = paste0('PPI clu', sub('mod', '', plot_coloc$mod)), cex.main = 1.3)

loc_gwas <- locus(data = gwas_sub, index_snp = lead_snp, 
                  flank = 2e5,
                  ens_db = "EnsDb.Hsapiens.v75")
locus_plot(loc_gwas, 
           #labels = c('index'), 
           xticks = 'none', 
           filter_gene_name = 'none',
           cex.axis = 1.2, cex.lab = 1.3, main = 'PLT', cex.main = 1.3)



