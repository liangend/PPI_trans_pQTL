library(data.table)
library(openxlsx)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(qvalue)
library(dplyr)
setwd('/project/xuanyao/jinghui')
source('software/manhattan_fun.R')
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
mod_annot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_annot.txt')
mod_source = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')

## snp annotation
# gtf = fread('pqtl/00_ref/hg38.knownGene.gtf.gz')
# gtf = gtf[,c(1,3,4,5)]
# gtf = gtf[gtf$V3 %in% c('CDS', 'exon', '3UTR', '5UTR'), ]
# pco_all$chr = paste0('chr', pco_all$chr)
# var_annot = c()
# for (i in 1:nrow(pco_all)) {
#   annot_ind = which(gtf$V1 == pco_all$chr[i] & gtf$V4 <= pco_all$bp_hg38[i] & gtf$V5 >= pco_all$bp_hg38[i])
#   if (length(annot_ind) == 0) {
#     var_annot[i] = 'non-coding'
#   } else {
#     var_annot[i] = paste(unique(gtf$V3[annot_ind]), collapse = ',')
#   }
#   if (i %% 100 == 0) {
#     print(i)
#   }
# }
# pco_all$snp_annot = var_annot

### number of PPI affected by uniq loci
pco_uniq = fread('pqtl/04_fdr/ukb/pco_mcl_uniq_loci.txt')
ggplot(pco_uniq, aes(x=n_mod_include)) + 
  geom_density() + 
  labs(title = '', y = "Density", x = "Number of PPI influenced") + 
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#### overall pco pqtl distribution and highlight novel signals
pco_manh = data.frame(SNP = pco_all$loci, P = pco_all$pco_p,
                      CHR = as.numeric(sub('chr', '', pco_all$chr)), 
                      BP = pco_all$bp_hg37, univar_max_logP = pco_all$univar_max_logP)
pco_manh$P[pco_manh$P < 1e-100] = 1e-100

chr_pos = readRDS('pqtl/00_ref/chromosome_location.rds')
pco_manh$chr_tot = chr_pos$tot[match(pco_manh$CHR, chr_pos$CHR)]
pco_manh$def_pos = pco_manh$BP + pco_manh$chr_tot

gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
pco_nov = pco_all[pco_all$is_nov, ]
uniq_gene_mod = unique(pco_all[, c('chr', 'mod', 'nearest_gene')])
hotspot_gene = as.data.frame(sort(table(uniq_gene_mod$nearest_gene), decreasing = T)[1:25])
colnames(hotspot_gene) = c('gene', 'n_mod')
hotspot_gene = cbind(hotspot_gene, gene_meta_hg37[match(hotspot_gene$gene, gene_meta_hg37$gene_name), 
                                                  c('chr', 'start')])
hotspot_gene$chr = as.numeric(sub('chr', '', hotspot_gene$chr))
hotspot_gene$chr_tot = chr_pos$tot[match(hotspot_gene$chr, chr_pos$CHR)]
hotspot_gene$def_pos = hotspot_gene$start + hotspot_gene$chr_tot
hotspot_gene = hotspot_gene[order(hotspot_gene$def_pos), ]
hotspot_gene$annot = paste0(hotspot_gene$gene, ' (', hotspot_gene$n_mod, ')')
rownames(hotspot_gene) = 1:nrow(hotspot_gene)
hotspot_gene_annot = hotspot_gene
hotspot_gene_annot$annot[c(4,7,9,18,20)] = paste0(hotspot_gene_annot$annot[c(4,7,9,18,20)], 
                                                    ' \n ', 
                                                    hotspot_gene_annot$annot[c(5,8,10,19,21)])
hotspot_gene_annot = hotspot_gene_annot[-c(5,8,10,19,21), ]
hotspot_gene_annot$annot[11] = paste(hotspot_gene_annot$annot[9:12], collapse = ' \n ')
hotspot_gene_annot = hotspot_gene_annot[-c(9,10,12), ]

## trans-eQTL hotspot
# eqtlgen_trans = fread('gtex/00_ref/eQTLGen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz')
# ## significant loci prot pairs
# loci_window = 250000
# eqtlgen_trans_loci = c()
# sig_gene = unique(eqtlgen_trans$GeneSymbol)
# loci_i = 1
# for (i in sig_gene) {
#   eqtlgen_i = eqtlgen_trans[eqtlgen_trans$GeneSymbol == i, ]
#   while (nrow(eqtlgen_i > 0)) {
#     lead_sig = which.min(eqtlgen_i$Pvalue)
#     lead_snp = eqtlgen_i[lead_sig, ]
#     lead_snp$loci = loci_i
#     lead_snp$loci_start = lead_snp$SNPPos - loci_window
#     lead_snp$loci_end = lead_snp$SNPPos + loci_window
#     row_rm = which(eqtlgen_i$SNPChr == lead_snp$SNPChr & eqtlgen_i$SNPPos >= lead_snp$loci_start &
#                      eqtlgen_i$SNPPos <= lead_snp$loci_end)
#     eqtlgen_trans_loci = rbind(eqtlgen_trans_loci, lead_snp)
#     loci_i = loci_i + 1
#     eqtlgen_i = eqtlgen_i[-row_rm, ]
#   }
#   print(i)
# }
# eqtlgen_trans_loci = eqtlgen_trans_loci[order(eqtlgen_trans_loci$GeneSymbol, 
#                                               eqtlgen_trans_loci$SNPChr, 
#                                               eqtlgen_trans_loci$SNPPos), ]
# ## merge overlapping windows
# eqtlgen_trans_new = eqtlgen_trans_loci[1, ]
# for (i in 2:nrow(eqtlgen_trans_loci)) {
#   chr_i = eqtlgen_trans_loci$SNPChr[i]
#   start_i = eqtlgen_trans_loci$loci_start[i]
#   gene_i = eqtlgen_trans_loci$GeneSymbol[i]
#   bp19_i = eqtlgen_trans_loci$SNPPos[i]
#   chr_prev = eqtlgen_trans_new$SNPChr[nrow(eqtlgen_trans_new)]
#   end_prev = eqtlgen_trans_new$loci_end[nrow(eqtlgen_trans_new)]
#   gene_prev = eqtlgen_trans_new$GeneSymbol[nrow(eqtlgen_trans_new)]
#   bp19_prev = eqtlgen_trans_new$SNPPos[nrow(eqtlgen_trans_new)]
#   # merge MHC regions
#   if (gene_prev == gene_i & chr_i == 6 & chr_prev == 6 & bp19_i > 28477797 & bp19_i < 33448354 &
#       bp19_prev > 28477797 & bp19_prev < 33448354) {
#     # merge two overlapping loci
#     add_pool_i = rbind(eqtlgen_trans_new[nrow(eqtlgen_trans_new), ], eqtlgen_trans_loci[i, ])
#     add_loci_i = add_pool_i[which.min(add_pool_i$Pvalue), ]
#     add_loci_i$loci_start = min(add_pool_i$loci_start)
#     add_loci_i$loci_end = max(add_pool_i$loci_end)
#     eqtlgen_trans_new[nrow(eqtlgen_trans_new), ] = add_loci_i
#     next
#   }
#   # keep non-overlapping windows
#   if (gene_prev != gene_i | chr_i != chr_prev | 
#       (gene_prev == gene_i & chr_i == chr_prev & start_i > end_prev)) {
#     # add the loci directly if two loci do not overlap
#     eqtlgen_trans_new = rbind(eqtlgen_trans_new, eqtlgen_trans_loci[i, ])
#   } else {
#     # merge two overlapping loci
#     add_pool_i = rbind(eqtlgen_trans_new[nrow(eqtlgen_trans_new), ], eqtlgen_trans_loci[i, ])
#     add_loci_i = add_pool_i[which.min(add_pool_i$Pvalue), ]
#     add_loci_i$loci_start = min(add_pool_i$loci_start)
#     add_loci_i$loci_end = max(add_pool_i$loci_end)
#     eqtlgen_trans_new[nrow(eqtlgen_trans_new), ] = add_loci_i
#   }
#   if (i %% 100 == 0) {
#     print(i)
#   }
# }
# 
# gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
# gene_meta_sub = gene_meta_hg37[gene_meta_hg37$gene_type == 'protein_coding', ]
# gene_meta_sub = gene_meta_sub[!grepl('gene_status', gene_meta_sub$gene_name), ]
# gene_meta_sub = gene_meta_sub[!grepl('ENSG', gene_meta_sub$gene_name), ]
# eqtlgen_near_gene = c()
# for (i in 1:nrow(eqtlgen_trans_new)) {
#   chr_i = paste0('chr', eqtlgen_trans_new$SNPChr[i])
#   gene_meta_i = gene_meta_sub[gene_meta_sub$chr == chr_i, ]
#   eqtlgen_near_gene[i] = gene_meta_i$gene_name[which.min(abs(gene_meta_i$start - eqtlgen_trans_new$SNPPos[i]))]
# }
# eqtlgen_trans_new$nearest_gene = eqtlgen_near_gene
# fwrite(eqtlgen_trans_new, 'gtex/00_ref/eQTLGen/trans_500k_loci.txt', sep = '\t')

## DGN hotspots
# dgn_qtl_intra = read.xlsx('pqtl/00_ref/dgn_qtl.xlsx', sheet = 4)
# dgn_qtl_inter = read.xlsx('pqtl/00_ref/dgn_qtl.xlsx', sheet = 6)
# dgn_qtl_inter$LOG_PVAL = as.numeric(dgn_qtl_inter$LOG_PVAL)
# dgn_qtl_inter_sub = dgn_qtl_inter %>%
#   group_by(GENE_NAME) %>%
#   slice(which.max(LOG_PVAL))
# dgn_trans = rbind(dgn_qtl_intra, dgn_qtl_inter_sub)
# gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
# gene_meta_sub = gene_meta_hg37[gene_meta_hg37$gene_type == 'protein_coding', ]
# gene_meta_sub = gene_meta_sub[!grepl('gene_status', gene_meta_sub$gene_name), ]
# gene_meta_sub = gene_meta_sub[!grepl('ENSG', gene_meta_sub$gene_name), ]
# dgn_near_gene = c()
# for (i in 1:nrow(dgn_trans)) {
#   chr_i = paste0('chr', dgn_trans$SNP_CHROM[i])
#   gene_meta_i = gene_meta_sub[gene_meta_sub$chr == chr_i, ]
#   dgn_near_gene[i] = gene_meta_i$gene_name[which.min(abs(gene_meta_i$start - dgn_trans$SNP_POS[i]))]
# }
# dgn_trans$nearest_gene = dgn_near_gene
# fwrite(dgn_trans, 'gtex/00_ref/dgn_trans_qtl.txt', sep = '\t')

eqtlgen_trans_snp = fread('gtex/00_ref/eQTLGen/trans_snp.txt.gz')
eqtlgen_trans_snp = eqtlgen_trans_snp[eqtlgen_trans_snp$SNPChr > 0, ]
prot_gene_hg37 = gene_meta_hg37[gene_meta_hg37$gene_type == 'protein_coding', ]
prot_gene_hg37 = prot_gene_hg37[!grepl('gene_status', prot_gene_hg37$gene_name), ]
prot_gene_hg37 = prot_gene_hg37[!grepl('ENSG', prot_gene_hg37$gene_name), ]
eqtlgen_near_gene = c()
for (i in 1:nrow(eqtlgen_trans_snp)) {
  chr_i = paste0('chr', eqtlgen_trans_snp$SNPChr[i])
  bp_i = eqtlgen_trans_snp$SNPPos[i]
  prot_gene_sub = prot_gene_hg37[prot_gene_hg37$chr == chr_i, ]
  eqtlgen_near_gene[i] = prot_gene_sub$gene_name[which.min(abs(bp_i - prot_gene_sub$start))]
}
eqtlgen_trans_snp$nearest_gene = eqtlgen_near_gene
sum(hotspot_gene_annot$gene %in% eqtlgen_near_gene)
eqtlgen_trans_new = fread('gtex/00_ref/eQTLGen/trans_500k_loci.txt')
trans_eqtl_hotspot = names(which(table(eqtlgen_trans_new$nearest_gene) > 10))
hotspot_gene_annot$is_eqtlgen_gene = (hotspot_gene_annot$gene %in% eqtlgen_trans_snp$nearest_gene)
hotspot_gene_annot$is_eqtlgen_hot = (hotspot_gene_annot$gene %in% trans_eqtl_hotspot)

eqtl_gene_annot = hotspot_gene_annot[hotspot_gene_annot$is_eqtlgen_hot, ]
pqtl_gene_annot = hotspot_gene_annot[!hotspot_gene_annot$is_eqtlgen_hot, ]
ggplot() +
  geom_point(data = pco_manh, aes(x = def_pos, y = -log10(P), color = factor(CHR))) +
  scale_color_manual(values = rep(c('lightblue','grey'), 11), guide = "none") +
  geom_point(data = pco_manh[pco_manh$univar_max_logP <  -log10(5e-8/2922), ], 
             aes(x = def_pos, y = -log10(P)), color = 'darkred') +
  labs(title = '', y = "-log10(p)", x = "Chromosome") + 
  scale_x_continuous(limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
                     label = chr_pos$CHR, breaks = chr_pos$center, expand = c(0, 0)) +
  annotate(geom="text", x=eqtl_gene_annot$def_pos, y=105, 
           label=eqtl_gene_annot$annot, size = 2.5, angle = 90, hjust = 0, color = 'black') + 
  annotate(geom="text", x=pqtl_gene_annot$def_pos, y=105, 
           label=pqtl_gene_annot$annot, size = 2.5, angle = 90, hjust = 0, color = 'red') + 
  coord_cartesian(ylim = c(11, 100), clip = "off") + 
  # geom_vline(xintercept = hotspot_gene$def_pos, 
  #            linetype="dotted", linewidth = 0.3) +
  theme(text = element_text(size=15, colour = "black"), 
        plot.margin = unit(c(1.5,1,0.5,0.5), "cm"),
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = 'black', size = 11),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

##### loci pleiotropy
univar_z = strsplit(pco_all$univar_z, ',', fixed = T)
univar_z = sapply(univar_z, as.numeric)
n_prot = sapply(strsplit(pco_all$univar_trans_prot, ',', fixed = T), length)

### some univar_z have more z values than number of proteins due to different SNPs at the same genomic position
for (i in 1:length(n_prot)) {
  univar_z[[i]] = univar_z[[i]][1:n_prot[i]]
}
univar_p = sapply(univar_z, function(x){pnorm(abs(x), lower.tail = F)*2})
pco_pleio = sapply(univar_p, function(x){sum(x < 0.05/length(x))})

pco_all$n_univar_p_bonf = pco_pleio
pco_all$pleio = (pco_all$n_univar_p_bonf >= 4 | (pco_all$n_univar_p_bonf > 1 & 
                                                 pco_all$n_univar_p_bonf < 4 & 
                                                 pco_all$n_univar_p_bonf / pco_all$n_prot > 0.2))

boot_mean = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, mean)
  return(aggre_mean[, 2])
}
pleio_boot = boot(pco_all, boot_mean, 1000, formula = pleio ~ is_nov)
pleio_plot = data.frame(pleio = apply(pleio_boot$t, 2, mean),
                        quant_lower = apply(pleio_boot$t, 2, function(x){quantile(x, 0.025)}), 
                        quant_higher = apply(pleio_boot$t, 2, function(x){quantile(x, 0.975)}),
                        is_nov = c('Shared', 'Novel'))
ggplot(pleio_plot, aes(x = is_nov, y = pleio, 
                       fill = is_nov)) + 
  geom_bar(stat="identity", width = 0.5) +
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.3) +
  labs(x = "", y = "Proportion of loci affecting PPI", color = '') +
  coord_cartesian(ylim = c(0, 0.6)) + 
  scale_fill_manual(values=c('darkred', 'lightblue')) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

# pleiotropy example
pco_ex = pco_all[pco_all$loci == 23656, ]
univar_z = as.numeric(unlist(strsplit(pco_ex$univar_z, ',', fixed = T)))
univar_logp = -((pnorm(abs(univar_z), lower.tail = F, log.p = T) + log(2)) * log10(exp(1)))
prot_name = unlist(strsplit(pco_ex$univar_trans_prot, ',', fixed = T))
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
prot_chr = ukb_prot$chr[match(prot_name, ukb_prot$gene_name)]
logp_tab = data.frame(logp = c(univar_logp, -log10(pco_ex$pco_p)), 
                      target = c(paste0(prot_name,' (chr', prot_chr, ')'), 
                                 paste0('PPI ', sub('mod', 'clu', pco_ex$mod))), 
                      group = c(rep('Univariate', length(univar_z)), 'Trans-PCO'))
logp_tab$target = factor(logp_tab$target, levels = c(logp_tab$target[order(univar_logp, decreasing = F)],
                                                     paste0('PPI ', sub('mod', 'clu', pco_ex$mod))))
ggplot(logp_tab, aes(x = target, y = logp, fill = group)) + 
  geom_bar(stat="identity", width = 0.7) +
  coord_flip() +
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
 
# near_gene = c()
# for (i in 1:nrow(pco_all)) {
#   chr_i = paste0('chr', pco_all$chr[i])
#   gene_meta_i = gene_meta_sub[gene_meta_sub$chr == chr_i, ]
#   near_gene[i] = gene_meta_i$gene_name[which.min(abs(gene_meta_i$start - pco_all$bp_hg37[i]))]
# }
# pco_all$nearest_gene = near_gene
# pco_all$nearest_gene_is_tf = (pco_all$nearest_gene %in% tf_gene_id$gene_name)

# near_ppi_gene = c()
# near_ppi_dist = c()
# for (i in 1:nrow(pco_all)) {
#   chr_i = paste0('chr', pco_all$chr[i])
#   prot_i = unlist(strsplit(pco_all$prot_all[i], ',', fixed = T))
#   chr_prot_i = gene_meta_sub$chr[match(prot_i, gene_meta_sub$gene_name)]
#   prot_i_sub = prot_i[which(chr_prot_i == chr_i)]
#   if (length(prot_i_sub) < 1) {
#     near_ppi_gene[i] = NA
#     near_ppi_dist[i] = NA
#   } else {
#     dist_i = abs(gene_meta_sub$start[match(prot_i_sub, gene_meta_sub$gene_name)] - pco_all$bp_hg37[i])
#     near_ppi_gene[i] = prot_i_sub[which.min(dist_i)]
#     near_ppi_dist[i] = min(dist_i)
#   }
# }
# pco_all$near_ppi_gene = near_ppi_gene
# pco_all$near_ppi_dist = near_ppi_dist

# TF enrichment
gene_meta_sub$gene_length = gene_meta_sub$end - gene_meta_sub$start
background_tf = matrix(rep(F, nrow(pco_all) * 1000), ncol = 1000)
for (i in 1:nrow(pco_all)) {
  gene_i = pco_all$nearest_gene[i]
  gene_length_i = gene_meta_sub$gene_length[match(gene_i, gene_meta_sub$gene_name)]
  gene_meta_i = gene_meta_sub[gene_meta_sub$gene_length >= 0.9*gene_length_i & 
                                gene_meta_sub$gene_length <= 1.1*gene_length_i, ]
  background_gene_i = sample(gene_meta_i$gene_name, 1000, replace = T)
  background_tf[i, ] = background_gene_i %in% tf_gene_id$gene_name
  if (i %% 1000 == 0) {
    print(i)
  }
}
background_tf = apply(background_tf, 2, sum)

pco_nov = pco_all[pco_all$is_nov, ]
tf_enrich = sum(pco_all$nearest_gene_is_tf) / background_tf
tf_enrich_nov = sum(pco_nov$nearest_gene_is_tf) / nrow(pco_nov) / 
  (background_tf / nrow(pco_all))

#### pco pqtl enrichment in PPI
## calculate PPI based on if the nearest gene is in PPI cluster
mod_prot = strsplit(pco_all$prot_all, ',', fixed = T)
# background null enrichment in PPI
# n_prot = sapply(mod_prot, length)
# prot_all = unique(unlist(mod_prot))
# n_ite = 1000
# n_ppi = c()
# for (i in 1:n_ite) {
#   is_ppi = c()
#   sample_i = sapply(n_prot, function(x){sample(prot_all, x)})
#   for (j in 1:nrow(pco_all)) {
#     is_ppi[j] = pco_all$nearest_gene[j] %in% sample_i[[j]]
#   }
#   n_ppi[i] = sum(is_ppi)
#   print(i)
# }
# fwrite(as.data.frame(n_ppi), 'pqtl/11_prot_complex/backgroud_ppi_mcl_pqtl.txt')
backgroud_ppi = fread('pqtl/11_prot_complex/backgroud_ppi_mcl_pqtl.txt')
null_ppi = backgroud_ppi$n_ppi / nrow(pco_all)

is_ppi = c()
for (j in 1:nrow(pco_all)) {
  is_ppi[j] = pco_all$nearest_gene[j] %in% mod_prot[[j]]
}
obs_ppi = sum(is_ppi)
ppi_enrich_all = obs_ppi / nrow(pco_all) / null_ppi
# pco_all$nearest_gene_is_ppi = is_ppi

ppi_enrich = sum(pco_all$nearest_gene_is_ppi) / nrow(pco_all) / 
  null_ppi
ppi_enrich_nov = sum(pco_nov$nearest_gene_is_ppi) / 
  nrow(pco_nov) / null_ppi

## calculate PPI based on if the nearest gene is in PPI with any of the PPI subunits
ppi_list = readRDS('pqtl/11_prot_complex/ppi_list.rds')
univar_prot = strsplit(pco_all$univar_trans_prot, ',', fixed = T)
# background null enrichment in PPI
n_ite = 1000
n_ppi = c()
for (i in 1:n_ite) {
  is_ppi = c()
  sample_i = sample(pco_all$nearest_gene)
  prot_pair_i = sapply(1:length(univar_prot), function(x){paste0(univar_prot[[x]], ',', 
                                                                 sample_i[x])})
  ppi_i = sapply(prot_pair_i, function(x){sum(x %in% ppi_list)})
  n_ppi[i] = sum(is_ppi)
  print(i)
}
fwrite(as.data.frame(n_ppi), 'pqtl/11_prot_complex/backgroud_ppi_mcl_pqtl.txt')

prot_pair = sapply(1:length(univar_prot), function(x){paste0(univar_prot[[x]], ',', 
                                                             pco_all$nearest_gene[x])})
obs_ppi = sapply(prot_pair, function(x){sum(x %in% ppi_list)})


# plot ppi and tf enrichment together
ppi_tf = data.frame(enrich = c(mean(tf_enrich_nov), mean(ppi_enrich_nov)),
                    ci_low = c(quantile(tf_enrich_nov, 0.025),
                               quantile(ppi_enrich_nov, 0.025)), 
                    ci_high = c(quantile(tf_enrich_nov, 0.975), 
                              quantile(ppi_enrich_nov, 0.975)),
                    group = c('TF', 'PPI'))

ggplot(ppi_tf, aes(x=group, y=enrich, fill = group)) +
  geom_bar(stat="identity", position=position_dodge(0.1), width = 0.5) +
  geom_errorbar(position=position_dodge(0.1), aes(ymin=ci_low, ymax=ci_high), width=.2) +
  labs(x = "", y = 'Enrichment', title = '', fill = '') + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  scale_fill_manual(values = brewer.pal(3, "Set1")[1:2]) + 
  # scale_fill_manual(values = brewer.pal(3, "Set1")[2:3]) + 
  # geom_segment(aes(x = 0.85, y = 3.9, xend = 1.15, yend = 3.9), color = "black", linewidth=0.22) +
  # annotate(geom="text", x=1, y=4.1, label="***", size = 5) + 
  # geom_segment(aes(x = 1.85, y = 2, xend = 2.15, yend = 2), color = "black", linewidth=0.22) +
  # annotate(geom="text", x=2, y=2.2, label="***", size = 5) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

### GO terms of nearby genes
mcl_near_gene_go = data.frame(go_term = c('GO:MF', 'GO:MF', 'GO:MF', 
                                          'GO:BP', 'GO:BP', 'GO:BP',
                                          'GO:CC', 'GO:CC', 'GO:CC'), 
                              name = c('protein binding', 'DNA-binding transcription activator',
                                       'cis-regulatory region binding', 
                                       'developmental process', 'response to stimulus', 
                                       'regulation of cell communication',
                                       'cytoplasm', 'endomembrane system', 'vesicle'),
                              neg_logP = -log10(c(1e-78, 5.8e-17, 5.8e-17,
                                                  6.9e-106, 2.6e-80, 2.1e-48,
                                                  6.3e-88, 5.2e-57,1.2e-37)))
mcl_near_gene_go$name = factor(mcl_near_gene_go$name, levels= rev(mcl_near_gene_go$name))
mcl_near_gene_go$go_term = factor(mcl_near_gene_go$go_term, levels=unique(mcl_near_gene_go$go_term))
mcl_near_gene_go$neg_logP[mcl_near_gene_go$neg_logP > 30] = 30

ggplot(data = mcl_near_gene_go, aes(x = name, y = neg_logP, fill = go_term)) +
  geom_bar(stat="identity", width = 0.7) + coord_flip() +
  ylim(c(0,30)) + 
  labs(x = '', y = '-log10(p)', fill = '', 
       title = '') +
  scale_fill_manual(values = brewer.pal(3, "Set2")) + 
  theme(text = element_text(size=13, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 13),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())










