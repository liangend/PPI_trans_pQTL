library(ggplot2)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(CMplot)
library(dplyr)
# library(tidyverse)
setwd('/project/xuanyao/jinghui')
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
pco_pleio = pco_all[pco_all$pleio, ]
hotspot_gene = names(which(table(pco_all$nearest_gene) >= 100))

pp4_all = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_all_gwas_5e8.txt')
pp4_sig = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_sig_gwas_5e8.txt')
pp4_all$pco_p = pco_all$pco_p[match(pp4_all$loci, pco_all$loci)]
# remove master regulator's effect
pp4_all$near_gene = pco_all$nearest_gene[match(pp4_all$loci, pco_all$loci)]
pp4_all = pp4_all[!pp4_all$near_gene %in% hotspot_gene, ]

#pp4_all = pp4_all[pp4_all$loci %in% pco_pleio$loci, ]

pp4_sig$near_gene = pco_all$nearest_gene[match(pp4_sig$loci, pco_all$loci)]
pp4_sig = pp4_sig[!pp4_sig$near_gene %in% hotspot_gene, ]

mod_go = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_annot.txt')
mod_go = mod_go[, c(1,4:6,10:12)]
mod_go$mod = paste0('mod', mod_go$mod)

all_trait = unique(pp4_sig$trait)
immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA', 'asthma',"T1D")
blood_trait = c("BASO","BASO-P","EO","EO-P","HCT",
                "HGB","HLSR","HLSR-P","IRF","LYMPH",
                "LYMPH-P","MCH","MCHC","MCV","MONO","MONO-P","MPV","MRV",
                "MSCV","NEUT","NEUT-P","PDW","PLT","PLTC","RBC","RBC-W",
                "RET","RET-P","WBC")
other_trait = c("AD","BMI","HDL","LDL","PD","T2D","WHR",
                "bipolar","height","hypertension",
                "schizo","stroke","weight")

mod_coloc = mod_go[, c('mod')]
mod_coloc = as.data.frame(mod_coloc)
for (i in 1:length(all_trait)) {
  pp4_all_i = pp4_all[pp4_all$trait == all_trait[i], ]
  if (all_trait[i] %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {  ## gwas with hg38 bp
    pp4_all_i$POS = pp4_all_i$bp_hg38
  } else {
    pp4_all_i$POS = pp4_all_i$bp_hg37
  }
  
  gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', all_trait[i], '.txt'))
  gwas_i = gwas_i[order(gwas_i$P)[1:50], ]
  
  keep_pp4_i = c()
  for (j in 1:nrow(gwas_i)) {
    keep_j = which(pp4_all_i$chr == gwas_i$CHR[j] &
                     pp4_all_i$POS >= gwas_i$loci_start[j] &
                     pp4_all_i$POS <= gwas_i$loci_end[j])
    keep_pp4_i = c(keep_pp4_i, keep_j)
  }
  pp4_all_sub = pp4_all_i[keep_pp4_i, ]
  pp4_max = aggregate(pp4 ~ mod, data = pp4_all_sub, max)
  
  # pp4_max = pp4_all_i %>%
  #   group_by(mod) %>%
  #   slice(order(pco_p)[1:10]) %>%
  #   slice(which.max(pp4)) %>%
  #   select(pp4, mod)
  
  mod_coloc = cbind(mod_coloc, pp4_max$pp4[match(mod_coloc$mod, pp4_max$mod)])
}
colnames(mod_coloc)[2:51] = all_trait
mod_coloc[is.na(mod_coloc)] = 0

mod_coloc_sub = mod_coloc[, -1]
mod_coloc_plot = t(as.matrix(mod_coloc_sub))
colnames(mod_coloc_plot) = mod_coloc$mod
coloc_sum = colSums(mod_coloc_plot)
mod_coloc_plot = mod_coloc_plot[, -which(coloc_sum == 0)]

immune_trait0 = c('CD', 'IBD', 'MS', 'UC', 'SLE', 'RA', 'Asthma',"T1D")
other_trait0 = c("AD","BMI","HDL","LDL","PD","T2D","WHR",
                "Bipolar","Height","Hypertension",
                "Schizophrenia","Stroke","Weight")
trait_annot = data.frame(row.names = c(immune_trait0, blood_trait, other_trait0),
                         trait = c(rep('Autoimmune', length(immune_trait0)),
                                   rep('Blood', length(blood_trait)),
                                   rep('Other', length(other_trait0))))
# get the group of mod after plotting heatmap (line127-130)
mod_annot = data.frame(group = c(rep('group1',length(mod1)), rep('group2',length(mod2)),
                                 rep('group3',length(mod3)), rep('group4',length(mod4))))
rownames(mod_annot) = c(mod1, mod2, mod3, mod4)
trait_annot$trait = factor(trait_annot$trait, levels = c('Blood', 'Autoimmune', 'Other'))
trait_colors = list(trait = c('Blood' = 'lightpink1',
                              'Autoimmune' = 'steelblue1',
                              'Other' = 'grey'))

mod_coloc_plot0 = mod_coloc_plot[c(blood_trait, immune_trait, other_trait), ]
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'lupus')] = 'SLE'
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'asthma')] = 'Asthma'
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'bipolar')] = 'Bipolar'
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'height')] = 'Height'
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'hypertension')] = 'Hypertension'
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'schizo')] = 'Schizophrenia'
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'stroke')] = 'Stroke'
rownames(mod_coloc_plot0)[which(rownames(mod_coloc_plot0) == 'weight')] = 'Weight'
mod_coloc_plot0 = mod_coloc_plot0[rev(as.character(coloc_tab$trait)), ] # trait order from 42_gwas_coloc.R

half_blue = colorRampPalette(c("white","#2A788EFF"))(100)[10]
all_trait_heat = pheatmap(mod_coloc_plot0, 
        col = c(colorRampPalette(c("white",half_blue))(50), 
                colorRampPalette(c(half_blue, '#2A788EFF'))(50)), 
        legend = T, show_colnames = F, 
        cluster_rows = F, 
        # treeheight_col = 0, 
        # annotation_col = mod_annot, 
        # annotation_row = trait_annot, 
        annotation_colors = trait_colors,
        fontsize = 8)

# mod1 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[1:75]]
# mod2 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[76:426]]
# mod3 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[427:629]]
# mod4 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[630:993]]

mod1 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[1:62]]
mod2 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[63:116]]
mod3 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[117:831]]
mod4 = colnames(mod_coloc_plot)[all_trait_heat$tree_col$order[832:983]]

mod_prot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')
mod_prot$mod = paste0('mod', mod_prot$mod)
prot1 = data.frame(prot = mod_prot$prot[which(mod_prot$mod %in% mod1)])
prot2 = data.frame(prot = mod_prot$prot[which(mod_prot$mod %in% mod2)])
prot3 = data.frame(prot = mod_prot$prot[which(mod_prot$mod %in% mod3)])
prot4 = data.frame(prot = mod_prot$prot[which(mod_prot$mod %in% mod4)])
fwrite(prot4, 'pqtl/prot_list2/prot4.txt', sep = '\t', col.names = F)


## heatmap for single trait
immune_show = c('RA', 'lupus', 'IBD')
gene_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_hg37 = gene_hg37[gene_hg37$gene_type == 'protein_coding', ]
gene_hg37 = gene_hg37[!grepl('gene_status', gene_hg37$gene_name), ]
gene_hg37 = gene_hg37[!grepl('ENSG', gene_hg37$gene_name), ]

gene_hg38 = fread('gtex/00_ref/gencode.v26.GRCh38.genes.gtf')
gene_hg38 = gene_hg38[gene_hg38$gene_type == 'protein_coding', ]

trait_i = 'IBD'
gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', 
                      trait_i,'.txt'))
gwas_i = gwas_i[gwas_i$P < 5e-8, ]

near_gene = c()
for (i in 1:nrow(gwas_i)) {
  chr_i = paste0('chr', gwas_i$CHR[i])
  pos_i = gwas_i$POS[i]
  if (trait_i %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {  ## gwas with hg38 bp
    gene_meta_i = gene_hg38[gene_hg38$chr == chr_i, ]
  } else {
    gene_meta_i = gene_hg37[gene_hg37$chr == chr_i, ]
  }
  near_gene[i] = gene_hg37$gene_name[which.min(abs(gene_meta_i$start - pos_i))]
}
gwas_i$near_gene = near_gene
pp4_sig$pleio = pco_all$pleio[match(pp4_sig$loci, pco_all$loci)]
pp4_sig = pp4_sig[pp4_sig$pleio, ]
pp4_sig_i = pp4_sig[pp4_sig$trait == trait_i, ]
trait_coloc_i = matrix(0, ncol = nrow(mod_go), nrow = nrow(gwas_i))
for (n in 1:nrow(gwas_i)) {
  gwas_coloc_n = gwas_i[n, ]
  chr_n = gwas_coloc_n$CHR
  start_n = gwas_coloc_n$loci_start
  end_n = gwas_coloc_n$loci_end
  if (trait_i %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {
    mod_n = pp4_sig_i$mod[which(pp4_sig_i$chr == chr_n & 
                                  pp4_sig_i$bp_hg38 >= start_n & pp4_sig_i$bp_hg38 <= end_n)]
  } else {
    mod_n = pp4_sig_i$mod[which(pp4_sig_i$chr == chr_n & 
                                  pp4_sig_i$bp_hg37 >= start_n & pp4_sig_i$bp_hg37 <= end_n)]
  }
  mod_n = as.numeric(sub('mod', '', mod_n))
  if (length(mod_n) > 0) {
    trait_coloc_i[n, mod_n] = 1
  }
}
colnames(trait_coloc_i) = paste0('clu', 1:ncol(trait_coloc_i))
rownames(trait_coloc_i) = gwas_i$near_gene
col_sum_i = apply(trait_coloc_i, 2, sum)
which.max(col_sum_i)
trait_coloc_plot = trait_coloc_i[, which(col_sum_i > 1)]
row_sum_i = apply(trait_coloc_plot, 1, sum)
trait_coloc_plot = trait_coloc_plot[which(row_sum_i > 0), ]

mod_annot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/3_immune_mod_annot.txt')
mod_annot = mod_annot[!mod_annot$annot %in% c('Inflammatory response',
                                             'Protein metabolism'), ]
ann_colors = list(group = c('T cell regulation' = brewer.pal(12, "Paired")[1],
                            'B cell regulation' = brewer.pal(12, "Paired")[2],
                            'Leukocyte regulation' = brewer.pal(12, "Paired")[3],
                            'CD40 signaling pathway' = brewer.pal(12, "Paired")[4],
                            'Cytokine response' = brewer.pal(12, "Paired")[5],
                            'Apoptosis' = brewer.pal(12, "Paired")[6],
                            'Cell adhesion' = brewer.pal(12, "Paired")[7],
                            #'Inflammatory response'= brewer.pal(12, "Paired")[8],
                            'Immune response' = brewer.pal(12, "Paired")[9],
                           'Response to stimulus'= brewer.pal(12, "Paired")[10],
                           'Tissue development' = brewer.pal(12, "Paired")[11],
                           #'Protein metabolism'= brewer.pal(12, "Paired")[12],
                           'Small molecule metabolism'= "black",
                           'Other' = "grey"))
mod_annot$annot = factor(mod_annot$annot, levels=names(ann_colors$group))
mod_annot = mod_annot[order(mod_annot$annot), ]
immune_annot = data.frame(group = mod_annot$annot)
rownames(immune_annot) = sub('mod', 'clu', mod_annot$mod)

trait_coloc_plot = trait_coloc_plot[, rownames(immune_annot)[rownames(immune_annot) %in% 
                                                               colnames(trait_coloc_plot)]]
pheatmap(trait_coloc_plot, cluster_cols = F, 
         col = colorRampPalette(c("white","red"))(100), 
         annotation_col = immune_annot, treeheight_row = 0, 
         annotation_colors = ann_colors,
         legend = F, show_colnames = T, main = trait_i)

##### specific example
gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)

pp4_all = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4.txt')
pp4_sig$near_ppi_gene = pco_all$near_ppi_gene[match(pp4_sig$loci, pco_all$loci)]
pp4_sig$near_ppi_dist = pco_all$near_ppi_dist[match(pp4_sig$loci, pco_all$loci)]

pp4_sig_i = pp4_sig[pp4_sig$trait == 'IBD', ]
pp4_sig_i$file = pp4_all$file[match(pp4_sig_i$loci, pp4_all$loci)]
pp4_sig_i$is_nov = pco_all$is_nov[match(pp4_sig_i$loci, pco_all$loci)]

gwas_locus = fread('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/IBD_KM.tsv.gz')
gwas_locus$P = as.numeric(gwas_locus$P)

plot_coloc = pp4_sig_i[pp4_sig_i$mod == 'mod1035', ] # SLE
plot_coloc = pp4_sig_i[pp4_sig_i$mod == 'mod964', ] # IBD
plot_coloc = pp4_sig_i[pp4_sig_i$mod == 'mod960', ] # RA

plot_coloc = plot_coloc[order(plot_coloc$chr, plot_coloc$bp_hg37), ]
# pco p file
pco_locus = c()
for (i in 1:nrow(plot_coloc)) {
  pco_locus_i = fread(paste0('pqtl/08_gwas_coloc/pco_mcl_loci/', plot_coloc$file[i]))
  pco_locus_i$bp_hg37 = as.numeric(sapply(strsplit(pco_locus_i$ID, ':', fixed = T), '[', 2))
  pco_locus_i = pco_locus_i[, c(1,8,6,7)]
  colnames(pco_locus_i) = c('chrom', 'pos', 'p', 'rsid')
  pco_locus_i$locus = i
  pco_locus = rbind(pco_locus, pco_locus_i)
}
# pco_locus = fread(paste0('pqtl/08_gwas_coloc/pco_mcl_loci/', plot_coloc$file))
# pco_locus$bp_hg37 = as.numeric(sapply(strsplit(pco_locus$ID, ':', fixed = T), '[', 2))
# pco_locus = pco_locus[, c(1,8,6,7)]
# colnames(pco_locus) = c('chrom', 'pos', 'p', 'rsid')
common_snp = intersect(pco_locus$rsid, gwas_locus$SNP)
pco_locus = pco_locus[pco_locus$rsid %in% common_snp, ]
pco_locus = pco_locus[, c('rsid', 'chrom', 'locus', 'pos', 'p')]

# gwas file
gwas_sub = gwas_locus[gwas_locus$SNP %in% common_snp, ]
gwas_sub = gwas_sub[, c('CHR', 'POS', 'P', 'SNP')]
colnames(gwas_sub) = c('chrom', 'pos', 'p', 'rsid')
gwas_sub$pos = pco_locus$pos[match(gwas_sub$rsid, pco_locus$rsid)]

pco_locus$gwas_p = gwas_sub$p[match(pco_locus$rsid, gwas_sub$rsid)]
pco_minp = aggregate(p ~ locus, data = pco_locus, min)
pco_minp$chrom = plot_coloc$chr
lead_snp = c()
lead_pos = c()
cm_plot_dat = c()
window_size = 2e5
for (i in 1:nrow(pco_minp)) {
  lead_snp_i = pco_locus$rsid[which(pco_locus$chrom == pco_minp$chrom[i] &
                                     pco_locus$p == pco_minp$p[i])]
  lead_pos_i = pco_locus$pos[which(pco_locus$chrom == pco_minp$chrom[i] &
                                     pco_locus$p == pco_minp$p[i])]
  pco_i = pco_locus[pco_locus$chrom == pco_minp$chrom[i] & 
                      pco_locus$pos > lead_pos_i - window_size &
                      pco_locus$pos < lead_pos_i + window_size, ]
  pco_i$pos = pco_i$pos - min(pco_i$pos) + 1
  cm_plot_dat = rbind(cm_plot_dat, pco_i)
  lead_snp = c(lead_snp, lead_snp_i)
  lead_pos = c(lead_pos, lead_pos_i)
}
lead_pos = round(lead_pos/1e6, 1)
cm_plot_dat$p = -log10(cm_plot_dat$p)
cm_plot_dat$gwas_p = -log10(cm_plot_dat$gwas_p)
# correct for small p values of SLE
# cm_plot_dat$p[which(cm_plot_dat$chrom == 6)] =
#   cm_plot_dat$p[which(cm_plot_dat$chrom == 6)]/10
# cm_plot_dat$gwas_p[which(cm_plot_dat$chrom == 6)] =
#   cm_plot_dat$gwas_p[which(cm_plot_dat$chrom == 6)]/10
CMplot(cm_plot_dat[, -2], plot.type="c", LOG10 = F,
       # chr.labels=paste0("chr", gwas_minp$chrom, ':', 
       #                   lead_pos - window_size/1e6, '-', 
       #                   lead_pos + window_size/1e6, ' Mb'),
       chr.labels = c(rep('', length(lead_pos))),
       r=8, band = 2, cir.axis=T, cir.axis.grid=T, 
       axis.cex = 1.5, H = 4,
       highlight = lead_snp, highlight.pch = 15, highlight.text = 'lead',
       outward=T, cir.axis.col="black", cir.chr.h=2, cir.band = 3,
       chr.den.col=NULL, file.output=F)
paste0("chr", pco_minp$chrom, ':', lead_pos - window_size/1e6, '-', 
       lead_pos + window_size/1e6, ' Mb')

### PPI clusters with multiple trans-pQTL coloc with the same trait
ppi_converge = c()
disease_trait = c(immune_trait, 'AD', 'PD', 'T2D', 'bipolar', 'hypertension',
                  'schizo', 'stroke')
pp4_sig_sub = pp4_sig[pp4_sig$pleio, ]
for (i in disease_trait) {
  pp4_sig_i = pp4_sig_sub[pp4_sig_sub$trait == i, ]
  n_mod = as.data.frame(table(pp4_sig_i$mod))
  colnames(n_mod) = c('mod', 'n')
  n_mod = n_mod[n_mod$n > 1, ]
  if (nrow(n_mod) > 0) {
    n_mod$trait = i
    ppi_converge = rbind(ppi_converge, n_mod)
  }
}
ppi_converge$mod = as.character(ppi_converge$mod)
ppi_converge = ppi_converge[ppi_converge$n > 3, ]
sum(ppi_converge$mod %in% mod_annot$mod)
ppi_converge = cbind(ppi_converge, mod_go[match(ppi_converge$mod, mod_go$mod), 2:7])
ppi_converge$annot = mod_annot$annot[match(ppi_converge$mod, mod_annot$mod)]
ukb_prot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_annot.txt')
ukb_prot$mod = paste0('mod', ukb_prot$mod)
ppi_converge$prot = ukb_prot$prot[match(ppi_converge$mod, ukb_prot$mod)]
fwrite(ppi_converge, 'pqtl/ppi_converge.txt', sep = '\t')
ppi_converge = fread('pqtl/ppi_converge.txt')


# eqtlgen file
eqtlgen_coloc = fread('pqtl/16_pco_cis_coloc/mcl_eqtlgen_eqtl_coloc.txt')
colnames(eqtlgen_coloc) = c('gene', 'loci', 'pp4')
eqtlgen_coloc = eqtlgen_coloc[eqtlgen_coloc$pp4 > 0.75, ]

gene_target = eqtlgen_coloc$gene[which(eqtlgen_coloc$loci == plot_coloc$loci)]
eqtlgen_locus = fread(paste0('gtex/00_ref/eQTLGen/cis_chr', plot_coloc$chr, '.txt.gz'))
eqtlgen_locus = eqtlgen_locus[eqtlgen_locus$SNPPos  >= min(pco_locus$pos) & 
                                eqtlgen_locus$SNPPos <= max(pco_locus$pos), ]
eqtlgen_locus = eqtlgen_locus[eqtlgen_locus$GeneSymbol %in% gene_target, ]
eqtlgen_locus = eqtlgen_locus[,c(1:3,6)]
colnames(eqtlgen_locus) = c('p', 'chrom', 'pos', 'gene')
eqtlgen_locus$rsid = pco_locus$rsid[match(eqtlgen_locus$pos, pco_locus$pos)]
eqtlgen_locus = eqtlgen_locus[eqtlgen_locus$rsid %in% common_snp, ]


library(locuszoomr)
library(EnsDb.Hsapiens.v75)
lead_snp = gwas_sub$rsid[which.min(gwas_sub$p)]
lead_pos = gwas_sub$pos[which.min(gwas_sub$p)]
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
           main = paste0('trans-PCO clu', sub('mod', '', plot_coloc$mod)), cex.main = 1.3)

loc_gwas <- locus(data = gwas_sub, index_snp = lead_snp, 
                 flank = 2e5,
                 ens_db = "EnsDb.Hsapiens.v75")
locus_plot(loc_gwas, 
           #labels = c('index'), 
           xticks = 'none', 
           filter_gene_name = 'none',
           cex.axis = 1.2, cex.lab = 1.3, main = 'SLE', cex.main = 1.3)

# multi_layout(nrow = 2,
#              plots = {
#                locus_plot(loc_pco, use_layout = FALSE, labels = c('index'), main = 'PCO', 
#                           xticks = 'none', filter_gene_name = 'none')
#                locus_plot(loc_gwas, use_layout = FALSE, labels = c('index'), main = 'IBD GWAS')
#              })


## GO term of the target protein module
mod_go1 = data.frame(name = c('detoxification', 'hydrogen peroxide \n metabolism', 
                             'response to superoxide', 'response to oxygen radical',
                             'hydrogen peroxide catabolism'),
                    neg_logP = c(10.9, 8.4, 6.8, 6.8, 5.6))
mod_go1$name = factor(mod_go1$name, levels= mod_go1$name)
ggplot(data = mod_go1, aes(x = name, y = neg_logP)) +
  geom_bar(stat="identity", width = 0.7, fill = brewer.pal(4, 'Accent')[1]) + 
  labs(x = '', y = '', title = 'group 1') +
  # scale_fill_manual(values = brewer.pal(4, "Set1")[1:4]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15, angle = 60, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

mod_go2 = data.frame(name = c('G protein-coupled \n receptor pathway', 
                              'nucleoside \n monophosphate metabolism',
                              'ribonucleoside \n monophosphate metabolism',
                              'pyrimidine \n nucleotide metabolism', 'AMP metabolism'),
                    neg_logP = c(15.4, 10.7, 8.2, 7.9, 5.2))
mod_go2$name = factor(mod_go2$name, levels= mod_go2$name)
ggplot(data = mod_go2, aes(x = name, y = neg_logP)) +
  geom_bar(stat="identity", width = 0.7, fill = brewer.pal(4, 'Accent')[3]) + 
  labs(x = '', y = '', title = 'group 2') +
  # scale_fill_manual(values = brewer.pal(4, "Set1")[1:4]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15, angle = 60, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

mod_go3 = data.frame(name = c('intracellular \n protein transport', 'regulation of cell cycle',
                              'cell cycle', 'muscle development', 'mitotic cell cycle'),
                    neg_logP = c(10.2, 10.2, 9.5, 9.4, 8.5))
mod_go3$name = factor(mod_go3$name, levels= mod_go3$name)
ggplot(data = mod_go3, aes(x = name, y = neg_logP)) +
  geom_bar(stat="identity", width = 0.7, fill = brewer.pal(4, 'Accent')[2]) + 
  labs(x = '', y = '', title = 'group 3') +
  # scale_fill_manual(values = brewer.pal(4, "Set1")[1:4]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15, angle = 60, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

mod_go4 = data.frame(name = c('cellular \n defense response', 'positive regulation \n of interleukin production',
                              'regulation of \n interleukin production', 
                             'interleukin production', 'complement activation'),
                    neg_logP = c(5.8, 4.9, 4.8, 4.8, 4.8))
mod_go4$name = factor(mod_go4$name, levels= mod_go4$name)
ggplot(data = mod_go4, aes(x = name, y = neg_logP)) +
  geom_bar(stat="identity", width = 0.7, fill = brewer.pal(4, 'Set2')[4]) + 
  labs(x = '', y = '', title = 'group 4') +
  # scale_fill_manual(values = brewer.pal(4, "Set1")[1:4]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15, angle = 60, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="none")


