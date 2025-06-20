library(ggplot2)
library(data.table)
library(openxlsx)
library(RColorBrewer)
library(pheatmap)
library(CMplot)
setwd('/project/xuanyao/jinghui')
pco_all0 = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')

rm_loci_all = fread('pqtl/25_cell_prop_expr/rm_loci_spearman.txt')
rm_loci = rm_loci_all$loci
pco_all = pco_all0[!pco_all0$loci %in% rm_loci, ]

### PPI trans-pQTL coloc with GWAS
pp4_all = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_all_gwas_5e8.txt')
pp4_sig = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_sig_gwas_5e8.txt')
pp4_all = pp4_all[pp4_all$loci %in% pco_all$loci, ]
pp4_sig = pp4_sig[pp4_sig$loci %in% pco_all$loci, ]

pp4_all$pco_p = pco_all$pco_p[match(pp4_all$loci, pco_all$loci)]

all_trait = unique(pp4_sig$trait)

### GWAS coloc with PPI trans-pQTL
pp4_gwas = fread('pqtl/08_gwas_coloc/gwas_coloc_mcl/pp4_all_gwas_5e8.txt')
pp4_gwas = pp4_gwas[pp4_gwas$min_pco_p < 5e-8, ]

# remove trans-pQTLs due to cell-composition eff
mcl_pqtl_rm = pco_all0[!(pco_all0$loci %in% pco_all$loci), ]
pp4_gwas_keep = as.data.frame(pp4_gwas)
pp4_gwas_keep$mcl_mod = paste0('mod', pp4_gwas_keep$mcl_mod)

mod_go = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_annot.txt')
mod_go$mod = paste0('mod', mod_go$mod)
# pp4_sig = cbind(pp4_sig, mod_go[match(pp4_sig$mod, mod_go$mod), -1])
# fwrite(pp4_sig, 'pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_sig_rm_blood_prop.txt', sep = '\t')

mod_go = mod_go[, c(1,4:6,10:12)]
mod_coloc = mod_go[, c('mod')]
mod_coloc = as.data.frame(mod_coloc)
n_loci = c()
n_coloc_mcl = c()
keep_trait = c()
## find the most significant coloc for each mod-trait pair
for (i in 1:length(all_trait)) {
  pp4_all_i = pp4_all[pp4_all$trait == all_trait[i], ]
  if (all_trait[i] %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {  ## gwas with hg38 bp
    pp4_all_i$POS = pp4_all_i$bp_hg38
  } else {
    pp4_all_i$POS = pp4_all_i$bp_hg37
  }
  
  gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', all_trait[i], '.txt'))
  gwas_sig_i = gwas_i[gwas_i$P < 5e-8, ]
  gwas_i = gwas_i[order(gwas_i$P)[1:200], ]  ## keep top GWAS loci only
  
  if (nrow(gwas_sig_i) >= 10) {
    keep_trait = c(keep_trait, all_trait[i])
    n_loci = c(n_loci, nrow(gwas_sig_i))
    keep_pp4_i = c()
    for (j in 1:nrow(gwas_i)) {
      keep_j = which(pp4_all_i$chr == gwas_i$CHR[j] &
                       pp4_all_i$POS >= gwas_i$loci_start[j] &
                       pp4_all_i$POS <= gwas_i$loci_end[j])
      keep_pp4_i = c(keep_pp4_i, keep_j)
    }
    pp4_all_sub = pp4_all_i[keep_pp4_i, ]
    pp4_max = aggregate(pp4 ~ mod, data = pp4_all_sub, max)
    
    pp4_gwas_i = pp4_gwas_keep[pp4_gwas_keep$trait == all_trait[i], ] 
    n_coloc_mcl = c(n_coloc_mcl, length(unique(pp4_gwas_i$gwas_loci)))
    mod_coloc = cbind(mod_coloc, pp4_max$pp4[match(mod_coloc$mod, pp4_max$mod)])
  }
  
}
coloc_tab = data.frame(trait = keep_trait, n_loci = n_loci, 
                       n_coloc = n_coloc_mcl, coloc_prop = n_coloc_mcl / n_loci)

white_blood_trait = c("BASO","BASO-P","EO","EO-P",
                      "LYMPH", "LYMPH-P","MONO","MONO-P",
                      "NEUT","NEUT-P", "WBC")
red_blood_trait = c("RBC","HCT","HGB","MCV","RBC-W","MCH","MCHC",
                    "RET","RET-P","MSCV","HLSR","HLSR-P","MRV","IRF")
plt_trait = c("PLT","PLTC","MPV","PDW")
immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA', 'asthma', "T1D")
cardio_trait = c("HDL","LDL","T2D","hypertension","stroke")
anthro_trait = c("BMI","WHR", "height", "weight")
nerv_trait = c("AD", "PD", "schizo")

coloc_tab$trait_cat = 'White blood cell'
coloc_tab$trait_cat[which(coloc_tab$trait %in% red_blood_trait)] = 'Red blood cell'
coloc_tab$trait_cat[which(coloc_tab$trait %in% plt_trait)] = 'Platelet'
coloc_tab$trait_cat[which(coloc_tab$trait %in% immune_trait)] = 'Immune'
coloc_tab$trait_cat[which(coloc_tab$trait %in% anthro_trait)] = 'Anthropometric'
coloc_tab$trait_cat[which(coloc_tab$trait %in% cardio_trait)] = 'Cardio'
coloc_tab$trait_cat[which(coloc_tab$trait %in% nerv_trait)] = 'Nerv'

coloc_tab$trait_cat = factor(coloc_tab$trait_cat, 
                             levels = c('Nerv', 'Anthropometric', 'Cardio', 'Immune', 
                                        'Platelet', 'White blood cell', 'Red blood cell'))
# rename the trait
coloc_tab$trait[coloc_tab$trait == 'asthma'] = 'Asthma'
coloc_tab$trait[coloc_tab$trait == 'lupus'] = 'SLE'
coloc_tab$trait[coloc_tab$trait == 'hypertension'] = 'Hypertension'
coloc_tab$trait[coloc_tab$trait == 'schizo'] = 'SCZ'
coloc_tab$trait[coloc_tab$trait == 'height'] = 'Height'
coloc_tab$trait[coloc_tab$trait == 'weight'] = 'Weight'

coloc_tab = coloc_tab[order(coloc_tab$trait_cat, coloc_tab$coloc_prop), ]
coloc_tab$trait = factor(coloc_tab$trait, levels = coloc_tab$trait)
rm_trait = c('WBC', 'PLT', 'LYMPH', 'NEUT',
             'PLTC', 'LYMPH-P', 'NEUT-P', "BASO-P", "EO-P", "MONO-P",
             "RBC","HCT","HGB", 
             "MCV","RBC-W","MCH","MCHC",
             "RET","RET-P","HLSR-P","IRF")

coloc_tab = coloc_tab[!coloc_tab$trait %in% rm_trait, ]
ggplot(data = coloc_tab, 
       aes(x = trait, y = coloc_prop, fill = trait_cat)) +
  coord_flip() +
  geom_bar(stat="identity", width = 0.8) +
  labs(x = '', y = '', title = '') +
  scale_fill_manual(values = brewer.pal(8, "Set3")[-2]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 14),
        axis.text.y = element_text(colour = 'black', size = 14),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

heatmap_trait = keep_trait
heatmap_trait[(length(heatmap_trait) - 5):length(heatmap_trait)] = c('Asthma', 'Height', 'Hypertension', 
                         'SLE', 'SCZ', 'Weight')
colnames(mod_coloc)[2:ncol(mod_coloc)] = heatmap_trait
mod_coloc[is.na(mod_coloc)] = 0

mod_coloc_sub = mod_coloc[, -1]
mod_coloc_plot = t(as.matrix(mod_coloc_sub))
colnames(mod_coloc_plot) = mod_coloc$mod

mod_coloc_plot0 = mod_coloc_plot[rev(as.character(coloc_tab$trait)), ]
coloc_sum = colSums(mod_coloc_plot0)
mod_coloc_plot0 = mod_coloc_plot0[, -which(coloc_sum == 0)]

## group PPIs into n clusters
n_clust = 3
ppi_clust = hclust(dist(t(mod_coloc_plot0)))
mod_annot = data.frame(group = cutree(ppi_clust, k = n_clust))
mod_color = list(group = c(brewer.pal(n_clust, "Dark2")))

all_trait_heat = pheatmap(mod_coloc_plot0, 
                          col = colorRampPalette(c("white","#2A788EFF"))(100), 
                          legend = T, show_colnames = F, cluster_rows = F, 
                          # treeheight_col = 0, 
                          #annotation_col = mod_annot, 
                          # annotation_row = trait_annot, 
                          annotation_colors = mod_color,
                          fontsize = 8)




