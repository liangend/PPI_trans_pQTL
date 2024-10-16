library(data.table)
library(ggplot2)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')

## trans-PCO pp4
# all_trait = list.files('pqtl/08_gwas_coloc/gwas_coloc_mcl/')
# pp4_mcl = c()
# for (i in all_trait) {
#   pp4_file = list.files(paste0('pqtl/08_gwas_coloc/gwas_coloc_mcl/', i))
#   pp4_all_i = c()
#   for (j in pp4_file) {
#     pp4_j = fread(paste0('pqtl/08_gwas_coloc/gwas_coloc_mcl/', i, '/', j))
#     pp4_all_i = rbind(pp4_all_i, pp4_j)
#   }
#   pp4_all_i$trait = i
#   pp4_mcl = rbind(pp4_mcl, pp4_all_i)
#   print(i)
# }
# 
# pp4_mcl = pp4_mcl[pp4_mcl$pp4 > 0.75, ]
# fwrite(pp4_mcl, 'pqtl/08_gwas_coloc/gwas_coloc_mcl/pp4_sig_all.txt', sep = '\t')

## UKB pqtl pp4
# pp4_pqtl = c()
# for (i in all_trait) {
#   pp4_file = list.files(paste0('pqtl/08_gwas_coloc/gwas_coloc_pqtl/', i))
#   pp4_all_i = c()
#   for (j in pp4_file) {
#     pp4_j = fread(paste0('pqtl/08_gwas_coloc/gwas_coloc_pqtl/', i, '/', j))
#     pp4_all_i = rbind(pp4_all_i, pp4_j)
#   }
#   pp4_all_i$trait = i
#   gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/',
#                         i,'.txt'))
#   pp4_all_i$gwas_chr = gwas_i$CHR[match(pp4_all_i$gwas_loci, gwas_i$loci)]
#   pp4_pqtl = rbind(pp4_pqtl, pp4_all_i)
#   print(i)
# }
# pp4_pqtl = pp4_pqtl[pp4_pqtl$pp4 > 0.75, ]
# fwrite(pp4_pqtl, 'pqtl/08_gwas_coloc/gwas_coloc_pqtl/pp4_sig_all.txt', sep = '\t')

pp4_mcl = fread('pqtl/08_gwas_coloc/gwas_coloc_mcl/pp4_sig_all.txt')
pp4_pqtl = fread('pqtl/08_gwas_coloc/gwas_coloc_pqtl/pp4_sig_all.txt')

all_trait = list.files('pqtl/08_gwas_coloc/gwas_coloc_mcl/')
all_trait = all_trait[!grepl('.txt', all_trait)]
immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA', 'asthma', "T1D")
blood_trait = c("BASO","BASO-P","EO","EO-P","HCT",
                "HGB","HLSR","HLSR-P","IRF","LYMPH",
                "LYMPH-P","MCH","MCHC","MCV","MONO","MONO-P","MPV","MRV",
                "MSCV","NEUT","NEUT-P","PDW","PLT","PLTC","RBC","RBC-W",
                "RET","RET-P","WBC")
other_trait = c("AD","BMI","HDL","LDL","PD","T2D","WHR",
                "bipolar","height","hypertension",
                "schizo","stroke","weight")

mcl_prot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')
pp4_pqtl_trans = pp4_pqtl[pp4_pqtl$prot %in% mcl_prot$file, ]
pp4_pqtl_trans$prot_chr = mcl_prot$chr[match(pp4_pqtl_trans$prot, mcl_prot$file)]
pp4_pqtl_trans = pp4_pqtl_trans[pp4_pqtl_trans$gwas_chr != pp4_pqtl_trans$prot_chr, ]

n_loci = c()
n_coloc_trans_pqtl = c()
n_coloc_mcl = c()
n_coloc_both = c()

for (i in 1:length(all_trait)) {
  gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', 
                        all_trait[i],'.txt'))
  gwas_i = gwas_i[gwas_i$P < 5e-8, ]  ## keep only genome-wide significant loci
  n_loci[i] = nrow(gwas_i)
  
  pp4_trans_pqtl_i = pp4_pqtl_trans[pp4_pqtl_trans$trait == all_trait[i], ]
  pp4_mcl_i = pp4_mcl[pp4_mcl$trait == all_trait[i], ]
  
  trans_pqtl_coloc_i = unique(pp4_trans_pqtl_i$gwas_loci)
  mcl_coloc_i = unique(pp4_mcl_i$gwas_loci)
  both_coloc_i = intersect(trans_pqtl_coloc_i, mcl_coloc_i)
  
  n_coloc_both[i] = sum(both_coloc_i %in% gwas_i$loci)
  n_coloc_trans_pqtl[i] = sum(setdiff(trans_pqtl_coloc_i, both_coloc_i) %in% gwas_i$loci)
  n_coloc_mcl[i] = sum(setdiff(mcl_coloc_i, both_coloc_i) %in% gwas_i$loci)
}

# coloc_tab = data.frame(trait = rep(all_trait, 4), 
#                        n_coloc = c(n_coloc_trans_pqtl, n_coloc_both, n_coloc_mcl, 
#                                    n_loci-n_coloc_both-n_coloc_trans_pqtl-n_coloc_mcl),
#                        coloc = rep(c('univar', 'both', 'mcl', 'none'),
#                                    each = length(all_trait)))
coloc_tab = data.frame(trait = all_trait, 
                       coloc_prop = (n_coloc_both + n_coloc_mcl) / n_loci)

immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA', 'asthma', "T1D")
white_blood_trait = c("BASO","BASO-P","EO","EO-P",
                "LYMPH", "LYMPH-P","MONO","MONO-P",
                "NEUT","NEUT-P", "WBC")
red_blood_trait = c("RBC","HCT","HGB","MCV","RBC-W","MCH","MCHC",
                    "RET","RET-P","MSCV","HLSR","HLSR-P","MRV","IRF")
plt_trait = c("PLT","PLTC","MPV","PDW")
cardio_trait = c("HDL","LDL","T2D","hypertension","stroke")

anthro_trait = c("BMI","WHR", "height", "weight")
nerv_trait = c("AD", "PD", "bipolar", "schizo")

coloc_tab$trait_cat = 'White blood cell'
coloc_tab$trait_cat[which(coloc_tab$trait %in% immune_trait)] = 'Immune'
coloc_tab$trait_cat[which(coloc_tab$trait %in% anthro_trait)] = 'Anthropometric'
coloc_tab$trait_cat[which(coloc_tab$trait %in% red_blood_trait)] = 'Red blood cell'
coloc_tab$trait_cat[which(coloc_tab$trait %in% plt_trait)] = 'Platelet'
coloc_tab$trait_cat[which(coloc_tab$trait %in% cardio_trait)] = 'Cardio'
coloc_tab$trait_cat[which(coloc_tab$trait %in% nerv_trait)] = 'Nerv'

coloc_tab$trait_cat = factor(coloc_tab$trait_cat, 
                             levels = c('Nerv', 'Anthropometric', 'Cardio', 'Immune', 
                                        'Platelet', 'White blood cell', 'Red blood cell'))

coloc_tab$trait[coloc_tab$trait == 'asthma'] = 'Asthma'
coloc_tab$trait[coloc_tab$trait == 'lupus'] = 'SLE'
coloc_tab$trait[coloc_tab$trait == 'hypertension'] = 'Hypertension'
coloc_tab$trait[coloc_tab$trait == 'schizo'] = 'Schizophrenia'
coloc_tab$trait[coloc_tab$trait == 'stroke'] = 'Stroke'
coloc_tab$trait[coloc_tab$trait == 'bipolar'] = 'Bipolar'
coloc_tab$trait[coloc_tab$trait == 'height'] = 'Height'
coloc_tab$trait[coloc_tab$trait == 'weight'] = 'Weight'

coloc_tab = coloc_tab[order(coloc_tab$trait_cat, coloc_tab$coloc_prop), ]
coloc_tab$trait = factor(coloc_tab$trait, levels = coloc_tab$trait)

ggplot(data = coloc_tab, 
       aes(x = trait, y = coloc_prop, fill = trait_cat)) +
  coord_flip() +
  geom_bar(stat="identity", width = 0.8) +
  labs(x = '', y = 'Proportion of coloc loci', title = '') +
  scale_fill_manual(values = brewer.pal(8, "Set3")[-2]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


## immune traits
coloc_immune = coloc_tab[coloc_tab$trait_cat == 'immune' & 
                           coloc_tab$coloc %in% c('both', 'mcl'), ]
ggplot(data = coloc_immune, 
       aes(x = reorder(trait, -n_coloc), 
           y = n_coloc, 
           fill = factor(coloc, levels = c('mcl','both')))) +
  geom_bar(stat="identity", width = 0.8) + # coord_flip() +
  labs(x = '', y = 'Number of coloc loci', fill = 'trans coloc', 
       title = 'Autoimmune diseases') +
  scale_fill_manual(values = brewer.pal(4, "Blues")[c(4,2)]) + 
  # scale_fill_brewer(palette = "Paired") +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 13),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

# coloc legend
coloc_immune$coloc[coloc_immune$coloc == 'mcl'] = 'Trans-PCO specific'
coloc_immune$coloc[coloc_immune$coloc == 'both'] = 'Shared'
ggplot(data = coloc_immune, 
       aes(x = trait, y = n_coloc, 
           fill = factor(coloc, levels = c('Trans-PCO specific', 'Shared')))) +
  geom_bar(stat="identity", width = 0.8) +
  labs(x = '', y = 'Number of coloc loci', fill = '', 
       title = 'Autoimmune diseases') +
  scale_fill_manual(values = brewer.pal(4, "Greys")[c(4,2)]) + 
  # scale_fill_brewer(palette = "Paired") +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


### blood traits
coloc_blood = coloc_tab[coloc_tab$trait_cat == 'blood' & 
                          coloc_tab$coloc %in% c('both', 'mcl'), ]
ggplot(data = coloc_blood, 
       aes(x = reorder(trait, -n_coloc), y = n_coloc, 
           fill = factor(coloc, levels = c('mcl', 'both')))) +
  geom_bar(stat="identity", width = 0.8) + # coord_flip() +
  labs(x = '', y = 'Number of coloc loci', fill = '', 
       title = 'Blood traits') +
  scale_fill_manual(values = brewer.pal(4, "Reds")[c(4,2)]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 13),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

### other traits
coloc_other = coloc_tab[coloc_tab$trait_cat == 'other' & 
                          coloc_tab$coloc %in% c('both', 'mcl'), ]
ggplot(data = coloc_other, 
       aes(x = reorder(trait, -n_coloc), y = n_coloc, 
           fill = factor(coloc, levels = c('mcl', 'both')))) +
  geom_bar(stat="identity", width = 0.8) + # coord_flip() +
  labs(x = '', y = 'Number of coloc loci', fill = '', 
       title = 'Other') +
  scale_fill_manual(values = brewer.pal(4, "Greens")[c(4,2)]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = 'black', size = 13),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


