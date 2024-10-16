library(data.table)
library(ggplot2)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')

## sqtl coloc
sqtl_pp4_dgn = fread('pqtl/16_pco_cis_coloc/mcl_dgn_sqtl_coloc.txt')
colnames(sqtl_pp4_dgn) = c('intron', 'pco_loci', 'pp4')
sqtl_coloc_dgn = sqtl_pp4_dgn[sqtl_pp4_dgn$pp4 > 0.75, ]
is_sqtl_dgn = rep(F, nrow(pco_all))
is_sqtl_dgn[which(pco_all$loci %in% unique(sqtl_coloc_dgn$pco_loci))] = T

sqtl_pp4_interval = fread('pqtl/16_pco_cis_coloc/mcl_interval_sqtl_coloc.txt')
colnames(sqtl_pp4_interval) = c('intron', 'pco_loci', 'pp4')
sqtl_coloc_interval = sqtl_pp4_interval[sqtl_pp4_interval$pp4 > 0.75, ]
is_sqtl_interval = rep(F, nrow(pco_all))
is_sqtl_interval[which(pco_all$loci %in% unique(sqtl_coloc_interval$pco_loci))] = T

## eqtl coloc
eqtl_pp4_gen = fread('pqtl/16_pco_cis_coloc/mcl_eqtlgen_eqtl_coloc.txt')
colnames(eqtl_pp4_gen) = c('gene', 'pco_loci', 'pp4')
eqtl_coloc_gen = eqtl_pp4_gen[eqtl_pp4_gen$pp4 > 0.75, ]
is_eqtl_gen = rep(F, nrow(pco_all))
is_eqtl_gen[which(pco_all$loci %in% unique(eqtl_coloc_gen$pco_loci))] = T

eqtl_pp4_interval = fread('pqtl/16_pco_cis_coloc/mcl_interval_eqtl_coloc.txt')
colnames(eqtl_pp4_interval) = c('gene', 'pco_loci', 'pp4')
eqtl_coloc_interval = eqtl_pp4_interval[eqtl_pp4_interval$pp4 > 0.75, ]
is_eqtl_interval = rep(F, nrow(pco_all))
is_eqtl_interval[which(pco_all$loci %in% unique(eqtl_coloc_interval$pco_loci))] = T

## cis-pqtl coloc
pqtl_pp4_all = fread('pqtl/16_pco_cis_coloc/mcl_cis_pqtl_coloc.txt')
colnames(pqtl_pp4_all) = c('prot', 'pco_loci', 'pp4')
pqtl_coloc = pqtl_pp4_all[pqtl_pp4_all$pp4 > 0.75, ]
length(unique(pqtl_coloc$pco_loci))
is_pqtl = rep(F, nrow(pco_all))
is_pqtl[which(pco_all$loci %in% unique(pqtl_coloc$pco_loci))] = T

pco_all$is_sqtl = (is_sqtl_interval | is_sqtl_dgn)
pco_all$is_eqtl = (is_eqtl_interval | is_eqtl_gen)
pco_all$is_pqtl = is_pqtl
# fwrite(pco_all, 'pqtl/04_fdr/ukb/pco_mcl_meta.txt', sep = '\t')

loci_mech = c()
for (i in 1:nrow(pco_all)) {
  is_eqtl_i = pco_all$is_eqtl[i]
  is_sqtl_i = pco_all$is_sqtl[i]
  is_pqtl_i = pco_all$is_pqtl[i]
  if (is_eqtl_i & is_sqtl_i & is_pqtl_i) {
    loci_mech[i] = 'eQTL & sQTL & pQTL'
  } else if (is_eqtl_i & !is_sqtl_i & !is_pqtl_i) {
    loci_mech[i] = 'eQTL only'
  } else if (is_sqtl_i & !is_eqtl_i & !is_pqtl_i) {
    loci_mech[i] = 'sQTL only'
  } else if (is_pqtl_i & !is_eqtl_i & !is_sqtl_i) {
    loci_mech[i] = 'pQTL only'
  } else if (is_eqtl_i & is_sqtl_i & !is_pqtl_i) {
    loci_mech[i] = 'eQTL & sQTL'
  } else if (is_pqtl_i & is_eqtl_i & !is_sqtl_i) {
    loci_mech[i] = 'eQTL & pQTL'
  } else if (is_pqtl_i & is_sqtl_i & !is_eqtl_i) {
    loci_mech[i] = 'pQTL & sQTL'
  } else {
    loci_mech[i] = 'none'
  }
}
pco_all$loci_mech = loci_mech

## define novel loci
pco_uniq_loci = fread('pqtl/04_fdr/ukb/pco_uniq_loci_nonoverlap.txt')
pco_all$is_nov = (pco_all$univar_max_logP < -log10(5e-8/2922))

is_nov_uniq = c()
is_sqtl_uniq = c()
is_eqtl_uniq = c()
is_pqtl_uniq = c()
for (i in 1:nrow(pco_uniq_loci)) {
  chr_i = pco_uniq_loci$chr[i]
  start_i = pco_uniq_loci$loci_start[i]
  end_i = pco_uniq_loci$loci_end[i]
  
  loci_i = which(pco_all$chr == chr_i & 
                   pco_all$loci_start >= start_i &
                   pco_all$loci_end <= end_i)
  is_nov_uniq[i] = (sum(pco_all$is_nov[loci_i]) > 0)
  is_sqtl_uniq[i] = (sum(pco_all$is_sqtl[loci_i]) > 0)
  is_eqtl_uniq[i] = (sum(pco_all$is_eqtl[loci_i]) > 0)
  is_pqtl_uniq[i] = (sum(pco_all$is_pqtl[loci_i]) > 0)
}
pco_uniq_loci$is_nov = is_nov_uniq
pco_uniq_loci$is_sqtl = is_sqtl_uniq
pco_uniq_loci$is_eqtl = is_eqtl_uniq
pco_uniq_loci$is_pqtl = is_pqtl_uniq

loci_mech_uniq = c()
for (i in 1:nrow(pco_uniq_loci)) {
  is_eqtl_i = pco_uniq_loci$is_eqtl[i]
  is_sqtl_i = pco_uniq_loci$is_sqtl[i]
  is_pqtl_i = pco_uniq_loci$is_pqtl[i]
  if (is_eqtl_i & is_sqtl_i & is_pqtl_i) {
    loci_mech_uniq[i] = 'all'
  } else if (is_eqtl_i & !is_sqtl_i & !is_pqtl_i) {
    loci_mech_uniq[i] = 'eQTL only'
  } else if (is_sqtl_i & !is_eqtl_i & !is_pqtl_i) {
    loci_mech_uniq[i] = 'sQTL only'
  } else if (is_pqtl_i & !is_eqtl_i & !is_sqtl_i) {
    loci_mech_uniq[i] = 'pQTL only'
  } else if (is_eqtl_i & is_sqtl_i & !is_pqtl_i) {
    loci_mech_uniq[i] = 'eQTL + sQTL'
  } else if (is_pqtl_i & is_eqtl_i & !is_sqtl_i) {
    loci_mech_uniq[i] = 'eQTL + pQTL'
  } else if (is_pqtl_i & is_sqtl_i & !is_eqtl_i) {
    loci_mech_uniq[i] = 'pQTL + sQTL'
  } else {
    loci_mech_uniq[i] = 'none'
  }
}
pco_uniq_loci$loci_mech = loci_mech_uniq
# fwrite(pco_uniq_loci, 'pqtl/04_fdr/ukb/pco_uniq_loci_nonoverlap.txt', sep = '\t')

## signal proportion
mech_tab = as.data.frame(table(pco_all$loci_mech))
colnames(mech_tab) = c('coloc', 'proportion')
mech_tab$proportion = mech_tab$proportion/ sum(mech_tab$proportion)
mech_tab$cat = ''
ggplot(mech_tab, aes(x = cat, y = proportion, 
                     fill = factor(coloc, levels = c('eQTL only', 'sQTL only', 'pQTL only', 
                                                     'eQTL & sQTL', 'eQTL & pQTL', 
                                                     'pQTL & sQTL', 'eQTL & sQTL & pQTL', 'none')))) + 
  geom_bar(stat="identity", color="black", width = 0.3) +
  scale_fill_manual(values = brewer.pal(8, "Dark2")) + 
  labs(x = "", fill = '', y = 'Proportion', title = 'PCO trans-pQTLs coloc with cis-QTL') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


mech_tab = as.data.frame(table(pco_all$loci_mech, pco_all$is_nov))
colnames(mech_tab) = c('coloc', 'is_nov', 'proportion')
mech_tab$proportion[1:8] = mech_tab$proportion[1:8]/sum(mech_tab$proportion[1:8])
mech_tab$proportion[9:16] = mech_tab$proportion[9:16]/sum(mech_tab$proportion[9:16])
ggplot(mech_tab, aes(x = is_nov, y = proportion, 
                     fill = factor(coloc, levels = c('eQTL only', 'sQTL only', 'pQTL only', 
                                                     'eQTL & sQTL', 'eQTL & pQTL', 
                                                     'pQTL & sQTL', 'eQTL & sQTL & pQTL', 'none')))) + 
  geom_bar(stat="identity", color="black", width = 0.3) +
  scale_fill_manual(values = brewer.pal(8, "Dark2")) + 
  labs(fill = '', y = 'Proportion', title = 'PCO trans-pQTLs coloc with cis-QTL') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


