library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(openxlsx)

setwd('/project/xuanyao/jinghui')
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
pco_all = as.data.frame(pco_all)
pco_pp4 = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4.txt')
pco_pp4 = as.data.frame(pco_pp4)

## keep coloc with gwas loci p < 5e-8
all_trait = colnames(pco_pp4)[-(1:2)]
pp4_sig = c()
for (i in all_trait) {
  gwas_loci_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', 
                             i, '.txt'))
  gwas_loci_i = gwas_loci_i[gwas_loci_i$P < 5e-8, ]
  pp4_i = pco_pp4[, c('loci', i)]
  colnames(pp4_i)[2] = 'pp4'
  # pp4_i = pp4_i[pp4_i$pp4 > 0.75, ]
  pp4_i = cbind(pp4_i, pco_all[match(pp4_i$loci, pco_all$loci), c('chr', 'bp_hg37', 'bp_hg38', 'mod')])
  pp4_keep = c()
  for (j in 1:nrow(pp4_i)) {
    chr_j = pp4_i$chr[j]
    if (i %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {  ## gwas with hg38 bp
      bp_j = pp4_i$bp_hg38[j]
    } else {
      bp_j = pp4_i$bp_hg37[j]
    }
    gwas_sub = gwas_loci_i[gwas_loci_i$CHR == chr_j, ]
    if_in_gwas = which(gwas_sub$loci_start <= bp_j & gwas_sub$loci_end >= bp_j)
    if (length(if_in_gwas) > 0) {
      pp4_keep = c(pp4_keep, j)
    }
  }
  pp4_i = pp4_i[pp4_keep, ]
  pp4_i$trait = i
  pp4_sig = rbind(pp4_sig, pp4_i)
  print(i)
}
fwrite(pp4_sig, 'pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_all_gwas_5e8.txt')


pco_pp4$n_coloc = apply(pco_pp4[,3:52], 1, function(x){sum(x > 0.75)})
pco_pp4$is_coloc = (pco_pp4$n_coloc > 0)
pco_pp4$is_nov = pco_all$is_nov[match(pco_pp4$loci, pco_all$loci)]
aggregate(is_coloc ~ is_nov, pco_pp4, mean)

pco_uniq_loci = fread('pqtl/04_fdr/ukb/pco_mcl_uniq_loci.txt')
pco_uniq_loci$n_coloc = pco_pp4$n_coloc[match(pco_uniq_loci$loci, pco_pp4$loci)]
pco_uniq_loci$is_nov = pco_all$is_nov[match(pco_uniq_loci$loci, pco_all$loci)]

## number of coloc per locus
n_coloc_tab = data.frame(n_coloc_per_locus = c(aggregate(n_coloc ~ is_nov, pco_pp4, mean)$n_coloc,
                                               aggregate(n_coloc ~ is_nov, pco_uniq_loci, mean)$n_coloc),
                         se = c(aggregate(n_coloc ~ is_nov, pco_pp4, sd)$n_coloc/
                                  sqrt(aggregate(n_coloc ~ is_nov, pco_pp4, sum)$n_coloc),
                                aggregate(n_coloc ~ is_nov, pco_uniq_loci, sd)$n_coloc/
                                  sqrt(aggregate(n_coloc ~ is_nov, pco_uniq_loci, sum)$n_coloc)),
                         is_nov = c(F, T, F, T),
                         group = rep(c('locus mod pair', 'unique locus'), each = 2))
ggplot(n_coloc_tab, aes(x = group, y = n_coloc_per_locus, fill = is_nov)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.4) +
  geom_errorbar(position=position_dodge(0.4), aes(ymin=n_coloc_per_locus-se, ymax=n_coloc_per_locus+se), width=.1) +
  labs(x = "", y = "# coloc per locus", fill = 'is novel', title = '') +
  scale_fill_brewer(palette="Paired") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## total number of coloc for each trait
n_coloc_all = c()
for (t in colnames(pco_pp4)[3:52]) {
  trait_i = pco_pp4[, c('is_nov', t)]
  trait_i = trait_i[trait_i[,2] > 0.75, ]
  n_coloc_all = c(n_coloc_all, unlist(table(trait_i$is_nov)))
}

n_coloc_all_tab = data.frame(n_coloc = n_coloc_all,
                                trait = rep(colnames(pco_pp4)[2:51], each = 2),
                                is_nov = names(n_coloc_all))

ggplot(n_coloc_all_tab, aes(x = reorder(trait, n_coloc), y = n_coloc, fill = factor(is_nov, levels = c(T, F)))) + 
  geom_bar(stat="identity", width = 0.5) +
  coord_flip() +
  labs(x = "", y = "# coloc", fill = 'is novel', title = '') +
  scale_fill_brewer(palette="Paired", direction=-1) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA', 'asthma')
pp4_immune = pco_pp4[, c('loci', 'is_nov', immune_trait)]
pp4_immune$n_coloc = apply(pp4_immune[,3:9], 1, function(x){sum(x > 0.75)})
pp4_immune = pp4_immune[pp4_immune$n_coloc > 0, ]


# find twas putative causal gene and module annotation for each immune trait
mcl_mod_annot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_annot.txt')
mcl_mod_annot$mod = paste0('mod', mcl_mod_annot$mod)
for (t in immune_trait) {
  trait_i = pp4_immune[, c('loci', t)]
  colnames(trait_i)[2] = 'pp4'
  trait_i = trait_i[trait_i$pp4 > 0.75, ]
  trait_i = cbind(trait_i, pco_all[match(trait_i$loci, pco_all$loci), -5])
  twas_i = read.xlsx(paste0('pqtl/08_gwas_coloc/TWAS_cGene/', t,'.xlsx'))
  twas_i = twas_i[twas_i$Tissue == 'Whole Blood', ]
  trait_i$nearest_gene_twas = (trait_i$nearest_gene %in% twas_i$Gene)
  ppi_prot_twas = c()
  for (i in 1:nrow(trait_i)) {
    prot_i = unlist(strsplit(trait_i$univar_trans_prot[i], ',', fixed = T))
    ppi_prot_twas[i] = paste(prot_i[prot_i %in% twas_i$Gene], collapse = ',')
  }
  trait_i$ppi_prot_twas = ppi_prot_twas
  trait_i = cbind(trait_i, mcl_mod_annot[match(trait_i$mod, mcl_mod_annot$mod) ,-(1:2)])
  fwrite(trait_i, paste0('pqtl/08_gwas_coloc/pco_mcl_coloc/pp4_immune/', 
                         t,'_full.txt'), sep = '\t')
  print(t)
}





