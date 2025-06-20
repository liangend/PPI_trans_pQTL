library(data.table)
library(ggplot2)
library(RColorBrewer)
library(simtrait)
library(gridExtra)
setwd('/project/xuanyao/jinghui')

make_qq_df <- function(pvals, label) {
  pvals <- sort(pvals)
  n <- length(pvals)
  expected <- -log10(ppoints(n))      # Expected under null
  observed <- -log10(pvals)           # Observed
  data.frame(expected, observed, group = label)
}

cell_prop_trait = c('RBC', 'WBC', 'PLT', 'LYMPH', 'NEUT')

## PPI pqtl corresponding to cell-composition GWAS
top_n = 200
pco_cell_comp_p = fread('pqtl/23_cell_prop_gwas_in_trans/pco_in_mvp_cell_comp.txt')
rm_loci = fread('pqtl/25_cell_prop_expr/rm_loci_spearman.txt')
rm_loci = rm_loci$loci
pco_cell_comp_p$is_rm = pco_cell_comp_p$loci %in% rm_loci

pco_all = pco_cell_comp_p[pco_cell_comp_p$trait == 'RBC', ]
rank_p_tab = aggregate(pco_p ~ snp_hg38, data = pco_all, min)
rank_p_tab = rank_p_tab[order(rank_p_tab$pco_p), ]
top_pco_snp = pco_all$loci[match(rank_p_tab$snp_hg38[1:top_n], pco_all$snp_hg38)]

pco_cell_comp_p_top = pco_cell_comp_p[pco_cell_comp_p$loci %in% top_pco_snp, ]
pco_cell_comp_p_top = na.omit(pco_cell_comp_p_top)
pco_infl = aggregate(trait_p ~ trait, pco_cell_comp_p_top, pval_infl)
pco_infl$gwas = 'PPI trans-pQTLs'

pco_cell_comp_p_top$trait_p[pco_cell_comp_p_top$trait_p < 1e-100] = 1e-100
qq_pco = list()
for (i in 1:length(cell_prop_trait)) {
  qq_pco[[i]] = make_qq_df(pco_cell_comp_p_top$trait_p[which(pco_cell_comp_p_top$trait == cell_prop_trait[i])], 
                           'PPI trans-pQTLs')
}

# PPI pQTL after removing cell-composition eff
pco_cell_comp_p_sub = pco_cell_comp_p[!pco_cell_comp_p$is_rm, ]
pco_sub = pco_cell_comp_p_sub[pco_cell_comp_p_sub$trait == 'RBC', ]
rank_p_sub_tab = aggregate(pco_p ~ snp_hg38, data = pco_sub, min)
rank_p_sub_tab = rank_p_sub_tab[order(rank_p_sub_tab$pco_p), ]
top_pco_sub_snp = pco_sub$loci[match(rank_p_sub_tab$snp_hg38[1:top_n], pco_sub$snp_hg38)]

pco_cell_comp_p_sub = pco_cell_comp_p_sub[pco_cell_comp_p_sub$loci %in% top_pco_sub_snp, ]
pco_cell_comp_p_sub = na.omit(pco_cell_comp_p_sub)
pco_sub_infl = aggregate(trait_p ~ trait, pco_cell_comp_p_sub, pval_infl)
pco_sub_infl$gwas = 'Filtered PPI trans-pQTLs'

pco_cell_comp_p_sub$trait_p[pco_cell_comp_p_sub$trait_p < 1e-100] = 1e-100
qq_pco_sub = list()
for (i in 1:length(cell_prop_trait)) {
  qq_pco_sub[[i]] = make_qq_df(pco_cell_comp_p_sub$trait_p[which(pco_cell_comp_p_sub$trait == cell_prop_trait[i])], 
                           'Filtered PPI trans-pQTLs')
}

## gwas corresponding to cell-composition GWAS
gwas_list = c("AD",  "BMI", "CD", "HDL",  "IBD", "LDL",
              "MS", "PD", "RA", "T1D",  "T2D", "UC",
              "Asthma", "Height", "Hypertension", "SLE",
              "Schizo", "Weight", 'MONO', 'BASO', 'EO', 
              'RBC', 'WBC', 'PLT', 'LYMPH', 'NEUT')
gwas_infl = c()
gwas_qq = list()
n_loci = c()
for (i in 1:length(gwas_list)) {
  gwas_i = fread(paste0('pqtl/23_cell_prop_gwas_in_trans/gwas_in_mvp_cell_comp/',
                        gwas_list[i], '.txt'))
  n_loci[i] = nrow(gwas_i)/5
  gwas_sub = gwas_i[gwas_i$trait == 'RBC', ]
  top_snp = gwas_sub$SNP[order(gwas_sub$P)[1:min(top_n, nrow(gwas_sub))]]
  gwas_i = gwas_i[gwas_i$SNP %in% top_snp, ]
  gwas_i = na.omit(gwas_i)
  gwas_i$trait_p[gwas_i$trait_p < 1e-100] = 1e-100
  gwas_infl_i = aggregate(trait_p ~ trait, gwas_i, pval_infl)
  gwas_infl_i$gwas = gwas_list[i]
  gwas_qq_i = list()
  for (j in 1:length(cell_prop_trait)) {
    gwas_qq_i[[j]] = make_qq_df(gwas_i$trait_p[which(gwas_i$trait == cell_prop_trait[j])], 
                                 gwas_list[i])
  }
  gwas_qq[[i]] = gwas_qq_i
  gwas_infl = rbind(gwas_infl, gwas_infl_i)
}

qq_plt = list()
for (i in 1:5) {
  cell_i = cell_prop_trait[i]
  qq_plt_i = rbind(qq_pco[[i]], qq_pco_sub[[i]], 
                   gwas_qq[[which(gwas_list == cell_i)]][[i]], 
                   gwas_qq[[which(gwas_list == 'HDL')]][[i]], 
                   gwas_qq[[which(gwas_list == 'IBD')]][[i]], 
                   gwas_qq[[which(gwas_list == 'BMI')]][[i]],
                   gwas_qq[[which(gwas_list == 'Schizo')]][[i]])
  qq_plt_i$group = factor(qq_plt_i$group, levels = unique(qq_plt_i$group))
  qq_plt[[i]] = ggplot(qq_plt_i, aes(x = expected, y = observed, color = group)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(x = expression(paste("Expected -",log[10],"(p)")), 
         y = expression(paste("Observed -",log[10],"(p)")), 
         color = '', 
         title = paste0('Q-Q plot of ', cell_prop_trait[i], ' P values')) +
    ylim(c(0,20)) + xlim(c(0,1)) + 
    scale_color_manual(values = c(brewer.pal(8, "Dark2")[1:2], brewer.pal(9, "Greys")[3:8])) + 
    theme(text = element_text(size=16, colour = 'black'),
          axis.text.x = element_text(colour = 'black', size = 15),
          axis.text.y = element_text(colour = 'black', size = 15),
          #axis.ticks.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}
qq_plt[[1]]
grid.arrange(qq_plt[[2]],qq_plt[[3]],
             qq_plt[[4]],qq_plt[[5]],nrow = 2)

