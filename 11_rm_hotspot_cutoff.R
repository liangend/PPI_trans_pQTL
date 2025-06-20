library(ggplot2)
library(data.table)
library(RColorBrewer)
library(dplyr)
setwd('/project/xuanyao/jinghui')

pco_all0 = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
cell_prop_trait = c('RBC', 'WBC', 'PLT', 'LYMPH', 'NEUT')

gene_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_hg37 = gene_hg37[gene_hg37$gene_type == 'protein_coding', ]
gene_hg37 = gene_hg37[!grepl('gene_status', gene_hg37$gene_name), ]
gene_hg37 = gene_hg37[!grepl('ENSG', gene_hg37$gene_name), ]
## genes close to GWAS loci of blood cell compositions
cell_prop_gene = vector("list", 5)
for (i in 1:length(cell_prop_trait)) {
  gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', 
                        cell_prop_trait[i], '.txt'))
  gwas_i = gwas_i[gwas_i$P < 5e-8, ]   # all gwas loci
  # gwas_i = gwas_i[order(gwas_i$P)[1:floor(0.1*nrow(gwas_i))], ]  # top 10% gwas loci
  for (j in 1:nrow(gwas_i)) {
    chr_j = paste0('chr', gwas_i$CHR[j])
    pos_j = gwas_i$POS[j]
    gene_meta_j = gene_hg37[gene_hg37$chr == chr_j, ]
    cell_prop_gene[[i]][j] = gene_meta_j$gene_name[which.min(abs(gene_meta_j$start - pos_j))]
  }
}
cell_prop_gene_all = unique(unlist(cell_prop_gene))

## CDF of trans-pQLTs for number of PPIs each locus affect
pqtl_hot_gene = as.data.frame(table(pco_all0$nearest_gene))

summ_hotspot = function(pqtl_hot_gene, genes, title, color, xintercept){
  ## if a trans-pQTL overlaps with cell composition gwas
  pqtl_hot_gene$group = pqtl_hot_gene$Var1 %in% genes
  pqtl_hot_gene1 = pqtl_hot_gene[pqtl_hot_gene$group, ]
  pqtl_hot_gene2 = pqtl_hot_gene[!pqtl_hot_gene$group, ]
  pqtl_hot_gene_count1 = as.data.frame(table(pqtl_hot_gene1$Freq))
  pqtl_hot_gene_count1$Var1 = as.numeric(as.character(pqtl_hot_gene_count1$Var1))
  pqtl_hot_gene_count1$n_qtl = pqtl_hot_gene_count1$Var1 * pqtl_hot_gene_count1$Freq
  total_qtl1 = c()
  for (i in 1:nrow(pqtl_hot_gene_count1)) {
    total_qtl1[i] = sum(pqtl_hot_gene_count1$n_qtl[1:i])
  }
  pqtl_hot_gene_count1$total_qtl = total_qtl1/sum(pqtl_hot_gene_count1$n_qtl)
  pqtl_hot_gene_count1$is_cell = T
  pqtl_hot_gene_count2 = as.data.frame(table(pqtl_hot_gene2$Freq))
  pqtl_hot_gene_count2$Var1 = as.numeric(as.character(pqtl_hot_gene_count2$Var1))
  pqtl_hot_gene_count2$n_qtl = pqtl_hot_gene_count2$Var1 * pqtl_hot_gene_count2$Freq
  total_qtl2 = c()
  for (i in 1:nrow(pqtl_hot_gene_count2)) {
    total_qtl2[i] = sum(pqtl_hot_gene_count2$n_qtl[1:i])
  }
  pqtl_hot_gene_count2$total_qtl = total_qtl2/sum(pqtl_hot_gene_count2$n_qtl)
  pqtl_hot_gene_count2$is_cell = F
  pqtl_hot_gene_count = rbind(pqtl_hot_gene_count1, pqtl_hot_gene_count2)
  
  ppi_cutoff = min(pqtl_hot_gene_count2$Var1[which(pqtl_hot_gene_count2$total_qtl > 0.88)])
  p = ggplot(pqtl_hot_gene_count, aes(x = Var1, y = total_qtl, color = is_cell)) +
    geom_point() +
    labs(x = 'Number of PPIs that independent loci affect', y = 'Proportion of PPI trans-pQTLs', 
         title = title,
         color = color) +
    geom_hline(yintercept = 0.88, linetype = 2) + 
    geom_vline(xintercept = xintercept, linetype = 2) + 
    scale_color_manual(values = c("black", "red")) + 
    theme(text = element_text(size=13, colour = 'black'),
          axis.text.x = element_text(colour = 'black', size = 12),
          axis.text.y = element_text(colour = 'black', size = 12),
          #axis.ticks.y = element_blank(),
          axis.line = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(p)
}
summ_hotspot(pqtl_hot_gene, cell_prop_gene_all, '', 'Overlap with cell composition GWAS', 50)


