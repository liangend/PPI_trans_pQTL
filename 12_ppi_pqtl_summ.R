library(data.table)
library(openxlsx)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(qvalue)
setwd('/project/xuanyao/jinghui')
pco_all0 = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')

####### remove trans-pQTLs that are possibly cell-type-prop effects
#### 1. remove pQTLs that have correlated univar beta and gene expr in each cell type
## immune cell gene expr
immune_cell_expr = fread('pqtl/25_cell_prop_expr/rna_immune_cell.tsv')
immune_cell_expr_wide = reshape(immune_cell_expr[,2:4], 
                                idvar = "Gene name", timevar = "Immune cell", 
                                direction = "wide")
colnames(immune_cell_expr_wide)[2:20] = sub('TPM.', '', 
                                            colnames(immune_cell_expr_wide)[2:20])
cell_expr = immune_cell_expr_wide[, -20]
cell_expr = as.data.frame(cell_expr)
keep_cell_type = c('neutrophil', 'classical monocyte', 'eosinophil', 'memory CD4 T-cell', 
                   'memory CD8 T-cell', 'naive CD4 T-cell', 'NK-cell')
cell_expr = cell_expr[, c('Gene name', keep_cell_type)]

## red blood cell gene expr
rbc_gene_expr = fread('pqtl/25_cell_prop_expr/rna_RBC_TPM.txt')

## platelet gene expr
gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.gtf')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
plt_gene_expr = fread('pqtl/25_cell_prop_expr/GSE68086_TEP_data_matrix.txt')
plt_gene_expr_sub = plt_gene_expr[, c(1,135:143)]
plt_gene_expr_sub$mean = apply(plt_gene_expr_sub[,-1], 1, mean)
plt_gene_expr_sub$gene = gene_meta$gene_name[match(plt_gene_expr_sub$V1, 
                                                   gene_meta$gene_id)]

## GTEx liver gene expr 
liver_gene_expr = fread('pqtl/25_cell_prop_expr/rna_liver.txt')

cell_expr$rbc = rbc_gene_expr$tpm[match(cell_expr$`Gene name`, rbc_gene_expr$gene_name)]
cell_expr$plt = plt_gene_expr_sub$mean[match(cell_expr$`Gene name`, plt_gene_expr_sub$gene)]
cell_expr$liver = liver_gene_expr$expr[match(cell_expr$`Gene name`, liver_gene_expr$gene_name)]

## standardize gene expr
library(limma)
library(RNOmni)
cell_expr_std = as.data.frame(cell_expr)
cell_expr_std = na.omit(cell_expr_std)
gene_rowsum = apply(cell_expr_std[,-1], 1, sum)
cell_expr_std = cell_expr_std[which(gene_rowsum > 0), ]

cell_expr_std2 = normalizeQuantiles(cell_expr_std[,-1])
cell_expr_std2 = apply(cell_expr_std2, 1, RankNorm)
cell_expr_std2 = t(cell_expr_std2)
cell_expr_std2 = as.data.frame(cell_expr_std2)
cell_expr_std2 = cbind(gene = cell_expr_std$`Gene name`, cell_expr_std2)
cell_expr_std2 = cell_expr_std2[cell_expr_std2$gene %in% ukb_prot$gene_name, ]

## spearman's cor between SNP eff on proteins and gene expr in each cell type
beta_all = fread('pqtl/25_cell_prop_expr/univar_eff_by_chr/beta_all.txt')
p_all = fread('pqtl/25_cell_prop_expr/univar_eff_by_chr/p_all.txt')

log10P_cutoff = 5 # cutoff to keep the univariate effect
loci_cor = c()
loci_cor_p = c()
loci_cor_cell = c() # which cell type has the most significant cor

for (i in 1:nrow(beta_all)) {
  loci_beta = unlist(beta_all[i, -1])
  loci_p = unlist(p_all[i, -1])
  beta_sub = loci_beta[which(loci_p > log10P_cutoff)]
  overlap_gene = intersect(names(beta_sub), cell_expr_std2$gene)
  if (length(overlap_gene) >= 10) {
    cell_expr_std_sub = cell_expr_std2[match(overlap_gene, cell_expr_std2$gene), ]
    beta_sub = beta_sub[overlap_gene]

    cor_i = apply(cell_expr_std_sub[, c(-1, -ncol(cell_expr_std_sub))], 2, 
                  function(x){cor(x, beta_sub, method = 'spearman', use = 'na.or.complete')})
    cor_p_i = apply(cell_expr_std_sub[, c(-1, -ncol(cell_expr_std_sub))], 2, 
                    function(x){cor.test(x, beta_sub, method = 'spearman')$p.value})
    
    index_i = which.min(cor_p_i)[1]
    loci_cor[i] = cor_i[index_i]
    loci_cor_p[i] = cor_p_i[index_i]
    loci_cor_cell[i] = names(cor_i)[index_i]
    
  } else{
    loci_cor[i] = NA
    loci_cor_p[i] = NA
    loci_cor_cell[i] = NA
  }
}

loci_eff_tab = data.frame(snp = beta_all$snp_hg38,
                          loci_cor = loci_cor,
                          loci_cor_p = loci_cor_p,
                          loci_cor_cell = loci_cor_cell)

pco_all0$snp_hg38 = paste0(sub('chr', '', pco_all0$chr), ':', 
                           pco_all0$bp_hg38)
pco_all0 = cbind(pco_all0, loci_eff_tab[match(pco_all0$snp_hg38, 
                                              loci_eff_tab$snp), -1])

#### 2. remove pQTLs that are trans-pQTL hotspots
hotspot_gene = names(which(table(pco_all0$nearest_gene) > 50))
rm_loci = pco_all0$loci[which(pco_all0$loci_cor_p < 0.001 | 
                                (pco_all0$nearest_gene %in% hotspot_gene))]
# fwrite(data.frame(loci = rm_loci), 'pqtl/25_cell_prop_expr/rm_loci_spearman.txt')

rm_loci = fread('pqtl/25_cell_prop_expr/rm_loci_spearman.txt')
rm_loci = rm_loci$loci
# pco_all0$is_rm = pco_all0$loci %in% rm_loci
# fwrite(pco_all0, 'pqtl/supp_file/pco_qtl.txt', sep = '\t')

pco_all = pco_all0[!pco_all0$loci %in% rm_loci, ]

## check significant cor enrichment near cell-composition GWAS
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
  gwas_i = gwas_i[gwas_i$P < 5e-8, ]
  for (j in 1:nrow(gwas_i)) {
    chr_j = paste0('chr', gwas_i$CHR[j])
    pos_j = gwas_i$POS[j]
    gene_meta_j = gene_hg37[gene_hg37$chr == chr_j, ]
    cell_prop_gene[[i]][j] = gene_meta_j$gene_name[which.min(abs(gene_meta_j$start - pos_j))]
  }
}
cell_prop_gene_all = unique(unlist(cell_prop_gene))
pco_all0$near_cell_prop_gwas = pco_all0$nearest_gene %in% cell_prop_gene_all
pco_all0$rm_cor = pco_all0$loci_cor_p < 0.001
pco_all0$rm_cor[is.na(pco_all0$rm_cor)] = F

table(pco_all0$near_cell_prop_gwas, pco_all0$rm_cor)

## analysis using trans-pQTLs after filtering cell-composition effects
pco_all = pco_all0[!pco_all0$loci %in% rm_loci, ]

mod_annot = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_annot.txt')
mod_source = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')

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


##### loci pleiotropy
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
  labs(x = "", y = "Proportion of loci affecting multiple proteins", color = '') +
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
backgroud_ppi = fread('pqtl/11_prot_complex/backgroud_ppi_mcl_pqtl_rm_cell_prop.txt')
null_ppi = backgroud_ppi$n_ppi / nrow(pco_all)

is_ppi = c()
for (j in 1:nrow(pco_all)) {
  is_ppi[j] = pco_all$nearest_gene[j] %in% mod_prot[[j]]
}
obs_ppi = sum(is_ppi)
ppi_enrich_all = obs_ppi / nrow(pco_all) / null_ppi

ppi_enrich = sum(pco_all$nearest_gene_is_ppi) / nrow(pco_all) / 
  null_ppi
ppi_enrich_nov = sum(pco_nov$nearest_gene_is_ppi) / 
  nrow(pco_nov) / null_ppi


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
  theme(text = element_text(size=18, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 18),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')











