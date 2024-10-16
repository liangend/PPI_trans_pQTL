library(data.table)
library(ggplot2)
library(ggrepel)
setwd('/project/xuanyao/jinghui')
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
ukb_qtl = fread('pqtl/04_fdr/ukb/ukb_loci_bonP_500kb_nonoverlap.txt')
ukb_qtl = ukb_qtl[ukb_qtl$CHROM < 23, ]
ukb_qtl$CHROM = paste0('chr', ukb_qtl$CHROM)
ukb_qtl = na.omit(ukb_qtl)
gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta_sub = gene_meta_hg37[gene_meta_hg37$gene_type == 'protein_coding', ]
gene_meta_sub = gene_meta_sub[!grepl('gene_status', gene_meta_sub$gene_name), ]
gene_meta_sub = gene_meta_sub[!grepl('ENSG', gene_meta_sub$gene_name), ]

near_gene = c()
for (i in 1:nrow(ukb_qtl)) {
  chr_i = ukb_qtl$CHROM[i]
  gene_sub = gene_meta_sub[gene_meta_sub$chr == chr_i, ]
  near_gene[i] = gene_sub$gene_name[which.min(abs(gene_sub$start - ukb_qtl$bp_hg19[i]))]
}
ukb_qtl$near_gene = near_gene
ukb_trans = ukb_qtl[ukb_qtl$cis_trans == 'trans', ]
fwrite(as.data.frame(unique(ukb_trans$near_gene)), 'pqtl/ukb_trans_near_gene.txt', col.names = F)

## GO term of near genes of PCO and univaritate trans-pQTLs
mod_go = data.frame(go_term = c('GO:MF', 'GO:MF',
                                'GO:BP', 'GO:BP',
                                'GO:CC', 'GO:CC',
                                'KEGG', 'KEGG'), 
                    name = c('coreceptor activity', 'virus receptor activity', 
                             'regulation of lymphocyte proliferation', 'T cell proliferation',
                             'external side of plasma membrane', 'protein complex involved in cell adhesion',
                             'Autoimmune thyroid disease', 'Rheumatoid arthritis'),
                    neg_logP = -log10(c(1.6e-6, 5.8e-4,
                                        1.1e-6, 1.1e-6,
                                        1.8e-6, 1.8e-6,
                                        1.7e-8, 9.2e-8)))

mod_go$name = factor(mod_go$name, levels= rev(mod_go$name))
mod_go$go_term = factor(mod_go$go_term, levels=unique(mod_go$go_term))

ggplot(data = mod_go, aes(x = name, y = neg_logP, fill = go_term)) +
  geom_bar(stat="identity", width = 0.7) + coord_flip() +
  labs(x = '', y = '-log10(p)', fill = '', 
       title = paste0('Function enrichment of module ', sub('mod', '', plot_coloc$mod))) +
  scale_fill_manual(values = brewer.pal(4, "Set1")[1:4]) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="none")


### nearest genes of trans-pQTLs
ukb_trans_near_gene = as.data.frame(sort(table(ukb_trans$near_gene), decreasing = T)[1:30])
colnames(ukb_trans_near_gene) = c('gene', 'n_univar_qtl')
pco_near_gene = as.data.frame(sort(table(pco_all$nearest_gene), decreasing = T)[1:30])
colnames(pco_near_gene) = c('gene', 'n_pco_qtl')

near_gene_tab = merge(ukb_trans_near_gene, pco_near_gene, by = 'gene')
near_gene_tab$chr = gene_meta_sub$chr[match(near_gene_tab$gene, gene_meta_sub$gene_name)]
near_gene_tab$gene = paste0(near_gene_tab$chr, ': ', near_gene_tab$gene)
ggplot(near_gene_tab, aes(x = n_univar_qtl, y = n_pco_qtl))+
  geom_point(color = 'steelblue') + geom_text_repel(aes(label = gene), size = 3.5) +
  labs(x = "# univariate trans-pQTL", y = "# PCO trans-pQTL", title = '# pQTLs near gene') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### nearest genes of coloc trans-pQTLs with immune diseases
immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA')
pp4_univar = fread('pqtl/08_gwas_coloc/ukb_coloc/ukb_nonoverlap_pp4.txt')
pp4_univar = pp4_univar[pp4_univar$cis_trans == 'trans', ]
pp4_univar = as.data.frame(pp4_univar)
pp4_univar$nearest_gene = ukb_qtl$near_gene[match(pp4_univar$loci, ukb_qtl$loci)]
all_trait = colnames(pp4_univar)[6:55]
immune_univar = c()
for (i in immune_trait) {
  immune_i = pp4_univar[, c('loci', 'chr', 'prot', 'nearest_gene', i)]
  colnames(immune_i)[5] = 'pp4'
  immune_i = immune_i[immune_i$pp4 > 0.75, ]
  immune_i$trait = i
  immune_univar = rbind(immune_univar, immune_i)
}

pp4_mcl = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/pco_mcl_pp4.txt')
pp4_mcl$nearest_gene = pco_all$nearest_gene[match(pp4_mcl$loci, pco_all$loci)]
pp4_mcl = as.data.frame(pp4_mcl)
immune_mcl = c()
for (i in immune_trait) {
  immune_i = pp4_mcl[, c('loci', 'nearest_gene', i)]
  colnames(immune_i)[3] = 'pp4'
  immune_i = immune_i[immune_i$pp4 > 0.75, ]
  immune_i$trait = i
  immune_mcl = rbind(immune_mcl, immune_i)
}
immune_mcl$is_nov = pco_all$is_nov[match(immune_mcl$loci, pco_all$loci)]

immune_univar_near_gene = as.data.frame(sort(table(immune_univar$nearest_gene), decreasing = T)[1:30])
colnames(immune_univar_near_gene) = c('gene', 'n_univar_coloc')
immune_pco_near_gene = as.data.frame(sort(table(immune_mcl$nearest_gene), decreasing = T)[1:30])
colnames(immune_pco_near_gene) = c('gene', 'n_pco_coloc')

near_gene_coloc_tab = merge(immune_univar_near_gene, immune_pco_near_gene, by = 'gene')
near_gene_coloc_tab$chr = gene_meta_sub$chr[match(near_gene_coloc_tab$gene, gene_meta_sub$gene_name)]
near_gene_coloc_tab$gene = paste0(near_gene_coloc_tab$chr, ': ', near_gene_coloc_tab$gene)

ggplot(near_gene_coloc_tab, aes(x = n_univar_coloc, y = n_pco_coloc))+
  geom_point(color = 'steelblue') + geom_text_repel(aes(label = gene), size = 3.5) +
  labs(x = "# univariate trans-pQTL", y = "# PCO trans-pQTL", title = '# immune coloc near gene') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



