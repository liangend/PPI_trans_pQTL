setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(RColorBrewer)
### Genome GTF
gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta = gene_meta[gene_meta$gene_type == 'protein_coding', ]
gene_meta = gene_meta[!grepl('gene_status', gene_meta$gene_name), ]
gene_meta = gene_meta[!grepl('ENSG', gene_meta$gene_name), ]
gene_meta$chr = as.numeric(sub('chr', '', gene_meta$chr))

ukb_qtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
ukb_qtl_at_interface = fread('pqtl/18_ppi_interface/ukb_qtl_at_interface.bed')
ukb_qtl$n_interface = ukb_qtl_at_interface$V4
ukb_trans = ukb_qtl[ukb_qtl$cis_trans == 'trans', ]
ukb_trans$near_gene_tss = gene_meta$start[match(ukb_trans$nearest_gene, gene_meta$gene_name)]
ukb_trans$near_gene_dist = abs(ukb_trans$near_gene_tss - ukb_trans$bp_hg19)
ukb_trans$maf = pmin(ukb_trans$A1FREQ, 1-ukb_trans$A1FREQ)
ukb_trans = ukb_trans[ukb_trans$coding == 'coding', ]


dgn_trans = fread('pqtl/18_ppi_interface/dgn_trans_org.txt')
dgn_trans = dgn_trans[dgn_trans$coding == 'coding', ]

ukb_rand = fread('pqtl/18_ppi_interface/ukb_background.txt')

## DGN enrichment
dgn_background = matrix(rep(0, nrow(dgn_trans) * 1000), nrow = 1000)
set.seed(101)
for (j in 1:nrow(dgn_trans)) {
  maf_j = dgn_trans$maf[j]
  dist_j = dgn_trans$near_gene_dist[j]
  ukb_sub = ukb_rand[which(abs(ukb_rand$maf - maf_j)/maf_j < 0.1 &
                             abs(ukb_rand$near_gene_dist - dist_j)/dist_j < 0.1), ]
  dgn_background[, j] = sample(ukb_sub$n_interface, 1000, replace = T)
}
dgn_background = (dgn_background > 0)
dng_background_is_interface = rowSums(dgn_background)
dgn_interface_enrich = sum(dgn_trans$n_interface > 0) / dng_background_is_interface
dgn_interface_enrich = dgn_interface_enrich[dgn_interface_enrich < 100]

## ukb enrichment
ukb_background = matrix(rep(0, nrow(ukb_trans) * 1000), nrow = 1000)
set.seed(101)
for (j in 1:nrow(ukb_trans)) {
  maf_j = ukb_trans$maf[j]
  dist_j = ukb_trans$near_gene_dist[j]
  ukb_sub = ukb_rand[which(abs(ukb_rand$maf - maf_j)/maf_j < 0.1 &
                             abs(ukb_rand$near_gene_dist - dist_j)/dist_j < 0.1), ]
  if (nrow(ukb_sub) == 0) {
    ukb_sub = ukb_rand[which(abs(ukb_rand$maf - maf_j)/maf_j < 0.2 &
                               abs(ukb_rand$near_gene_dist - dist_j)/dist_j < 0.2), ]
  }
  ukb_background[, j] = sample(ukb_sub$n_interface, 1000, replace = T)
}
ukb_background = (ukb_background > 0)
ukb_background_is_interface = rowSums(ukb_background)
ukb_interface_enrich = sum(ukb_trans$n_interface > 0) / ukb_background_is_interface


interface_enrich = data.frame(enrich = c(mean(dgn_interface_enrich), mean(ukb_interface_enrich)),
                    ci_low = c(quantile(dgn_interface_enrich, 0.025),
                               quantile(ukb_interface_enrich, 0.025)), 
                    ci_high = c(quantile(dgn_interface_enrich, 0.975), 
                                quantile(ukb_interface_enrich, 0.975)),
                    group = c('trans-eQTL', 'trans-pQTL'))

ggplot(interface_enrich, aes(x=group, y=enrich, fill = group)) +
  geom_bar(stat="identity", position=position_dodge(0.1), width = 0.5) +
  geom_errorbar(position=position_dodge(0.1), aes(ymin=ci_low, ymax=ci_high), width=.2) +
  labs(x = "", y = 'Enrichment at PPI interfaces', title = '', fill = '') + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  scale_fill_manual(values = brewer.pal(5,"Dark2")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')







