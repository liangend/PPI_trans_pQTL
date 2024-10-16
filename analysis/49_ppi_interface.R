setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(RColorBrewer)
#### PPI interface data orgnization
# ppi_interface = fread('pqtl/18_ppi_interface/all.bed.gz', sep = '|')
# colnames(ppi_interface) = 'info'
# prot_index = grep('track', ppi_interface$info)
# prot_info = ppi_interface$info[prot_index]
# 
# interface_coor = strsplit(ppi_interface$info, '\t', fixed = T)
# ppi_interface$chr = sapply(interface_coor, '[', 1)
# ppi_interface$start = sapply(interface_coor, '[', 2)
# ppi_interface$end = sapply(interface_coor, '[', 3)
# 
# prot_info = strsplit(prot_info, ' ', fixed = T)
# prot_pair = sapply(prot_info, '[', 6)
# prot1 = sapply(strsplit(prot_pair, '_', fixed = T), '[', 1)
# prot2 = sapply(strsplit(prot_pair, '_', fixed = T), '[', 2)
# ppi_source = sapply(prot_info, '[', 8)
# 
# org_interface = as.data.frame(ppi_interface)
# org_interface = org_interface[-prot_index, -1]
# 
# prot_index = c(prot_index, nrow(ppi_interface)+1)
# rep_time = prot_index[-1] - prot_index[-length(prot_index)] - 1
# org_interface$prot1 = rep(prot1, times = rep_time)
# org_interface$prot2 = rep(prot2, times = rep_time)
# org_interface$ppi_source = rep(ppi_source, times = rep_time)
# 
# fwrite(org_interface, 'pqtl/18_ppi_interface/insider_ppi_interface.txt', sep = '\t')


# ukb_qtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
# ukb_qtl_annot = fread('pqtl/17_snpEff/ukb_all_annot.vcf')
# annot_coding = function(x){
#   x = cbind(x, data.frame(annot = sapply(strsplit(x$INFO, '|', fixed = T), '[', 2), 
#                           coding = 'noncoding'))
#   x$coding[which(grepl('missense_variant', x$annot) | 
#                    grepl('synonymous_variant', x$annot) |
#                    grepl('start_lost', x$annot) |
#                    grepl('start_gained', x$annot) |
#                    grepl('stop_lost', x$annot) |
#                    grepl('stop_gained', x$annot) |
#                    grepl('inframe', x$annot) |
#                    grepl('frameshift', x$annot))] = 'coding'
#   return(x)
# }
# ukb_qtl_annot = annot_coding(ukb_qtl_annot)
# 
# ukb_qtl$annot = ukb_qtl_annot$annot
# ukb_qtl$coding = ukb_qtl_annot$coding
# fwrite(ukb_qtl, 'pqtl/04_fdr/ukb/ukb_500kb_loci.txt', sep = '\t')

### make bed file and use bedtools to find the overlap between qtl and ppi interface
# ukb_qtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
# ukb_qtl$CHROM = paste0('chr', ukb_qtl$CHROM)
# 
# ukb_qtl_bed = data.frame(chr = ukb_qtl$CHROM, start = ukb_qtl$bp_hg38, 
#                          end = ukb_qtl$bp_hg38)
# 
# ppi_interface = fread('pqtl/18_ppi_interface/insider_ppi_interface.txt.gz')
# ppi_interface_bed = ppi_interface[, 1:3]
# 
# dgn_trans = fread('pqtl/17_snpEff/dgn_trans_snp_annot.vcf')
# dgn_trans_bed = data.frame(chr = dgn_trans$`#CHROM`, start = dgn_trans$POS, 
#                            end = dgn_trans$POS)
# dgn_trans_bed$chr = paste0('chr', dgn_trans_bed$chr)

# ukb_rand = fread('pqtl/00_ref/ukb_geno/rand_sub.txt')
# ukb_rand_bed = data.frame(chr = paste0('chr', ukb_rand$CHROM), start = ukb_rand$bp_hg38, 
#                           end = ukb_rand$bp_hg38)
#   
# fwrite(ukb_qtl_bed, 'pqtl/18_ppi_interface/ukb_qtl.bed', sep = '\t', col.names = F)
# fwrite(ppi_interface_bed, 'pqtl/18_ppi_interface/ppi_interface.bed', sep = '\t', col.names = F)
# fwrite(dgn_trans_bed, 'pqtl/18_ppi_interface/dgn_trans_hg19.bed', sep = '\t', col.names = F)
# fwrite(ukb_rand_bed, 'pqtl/18_ppi_interface/ukb_rand.bed', sep = '\t', col.names = F)


#### dgn file information completion
# dgn_trans = fread('pqtl/17_snpEff/dgn_trans_snp_annot.vcf')
# dgn_trans = annot_coding(dgn_trans)
# dgn_trans_at_interface = fread('pqtl/18_ppi_interface/dgn_trans_at_interface.bed')
# dgn_trans$n_interface = dgn_trans_at_interface$V4
# dgn_trans$maf = dgn_trans_new$af
# colnames(dgn_trans)[1] = 'chr'
# dgn_trans = dgn_trans[, c(1:2,4:5,9:12)]
# near_gene = c()
# near_gene_dist = c()
# for (i in 1:nrow(dgn_trans)) {
#   chr_i = dgn_trans$chr[i]
#   bp_i = dgn_trans$POS[i]
#   gene_sub = gene_meta[gene_meta$chr == chr_i, ]
#   near_gene[i] = gene_sub$gene_name[which.min(abs(gene_sub$start - bp_i))]
#   near_gene_dist[i] = min(abs(gene_sub$start - bp_i))
# }
# dgn_trans$near_gene = near_gene
# dgn_trans$near_gene_dist = near_gene_dist
# fwrite(dgn_trans, 'pqtl/18_ppi_interface/dgn_trans_org.txt', sep = '\t')


### UKB random background of coding var at PPI interface
# ukb_rand = fread('pqtl/17_snpEff/ukb_rand_annot.vcf')
# ukb_rand = annot_coding(ukb_rand)
# ukb_rand = ukb_rand[, -8]
# 
# ukb_rand_bed = fread('pqtl/00_ref/ukb_geno/rand_sub.txt')
# ukb_rand$maf = pmin(ukb_rand_bed$A1FREQ, 1-ukb_rand_bed$A1FREQ)
# ukb_at_interface = fread('pqtl/18_ppi_interface/ukb_rand_at_interface.bed')
# ukb_rand$n_interface = ukb_at_interface$V4
# ukb_rand = ukb_rand[, c(1,2,8:11)]
# ukb_rand = ukb_rand[ukb_rand$coding == 'coding', ]
# colnames(ukb_rand)[1] = 'chr'
# near_gene = c()
# near_gene_dist = c()
# for (i in 1:nrow(ukb_rand)) {
#   chr_i = ukb_rand$chr[i]
#   bp_i = ukb_rand$POS[i]
#   gene_sub = gene_meta[gene_meta$chr == chr_i, ]
#   near_gene[i] = gene_sub$gene_name[which.min(abs(gene_sub$start - bp_i))]
#   near_gene_dist[i] = min(abs(gene_sub$start - bp_i))
# }
# ukb_rand$near_gene = near_gene
# ukb_rand$near_gene_dist = near_gene_dist
# fwrite(ukb_rand, 'pqtl/18_ppi_interface/ukb_background.txt', sep = '\t')


annot_coding = function(x){
  x = cbind(x, data.frame(annot = sapply(strsplit(x$INFO, '|', fixed = T), '[', 2),
                          coding = 'noncoding'))
  x$coding[which(grepl('missense_variant', x$annot) |
                   grepl('synonymous_variant', x$annot) |
                   grepl('start_lost', x$annot) |
                   grepl('start_gained', x$annot) |
                   grepl('stop_lost', x$annot) |
                   grepl('stop_gained', x$annot) |
                   grepl('inframe', x$annot) |
                   grepl('frameshift', x$annot))] = 'coding'
  return(x)
}

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
  scale_fill_manual(values = brewer.pal(5,"Set1")[c(3,5)]) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')







