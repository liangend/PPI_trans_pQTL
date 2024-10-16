library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
ukb_pqtl = read.xlsx('pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                     startRow = 5, cols = c(1:11, 20), sheet = 10)
dgn_beta_ukb_pqtl = fread('pqtl/12_beta_across_two_data/dgn_beta_ukb_pqtl_all.txt')
eqtlgen_z_ukb_pqtl = fread('pqtl/12_beta_across_two_data/eqtlgen_beta_ukb_pqtl_all.txt')

colnames(dgn_beta_ukb_pqtl) = c(colnames(ukb_pqtl)[1:11], 'cis_trans', 'corr_beta', 'corr_se', 
                                'corr_p', 'extr_var_id', 'extr_beta', 'extr_se', 'extr_p')
# hg38_gtf = fread('pqtl/00_ref/hg38.knownGene.gtf.gz')
# hg38_gtf = hg38_gtf[, c(1,4,5,3)]
# hg38_gtf = unique(hg38_gtf)
# hg38_gtf_3utr = hg38_gtf[hg38_gtf$V3 == '3UTR', ]
# hg38_gtf_5utr = hg38_gtf[hg38_gtf$V3 == '5UTR', ]
# hg38_gtf_3utr = hg38_gtf_3utr[order(hg38_gtf_3utr$V1, hg38_gtf_3utr$V4, hg38_gtf_3utr$V5), ]
# hg38_gtf_5utr = hg38_gtf_5utr[order(hg38_gtf_5utr$V1, hg38_gtf_5utr$V4, hg38_gtf_5utr$V4), ]
# 
# fwrite(hg38_gtf, 'pqtl/13_UTR_enrich/hg38_gtf.bed', sep = '\t', col.names = F)
# fwrite(hg38_gtf_3utr, 'pqtl/13_UTR_enrich/hg38_3utr.bed', sep = '\t', col.names = F)
# fwrite(hg38_gtf_5utr, 'pqtl/13_UTR_enrich/hg38_5utr.bed', sep = '\t', col.names = F)

## 3UTR and 5UTR regions after merging
hg38_gtf_3utr = fread('pqtl/13_UTR_enrich/hg38_3utr_merge.bed')
hg38_gtf_5utr = fread('pqtl/13_UTR_enrich/hg38_5utr_merge.bed')
prop_3utr = sum(hg38_gtf_3utr$V3 - hg38_gtf_3utr$V2 + 1) / 3088269832 # proportion of 3'UTR in whole genome
prop_5utr = sum(hg38_gtf_5utr$V3 - hg38_gtf_5utr$V2 + 1) / 3088269832 # proportion of 5'UTR in whole genome

## all UKB pQTLs
pqtl_bed = data.frame(chr = ukb_pqtl$CHROM, start = ukb_pqtl$`GENPOS.(hg38)` - 1,
                      end = ukb_pqtl$`GENPOS.(hg38)`, cis_trans = ukb_pqtl$`cis/trans`)
pqtl_bed = unique(pqtl_bed)
pqtl_bed$chr = paste0('chr', pqtl_bed$chr)
#fwrite(pqtl_bed, 'pqtl/13_UTR_enrich/ukb_pqtl_hg38.bed', sep = '\t', col.names = F)

## pQTL with eQTL information
dgn_beta_ukb_pqtl = dgn_beta_ukb_pqtl[dgn_beta_ukb_pqtl$corr_se > 0, ]

## Loci are both cis-pQTL and cis-eQTL (based on DGN)
ukb_pqtl_w_eqtl = dgn_beta_ukb_pqtl[dgn_beta_ukb_pqtl$corr_p < 1e-5 & 
                                      dgn_beta_ukb_pqtl$cis_trans == 'cis', ]
pqtl_w_eqtl_bed = data.frame(chr = ukb_pqtl_w_eqtl$CHROM, start = ukb_pqtl_w_eqtl$`GENPOS.(hg38)` - 1,
                             end = ukb_pqtl_w_eqtl$`GENPOS.(hg38)`, 
                             cis_trans = ukb_pqtl_w_eqtl$cis_trans)
pqtl_w_eqtl_bed = unique(pqtl_w_eqtl_bed)
pqtl_w_eqtl_bed$chr = paste0('chr', pqtl_w_eqtl_bed$chr)
#fwrite(pqtl_w_eqtl_bed, 'pqtl/13_UTR_enrich/ukb_pqtl_w_eqtl_hg38.bed', sep = '\t', col.names = F)

## Loci are cis-pQTL but not cis-eQTL (based on DGN)
ukb_pqtl_wo_eqtl = dgn_beta_ukb_pqtl[dgn_beta_ukb_pqtl$corr_p >= 1e-5 & 
                                      dgn_beta_ukb_pqtl$cis_trans == 'cis', ]
pqtl_wo_eqtl_bed = data.frame(chr = ukb_pqtl_wo_eqtl$CHROM, start = ukb_pqtl_wo_eqtl$`GENPOS.(hg38)` - 1,
                             end = ukb_pqtl_wo_eqtl$`GENPOS.(hg38)`, 
                             cis_trans = ukb_pqtl_wo_eqtl$cis_trans)
pqtl_wo_eqtl_bed = unique(pqtl_wo_eqtl_bed)
pqtl_wo_eqtl_bed$chr = paste0('chr', pqtl_wo_eqtl_bed$chr)
#fwrite(pqtl_wo_eqtl_bed, 'pqtl/13_UTR_enrich/ukb_pqtl_wo_eqtl_hg38.bed', sep = '\t', col.names = F)

## Loci are both cis-pQTL and cis-eQTL (based on eQTLgen)
ukb_pqtl_w_eqtl2 = eqtlgen_z_ukb_pqtl[eqtlgen_z_ukb_pqtl$eqtlgen_bon_pval < 0.05, ]
pqtl_w_eqtl_bed2 = data.frame(chr = ukb_pqtl_w_eqtl2$CHROM, start = ukb_pqtl_w_eqtl2$`GENPOS.(hg38)` - 1,
                             end = ukb_pqtl_w_eqtl2$`GENPOS.(hg38)`)
pqtl_w_eqtl_bed2 = unique(pqtl_w_eqtl_bed2)
pqtl_w_eqtl_bed2$chr = paste0('chr', pqtl_w_eqtl_bed2$chr)
#fwrite(pqtl_w_eqtl_bed2, 'pqtl/13_UTR_enrich/ukb_pqtl_w_eqtl2_hg38.bed', sep = '\t', col.names = F)

## Loci are cis-pQTL but not cis-eQTL (based on DGN)
ukb_pqtl_wo_eqtl2 = eqtlgen_z_ukb_pqtl[eqtlgen_z_ukb_pqtl$eqtlgen_bon_pval >= 0.05, ]
pqtl_wo_eqtl_bed2 = data.frame(chr = ukb_pqtl_wo_eqtl2$CHROM, start = ukb_pqtl_wo_eqtl2$`GENPOS.(hg38)` - 1,
                              end = ukb_pqtl_wo_eqtl2$`GENPOS.(hg38)`)
pqtl_wo_eqtl_bed2 = unique(pqtl_wo_eqtl_bed2)
pqtl_wo_eqtl_bed2$chr = paste0('chr', pqtl_wo_eqtl_bed2$chr)
#fwrite(pqtl_wo_eqtl_bed2, 'pqtl/13_UTR_enrich/ukb_pqtl_wo_eqtl2_hg38.bed', sep = '\t', col.names = F)

## Blood eQTL (based on GTEx)
gtex_eqtl = fread('gtex/00_ref/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz')
eqtl_bed = data.frame(chr = gtex_eqtl$chr,
                      start = gtex_eqtl$variant_pos - 1,
                      end = gtex_eqtl$variant_pos)
eqtl_bed = unique(eqtl_bed)
#fwrite(eqtl_bed, 'pqtl/13_UTR_enrich/gtex_eqtl.bed', sep = '\t', col.names = F)

## Blood eQTL (based on eQTLgen)
eqtlgen_annot = fread('/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
eqtlgen_annot_sub = eqtlgen_annot[, c(1,3,4,8)]
eqtlgen_annot_sub = eqtlgen_annot_sub[!duplicated(eqtlgen_annot_sub[, c(1,4)]), ]
lead_snp = aggregate(Pvalue ~ Gene, data = eqtlgen_annot, min)
lead_snp = merge(lead_snp, eqtlgen_annot_sub, by = c('Pvalue', 'Gene'))
eqtlgen_bed = data.frame(chr = lead_snp$SNPChr,
                         start = lead_snp$SNPPos - 1,
                         end = lead_snp$SNPPos)
eqtlgen_bed$chr = paste0('chr', eqtlgen_bed$chr)
#fwrite(eqtlgen_bed, 'pqtl/13_UTR_enrich/eqtlgen_eqtl.bed', sep = '\t', col.names = F)

## annotation of all pQTLs
pqtl_annot = fread('pqtl/13_UTR_enrich/ukb_pqtl_annot.bed')
pqtl_3UTR = pqtl_annot[pqtl_annot$V8 == '3UTR', 1:4]
all_pqtl_3utr = nrow(unique(pqtl_3UTR)) / nrow(pqtl_bed) / prop_3utr
cis_pqtl_3utr = nrow(unique(pqtl_3UTR[pqtl_3UTR$V4 == 'cis', ])) / 
  nrow(pqtl_bed[pqtl_bed$cis_trans == 'cis', ]) / prop_3utr
trans_pqtl_3utr = nrow(unique(pqtl_3UTR[pqtl_3UTR$V4 == 'trans', ])) / 
  nrow(pqtl_bed[pqtl_bed$cis_trans == 'trans', ]) / prop_3utr

## annotation of cis-pQTLs that are also cis-eQTLs (DGN)
pqtl_w_eqtl_annot = fread('pqtl/13_UTR_enrich/ukb_pqtl_w_eqtl_annot.bed')
pqtl_w_eqtl_3UTR = pqtl_w_eqtl_annot[pqtl_w_eqtl_annot$V8 == '3UTR', 1:4]
cis_pqtl_eqtl_3utr = nrow(unique(pqtl_w_eqtl_3UTR)) / nrow(pqtl_w_eqtl_bed) / prop_3utr

## annotation of cis-pQTLs that are not cis-eQTLs (DGN)
pqtl_wo_eqtl_annot = fread('pqtl/13_UTR_enrich/ukb_pqtl_wo_eqtl_annot.bed')
pqtl_wo_eqtl_3UTR = pqtl_wo_eqtl_annot[pqtl_wo_eqtl_annot$V8 == '3UTR', 1:4]
cis_pqtl_no_eqtl_3utr = nrow(unique(pqtl_wo_eqtl_3UTR)) / nrow(pqtl_wo_eqtl_bed) / prop_3utr

## annotation of cis-pQTLs that are also cis-eQTLs (eQTLgen)
pqtl_w_eqtl_annot2 = fread('pqtl/13_UTR_enrich/ukb_pqtl_w_eqtl2_annot.bed')
pqtl_w_eqtl2_3UTR = pqtl_w_eqtl_annot2[pqtl_w_eqtl_annot2$V7 == '3UTR', 1:3]
cis_pqtl_eqtl2_3utr = nrow(unique(pqtl_w_eqtl2_3UTR)) / nrow(pqtl_w_eqtl_bed2) / prop_3utr

## annotation of cis-pQTLs that are not cis-eQTLs (eQTLgen)
pqtl_wo_eqtl_annot2 = fread('pqtl/13_UTR_enrich/ukb_pqtl_wo_eqtl2_annot.bed')
pqtl_wo_eqtl2_3UTR = pqtl_wo_eqtl_annot2[pqtl_wo_eqtl_annot2$V7 == '3UTR', 1:3]
cis_pqtl_no_eqtl2_3utr = nrow(unique(pqtl_wo_eqtl2_3UTR)) / nrow(pqtl_wo_eqtl_bed2) / prop_3utr

## annotation of GTEx cis-eQTLs
eqtl_annot = fread('pqtl/13_UTR_enrich/gtex_eqtl_annot.bed')
eqtl_3UTR = eqtl_annot[eqtl_annot$V7 == '3UTR', 1:3]
cis_eqtl_3utr = nrow(unique(eqtl_3UTR)) / nrow(eqtl_bed) / prop_3utr

## annotation of eQTLgen cis-eQTLs
eqtl_annot2 = fread('pqtl/13_UTR_enrich/eqtlgen_eqtl_annot.bed')
eqtl2_3UTR = eqtl_annot2[eqtl_annot2$V7 == '3UTR', 1:3]
cis_eqtl2_3utr = nrow(unique(eqtl2_3UTR)) / nrow(eqtlgen_bed) / prop_3utr

dat_plot_dgn = data.frame(enrich = c(all_pqtl_3utr, cis_pqtl_3utr,
                                 trans_pqtl_3utr, cis_pqtl_eqtl_3utr, 
                                 cis_pqtl_no_eqtl_3utr,
                                 cis_eqtl_3utr),
                      n_snp = c(nrow(pqtl_bed), 
                                nrow(pqtl_bed[pqtl_bed$cis_trans == 'cis', ]),
                                nrow(pqtl_bed[pqtl_bed$cis_trans == 'trans', ]),
                                nrow(pqtl_w_eqtl_bed),
                                nrow(pqtl_wo_eqtl_bed),
                                nrow(eqtl_bed)), 
                      group = c('all pQTL', 'cis-pQTL', 'trans-pQTL',
                                'cis-pQTL & cis-eQTL', 'cis-pQTL not cis-eQTL',
                                'GTEx blood cis-eQTL'))
dat_plot_eqtlgen = data.frame(enrich = c(all_pqtl_3utr, cis_pqtl_3utr,
                                     trans_pqtl_3utr, cis_pqtl_eqtl2_3utr, 
                                     cis_pqtl_no_eqtl2_3utr,
                                     cis_eqtl2_3utr),
                          n_snp = c(nrow(pqtl_bed), 
                                    nrow(pqtl_bed[pqtl_bed$cis_trans == 'cis', ]),
                                    nrow(pqtl_bed[pqtl_bed$cis_trans == 'trans', ]),
                                    nrow(pqtl_w_eqtl_bed2),
                                    nrow(pqtl_wo_eqtl_bed2),
                                    nrow(eqtlgen_bed)), 
                          group = c('all pQTL', 'cis-pQTL', 'trans-pQTL',
                                    'cis-pQTL & cis-eQTL', 'cis-pQTL not cis-eQTL',
                                    'eQTLgen cis-eQTL'))

p1 = ggplot(dat_plot_dgn, aes(x=reorder(group, enrich), y=enrich, label = n_snp)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue', width = 0.5) +
  labs(title = "3'UTR enrich based on DGN", x = '', y = 'enrichment') + coord_flip() +
  geom_text(hjust=-0.2, size = 4.5) + ylim(c(0,7)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


p2 = ggplot(dat_plot_eqtlgen, aes(x=reorder(group, enrich), y=enrich, label = n_snp)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue', width = 0.5) +
  labs(title = "3'UTR enrich based on eQTLgen", x = '', y = 'enrichment') + coord_flip() +
  geom_text(hjust=-0.2, size = 4.5) + ylim(c(0,7)) + 
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, nrow = 2)

