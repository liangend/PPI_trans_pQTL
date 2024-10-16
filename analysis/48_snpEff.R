setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(boot)
dgn_cis_annot = fread('pqtl/17_snpEff/dgn_cis_annot.vcf')
dgn_trans_annot = fread('pqtl/17_snpEff/dgn_trans_snp_annot.vcf')
dgn_pco_annot = fread('pqtl/17_snpEff/dgn_trans_pco_annot.vcf')

ukb_all_annot = fread('pqtl/17_snpEff/ukb_all_annot.vcf')
ukb_pco_annot = fread('pqtl/17_snpEff/ukb_pco_annot.vcf')

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

dgn_cis_annot = annot_coding(dgn_cis_annot)
dgn_trans_annot = annot_coding(dgn_trans_annot)
dgn_pco_annot = annot_coding(dgn_pco_annot)

ukb_all_annot = annot_coding(ukb_all_annot)
ukb_cis_annot = ukb_all_annot[ukb_all_annot$FILTER == 'cis', ]
ukb_trans_annot = ukb_all_annot[ukb_all_annot$FILTER == 'trans', ]
ukb_pco_annot = annot_coding(ukb_pco_annot)


# pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
# pco_all$snp_annot = ukb_pco_annot$coding
# table(pco_all$snp_annot) / nrow(pco_all)
# table(pco_all$snp_annot, pco_all$is_nov)

## coding vs noncoding snps
snp_annot = rbind(as.data.frame(table(dgn_cis_annot$coding) / nrow(dgn_cis_annot)),
                  as.data.frame(table(dgn_trans_annot$coding) / nrow(dgn_trans_annot)),
                  as.data.frame(table(ukb_cis_annot$coding) / nrow(ukb_cis_annot)),
                  as.data.frame(table(ukb_trans_annot$coding) / nrow(ukb_trans_annot)))
colnames(snp_annot) = c('SNP', 'prop')
snp_annot$data = rep(c('cis-eQTL', 'trans-eQTL', 
                       'cis-pQTL', 'trans-pQTL'), each = 2)
ggplot(snp_annot, aes(x = data, y = prop, fill = SNP, label = round(prop, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  labs(title = 'Annotation of trans-SNPs', x = "", y = "Proportion", color = '') +
  geom_text(position=position_dodge(width = 0.8), vjust = -0.4) + 
  scale_fill_manual(values = brewer.pal(3, "Paired")[c(2,1)]) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

snp_annot_sub = snp_annot[snp_annot$data %in% c('trans-eQTL', 'trans-pQTL') &
                            snp_annot$SNP == 'coding', ]
chisq.test(rbind(c(15,558), c(3149, 19032)))
ggplot(snp_annot_sub, aes(x = data, y = prop)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5, fill = 'steelblue') +
  labs(title = 'Coding variants', x = "", y = "Proportion", color = '') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## detailed annotation
annot_group = function(x){
  annot2 = x$annot
  annot2[which(annot2 %in% c('intron_variant', 'intragenic_variant', 'splice_region_variant&intron_variant',
                             'splice_donor_variant&intron_variant', 'splice_region_variant',
                             'splice_donor_variant&splice_region_variant&intron_variant',
                             'splice_acceptor_variant&intron_variant',
                             'splice_region_variant&non_coding_transcript_exon_variant'))] = 'Intronic'
  annot2[which(annot2 == '3_prime_UTR_variant')] = "3' UTR"
  annot2[which(annot2 %in% c('5_prime_UTR_variant', 
                             '5_prime_UTR_premature_start_codon_gain_variant'))] = "5' UTR"
  annot2[which(annot2 %in% c('intergenic_region', 'sequence_feature'))] = 'Intergenic'
  annot2[which(annot2 %in% c('upstream_gene_variant', 'downstream_gene_variant',
                             'structural_interaction_variant', 'TF_binding_site_variant',
                             'non_coding_transcript_exon_variant'))] = 'Upstream & downstream gene'
  annot2[which(annot2 %in% c('synonymous_variant', 'splice_region_variant&synonymous_variant'))] = 'Synonymous'
  annot2[which(annot2 %in% c('missense_variant&splice_region_variant', 
                             'missense_variant', 'conservative_inframe_deletion',
                             'start_lost&splice_region_variant', 'stop_gained&splice_region_variant',
                             'disruptive_inframe_deletion', 'stop_gained', 'disruptive_inframe_insertion',
                             'frameshift_variant', 'stop_lost', 'start_lost'))] = 'Missense'
                   
  x = cbind(x, annot2 = annot2)
  return(x)
}
dgn_trans_annot = annot_group(dgn_trans_annot)
ukb_trans_annot = annot_group(ukb_trans_annot)
ukb_pco_annot = annot_group(ukb_pco_annot)

dgn_trans_annot$annot2 = as.factor(dgn_trans_annot$annot2)

boot_annot = function(data, i){
  df = data[i, ]
  table(df$annot2)/nrow(df)
}
dgn_trans_boot = boot(dgn_trans_annot, boot_annot, R = 1000)
ukb_trans_boot = boot(ukb_trans_annot, boot_annot, R = 1000)
annot_enrich = ukb_trans_boot$t / dgn_trans_boot$t
empirical_p = apply(annot_enrich, 2, function(x){mean(x > 1)})
test_p = 2 * ifelse(empirical_p > 0.5, 1 - empirical_p, empirical_p)
names(test_p) = c("3' UTR", "5' UTR", "Intergenic", "Intronic", 
                  "Missense", "Synonymous", "Upstream & downstream gene")
  
annot_enrich[annot_enrich > 100] = NA
annot_plot = data.frame(enrich = log2(apply(annot_enrich, 2, function(x){mean(x, na.rm = T)})),
                      quant_lower = log2(apply(annot_enrich, 2, function(x){quantile(x, 0.025, na.rm = T)})), 
                      quant_higher = log2(apply(annot_enrich, 2, function(x){quantile(x, 0.975, na.rm = T)})),
                      annot = c("3' UTR", "5' UTR", "Intergenic", "Intronic", 
                                "Missense", "Synonymous", "Upstream & downstream gene"))
ggplot(annot_plot, aes(x = reorder(annot, enrich), y = enrich)) + 
  geom_point(stat="identity", size = 2, color = 'black') +
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.2) +
  geom_hline(yintercept = 0, lty = 2) + 
  coord_flip() + 
  labs(title = '', x = "", y = "log2 fold change of trans-pQTL vs trans-eQTL", fill = '') +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## trans-pQTLs coloc with GWAS vs. trans-pQTL not coloc with GWAS
# UKB trans-pQTL
ukb_loci = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
ukb_trans = ukb_loci[ukb_loci$cis_trans == 'trans']
ukb_qtl_annot$loci = ukb_trans$loci
pp4_sig = fread('pqtl/08_gwas_coloc/ukb_coloc/ukb_pp4_gwas_5e8.txt')
pp4_trans_sig = pp4_sig[pp4_sig$cis_trans == 'trans', ]
ukb_qtl_annot$coloc = ukb_qtl_annot$loci %in% pp4_trans_sig$loci

# Trans-PCO trans-pQTL
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
ukb_pco_annot$loci = pco_all$loci
pp4_pco = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4_sig_gwas_5e8.txt')
ukb_pco_annot$coloc = ukb_pco_annot$loci %in% pp4_pco$loci

coloc_annot = rbind(as.data.frame(table(ukb_qtl_annot$coding, ukb_qtl_annot$coloc)/
                                    nrow(ukb_qtl_annot)),
                    as.data.frame(table(ukb_pco_annot$coding, ukb_pco_annot$coloc)/
                                    nrow(ukb_pco_annot)))
colnames(coloc_annot) = c('SNP', 'coloc', 'prop')
coloc_annot$data = rep(c('UKB', 'UKB PCO'), each = 4)
coloc_annot$prop_coloc = coloc_annot$prop
coloc_annot$prop_coloc[c(1,3)] = coloc_annot$prop_coloc[c(1,3)] / 
  sum(coloc_annot$prop_coloc[c(1,3)])
coloc_annot$prop_coloc[c(2,4)] = coloc_annot$prop_coloc[c(2,4)] / 
  sum(coloc_annot$prop_coloc[c(2,4)])
coloc_annot$prop_coloc[c(5,7)] = coloc_annot$prop_coloc[c(5,7)] / 
  sum(coloc_annot$prop_coloc[c(5,7)])
coloc_annot$prop_coloc[c(6,8)] = coloc_annot$prop_coloc[c(6,8)] / 
  sum(coloc_annot$prop_coloc[c(6,8)])

p1 = ggplot(coloc_annot[coloc_annot$data == 'UKB', ], 
       aes(x = SNP, y = prop_coloc, fill = coloc, label = round(prop_coloc, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  labs(title = 'Annotation of UKB trans-SNPs', x = "", y = "Proportion", color = '') +
  geom_text(position=position_dodge(width = 0.8), vjust = -0.4) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2 = ggplot(coloc_annot[coloc_annot$data == 'UKB PCO', ], 
       aes(x = SNP, y = prop_coloc, fill = coloc, label = round(prop_coloc, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  labs(title = 'Annotation of UKB PCO trans-SNPs', x = "", y = "Proportion", color = '') +
  geom_text(position=position_dodge(width = 0.8), vjust = -0.4) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)

