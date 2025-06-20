setwd('/project/xuanyao/jinghui')
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(boot)
library(openxlsx)
dgn_trans_annot = fread('pqtl/17_snpEff/dgn_trans_snp_annot.vcf')
dgn_trans_eqtl = fread('pqtl/00_ref/dgn_trans_eqtl.txt')
dgn_trans_annot$maf = dgn_trans_eqtl$af


ukb_all_annot = fread('pqtl/17_snpEff/ukb_all_annot.vcf')
ukb_pqtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
ukb_all_annot$maf = pmin(ukb_pqtl$A1FREQ, 1-ukb_pqtl$A1FREQ)

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

dgn_trans_annot = annot_coding(dgn_trans_annot)

ukb_all_annot = ukb_all_annot[ukb_all_annot$maf > 0.05, ]
ukb_all_annot = annot_coding(ukb_all_annot)
ukb_trans_annot = ukb_all_annot[ukb_all_annot$FILTER == 'trans', ]

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
  x = x[x$annot2 %in% c('Intronic', "3' UTR", "5' UTR", 'Intergenic', 
                        'Upstream & downstream gene', 'Synonymous', 'Missense')]
  return(x)
}
dgn_trans_annot = annot_group(dgn_trans_annot)
ukb_trans_annot = annot_group(ukb_trans_annot)

dgn_trans_annot$annot2 = as.factor(dgn_trans_annot$annot2)
ukb_trans_annot$annot2 = as.factor(ukb_trans_annot$annot2)

boot_annot = function(data, i){
  df = data[i, ]
  table(df$annot2)/nrow(df)
}
dgn_trans_boot = boot(dgn_trans_annot, boot_annot, R = 1000)
ukb_trans_boot = boot(ukb_trans_annot, boot_annot, R = 1000)

trans_annot_enrich = ukb_trans_boot$t / dgn_trans_boot$t
trans_empirical_p = apply(trans_annot_enrich, 2, function(x){mean(x > 1)})
trans_test_p = 2 * ifelse(trans_empirical_p > 0.5, 1 - trans_empirical_p, trans_empirical_p)
names(trans_test_p) = c("3' UTR", "5' UTR", "Intergenic", "Intronic", 
                  "Missense", "Synonymous", "Upstream & downstream gene")
  
trans_annot_enrich[trans_annot_enrich > 100] = NA
trans_annot_plot = data.frame(enrich = log2(apply(trans_annot_enrich, 2, function(x){mean(x, na.rm = T)})),
                      quant_lower = log2(apply(trans_annot_enrich, 2, function(x){quantile(x, 0.025, na.rm = T)})), 
                      quant_higher = log2(apply(trans_annot_enrich, 2, function(x){quantile(x, 0.975, na.rm = T)})),
                      annot = c("3' UTR", "5' UTR", "Intergenic", "Intronic", 
                                "Missense", "Synonymous", "Upstream & \ndownstream genes"),
                      group = 'trans')
ggplot(trans_annot_plot, aes(x = reorder(annot, enrich), y = enrich)) + 
  geom_point(stat="identity", size = 2, color = 'black') +
  geom_errorbar(aes(ymin=quant_lower, ymax=quant_higher), width=.2) +
  geom_hline(yintercept = 0, lty = 2) + 
  coord_flip() + 
  labs(title = '', x = "", y = "log2 fold enrichment of trans-pQTL vs trans-eQTL", fill = '') +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


