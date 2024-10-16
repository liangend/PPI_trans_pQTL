library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui/')
eqtlgen_trans = fread('pqtl/12_beta_across_two_data/ukb_beta_eqtlgen_trans_all.txt')
colnames(eqtlgen_trans) = c("eqtlgen_p", "eqtlgen_snp", "eqtlgen_chr", "eqtlgen_pos",
                            "eqtlgen_A1", "eqtlgen_A0", "eqtlgen_z", "eqtlgen_gene_id",
                            "eqtlgen_gene_name", "eqtlgen_gene_chr", "eqtlgen_gene_pos", 
                            "eqtlgen_NrCohorts", "eqtlgen_NrSamples", "eqtlgen_fdr",
                            "eqtlgen_bon_p", "ukb_A0", "ukb_A1", "ukb_A1_frq",        
                            "ukb_beta", "ukb_se", "ukb_logP")

# clump eQTLGen SNPs
uniq_gene = unique(eqtlgen_trans$eqtlgen_gene_id)
eqtlgen_trans_loci = c()
window_size = 250000
for (i in uniq_gene) {
  ukb_beta_sub = eqtlgen_trans[eqtlgen_trans$eqtlgen_gene_id == i, ]
  while(nrow(ukb_beta_sub > 0)) {
    lead_snp = which.min(ukb_beta_sub$eqtlgen_p)[1]
    lead_chr = ukb_beta_sub$eqtlgen_chr[lead_snp]
    lead_pos = ukb_beta_sub$eqtlgen_pos[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(ukb_beta_sub$eqtlgen_chr == lead_chr & 
                       ukb_beta_sub$eqtlgen_pos >= start_pos &
                       ukb_beta_sub$eqtlgen_pos <= end_pos)
    eqtlgen_trans_loci = rbind(eqtlgen_trans_loci, ukb_beta_sub[lead_snp, ])
    ukb_beta_sub = ukb_beta_sub[-rm_index, ]
  }
}
eqtlgen_trans_loci = na.omit(eqtlgen_trans_loci)


## tf enrichment
tf_gene_id = fread('pqtl/15_TF/TFs_Ensembl_v_1.01.txt', header = F)
hg37_gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
hg37_gene_meta$gene_id = sapply(strsplit(hg37_gene_meta$gene_id, '.', fixed = T), '[', 1)
hg37_gene_meta = hg37_gene_meta[hg37_gene_meta$gene_type == 'protein_coding', ]

eqtl_gene = c()
for (i in 1:nrow(eqtlgen_trans_loci)) {
  chr_i = eqtlgen_trans_loci$eqtlgen_chr[i]
  bp_i = eqtlgen_trans_loci$eqtlgen_pos[i]
  gene_sub = hg37_gene_meta[hg37_gene_meta$chr == paste0('chr', chr_i), ]
  eqtl_gene[i] = gene_sub$gene_id[which.min(abs(gene_sub$start - bp_i))]
}

eqtlgen_trans_loci$near_gene = eqtl_gene
eqtlgen_trans_loci$near_gene_name = hg37_gene_meta$gene_name[match(eqtlgen_trans_loci$near_gene, 
                                                                   hg37_gene_meta$gene_id)]
eqtlgen_trans_loci$is_tf = (eqtl_gene %in% tf_gene_id$V1)

boot_tf = function(trans_all, ukb_p_thre1, ukb_p_thre2, 
                   eqtl_p_thre1, eqtl_p_thre2, n_sample){
  trans_sub = trans_all[trans_all$eqtlgen_p >= eqtl_p_thre1 & 
                          trans_all$eqtlgen_p < eqtl_p_thre2 & 
                          trans_all$ukb_logP >= ukb_p_thre1 &
                          trans_all$ukb_logP < ukb_p_thre2, ]
  print(paste0('Number of cases: ', nrow(trans_sub)))
  print(paste0('Number of TF: ', sum(trans_sub$is_tf)))
  #baseline = trans_all[trans_all$ukb_logP > ukb_p_thre, ]
  n_qtl = nrow(trans_sub)
  enrich = c()
  for (i in 1:n_sample) {
    sample_i = trans_all[sample(1:nrow(trans_all), n_qtl), ]
    enrich[i] = sum(trans_sub$is_tf) / sum(sample_i$is_tf)
  }
  return(enrich)
}

se = function(x) {
  x = na.omit(x)
  x = x[x < 100]
  return(sd(x) / sqrt(length(x)))
}

## TF enrichment
tf_enrich1 = boot_tf(eqtlgen_trans_loci, 4, 500, 0, 1, 100)
tf_enrich2 = boot_tf(eqtlgen_trans_loci, 2, 4, 0, 1, 100)
tf_enrich3 = boot_tf(eqtlgen_trans_loci, -log10(0.05), 2, 0, 1, 100)
tf_enrich4 = boot_tf(eqtlgen_trans_loci, 1, -log10(0.05), 0, 1, 100)

tf_table = data.frame(enrich = c(mean(tf_enrich1), mean(tf_enrich2),
                                 mean(tf_enrich3), mean(tf_enrich4)),
                      se = c(se(tf_enrich1), se(tf_enrich2), 
                             se(tf_enrich3), se(tf_enrich4)),
                      ukb_p = c('q1', 'q2', 'q3', 'q4'))
tf_table$type = 'TF'
ggplot(tf_table, aes(x = factor(ukb_p, levels = c('q1', 'q2', 'q3', 'q4')),
                     y = enrich)) + 
  geom_point(size = 2) +
  # geom_hline(yintercept = enrich_table$enrich[4], lty = 2, color = "#00BFC4") + 
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=.05) +
  labs(x = "UKB p value", y = "enrich", color = '',
       title = 'Enrichment for nearby genes of trans-eQTLs in TF') + 
  scale_x_discrete(labels=c("< 0.001", "0.001 ~ 0.01", "0.01 ~ 0.05", "0.05 ~ 0.1")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### protein-protein interaction and complex
# BioPlex database
bioplex_mod = fread('pqtl/11_prot_complex/bioplex_mod.txt')
bioplex_mod$mod = paste0('bioplex_', bioplex_mod$mod)

# CORUM database
corum = read.xlsx('pqtl/11_prot_complex/CORUM_2022_09_12.xlsx')
corum = corum[corum$Organism == 'Human', ]
prot_list = strsplit(corum$`subunits(Gene.name)`, ';', fixed = T)
corum_size = sapply(prot_list, length)
corum_mod = data.frame(prot = unlist(prot_list), 
                       mod = rep(1:length(prot_list), corum_size))
corum_mod$mod = paste0('corum_', corum_mod$mod)

# STRING database
string_mod = fread('pqtl/11_prot_complex/string_mod.txt')
string_mod$mod = paste0('string_', string_mod$mod)
string_mod = string_mod[, c(3,2)]
colnames(string_mod)[1] = 'prot'
# HIPPIE database
hippie_mod = fread('pqtl/11_prot_complex/hippie_mod.txt')
hippie_mod$mod = paste0('hippie_', hippie_mod$mod)

# all ppi
ppi_all = rbind(bioplex_mod, corum_mod, string_mod, hippie_mod)
ppi_all = ppi_all[ppi_all$prot %in% c(eqtlgen_trans_loci$eqtlgen_gene_name, 
                                      eqtlgen_trans_loci$near_gene_name), ]
## PPI enrichment
ppi_overlap = apply(eqtlgen_trans_loci, 1, function(x){
  cis_mod = ppi_all$mod[which(ppi_all$prot == x[23])]
  trans_mod = ppi_all$mod[which(ppi_all$prot == x[9])]
  return(length(intersect(cis_mod, trans_mod)))
})
eqtlgen_trans_loci$ppi = ppi_overlap
sum(ppi_overlap > 0) / nrow(eqtlgen_trans_loci)

## background PPI
n_ite = 100
n_ppi = c()
trans_sample = eqtlgen_trans_loci[, c('eqtlgen_gene_name', 'near_gene_name')]
for (i in 1:n_ite) {
  trans_sample$eqtlgen_gene_name = sample(trans_sample$eqtlgen_gene_name)
  trans_sample$near_gene_name = sample(trans_sample$near_gene_name)
  ppi_sample = apply(trans_sample, 1, function(x){
    cis_mod = ppi_all$mod[which(ppi_all$prot == x[1])]
    trans_mod = ppi_all$mod[which(ppi_all$prot == x[2])]
    return(length(intersect(cis_mod, trans_mod)))
  })
  n_ppi[i] = sum(ppi_sample > 0)
  print(i)
}

null_ppi = n_ppi / nrow(trans_sample)

ppi_enrich1 = sum(eqtlgen_trans_loci[eqtlgen_trans_loci$ukb_logP > 4, ]$ppi > 0) / 
  sum(eqtlgen_trans_loci$ukb_logP > 4)
ppi_enrich2 = sum(eqtlgen_trans_loci[eqtlgen_trans_loci$ukb_logP > 2 & eqtlgen_trans_loci$ukb_logP < 4, ]$ppi > 0) / 
  sum(eqtlgen_trans_loci$ukb_logP > 2 & eqtlgen_trans_loci$ukb_logP < 4)
ppi_enrich3 = sum(eqtlgen_trans_loci[eqtlgen_trans_loci$ukb_logP > -log10(0.05) & eqtlgen_trans_loci$ukb_logP < 2, ]$ppi > 0) / 
  sum(eqtlgen_trans_loci$ukb_logP > -log10(0.05) & eqtlgen_trans_loci$ukb_logP < 2)
ppi_enrich4 = sum(eqtlgen_trans_loci[eqtlgen_trans_loci$ukb_logP > 1 & eqtlgen_trans_loci$ukb_logP < -log10(0.05), ]$ppi > 0) / 
  sum(eqtlgen_trans_loci$ukb_logP > 1 & eqtlgen_trans_loci$ukb_logP < -log10(0.05))

ppi_table = data.frame(enrich = c(mean(ppi_enrich1/null_ppi), mean(ppi_enrich2/null_ppi),
                                  mean(ppi_enrich3/null_ppi), mean(ppi_enrich4/null_ppi)),
                       se = c(se(ppi_enrich1/null_ppi), se(ppi_enrich2/null_ppi),
                              se(ppi_enrich3/null_ppi), se(ppi_enrich4/null_ppi)),
                       ukb_p = c('q1', 'q2', 'q3', 'q4'))
ppi_table$type = 'PPI'
ggplot(ppi_table, aes(x = factor(ukb_p, levels = c('q1', 'q2', 'q3', 'q4')),
                      y = enrich)) + 
  geom_point(size = 2) +
  # geom_hline(yintercept = enrich_table$enrich[4], lty = 2, color = "#00BFC4") + 
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=.05) +
  labs(x = "DGN p value", y = "enrich", color = '',
       title = 'Enrichment for nearby genes of trans-eQTLs in PPI') + 
  scale_x_discrete(labels=c("< 0.001", "0.001 ~ 0.01", "0.01 ~ 0.05", "0.05 ~ 0.1")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

enrich_table = rbind(tf_table, ppi_table)
ggplot(enrich_table, aes(x = factor(ukb_p, levels = c('q1', 'q2', 'q3', 'q4')),
                         y = enrich, color = type)) + 
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=.05) +
  labs(x = "UKB p value", y = "enrich", color = '',
       title = 'Enrichment for nearby genes of trans-eQTLs') + 
  scale_x_discrete(labels=c("< 0.0001", "0.0001 ~ 0.01", "0.01 ~ 0.05", "0.05 ~ 0.1")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
# library(scales)
# hue_pal()(2)




