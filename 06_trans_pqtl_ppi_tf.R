library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(boot)
setwd('/project/xuanyao/jinghui/')
ukb_sig = fread('pqtl/12_beta_across_two_data/dgn_beta_ukb_sig_all.txt')
colnames(ukb_sig) = c('ukb_chr', 'ukb_bp_hg38', 'ukb_bp_hg19', 'ukb_A0', 'ukb_A1',
                      'ukb_A1_frq', 'ukb_N', 'ukb_beta', 'ukb_se', 'ukb_logP',
                      'ukb_file', 'ukb_gene_chr', 'ukb_tss', 'cis_trans', 'ukb_gene_id',
                      'ukb_gene_name', 'ukb_loci', 'ukb_loci_start', 'ukb_loci_end', 
                      'dgn_beta', 'dgn_se', 'dgn_p')
ukb_sig = ukb_sig[ukb_sig$ukb_logP > -log10(5e-8/2923), ]

#### tf enrichment
tf_gene_id = fread('pqtl/15_TF/humantf_ccbr/TFs_Ensembl_v_1.01.txt', header = F)
hg38_gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.unadj.gtf')
hg38_gene_meta$gene_id = sapply(strsplit(hg38_gene_meta$gene_id, '.', fixed = T), '[', 1)
hg38_gene_meta = hg38_gene_meta[hg38_gene_meta$gene_type == 'protein_coding']
tf_gene_id$gene_name = hg38_gene_meta$gene_name[match(tf_gene_id$V1, hg38_gene_meta$gene_id)]

gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta_sub = gene_meta_hg37[gene_meta_hg37$gene_type == 'protein_coding', ]
gene_meta_sub = gene_meta_sub[!grepl('gene_status', gene_meta_sub$gene_name), ]
gene_meta_sub = gene_meta_sub[!grepl('ENSG', gene_meta_sub$gene_name), ]
gene_meta_sub$gene_length = gene_meta_sub$end - gene_meta_sub$start

## tf enrich of trans-pQTL
pqtl_gene = c()
for (i in 1:nrow(ukb_sig)) {
  chr_i = ukb_sig$ukb_chr[i]
  bp_i = ukb_sig$ukb_bp_hg19[i]
  gene_sub = gene_meta_sub[gene_meta_sub$chr == paste0('chr', chr_i), ]
  pqtl_gene[i] = gene_sub$gene_id[which.min(abs(gene_sub$start - bp_i))]
}

ukb_sig$near_gene_id = pqtl_gene
ukb_sig$near_gene_name = gene_meta_hg37$gene_name[match(ukb_sig$near_gene_id, 
                                                        gene_meta_hg37$gene_id)]
ukb_sig$is_tf = (ukb_sig$near_gene_name %in% tf_gene_id$gene_name)

ukb_trans_pqtl = ukb_sig[ukb_sig$cis_trans == 'trans', ]

background_tf_pqtl = matrix(rep(F, nrow(ukb_trans_pqtl) * 1000), ncol = 1000)
for (i in 1:nrow(ukb_trans_pqtl)) {
  gene_i = ukb_trans_pqtl$near_gene_name[i]
  gene_length_i = gene_meta_sub$gene_length[match(gene_i, gene_meta_sub$gene_name)]
  gene_meta_i = gene_meta_sub[gene_meta_sub$gene_length >= 0.9*gene_length_i & 
                                gene_meta_sub$gene_length <= 1.1*gene_length_i, ]
  background_gene_i = sample(gene_meta_i$gene_name, 1000, replace = T)
  background_tf_pqtl[i, ] = background_gene_i %in% tf_gene_id$gene_name
  if (i %% 1000 == 0) {
    print(i)
  }
}
background_tf_pqtl = apply(background_tf_pqtl, 2, sum)
obs_tf_pqtl = sum(ukb_trans_pqtl$is_tf)
tf_enrich_pqtl = obs_tf_pqtl / background_tf_pqtl

## tf enrich of trans-eQTL (DGN)
dgn_trans = fread('pqtl/00_ref/dgn_trans_eqtl.txt')
eqtl_gene = c()
for (i in 1:nrow(dgn_trans)) {
  chr_i = dgn_trans$var_chr[i]
  bp_i = dgn_trans$var_bp[i]
  gene_sub = gene_meta_sub[gene_meta_sub$chr == chr_i, ]
  eqtl_gene[i] = gene_sub$gene_id[which.min(abs(gene_sub$start - bp_i))]
}
dgn_trans$near_gene_id = eqtl_gene
dgn_trans$near_gene_name = gene_meta_hg37$gene_name[match(dgn_trans$near_gene_id, 
                                                          gene_meta_hg37$gene_id)]
dgn_trans$is_tf = (dgn_trans$near_gene_name %in% tf_gene_id$gene_name)

background_tf_eqtl = matrix(rep(F, nrow(dgn_trans) * 1000), ncol = 1000)
for (i in 1:nrow(dgn_trans)) {
  gene_i = dgn_trans$near_gene_name[i]
  gene_length_i = gene_meta_sub$gene_length[match(gene_i, gene_meta_sub$gene_name)]
  gene_meta_i = gene_meta_sub[gene_meta_sub$gene_length >= 0.9*gene_length_i & 
                                gene_meta_sub$gene_length <= 1.1*gene_length_i, ]
  background_gene_i = sample(gene_meta_i$gene_name, 1000, replace = T)
  background_tf_eqtl[i, ] = background_gene_i %in% tf_gene_id$gene_name
}
background_tf_eqtl = apply(background_tf_eqtl, 2, sum)
obs_tf_eqtl = sum(dgn_trans$is_tf)
tf_enrich_eqtl = obs_tf_eqtl / background_tf_eqtl

enrich_tf = function(trans_all, dgn_p_thre1, dgn_p_thre2){
  ukb_trans_more_eqtl = trans_all[trans_all$dgn_p >= dgn_p_thre1 & 
                                    trans_all$dgn_p < dgn_p_thre2, ]
  print(paste0('Number of cases: ', nrow(ukb_trans_more_eqtl)))
  print(paste0('Number of TF: ', sum(ukb_trans_more_eqtl$is_tf)))
  #baseline = trans_all[trans_all$ukb_logP > ukb_p_thre, ]
  n_qtl = nrow(ukb_trans_more_eqtl)
  enrich = sum(ukb_trans_more_eqtl$is_tf) / n_qtl / 
    (background_tf / nrow(ukb_trans_pqtl))
  return(enrich)
}

ukb_trans_sig = ukb_trans_pqtl[ukb_trans_pqtl$ukb_logP > -log10(5e-8/2922), ]
ukb_trans_sig = ukb_trans_sig[ukb_trans_sig$dgn_se != 0, ]
## TF enrichment
tf_enrich_dgn1 = enrich_tf(ukb_trans_sig, 0, 0.01)
tf_enrich_dgn2 = enrich_tf(ukb_trans_sig, 0.01, 0.05)
tf_enrich_dgn3 = enrich_tf(ukb_trans_sig, 0.05, 0.1)
tf_enrich_dgn4 = enrich_tf(ukb_trans_sig, 0.1, 1)

tf_table = data.frame(enrich = c(mean(tf_enrich_dgn1), mean(tf_enrich_dgn2),
                                 mean(tf_enrich_dgn3), mean(tf_enrich_dgn4)),
                      ci_low = c(quantile(tf_enrich_dgn1, 0.025), quantile(tf_enrich_dgn2, 0.025), 
                                 quantile(tf_enrich_dgn3, 0.025), quantile(tf_enrich_dgn4, 0.025)),
                      ci_high = c(quantile(tf_enrich_dgn1, 0.975), quantile(tf_enrich_dgn2, 0.975), 
                                 quantile(tf_enrich_dgn3, 0.975), quantile(tf_enrich_dgn4, 0.975)),
                      dgn_p = c('q1', 'q2', 'q3', 'q4'),
                      type = 'TF')

p1 = ggplot(tf_table, aes(x = factor(dgn_p, levels = c('q1', 'q2', 'q3', 'q4')),
                          y = enrich, group = type)) + 
  geom_point(size = 2, color = brewer.pal(3, "Set1")[2]) +
 # geom_hline(yintercept = enrich_table$enrich[4], lty = 2, color = "#00BFC4") + 
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.05, color = brewer.pal(3, "Set1")[2]) +
  geom_line(linetype = "dashed", color = brewer.pal(3, "Set1")[2]) +
  labs(x = "", y = "Enrichment", color = '',
       title = '') + 
  scale_x_discrete(labels=c('', '', '', '')) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### protein-protein interaction and complex
ppi_list = readRDS('pqtl/11_prot_complex/ppi_list.rds')
background_ppi = c()
for (i in 1:1000) {
  gene_pair_i = paste0(sample(c(ukb_trans_pqtl$near_gene_name, 
                                dgn_trans$near_gene_name)), ',', 
                       c(ukb_trans_pqtl$ukb_gene_name, dgn_trans$phenotype_id))
  background_ppi[i] = sum(gene_pair_i %in% ppi_list)
  if (i %% 100 == 0) {
    print(i)
  }
}

background_ppi_prop = background_ppi/nrow(ukb_trans_pqtl)

obs_ppi_pqtl = sum(paste0(ukb_trans_pqtl$near_gene_name, ',', ukb_trans_pqtl$ukb_gene_name) %in% 
                     ppi_list)
obs_ppi_eqtl = sum(paste0(dgn_trans$near_gene_name, ',', dgn_trans$phenotype_id) %in% 
                     ppi_list)

enrich_ppi = function(df, dgn_p_low, dgn_p_high) {
  df_sub = df[df$dgn_p > dgn_p_low & df$dgn_p < dgn_p_high, ]
  print(paste0('# pairs: ', nrow(df_sub)))
  gene_pair_sub = paste0(df_sub$near_gene_name, ',', df_sub$ukb_gene_name)
  ppi_prop = sum(gene_pair_sub %in% ppi_list) / nrow(df_sub)
  return(ppi_prop / background_ppi_prop)
}
ppi_enrich_dgn1 = enrich_ppi(ukb_trans_pqtl, 0, 0.01)
ppi_enrich_dgn2 = enrich_ppi(ukb_trans_pqtl, 0.01, 0.05)
ppi_enrich_dgn3 = enrich_ppi(ukb_trans_pqtl, 0.05, 0.1)
ppi_enrich_dgn4 = enrich_ppi(ukb_trans_pqtl, 0.1, 1)

ppi_table = data.frame(enrich = c(mean(ppi_enrich_dgn1), mean(ppi_enrich_dgn2),
                                  mean(ppi_enrich_dgn3), mean(ppi_enrich_dgn4)),
                       ci_low = c(quantile(ppi_enrich_dgn1, 0.025), quantile(ppi_enrich_dgn2, 0.025), 
                                  quantile(ppi_enrich_dgn3, 0.025), quantile(ppi_enrich_dgn4, 0.025)),
                       ci_high = c(quantile(ppi_enrich_dgn1, 0.975), quantile(ppi_enrich_dgn2, 0.975), 
                                   quantile(ppi_enrich_dgn3, 0.975), quantile(ppi_enrich_dgn4, 0.975)),
                       dgn_p = c('q1', 'q2', 'q3', 'q4'),
                       type = 'PPI')

p2 = ggplot(ppi_table, aes(x = factor(dgn_p, levels = c('q1', 'q2', 'q3', 'q4')),
                           y = enrich, group = type)) + 
  geom_point(size = 2, color = brewer.pal(3, "Set1")[1]) +
  # geom_hline(yintercept = enrich_table$enrich[4], lty = 2, color = "#00BFC4") + 
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.05, color = brewer.pal(3, "Set1")[1]) +
  geom_line(linetype = "dashed", color = brewer.pal(3, "Set1")[1]) +
  labs(x = "eQTL p value", y = "Enrichment",
       title = '') + 
  scale_x_discrete(labels=c("< 0.01", "0.01~0.05", "0.05~0.1", "> 0.1")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 2)

### plot tf and ppi together
enrich_all = data.frame(enrich = c(mean(obs_tf_pqtl/background_tf_pqtl), mean(obs_ppi_pqtl/background_ppi), 
                                   mean(obs_tf_eqtl/background_tf_eqtl), mean(obs_ppi_eqtl/background_ppi)),
                        ci_low = c(quantile(obs_tf_pqtl/background_tf_pqtl, 0.025), 
                                   quantile(obs_ppi_pqtl/background_ppi, 0.025),
                                   quantile(obs_tf_eqtl/background_tf_eqtl, 0.025), 
                                   quantile(obs_ppi_eqtl/background_ppi, 0.025)),
                        ci_high = c(quantile(obs_tf_pqtl/background_tf_pqtl, 0.975), 
                                    quantile(obs_ppi_pqtl/background_ppi, 0.975),
                                    quantile(obs_tf_eqtl/background_tf_eqtl, 0.975), 
                                    quantile(obs_ppi_eqtl/background_ppi, 0.975)),
                        type = c('TF', 'PPI', 'TF', 'PPI'),
                        trait = c('trans-pQTL', 'trans-pQTL', 'trans-eQTL', 'trans-eQTL'))
# p value of tf enrichment
tf_enrich_p = obs_tf/background_tf - 1
empirical_p_tf = mean(tf_enrich_p > 0) 
2*min(c(empirical_p_tf, 1-empirical_p_tf))

# p value of ppi enrichment
ppi_enrich_p = obs_ppi/background_ppi - 1
empirical_p_ppi = mean(ppi_enrich_p > 0) 
2*min(c(empirical_p_ppi, 1-empirical_p_ppi))

ggplot(enrich_all, aes(x = type, y = enrich, fill = trait)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  geom_errorbar(position=position_dodge(0.5), aes(ymin=ci_low, ymax=ci_high), width=.2) + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  labs(x = "", y = "Enrichment", fill = '',  title = '') + 
  scale_fill_manual(values = brewer.pal(3, "Dark2")[1:2]) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## TF and PPI enrichment after removing trans-pQTLs coloc with cell-composition GWAS
pp4_sig = fread('pqtl/08_gwas_coloc/ukb_coloc/ukb_pp4_gwas_5e8.txt')
cell_prop_trait = c('RBC', 'WBC', 'PLT', 'LYMPH', 'NEUT')
cell_comp_loci = unique(pp4_sig$loci[which(pp4_sig$trait %in% cell_prop_trait)])

ukb_trans_sub = ukb_trans_pqtl[!ukb_trans_pqtl$ukb_loci %in% cell_comp_loci, ]
obs_tf_sub = sum(ukb_trans_sub$is_tf)
obs_ppi_sub = sum(paste0(ukb_trans_sub$near_gene_name, ',', ukb_trans_sub$ukb_gene_name) %in%
                    ppi_list)

enrich_sub = data.frame(enrich = c(mean(obs_tf_sub/background_tf), mean(obs_ppi_sub/background_ppi)),
                        ci_low = c(quantile(obs_tf_sub/background_tf, 0.025), 
                                   quantile(obs_ppi_sub/background_ppi, 0.025)),
                        ci_high = c(quantile(obs_tf_sub/background_tf, 0.975), 
                                    quantile(obs_ppi_sub/background_ppi, 0.975)),
                        type = c('TF', 'PPI'))
ggplot(enrich_sub, aes(x = type, y = enrich, fill = type)) + 
  geom_bar(stat="identity", position=position_dodge(0.1), width = 0.5) +
  geom_errorbar(position=position_dodge(0.1), aes(ymin=ci_low, ymax=ci_high), width=.2) + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  labs(x = "", y = "Enrichment", fill = '',  title = '') + 
  scale_fill_manual(values = brewer.pal(3, "Set1")[1:2]) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')




