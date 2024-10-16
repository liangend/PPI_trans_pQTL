library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui/')
hg38_gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.unadj.gtf')
hg38_gene_meta$gene_id = sapply(strsplit(hg38_gene_meta$gene_id, '.', fixed = T), '[', 1)
hg38_gene_meta = hg38_gene_meta[hg38_gene_meta$gene_type == 'protein_coding', ]
ukb_pqtl = fread('pqtl/12_beta_across_two_data/dgn_beta_ukb_pqtl_all.txt')
near_gene = c()
for (i in 1:nrow(ukb_pqtl)) {
  chr_i = ukb_pqtl$ukb_chr[i]
  bp_i = ukb_pqtl$ukb_bp_hg38[i]
  gene_sub = hg38_gene_meta[hg38_gene_meta$chr == paste0('chr',chr_i), ]
  near_gene[i] = gene_sub$gene_name[which.min(abs(gene_sub$start - bp_i))]
}
ukb_pqtl$near_gene = near_gene
ukb_pqtl = ukb_pqtl[ukb_pqtl$ukb_cis_trans == 'trans', ]

## protein-protein interaction and complex
hippie = fread('pqtl/11_prot_complex/hippie_current.txt')
hip_prot1 = sapply(strsplit(hippie$V1, '_', fixed = T), '[', 1)
hip_prot2 = sapply(strsplit(hippie$V3, '_', fixed = T), '[', 1)
cis_trans_inter = c()
for (i in 1:nrow(ukb_pqtl)) {
  prot1 = ukb_pqtl$ukb_gene[i]
  prot2 = ukb_pqtl$near_gene[i]
  row_index = which((hip_prot1 == prot1 & hip_prot2 == prot2) |
                      (hip_prot1 == prot2 & hip_prot2 == prot1))
  if (length(row_index) > 0) {
    cis_trans_inter[i] = hippie$V5[row_index]
  } else {
    cis_trans_inter[i] = NA
  }
  if (i %% 100 == 0) {
    print(i)
  }
}
mean(cis_trans_inter[which(!is.na(cis_trans_inter))])
ukb_pqtl$hippie_inter = cis_trans_inter

string_clust_prot = fread('pqtl/11_prot_complex/string_cluster/9606.clusters.proteins.v12.0.txt.gz')
string_clust = fread('pqtl/11_prot_complex/string_cluster/9606.clusters.info.v12.0.txt.gz')
string_prot = fread('pqtl/11_prot_complex/string_cluster/9606.protein.info.v12.0.txt.gz')
big_clust = string_clust$cluster_id[which(string_clust$cluster_size > 500)]
string_clust_prot = string_clust_prot[!(string_clust_prot$cluster_id %in% big_clust), ]
string_clust_prot$prot_name = string_prot$preferred_name[match(string_clust_prot$protein_id, 
                                                               string_prot$`#string_protein_id`)]
string_overlap = apply(ukb_pqtl, 1, function(x){
  cis_clust = string_clust_prot$cluster_id[which(string_clust_prot$prot_name == x[27])]
  trans_clust = string_clust_prot$cluster_id[which(string_clust_prot$prot_name == x[9])]
  return(length(intersect(cis_clust, trans_clust)))
})
ukb_pqtl$string_overlap = string_overlap

corum = read.xlsx('pqtl/11_prot_complex/CORUM_2022_09_12.xlsx')
corum = corum[corum$Organism == 'Human', ]
prot_list = strsplit(corum$`subunits(Gene.name)`, ';', fixed = T)
corum_size = sapply(prot_list, length)
corum_table = data.frame(prot = unlist(prot_list), 
                         complex = rep(1:length(prot_list), corum_size))

corum_overlap = apply(ukb_pqtl, 1, function(x){
  cis_mod = corum_table$complex[which(corum_table$prot == x[27])]
  trans_mod = corum_table$complex[which(corum_table$prot == x[9])]
  return(length(intersect(cis_mod, trans_mod)))
})
ukb_pqtl$corum_overlap = corum_overlap


boot_tf = function(trans_all, p_thre, n_sample){
  ukb_trans_more_eqtl = trans_all[trans_all$dgn_p < p_thre, ]
  n_qtl = nrow(ukb_trans_more_eqtl)
  hippie_enrich = c()
  string_enrich = c()
  corum_enrich = c()
  for (i in 1:n_sample) {
    sample_i = trans_all[sample(1:nrow(trans_all), n_qtl), ]
    hippie_enrich[i] = sum(!is.na(ukb_trans_more_eqtl$hippie_inter)) / 
      sum(!is.na(sample_i$hippie_inter))
    string_enrich[i] = sum(ukb_trans_more_eqtl$string_overlap > 0) / 
      sum(sample_i$string_overlap > 0)
    corum_enrich[i] = sum(ukb_trans_more_eqtl$corum_overlap > 0) / 
      sum(sample_i$corum_overlap > 0)
  }
  enrich = cbind(hippie_enrich, string_enrich, corum_enrich)
  return(enrich)
}

enrich_02 = boot_tf(ukb_pqtl, 0.2, 1000)
enrich_01 = boot_tf(ukb_pqtl, 0.1, 1000)
enrich_005 = boot_tf(ukb_pqtl, 0.05, 1000)
enrich_001 = boot_tf(ukb_pqtl, 0.01, 1000)

string_table = data.frame(enrich = c(mean(enrich_02[,2]), mean(enrich_01[,2]), 
                                  mean(enrich_005[,2]), mean(enrich_001[which(enrich_001[,2] < 100),2])),
                       se = c(sd(enrich_02[,2])/sqrt(1000), sd(enrich_01[,2])/sqrt(1000),
                              sd(enrich_005[,2])/sqrt(1000), 
                              sd(enrich_001[which(enrich_001[,2] < 100),2])/sqrt(sum(enrich_001[,2] < 100))),
                       dgn_p = c('< 0.2', '< 0.1', '< 0.05', '< 0.01'))

p1 = ggplot(string_table, aes(x = factor(dgn_p, levels = c('< 0.2', '< 0.1', '< 0.05', '< 0.01')), y = enrich)) + 
  geom_point(position=position_dodge(0.3), size = 3, color = 'steelblue') +
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=.07, 
                position=position_dodge(0.3), color = 'steelblue') +
  labs(x = "DGN p value", y = "enrichment", title = 'Nearby genes and trans genes in STRING') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


hippie_table = data.frame(enrich = c(mean(enrich_02[enrich_02[,1] < 100,1]), 
                                     mean(enrich_01[enrich_01[,1] < 100,1]), 
                                     mean(enrich_005[enrich_005[,1] < 100,1])),
                          se = c(sd(enrich_02[enrich_02[,1] < 100,1])/sqrt(sum(enrich_02[,1] < 100)), 
                                 sd(enrich_01[enrich_01[,1] < 100,1])/sqrt(sum(enrich_01[,1] < 100)),
                                 sd(enrich_005[enrich_005[,1] < 100,1])/sqrt(sum(enrich_005[,1] < 100))),
                          dgn_p = c('< 0.2', '< 0.1', '< 0.05'))

p2 = ggplot(hippie_table, aes(x = factor(dgn_p, levels = c('< 0.2', '< 0.1', '< 0.05')), y = enrich)) + 
  geom_point(position=position_dodge(0.3), size = 3, color = 'steelblue') +
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=.07, 
                position=position_dodge(0.3), color = 'steelblue') +
  labs(x = "DGN p value", y = "enrichment", title = 'Nearby genes and trans genes in HIPPIE') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


corum_table = data.frame(enrich = c(mean(enrich_02[enrich_02[,3] < 100, 3]), 
                                     mean(enrich_01[enrich_01[,3] < 100, 3])),
                          se = c(sd(enrich_02[enrich_02[,3] < 100, 3])/sqrt(sum(enrich_02[,3] < 100)), 
                                 sd(enrich_01[enrich_01[,3] < 100, 3])/sqrt(sum(enrich_01[,3] < 100))),
                          dgn_p = c('< 0.2', '< 0.1'))

p3 = ggplot(corum_table, aes(x = factor(dgn_p, levels = c('< 0.2', '< 0.1')), y = enrich)) + 
  geom_point(position=position_dodge(0.3), size = 3, color = 'steelblue') +
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=.07, 
                position=position_dodge(0.3), color = 'steelblue') +
  labs(x = "DGN p value", y = "enrichment", title = 'Nearby genes and trans genes in CORUM') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1,p2,p3, nrow = 1)

