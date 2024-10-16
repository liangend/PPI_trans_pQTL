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
ukb_sig = ukb_sig[ukb_sig$ukb_logP > -log10(5e-8/2922), ]

##### TF enrichment
## TF based on TFLink
tf1 = fread('pqtl/15_TF/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz')
ukb_prot_list = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
gene_meta = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
ukb_prot_list$gene_name = gene_meta$gene_name[match(ukb_prot_list$gene_id, gene_meta$gene_id)]

tf1 = tf1[tf1$Name.Target %in% ukb_prot_list$gene_name, ]
tf1 = unique(tf1[, c('Name.TF', 'Name.Target')])
colnames(tf1) = c('TF', 'target')

## TF based on hTFtarget
tf2 = fread('pqtl/15_TF/tf-target-infomation.txt')
tf2 = tf2[tf2$tissue == 'blood', ]
tf2 = tf2[tf2$target %in% ukb_prot_list$gene_name, ]
tf2 = unique(tf2)
tf2 = tf2[, 1:2]

tf_all = rbind(tf1, tf2)
tf_all = unique(tf_all)
tf_list = paste0(tf_all$TF, ',', tf_all$target)


## nearby protein-coding genes of ukb pQTL
gene_meta_prot = gene_meta[gene_meta$gene_type == 'protein_coding', ]
gene_meta_prot = gene_meta_prot[!grepl('gene_status', gene_meta_prot$gene_name), ]
gene_meta_prot = gene_meta_prot[!grepl('ENSG', gene_meta_prot$gene_name), ]

pqtl_gene = c()
for (i in 1:nrow(ukb_sig)) {
  chr_i = ukb_sig$ukb_chr[i]
  bp_i = ukb_sig$ukb_bp_hg19[i]
  gene_sub = gene_meta_prot[gene_meta_prot$chr == paste0('chr', chr_i), ]
  pqtl_gene[i] = gene_sub$gene_id[which.min(abs(gene_sub$start - bp_i))]
}
ukb_sig$near_gene_id = pqtl_gene
ukb_sig$near_gene_name = gene_meta_prot$gene_name[match(ukb_sig$near_gene, 
                                                        gene_meta_prot$gene_id)]

background_n = fread('pqtl/11_prot_complex/backgroud_gene_pair_ppi_tf.txt')
## TF enrichment
ukb_trans_pqtl = ukb_sig[ukb_sig$cis_trans == 'trans', ]
# background_tf = c()
# for (i in 1:1000) {
#   gene_pair_i = paste0(sample(ukb_trans_pqtl$near_gene_name), ',', ukb_trans_pqtl$ukb_gene_name)
#   background_tf[i] = sum(gene_pair_i %in% tf_list)
#   if (i %% 100 == 0) {
#     print(i)
#   }
# }
background_tf = background_n$tf
background_tf_prop = background_tf/nrow(ukb_trans_pqtl)
obs_tf = sum(paste0(ukb_trans_pqtl$near_gene_name, ',', ukb_trans_pqtl$ukb_gene_name) %in% 
               tf_list)

ukb_trans_w_dgn = ukb_trans_pqtl[ukb_trans_pqtl$dgn_se > 0, ]

enrich_tf = function(df, dgn_p_low, dgn_p_high) {
  df_sub = df[df$dgn_p > dgn_p_low & df$dgn_p < dgn_p_high, ]
  print(paste0('# pairs: ', nrow(df_sub)))
  gene_pair_sub = paste0(df_sub$near_gene_name, ',', df_sub$ukb_gene_name)
  tf_prop = sum(gene_pair_sub %in% tf_list) / nrow(df_sub)
  return(tf_prop / background_tf_prop)
}
tf_enrich_dgn1 = enrich_tf(ukb_trans_w_dgn, 0, 0.01)
tf_enrich_dgn2 = enrich_tf(ukb_trans_w_dgn, 0.01, 0.05)
tf_enrich_dgn3 = enrich_tf(ukb_trans_w_dgn, 0.05, 0.1)
tf_enrich_dgn4 = enrich_tf(ukb_trans_w_dgn, 0.1, 1)
mean(tf_enrich_dgn1)
tf_table = data.frame(enrich = c(mean(tf_enrich_dgn1), mean(tf_enrich_dgn2),
                                 mean(tf_enrich_dgn3), mean(tf_enrich_dgn4)),
                      ci_low = c(quantile(tf_enrich_dgn1, 0.025), quantile(tf_enrich_dgn2, 0.025), 
                                 quantile(tf_enrich_dgn3, 0.025), quantile(tf_enrich_dgn4, 0.025)),
                      ci_high = c(quantile(tf_enrich_dgn1, 0.975), quantile(tf_enrich_dgn2, 0.975), 
                                  quantile(tf_enrich_dgn3, 0.975), quantile(tf_enrich_dgn4, 0.975)),
                      dgn_p = c('q1', 'q2', 'q3', 'q4'),
                      type = 'TF')

ggplot(tf_table, aes(x = factor(dgn_p, levels = c('q1', 'q2', 'q3', 'q4')),
                     y = enrich)) + 
  geom_point(size = 2) +
  # geom_hline(yintercept = enrich_table$enrich[4], lty = 2, color = "#00BFC4") + 
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.05) +
  labs(x = "DGN p value", y = "enrich", color = '',
       title = 'Enrichment for nearby genes of trans-pQTLs in TF') + 
  scale_x_discrete(labels=c("< 0.001", "0.001 ~ 0.01", "0.01 ~ 0.1", "> 0.1")) + 
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
bioplex = fread('pqtl/11_prot_complex/BioPlex_293T_Network_10K_Dec_2019.tsv')
bioplex = unique(bioplex[, c('SymbolA', 'SymbolB')])
colnames(bioplex) = c('prot1', 'prot2')

# CORUM database
corum = read.xlsx('pqtl/11_prot_complex/CORUM_2022_09_12.xlsx')
corum = corum[corum$Organism == 'Human', ]
corum_prot_list = strsplit(corum$`subunits(Gene.name)`, ';', fixed = T)
corum_n_prot = sapply(corum_prot_list, function(x){length(unique(x))})
corum_prot_list = corum_prot_list[which(corum_n_prot > 1)]
corum_ppi = sapply(corum_prot_list, function(x){combn(x, 2)})
corum_ppi = sapply(corum_ppi, t)
corum = as.data.frame(do.call(rbind, corum_ppi))
corum = unique(corum)
colnames(corum) = c('prot1', 'prot2')

# STRING database
string_prot = fread('pqtl/11_prot_complex/string_cluster/9606.protein.links.detailed.v12.0.txt.gz')
string_prot_info = fread('pqtl/11_prot_complex/string_cluster/9606.protein.info.v12.0.txt.gz')
## use only known interactions (experimental and database)
string_prot = string_prot[, c('protein1', 'protein2', 'experimental', 'database')]
string_prot$experimental = string_prot$experimental/1000
string_prot$database = string_prot$database/1000

## calculate the combined score (https://string-db.org/cgi/help?sessionId=b2PwxBVzniM9)
p = 0.041
s_exp_nop = (string_prot$experimental - p) / (1 - p)
s_db_nop = (string_prot$database - p) / (1 - p)
s_tot_nop = 1 - (1 - s_exp_nop) * (1 - s_db_nop)
s_tot = s_tot_nop + p * (1 - s_tot_nop)
string_prot$combined_score = s_tot
string_fine_ppi = string_prot[string_prot$combined_score > 0.75, ]
string = unique(string_fine_ppi[, 1:2])
colnames(string) = c('prot1', 'prot2')
string$prot1 = string_prot_info$preferred_name[match(string$prot1, 
                                                     string_prot_info$`#string_protein_id`)]
string$prot2 = string_prot_info$preferred_name[match(string$prot2, 
                                                     string_prot_info$`#string_protein_id`)]

# HIPPIE database
hippie = fread('pqtl/11_prot_complex/hippie_current.txt')
hippie = hippie[hippie$V5 >= 0.75, c(1,3)]
hippie$V1 = sapply(strsplit(hippie$V1, '_', fixed = T), '[', 1)
hippie$V3 = sapply(strsplit(hippie$V3, '_', fixed = T), '[', 1)
hippie = unique(na.omit(hippie))
hippie = hippie[hippie$V1 != hippie$V3, ]
colnames(hippie) = c('prot1', 'prot2')

# all ppi
ppi_all = rbind(bioplex, corum, string, hippie)
ppi_all = unique(ppi_all)
ppi_list = c(paste0(ppi_all$prot1, ',', ppi_all$prot2),
             paste0(ppi_all$prot2, ',', ppi_all$prot1))

# background_ppi = c()
# for (i in 1:1000) {
#   gene_pair_i = paste0(sample(ukb_trans_pqtl$near_gene_name), ',', ukb_trans_pqtl$ukb_gene_name)
#   background_ppi[i] = sum(gene_pair_i %in% ppi_list)
#   if (i %% 100 == 0) {
#     print(i)
#   }
# }
background_ppi = background_n$ppi
background_ppi_prop = background_ppi/nrow(ukb_trans_pqtl)

obs_ppi = sum(paste0(ukb_trans_pqtl$near_gene_name, ',', ukb_trans_pqtl$ukb_gene_name) %in% 
                ppi_list)

enrich_ppi = function(df, dgn_p_low, dgn_p_high) {
  df_sub = df[df$dgn_p > dgn_p_low & df$dgn_p < dgn_p_high, ]
  print(paste0('# pairs: ', nrow(df_sub)))
  gene_pair_sub = paste0(df_sub$near_gene_name, ',', df_sub$ukb_gene_name)
  ppi_prop = sum(gene_pair_sub %in% ppi_list) / nrow(df_sub)
  return(ppi_prop / background_ppi_prop)
}
ppi_enrich_dgn1 = enrich_ppi(ukb_trans_w_dgn, 0, 0.01)
ppi_enrich_dgn2 = enrich_ppi(ukb_trans_w_dgn, 0.01, 0.05)
ppi_enrich_dgn3 = enrich_ppi(ukb_trans_w_dgn, 0.05, 0.1)
ppi_enrich_dgn4 = enrich_ppi(ukb_trans_w_dgn, 0.1, 1)

ppi_table = data.frame(enrich = c(mean(ppi_enrich_dgn1), mean(ppi_enrich_dgn2),
                                 mean(ppi_enrich_dgn3), mean(ppi_enrich_dgn4)),
                      ci_low = c(quantile(ppi_enrich_dgn1, 0.025), quantile(ppi_enrich_dgn2, 0.025), 
                                 quantile(ppi_enrich_dgn3, 0.025), quantile(ppi_enrich_dgn4, 0.025)),
                      ci_high = c(quantile(ppi_enrich_dgn1, 0.975), quantile(ppi_enrich_dgn2, 0.975), 
                                  quantile(ppi_enrich_dgn3, 0.975), quantile(ppi_enrich_dgn4, 0.975)),
                      dgn_p = c('q1', 'q2', 'q3', 'q4'),
                      type = 'PPI')

ggplot(ppi_table, aes(x = factor(dgn_p, levels = c('q1', 'q2', 'q3', 'q4')),
                     y = enrich)) + 
  geom_point(size = 2) +
  # geom_hline(yintercept = enrich_table$enrich[4], lty = 2, color = "#00BFC4") + 
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.05) +
  labs(x = "DGN p value", y = "enrich", color = '',
       title = 'Enrichment for nearby genes of trans-pQTLs in ppi') + 
  scale_x_discrete(labels=c("< 0.001", "0.001 ~ 0.01", "0.01 ~ 0.1", "> 0.1")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### plot tf and ppi together
enrich_all = data.frame(enrich = c(mean(obs_tf/background_tf), mean(obs_ppi/background_ppi)),
                        ci_low = c(quantile(obs_tf/background_tf, 0.025), quantile(obs_ppi/background_ppi, 0.025)),
                        ci_high = c(quantile(obs_tf/background_tf, 0.975), quantile(obs_ppi/background_ppi, 0.975)),
                        type = c('TF', 'PPI'))
ggplot(enrich_all, aes(x = type, y = enrich, fill = type)) + 
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

enrich_by_dgn = rbind(tf_table, ppi_table)
ggplot(enrich_by_dgn, 
       aes(x = factor(dgn_p, levels = c('q1', 'q2', 'q3', 'q4')),
           y = enrich, color = type, group = type)) + 
  geom_point(size = 3) +
  geom_line(linetype = "dashed") + 
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.05) +
  labs(x = "DGN p value", y = "Enrichment", color = '',
       title = '') + 
  scale_x_discrete(labels=c("< 0.01", "0.01 ~ 0.05", "0.05 ~ 0.1", "> 0.1")) + 
  scale_color_manual(values = brewer.pal(3, "Set1")[1:2]) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

# background_n = data.frame(ppi = background_ppi,
#                           tf = background_tf)
# fwrite(background_n, 'pqtl/11_prot_complex/backgroud_gene_pair_ppi_tf.txt')
# library(scales)
# hue_pal()(2)




