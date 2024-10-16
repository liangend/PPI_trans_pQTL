library(data.table)
library(qvalue)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
ukb_p_csf_sig = fread('pqtl/12_beta_across_two_data/ukb_beta_csf_sig_all.txt')
colnames(ukb_p_csf_sig) = c("csf_snp", "csf_A1", "csf_A2", "csf_logP", "csf_rsid", 
                            "csf_file", "csf_chr", "csf_bp", "csf_loci", "csf_loci_start", 
                            "csf_loci_end", 'csf_gene_id', 'ukb_beta', 'ukb_se', 'ukb_logP')
gene_meta = fread('gtex/00_ref/genecode.GRCh38.gene.meta.unadj.gtf')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
ukb_p_csf_sig$csf_gene_chr = gene_meta$chr[match(ukb_p_csf_sig$csf_gene_id, gene_meta$gene_id)]
ukb_p_csf_sig$csf_tss = gene_meta$start[match(ukb_p_csf_sig$csf_gene_id, gene_meta$gene_id)]
ukb_p_csf_sig = na.omit(ukb_p_csf_sig)
ukb_p_csf_sig$csf_cis_trans = 'trans'
ukb_p_csf_sig$csf_cis_trans[which(ukb_p_csf_sig$csf_chr == ukb_p_csf_sig$csf_gene_chr &
                                    abs(ukb_p_csf_sig$csf_bp - ukb_p_csf_sig$csf_tss) < 1000000)] = 'cis'
ukb_p_csf_sig = ukb_p_csf_sig[ukb_p_csf_sig$ukb_logP != 1, ]
ukb_p_csf_sig = ukb_p_csf_sig[ukb_p_csf_sig$csf_logP < log(5e-8/7028), ]

csf_p_ukb_sig = fread('pqtl/12_beta_across_two_data/csf_beta_ukb_sig_all.txt')
colnames(csf_p_ukb_sig) = c("ukb_chr", "ukb_bp_hg38", "ukb_bp_hg19", "ukb_A0", "ukb_A1", 
                            "ukb_A1_freq", "ukb_N", "ukb_beta", "ukb_se", "ukb_logP", 
                            "ukb_file", "ukb_gene_chr", "ukb_tss", "ukb_cis_trans", 
                            "ukb_gene_id", "ukb_gene_name", "ukb_loci", "ukb_loci_start", 
                            "ukb_loci_end", "csf_p")
csf_p_ukb_sig = csf_p_ukb_sig[csf_p_ukb_sig$csf_p != 1, ]
csf_p_ukb_sig = csf_p_ukb_sig[csf_p_ukb_sig$ukb_logP > -log(5e-8/2837), ]

### comparing ukb cis-pqtls with csf cis-pqtls
csf_p_ukb_cis = csf_p_ukb_sig[csf_p_ukb_sig$ukb_cis_trans == 'cis', ]
ukb_p_csf_cis = ukb_p_csf_sig[ukb_p_csf_sig$csf_cis_trans == 'cis', ]

# cis-pQTL in ukb plasma replicated as cis-pQTL in csf
qval_csf_cis = qvalue(csf_p_ukb_cis$csf_p)
1 - qval_csf_cis$pi0

# cis-pQTL in csf replicated as cis-pQTL in ukb plasma
qval_ukb_cis = qvalue(10^-(ukb_p_csf_cis$ukb_logP))
1 - qval_ukb_cis$pi0

### comparing ukb trans-pqtls with csf trans-pqtls 
csf_p_ukb_trans = csf_p_ukb_sig[csf_p_ukb_sig$ukb_cis_trans == 'trans', ]
ukb_p_csf_trans = ukb_p_csf_sig[ukb_p_csf_sig$csf_cis_trans == 'trans', ]

# trans-pQTL replication
qval_ukb_trans = qvalue(10^-(ukb_p_csf_trans$ukb_logP))
1-qval_ukb_trans$pi0

qval_csf_trans = qvalue(csf_p_ukb_trans$csf_p)
1-qval_csf_trans$pi0

## cis trans pi1 summary
pi_table = data.frame(pi1 = c(1-qval_ukb_cis$pi0, 1-qval_ukb_trans$pi0, 
                              1-qval_csf_cis$pi0, 1-qval_csf_trans$pi0), 
                      rep_group = c('csf in ukb plasma',
                                    'csf in ukb plasma',
                                    'ukb plasma in csf',
                                    'ukb plasma in csf'),
                      cis_trans = c('cis', 'trans', 'cis', 'trans'))
ggplot(pi_table, aes(x = rep_group, y = pi1, fill = cis_trans, label = round(pi1, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.4) +
  geom_text(position=position_dodge(0.4), vjust=-0.2, size = 4.5) + 
  labs(x = "", fill = '', title = 'ukb vs csf (p < 1e-8)') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



