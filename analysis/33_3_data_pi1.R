library(data.table)
library(qvalue)
setwd('/project/xuanyao/jinghui')

### DGN eQTL replication in eQTLGen
# cis
eqtlgen_z_dgn_cis = fread('pqtl/12_beta_across_two_data/eqtlgen_z_dgn_cis.txt')
colnames(eqtlgen_z_dgn_cis) = c("dgn_var", "dgn_gene", "dgn_p", "dgn_beta",
                                "dgn_se", "dgn_af", 'dgn_chr', 'dgn_bp', 
                                "dgn_A0", "dgn_A1", "eqtlgen_A1", "eqtlgen_A0", "eqtlgen_z", 
                                "eqtlgen_p", 'eqtlgen_fdr', 'eqtlgen_BonP')
# clump SNPs
uniq_gene = unique(eqtlgen_z_dgn_cis$dgn_gene)
eqtlgen_z_dgn_cis_loci = c()
window_size = 250000
for (i in uniq_gene) {
  eqtlgen_sub = eqtlgen_z_dgn_cis[eqtlgen_z_dgn_cis$dgn_gene == i, ]
  while(nrow(eqtlgen_sub) > 0) {
    lead_snp = which.min(eqtlgen_sub$dgn_p)[1]
    lead_chr = eqtlgen_sub$dgn_chr[lead_snp]
    lead_pos = eqtlgen_sub$dgn_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(eqtlgen_sub$dgn_chr == lead_chr & 
                       eqtlgen_sub$dgn_bp >= start_pos &
                       eqtlgen_sub$dgn_bp <= end_pos)
    eqtlgen_z_dgn_cis_loci = rbind(eqtlgen_z_dgn_cis_loci, eqtlgen_sub[lead_snp, ])
    eqtlgen_sub = eqtlgen_sub[-rm_index, ]
  }
}
eqtlgen_z_dgn_cis_loci = na.omit(eqtlgen_z_dgn_cis_loci)
eqtlgen_z_dgn_cis_str = eqtlgen_z_dgn_cis_loci[eqtlgen_z_dgn_cis_loci$dgn_p < 5e-8/12585, ]
qval_eqtlgen_z_dgn_cis_str = qvalue(eqtlgen_z_dgn_cis_str$eqtlgen_p)
sum(eqtlgen_z_dgn_cis_str$eqtlgen_p < 0.05/nrow(eqtlgen_z_dgn_cis_str)) /
  nrow(eqtlgen_z_dgn_cis_str)

# trans
eqtlgen_z_dgn_trans = fread('pqtl/12_beta_across_two_data/eqtlgen_z_dgn_trans.txt')
colnames(eqtlgen_z_dgn_trans) = c("dgn_var", "dgn_gene", "dgn_p", "dgn_beta", "dgn_se", 
                                "dgn_af", "dgn_gene_chr", "dgn_gene_start", "dgn_gene_end", 
                                'dgn_chr', 'dgn_bp', "dgn_A0", "dgn_A1", 'dgn_gene_id',
                                 "eqtlgen_A1", "eqtlgen_A0", "eqtlgen_z", 
                                "eqtlgen_p", 'eqtlgen_fdr', 'eqtlgen_BonP')
# clump SNPs
uniq_gene = unique(eqtlgen_z_dgn_trans$dgn_gene)
eqtlgen_z_dgn_trans_loci = c()
window_size = 250000
for (i in uniq_gene) {
  eqtlgen_sub = eqtlgen_z_dgn_trans[eqtlgen_z_dgn_trans$dgn_gene == i, ]
  while(nrow(eqtlgen_sub) > 0) {
    lead_snp = which.min(eqtlgen_sub$dgn_p)[1]
    lead_chr = eqtlgen_sub$dgn_chr[lead_snp]
    lead_pos = eqtlgen_sub$dgn_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(eqtlgen_sub$dgn_chr == lead_chr & 
                       eqtlgen_sub$dgn_bp >= start_pos &
                       eqtlgen_sub$dgn_bp <= end_pos)
    eqtlgen_z_dgn_trans_loci = rbind(eqtlgen_z_dgn_trans_loci, eqtlgen_sub[lead_snp, ])
    eqtlgen_sub = eqtlgen_sub[-rm_index, ]
  }
}
eqtlgen_z_dgn_trans_loci = na.omit(eqtlgen_z_dgn_trans_loci)
eqtlgen_z_dgn_trans_str = eqtlgen_z_dgn_trans_loci[eqtlgen_z_dgn_trans_loci$dgn_p < 5e-8/12585, ]
qval_eqtlgen_z_dgn_trans_str = qvalue(eqtlgen_z_dgn_trans_str$eqtlgen_p)
sum(eqtlgen_z_dgn_trans_str$eqtlgen_p < 0.05/nrow(eqtlgen_z_dgn_trans_str)) /
  nrow(eqtlgen_z_dgn_trans_str)

### eQTLGen eQTL replication in DGN
# cis
dgn_beta_eqtlgen_cis = fread('pqtl/12_beta_across_two_data/dgn_beta_eqtlgen_cis.txt')
colnames(dgn_beta_eqtlgen_cis) = c("eqtlgen_p", "eqtlgen_snp", "eqtlgen_chr", "eqtlgen_bp",
                                "eqtlgen_A1", "eqtlgen_A0", 'eqtlgen_z', 'eqtlgen_gene_id', 
                                "eqtlgen_gene", "eqtlgen_gene_chr", "eqtlgen_gene_start", 
                                "eqtlgen_Nrcohort", "eqtlgen_Nrsample", 'eqtlgen_fdr', 'eqtlgen_BonP',
                                'dgn_beta', 'dgn_se', 'dgn_p')
# clump SNPs
uniq_gene = unique(dgn_beta_eqtlgen_cis$eqtlgen_gene)
dgn_beta_eqtlgen_cis_loci = c()
window_size = 250000
for (i in uniq_gene) {
  dgn_sub = dgn_beta_eqtlgen_cis[dgn_beta_eqtlgen_cis$eqtlgen_gene == i, ]
  while(nrow(dgn_sub) > 0) {
    lead_snp = which.min(dgn_sub$eqtlgen_p)[1]
    lead_chr = dgn_sub$eqtlgen_chr[lead_snp]
    lead_pos = dgn_sub$eqtlgen_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(dgn_sub$eqtlgen_chr == lead_chr & 
                       dgn_sub$eqtlgen_bp >= start_pos &
                       dgn_sub$eqtlgen_bp <= end_pos)
    dgn_beta_eqtlgen_cis_loci = rbind(dgn_beta_eqtlgen_cis_loci, dgn_sub[lead_snp, ])
    dgn_sub = dgn_sub[-rm_index, ]
  }
}
dgn_beta_eqtlgen_cis_loci = na.omit(dgn_beta_eqtlgen_cis_loci)
dgn_beta_eqtlgen_cis_str = dgn_beta_eqtlgen_cis_loci[dgn_beta_eqtlgen_cis_loci$eqtlgen_p < 5e-8/19942, ]
qval_dgn_beta_eqtlgen_cis_str = qvalue(dgn_beta_eqtlgen_cis_str$dgn_p)
1-qval_dgn_beta_eqtlgen_cis_str$pi0

# trans
dgn_beta_eqtlgen_trans = fread('pqtl/12_beta_across_two_data/dgn_beta_eqtlgen_trans.txt')
colnames(dgn_beta_eqtlgen_trans) = c("eqtlgen_p", "eqtlgen_snp", "eqtlgen_chr", "eqtlgen_bp",
                                   "eqtlgen_A1", "eqtlgen_A0", 'eqtlgen_z', 'eqtlgen_gene_id', 
                                   "eqtlgen_gene", "eqtlgen_gene_chr", "eqtlgen_gene_start", 
                                   "eqtlgen_Nrcohort", "eqtlgen_Nrsample", 'eqtlgen_fdr', 'eqtlgen_BonP',
                                   'dgn_beta', 'dgn_se', 'dgn_p')
# clump SNPs
uniq_gene = unique(dgn_beta_eqtlgen_trans$eqtlgen_gene)
dgn_beta_eqtlgen_trans_loci = c()
window_size = 250000
for (i in uniq_gene) {
  dgn_sub = dgn_beta_eqtlgen_trans[dgn_beta_eqtlgen_trans$eqtlgen_gene == i, ]
  while(nrow(dgn_sub) > 0) {
    lead_snp = which.min(dgn_sub$eqtlgen_p)[1]
    lead_chr = dgn_sub$eqtlgen_chr[lead_snp]
    lead_pos = dgn_sub$eqtlgen_bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(dgn_sub$eqtlgen_chr == lead_chr & 
                       dgn_sub$eqtlgen_bp >= start_pos &
                       dgn_sub$eqtlgen_bp <= end_pos)
    dgn_beta_eqtlgen_trans_loci = rbind(dgn_beta_eqtlgen_trans_loci, dgn_sub[lead_snp, ])
    dgn_sub = dgn_sub[-rm_index, ]
  }
}
dgn_beta_eqtlgen_trans_loci = na.omit(dgn_beta_eqtlgen_trans_loci)
dgn_beta_eqtlgen_trans_str = dgn_beta_eqtlgen_trans_loci[dgn_beta_eqtlgen_trans_loci$eqtlgen_p < 5e-8/19942, ]
qval_dgn_beta_eqtlgen_trans_str = qvalue(dgn_beta_eqtlgen_trans_str$dgn_p)
1-qval_dgn_beta_eqtlgen_trans_str$pi0

## plot ukb, dgn and eQTLGen together
library(corrplot)
pi_cis = rbind(c(1, 0.8, 0.69),
               c(0.84, 1, 0.94),
               c(0.86, 0.98, 1))
colnames(pi_cis) = c('UKB-PPP \n pQTL', 'DGN \n eQTL', 'eQTLGen \n eQTL')
rownames(pi_cis) = c('UKB-PPP \n pQTL', 'DGN \n eQTL', 'eQTLGen \n eQTL')
#par(mfrow=c(1,2))
corrplot(pi_cis, method = 'color', is.corr = FALSE, col.lim = c(0, 1), diag = F,
         addCoef.col = "black", tl.col = "black", tl.srt = 45, tl.cex = 1.1, cl.cex = 1, 
         col=colorRampPalette(c("white","dodgerblue3"))(100),
         mar=c(0,0,2,0), cex.main = 2, font.main= 4, tl.pos = 'd', number.cex = 1.3)

pi_trans = rbind(c(1, 0.77, 0.5),
               c(0.11, 1, 0.78),
               c(0.24, 1, 1))
colnames(pi_trans) = c('UKB-PPP \n pQTL', 'DGN \n eQTL', 'eQTLGen \n eQTL')
rownames(pi_trans) = c('UKB-PPP \n pQTL', 'DGN \n eQTL', 'eQTLGen \n eQTL')
corrplot(pi_trans, method = 'color', is.corr = FALSE,  col.lim = c(0, 1), diag = F,
         addCoef.col = "black", tl.col = "black", tl.srt = 45, tl.cex = 1.1, cl.cex = 1,
         col=colorRampPalette(c("white","dodgerblue3"))(100),
         mar=c(0,0,2,0), cex.main = 2, font.main= 4, tl.pos = 'd', number.cex = 1.3)




