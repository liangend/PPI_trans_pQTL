library(data.table)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
### UKB cis vs. DGN cis
dgn_beta = fread('pqtl/12_beta_across_two_data/dgn_beta_all.txt')

adj_dgn_beta = c()
for (i in 1:nrow(dgn_beta)) {
  ukb_A1 = dgn_beta$ukb_A1[i]
  ukb_A2 = dgn_beta$ukb_A2[i]
  dgn_A1 = dgn_beta$dgn_A1[i]
  dgn_A2 = dgn_beta$dgn_A2[i]
  if (ukb_A1 == dgn_A1 & ukb_A2 == dgn_A2) {
    adj_dgn_beta[i] = -dgn_beta$dgn_beta[i]
  } else if (ukb_A1 == dgn_A2 & ukb_A2 == dgn_A1) {
    adj_dgn_beta[i] = dgn_beta$dgn_beta[i]
  } else {
    adj_dgn_beta[i] = 0
  }
}
dgn_beta$adj_dgn_beta = adj_dgn_beta
sum(adj_dgn_beta == 0)

sum(sign(dgn_beta$beta) == sign(dgn_beta$adj_dgn_beta)) / nrow(dgn_beta)

dgn_beta = as.data.frame(dgn_beta)
dgn_cis_beta = dgn_beta[which(dgn_beta$cis_trans == 'cis'), ]
sum(sign(dgn_cis_beta$beta) == sign(dgn_cis_beta$adj_dgn_beta)) / nrow(dgn_cis_beta)
p_cis = pnorm(abs(dgn_cis_beta$dgn_beta/dgn_cis_beta$dgn_se), lower.tail = F) * 2
dgn_cis_beta_sig = dgn_cis_beta[which(p_cis < 1e-5 & dgn_cis_beta$neg_logP > 5), ]
sum(sign(dgn_cis_beta_sig$beta) == sign(dgn_cis_beta_sig$adj_dgn_beta)) / nrow(dgn_cis_beta_sig)
cor(dgn_cis_beta_sig$beta / dgn_cis_beta_sig$se, 
     dgn_cis_beta_sig$adj_dgn_beta / dgn_cis_beta_sig$dgn_se)
par(mfrow=c(1,2))
plot(dgn_cis_beta_sig$beta / dgn_cis_beta_sig$se, 
     dgn_cis_beta_sig$adj_dgn_beta / dgn_cis_beta_sig$dgn_se, 
     xlab = 'ukb prot cis z', ylab = 'DGN gene cis z')
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

### UKB cis vs. eQTLgen cis
ukb_small_cis_eqtlgen = fread('pqtl/12_beta_across_two_data/eqtlgen_z_ukb_cis_all.txt')
ukb_small_cis_eqtlgen_sig = ukb_small_cis_eqtlgen[ukb_small_cis_eqtlgen$FDR < 0.05, ]
adj_eqtlgen_z = ukb_small_cis_eqtlgen_sig$Zscore * 
  (ukb_small_cis_eqtlgen_sig$ALLELE0 == ukb_small_cis_eqtlgen_sig$OtherAllele &
     ukb_small_cis_eqtlgen_sig$ALLELE1 == ukb_small_cis_eqtlgen_sig$AssessedAllele) - 
  ukb_small_cis_eqtlgen_sig$Zscore * 
  (ukb_small_cis_eqtlgen_sig$ALLELE1 == ukb_small_cis_eqtlgen_sig$OtherAllele &
     ukb_small_cis_eqtlgen_sig$ALLELE0 == ukb_small_cis_eqtlgen_sig$AssessedAllele)
sum(sign(adj_eqtlgen_z) == sign(ukb_small_cis_eqtlgen_sig$BETA)) / nrow(ukb_small_cis_eqtlgen_sig)
cor(adj_eqtlgen_z, ukb_small_cis_eqtlgen_sig$BETA/ukb_small_cis_eqtlgen_sig$SE)
ukb_small_cis_eqtlgen_sig$adj_eqtlgen_z = adj_eqtlgen_z
uniq_gene = unique(ukb_small_cis_eqtlgen_sig$gene_id)
eqtlgen_gene = c()
for (i in uniq_gene) {
  ukb_small_sub = ukb_small_cis_eqtlgen_sig[ukb_small_cis_eqtlgen_sig$gene_id == i, ]
  eqtlgen_gene = rbind(eqtlgen_gene, ukb_small_sub[which.max(ukb_small_sub$LOG10P), ])
}
cor(eqtlgen_gene$BETA/eqtlgen_gene$SE, eqtlgen_gene$adj_eqtlgen_z)
plot(eqtlgen_gene$BETA/eqtlgen_gene$SE, eqtlgen_gene$adj_eqtlgen_z,
     xlab = 'ukb prot cis z', ylab = 'eQTLgen cis z')
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

sum(sign(eqtlgen_gene$BETA) == sign(eqtlgen_gene$adj_eqtlgen_z)) / nrow(eqtlgen_gene)


