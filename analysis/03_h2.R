library(data.table)
h2_list = list.files('/project/xuanyao/jinghui/pqtl/05_h2/02_h2_53annot/cis_snp_wo_intercept/', pattern = 'results')

coeff_all = c()
nSNP_all = c()
i = h2_list[1]
for (i in h2_list) {
  h2_file = fread(paste0('/project/xuanyao/jinghui/pqtl/05_h2/02_h2_53annot/all_snp_wo_intercept/', i))
  coeff_all = rbind(coeff_all, h2_file$Coefficient)
  nSNP_all = rbind(nSNP_all, h2_file$M_annot)
}
coeff_all_ave = apply(coeff_all, 2, mean)
h2_all = nSNP_all %*% coeff_all_ave

coeff_cis = c()
nSNP_cis = c()
for (i in h2_list) {
  h2_file = fread(paste0('/project/xuanyao/jinghui/pqtl/05_h2/02_h2_53annot/cis_chr_snp_wo_intercept/', i))
  coeff_cis = rbind(coeff_cis, h2_file$Coefficient)
  nSNP_cis = rbind(nSNP_cis, h2_file$M_annot)
}
coeff_cis_ave = apply(coeff_cis, 2, mean)
h2_cis = nSNP_cis %*% coeff_cis_ave

coeff_trans = c()
nSNP_trans = c()
for (i in h2_list) {
  h2_file = fread(paste0('/project/xuanyao/jinghui/pqtl/05_h2/02_h2_53annot/trans_snp_wo_intercept/', i))
  coeff_trans = rbind(coeff_trans, h2_file$Coefficient)
  nSNP_trans = rbind(nSNP_trans, h2_file$M_annot)
}
coeff_trans_ave = apply(coeff_trans, 2, mean)
h2_trans = nSNP_trans %*% coeff_trans_ave

boxplot(h2_all[,1], h2_cis[,1], h2_trans[,1], ylim = c(0,0.15))
mean(h2_cis[,1]/h2_all[,1]) + 
mean(h2_trans[,1]/h2_all[,1])
mean(h2_all[,1])
mean(h2_cis[,1])
mean(h2_trans[,1])

library(ggplot2)
h2_plot = data.frame(h2 = h2_all, cat = 'all')
ggplot(h2_plot, aes(x=cat, y=h2)) + 
  geom_boxplot(width = 0.4) +
  labs(x = "", y = "h2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




