library(data.table)
setwd('/project/xuanyao/jinghui/pqtl')
pco_sig_loci = fread('04_fdr/ukb/pco_univar_comp.txt')
prot_mod = fread('11_prot_complex/prot_module_corum_2022.txt')
prot_mod$mod = paste0('mod', prot_mod$mod)

univar_gene = c()
univar_bp = c()
univar_p = c()

for (i in 1:nrow(pco_sig_loci)) {
  mod_i = pco_sig_loci$mod[i]
  chr_i = paste0('chr', pco_sig_loci$chr[i])
  start_i = pco_sig_loci$loci_start[i]
  end_i = pco_sig_loci$loci_end[i]
  
  z_mod_i = fread(paste0('/project/xuanyao/jinghui/pqtl/01_zmat/ukb/',
                         mod_i, '.', chr_i, '.txt.gz'))
  z_bp_i = sapply(strsplit(z_mod_i$V1, ':', fixed = T), '[', 2)
  z_bp_i = as.numeric(z_bp_i)
  
  z_sub_i = z_mod_i[which(z_bp_i >= start_i & z_bp_i <= end_i), ]
  z_sub_i = as.data.frame(z_sub_i)
  z_sub_i[, -1] = abs(z_sub_i[, -1])
  minp_index = which(z_sub_i == max(abs(z_sub_i[,-1])), arr.ind = T)[1,]
  
  univar_gene[i] = colnames(z_sub_i)[minp_index[2]]
  univar_bp[i] = as.numeric(strsplit(z_sub_i$V1[minp_index[1]], ':', fixed = T)[[1]][2])
  univar_p[i] = pnorm(max(abs(z_sub_i[,-1]))[1], lower.tail = F) * 2
  print(i)
}
pco_sig_loci$univar_bp = univar_bp
pco_sig_loci$univar_gene = univar_gene
pco_sig_loci$univar_p = univar_p
fwrite(pco_sig_loci, '04_fdr/ukb/pco_univar_comp_result.txt', sep = '\t')






