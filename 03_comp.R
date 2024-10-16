library(data.table)
univar_trans_minp = fread('/project/xuanyao/jinghui/pqtl/04_fdr/univar_pco_comp.txt')
prot_meta = fread('/project/xuanyao/jinghui/pqtl/00_ref/prot_meta.txt')
pco_mod = c()
pco_minp = c()
pco_snp = c()
window_size = 100000

for (i in 1:nrow(univar_trans_minp)) {
  targe_i = univar_trans_minp$V4[i]
  mod_i = prot_meta$module[which(prot_meta$file == targe_i)]
  chr_i = univar_trans_minp$V1[i]
  bp_center = univar_trans_minp$V2[i]
  bp_start = bp_center - window_size
  bp_end = bp_center + window_size
  min_p = c()
  for (j in mod_i) {
    pco_p_j = readRDS(paste0('/project/xuanyao/jinghui/pqtl/03_p/',
                           'p.mod', j, '.', chr_i, '.rds'))
    pco_p_bp = as.numeric(sapply(strsplit(names(pco_p_j), ':', fixed = T), '[', 2))
    pco_p_loci = pco_p_j[which(pco_p_bp >= bp_start & pco_p_bp <= bp_end)]
    min_p_j = which.min(pco_p_loci)[1]
    min_p = c(min_p, pco_p_loci[min_p_j])
  }
  minp_uniq = which.min(min_p)[1]
  pco_mod[i] = mod_i[minp_uniq]
  pco_minp[i] = min_p[minp_uniq]
  pco_snp[i] = names(min_p)[minp_uniq]
  print(i)
}
univar_trans_minp$pco_minp = pco_minp
univar_trans_minp$pco_mod = pco_mod
univar_trans_minp$pco_snp = pco_snp
fwrite(univar_trans_minp, '/project/xuanyao/jinghui/pqtl/04_fdr/univar_pco_comp.txt', sep = '\t')

