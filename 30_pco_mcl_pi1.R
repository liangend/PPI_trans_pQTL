library(data.table)
setwd('/project/xuanyao/jinghui')
args = commandArgs(trailingOnly = T)
chr_i = as.numeric(args[1])

pco_loci = fread('pqtl/04_fdr/ukb/pco_mcl_loci.txt')
pco_loci = as.data.frame(pco_loci)

prot_mod = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')
prot_mod$mod = paste0('mod', prot_mod$mod)

loci_i = pco_loci[pco_loci$chr == chr_i, ]
mod_i = unique(loci_i$mod)
print(paste0('# mod: ', length(mod_i)))
for (j in 1:length(mod_i)) {
  mod_j = mod_i[j]
  loci_j = loci_i[loci_i$mod == mod_j, ]
  
  univar_prot = c()
  univar_z = c()
  z_mat_j = fread(paste0('pqtl/01_zmat/ukb/', mod_j, '.chr', chr_i, '.txt.gz'))
  bp_j = as.numeric(sapply(strsplit(z_mat_j$V1, ':', fixed = T), '[', 2))
  for (k in 1:nrow(loci_j)) {
    bp_k = loci_j$bp_hg37[k]
    z_mat_sub = z_mat_j[bp_j == bp_k, -1]
    univar_prot[k] = paste(colnames(z_mat_sub), collapse = ',')
    univar_z[k] = paste(unlist(z_mat_sub), collapse = ',')
  }
  loci_j$prot = univar_prot
  loci_j$univar_z = univar_z
  fwrite(loci_j, paste0('pqtl/04_fdr/ukb/pco_univar_comp/chr', chr_i, mod_j, '.txt'), 
         sep = '\t', col.names = F)
  print(j)
  gc()
}


