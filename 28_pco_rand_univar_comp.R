library(data.table)
setwd('/project/xuanyao/jinghui')
args = commandArgs(trailingOnly = T)
chr_i = as.numeric(args[1])

pco_loci = fread('pqtl/04_fdr/ukb/pco_random_loci.txt')
pco_loci = as.data.frame(pco_loci)

prot_mod = fread('pqtl/11_prot_complex/random_mod.txt')
prot_mod$mod = paste0('mod', prot_mod$mod)

loci_i = pco_loci[pco_loci$chr == chr_i, ]
mod_i = unique(loci_i$mod)
print(paste0('# mod: ', length(mod_i)))
for (j in 1:length(mod_i)) {
  mod_j = mod_i[j]
  loci_j = loci_i[loci_i$mod == mod_j, ]
  
  # prot_i = prot_mod[prot_mod$mod == mod_j, ]
  # prot_i = prot_i[prot_i$chr != chr_i, ]
  
  # minp_all = c()
  # for (k in 1:nrow(prot_i)) {
  #   file_k = list.files(paste0('pqtl/UKB_PPP/', prot_i$file[k]))
  #   file_k = file_k[grep(paste0('chr', chr_i, '_'), file_k)]
  #   summ_stat_k = fread(paste0('pqtl/UKB_PPP/', prot_i$file[k], '/', file_k))
  #   summ_stat_k$bp_hg37 = as.numeric(sapply(strsplit(summ_stat_k$ID, ':', fixed = T), '[', 2))
  #   minp_k = c()
  #   for (n in 1:nrow(loci_j)) {
  #     start_n = loci_j$loci_start[n]
  #     end_n = loci_j$loci_end[n]
  #     minp_k[n] = max(summ_stat_k$LOG10P[which(summ_stat_k$bp_hg37 >= start_n & 
  #                                                summ_stat_k$bp_hg37 <= end_n)])
  #   }
  #   minp_all = rbind(minp_all, minp_k)
  # }
  # max_univar_logP = apply(minp_all, 2, max)
  
  max_univar_logP = c()
  z_mat_j = fread(paste0('pqtl/01_zmat/ukb_random/', mod_j, '.chr', chr_i, '.txt.gz'))
  bp_j = as.numeric(sapply(strsplit(z_mat_j$V1, ':', fixed = T), '[', 2))
  ## remove proteins on the same chr
  prot_mod_j = prot_mod[prot_mod$mod == mod_j, ]
  prot_j = prot_mod_j$prot[prot_mod_j$chr != chr_i]
  z_mat_j = as.data.frame(z_mat_j)
  z_mat_j = z_mat_j[, c('V1', prot_j)]
  
  for (k in 1:nrow(loci_j)) {
    start_k = loci_j$loci_start[k]
    end_k = loci_j$loci_end[k]
    z_mat_sub = z_mat_j[which(bp_j >= start_k & bp_j <= end_k), -1]
    log_p = pnorm(max(abs(z_mat_sub)), lower.tail = F, log.p = T) + log(2) # transfer p to ln(P)
    max_univar_logP[k] = -log_p * log10(exp(1)) # transfer ln(P) to log10(P)
  }
  loci_j$max_univar_logP = max_univar_logP
  fwrite(loci_j, paste0('pqtl/04_fdr/ukb/pco_univar_comp/chr', chr_i, mod_j, '.txt'), 
         sep = '\t', col.names = F)
  print(j)
  gc()
}


