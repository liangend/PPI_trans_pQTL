library(qvalue)
library(data.table)
trans_pqtl = fread('/project/xuanyao/jinghui/pqtl/04_fdr/sig_pair_all.txt')
pi1 = c()
for (i in 1:nrow(trans_pqtl)) {
  mod_i = trans_pqtl$mod[i]
  chr_i = trans_pqtl$chr[i]
  if (i == 1) {
    zmat_i = fread(paste0('/project/xuanyao/jinghui/pqtl/01_zmat/', mod_i, '.chr', chr_i, '.txt.gz'))
  } else {
    mod_i1 = trans_pqtl$mod[i-1]
    chr_i1 = trans_pqtl$chr[i-1]
    if (mod_i != mod_i1 | chr_i != chr_i1) {
      zmat_i = fread(paste0('/project/xuanyao/jinghui/pqtl/01_zmat/', mod_i, '.chr', chr_i, '.txt.gz'))
    }
  }
  z_i = unname(unlist(zmat_i[which(zmat_i$V1 == trans_pqtl$snp[i]), -1]))
  qval_i = try(qvalue(p = pnorm(abs(z_i), lower.tail = F) * 2), silent = T)  # pi0 estimate may not be available for small
  if (length(qval_i) == 8) {
    pi1[i] = 1 - qval_i$pi0
  } else {
    pi1[i] = NA
  }
  
  print(i)
}

trans_pqtl$pi1 = pi1
fwrite(trans_pqtl, '/project/xuanyao/jinghui/pqtl/04_fdr/sig_pair_all.txt', sep = '\t')

