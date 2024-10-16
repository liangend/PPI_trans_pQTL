library(data.table)
setwd('/project/xuanyao/jinghui')
args = commandArgs(trailingOnly = T)
trait_index = as.numeric(args[1])

gwas_trait = list.files('gtex/10_GWAS_coloc/GWAS_sum_stats/quant_trait')
gwas_trait = gwas_trait[grep('.gz', gwas_trait)]

gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/quant_trait/', gwas_trait[trait_index]))
gwas_i$POS = as.numeric(gwas_i$POS)
gwas_i$P = as.numeric(gwas_i$P)
gwas_i$CHR = as.numeric(gwas_i$CHR)
gwas_i = gwas_i[gwas_i$CHR %in% 1:22, ]
## Group SNPs into windows centered by leading SNP and flanking window with size = loci_window
loci_window = 250000
gwas_sig_loci = c()
loci_i = 1
gwas_i_sub = gwas_i
gwas_i_sub = gwas_i_sub[gwas_i_sub$P < 1e-5, ]
while (nrow(gwas_i_sub > 0)) {
  lead_sig = which.min(gwas_i_sub$P)
  lead_snp = gwas_i_sub[lead_sig, ]
  lead_snp$loci = loci_i
  lead_snp$loci_start = lead_snp$POS - loci_window
  lead_snp$loci_end = lead_snp$POS + loci_window
  row_rm = which(gwas_i_sub$CHR == lead_snp$CHR & gwas_i_sub$POS >= lead_snp$loci_start &
                   gwas_i_sub$POS <= lead_snp$loci_end)
  gwas_sig_loci = rbind(gwas_sig_loci, lead_snp)
  loci_i = loci_i + 1
  gwas_i_sub = gwas_i_sub[-row_rm, ]
}

gwas_sig_loci = gwas_sig_loci[order(gwas_sig_loci$CHR, gwas_sig_loci$POS), ]
# gwas_sig_loci = na.omit(gwas_sig_loci)

## merge overlapping windows
gwas_sig_loci_new = gwas_sig_loci[1, ]
for (i in 2:nrow(gwas_sig_loci)) {
  chr_i = gwas_sig_loci$CHR[i]
  start_i = gwas_sig_loci$loci_start[i]
  bp37_i = gwas_sig_loci$POS[i]
  
  chr_prev = gwas_sig_loci_new$CHR[nrow(gwas_sig_loci_new)]
  end_prev = gwas_sig_loci_new$loci_end[nrow(gwas_sig_loci_new)]
  bp37_prev = gwas_sig_loci_new$POS[nrow(gwas_sig_loci_new)]
  # merge MHC regions
  if (chr_i == 6 & chr_prev == 6 & bp37_i > 28477797 & bp37_i < 33448354 &
      bp37_prev > 28477797 & bp37_prev < 33448354) {
    # merge two overlapping loci
    add_pool_i = rbind(gwas_sig_loci_new[nrow(gwas_sig_loci_new), ], gwas_sig_loci[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$P), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    gwas_sig_loci_new[nrow(gwas_sig_loci_new), ] = add_loci_i
    next
  }
  # keep non-overlapping windows
  if (chr_i != chr_prev | 
      (chr_i == chr_prev & start_i > end_prev)) {
    # add the loci directly if two loci do not overlap
    gwas_sig_loci_new = rbind(gwas_sig_loci_new, gwas_sig_loci[i, ])
  } else {
    # merge two overlapping loci
    add_pool_i = rbind(gwas_sig_loci_new[nrow(gwas_sig_loci_new), ], gwas_sig_loci[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$P), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    gwas_sig_loci_new[nrow(gwas_sig_loci_new), ] = add_loci_i
  }
}
trait_i = unlist(strsplit(gwas_trait[trait_index], '_', fixed = T))[1]
fwrite(gwas_sig_loci_new, paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/quant_trait/gwas_loci_500kb/', 
                                 trait_i, '.txt'), sep = '\t')

if (!dir.exists(paste0('gtex/10_GWAS_coloc/gwas_loci/quant_trait/', trait_i))) {
  dir.create(paste0('gtex/10_GWAS_coloc/gwas_loci/quant_trait/', trait_i))
}

for (i in 1:nrow(gwas_sig_loci_new)) {
  chr_i = as.numeric(gwas_sig_loci_new$CHR[i])
  start_i = gwas_sig_loci_new$loci_start[i]
  end_i = gwas_sig_loci_new$loci_end[i]
  
  loci_i = gwas_i[which(gwas_i$CHR == chr_i & gwas_i$POS >= start_i & gwas_i$POS <= end_i), ]

  fwrite(loci_i, paste0('gtex/10_GWAS_coloc/gwas_loci/quant_trait/', trait_i, '/chr', chr_i, '_', 
                        gwas_sig_loci_new$loci[i], '.txt.gz'), sep = '\t')
  print(i)
}


