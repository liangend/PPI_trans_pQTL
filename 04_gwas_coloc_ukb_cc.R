library(coloc)
library(data.table)
setwd('/project/xuanyao/jinghui')

args = commandArgs(trailingOnly = T)
trait_index = as.numeric(args[1])

gwas_files = list.files('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/',
                        pattern = '.gz')
g = gwas_files[trait_index]

gwas = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/', g))
gwas = gwas[gwas$Freq_Tested_Allele > 0 & gwas$Freq_Tested_Allele < 1, ]
if (!is.null(gwas$SE)) {
  gwas = gwas[gwas$SE > 0, ]
}

gwas_trait = strsplit(g, '_', fixed = T)[[1]][1]
print(gwas_trait)
gwas$CHR = as.numeric(gwas$CHR)
gwas$P = as.numeric(gwas$P)

if (!dir.exists(paste0('pqtl/08_gwas_coloc/ukb_coloc/',
                       gwas_trait))) {
  dir.create(paste0('pqtl/08_gwas_coloc/ukb_coloc/',
                    gwas_trait))
}

loci_file = list.files('pqtl/08_gwas_coloc/ukb_qtl_loci/')
pp4_all = c()
for (i in 1:length(loci_file)) {
  loci = fread(paste0('pqtl/08_gwas_coloc/ukb_qtl_loci/', loci_file[i]))
  # loci = na.omit(loci)
  loci$BETA = as.numeric(loci$BETA)
  loci$SE = as.numeric(loci$SE)
  loci$hg_37 = as.numeric(sapply(strsplit(loci$ID, ':', fixed = T), '[', 2))
  chr_i = loci$CHROM[1]
  
  gwas_sub = gwas[gwas$CHR == chr_i, ]
  gwas_sub = na.omit(gwas_sub)
  if (gwas_trait %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {  ## gwas with hg38 bp
    common_rsid = intersect(gwas_sub$POS, loci$GENPOS)
    if (length(common_rsid) > 0) {
      loci = loci[match(common_rsid, loci$GENPOS), ]
      gwas_sub = gwas_sub[match(common_rsid, gwas_sub$POS), ]
      
      gwas_frq = gwas_sub$Freq_Tested_Allele
      gwas_maf = pmin(gwas_frq, 1 - gwas_frq)
      
      pqtl_frq = loci$A1FREQ
      pqtl_maf = pmin(pqtl_frq, 1 - pqtl_frq)
      
      gwas_sub$P[gwas_sub$P == 0] = 1e-310
      loci$P = 10^-(loci$LOG10P)
      loci$P[loci$P == 0] = 1e-310
      
      gwas_region = list('type' = 'cc', 
                         'snp' = gwas_sub$POS, 
                         'MAF' = gwas_maf,
                         'N' = gwas_sub$N,
                         "pvalues" = gwas_sub$P,
                         "s" = gwas_sub$N_case / gwas_sub$N)
      
      pqtl_region = list('pvalues' = loci$P,
                         'N' = loci$N,
                         'MAF' = pqtl_maf,
                         'type' = 'quant',
                         'snp' = loci$GENPOS)
      
      coloc_res = coloc.abf(pqtl_region, gwas_region)
      pp4_res = coloc_res$summary["PP.H4.abf"]
      
    } else {
      pp4_res = 0
    }
  } else {
    common_rsid = intersect(gwas_sub$POS, loci$hg_37)
    if (length(common_rsid) > 0) {
      loci = loci[match(common_rsid, loci$hg_37), ]
      gwas_sub = gwas_sub[match(common_rsid, gwas_sub$POS), ]
      
      gwas_frq = gwas_sub$Freq_Tested_Allele
      gwas_maf = pmin(gwas_frq, 1 - gwas_frq)
      
      pqtl_frq = loci$A1FREQ
      pqtl_maf = pmin(pqtl_frq, 1 - pqtl_frq)
      
      gwas_sub$P[gwas_sub$P == 0] = 1e-310
      loci$P = 10^-(loci$LOG10P)
      loci$P[loci$P == 0] = 1e-310
      
      gwas_region = list('type' = 'cc', 
                         'snp' = gwas_sub$POS, 
                         'MAF' = gwas_maf,
                         'N' = gwas_sub$N,
                         "pvalues" = gwas_sub$P,
                         "s" = gwas_sub$N_case / gwas_sub$N)
      
      pqtl_region = list('pvalues' = loci$P,
                         'N' = loci$N,
                         'MAF' = pqtl_maf,
                         'type' = 'quant',
                         'snp' = loci$hg_37)
      
      coloc_res = coloc.abf(pqtl_region, gwas_region)
      pp4_res = coloc_res$summary["PP.H4.abf"]
      
    } else {
      pp4_res = 0
    }
  }
  pp4_all[i] = pp4_res
  print(i)
}

pp4_table = data.frame(file = loci_file, pp4 = pp4_all)
fwrite(pp4_table, paste0('pqtl/08_gwas_coloc/ukb_coloc/', 
                         gwas_trait, '/pp4_all.txt'), sep = '\t')




