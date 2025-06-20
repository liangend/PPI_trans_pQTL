library(coloc)
library(data.table)
args = commandArgs(trailingOnly = T)
trait_index = as.numeric(args[1])
setwd('/project/xuanyao/jinghui')
gwas_files = list.files('gtex/10_GWAS_coloc/GWAS_sum_stats/quant_trait/',
                        pattern = '.gz')
g = gwas_files[trait_index]

gwas = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/quant_trait/', g))
gwas = gwas[gwas$Freq_Tested_Allele > 0 & gwas$Freq_Tested_Allele < 1, ]
gwas = gwas[gwas$SE > 0, ]
gwas_trait = strsplit(g, '_', fixed = T)[[1]][1]
gwas$CHR = as.numeric(gwas$CHR)
#gwas$SNP = sapply(strsplit(gwas$SNP, ':', fixed = T), '[', 1) # correct some gwas files with different snp id format

if (!dir.exists(paste0('pqtl/08_gwas_coloc/pco_mcl_coloc/',
                       gwas_trait))) {
  dir.create(paste0('pqtl/08_gwas_coloc/pco_mcl_coloc/',
                    gwas_trait))
}

loci_file = list.files('pqtl/08_gwas_coloc/pco_mcl_loci/')
pp4_all = c()
for (i in 1:length(loci_file)) {
  loci = fread(paste0('pqtl/08_gwas_coloc/pco_mcl_loci/', loci_file[i]))
  # loci = na.omit(loci)
  loci$P = as.numeric(loci$P)
  chr_i = loci$CHROM[1]
  gwas_sub = gwas[gwas$CHR == chr_i, ]
  
  common_rsid = intersect(loci$bp_hg37, gwas_sub$POS)
  if (length(common_rsid) > 0) {
    loci = loci[match(common_rsid, loci$bp_hg37), ]
    gwas_sub = gwas_sub[match(common_rsid, gwas_sub$POS), ]
    
    gwas_frq = gwas_sub$Freq_Tested_Allele
    gwas_maf = pmin(gwas_frq, 1 - gwas_frq)
    
    pqtl_frq = loci$A1FREQ
    pqtl_maf = pmin(pqtl_frq, 1 - pqtl_frq)
    
    ### coloc has an error when all the p values are 0
    if (sum(loci$P) == 0) {
      loci$P = pmax(loci$P, 1e-300)
    }
    
    gwas_region = list('beta' = gwas_sub$BETA, 
                       'varbeta' = (gwas_sub$SE)^2,
                       'type' = 'quant', 
                       'snp' = gwas_sub$POS, 
                       'MAF' = gwas_maf,
                       'N' = gwas_sub$N)
    pqtl_region = list('pvalues' = loci$P,
                       'type' = 'quant',
                       'snp' = loci$bp_hg37,
                       'MAF' = pqtl_maf,
                       'N' = loci$N)
    
    coloc_res = coloc.abf(pqtl_region, gwas_region)
    pp4_res = coloc_res$summary["PP.H4.abf"]
      
    # coloc_table = coloc_res$results
    # coloc_save = data.frame(snp = coloc_table$snp, pp4 = coloc_table$SNP.PP.H4)
    # coloc_save$gwas_p = gwas_sub$P[match(coloc_table$snp, gwas_sub$SNP)]
    # coloc_save$qtl_p = loci$P[match(coloc_table$snp, loci$snp)]
  } else {
    pp4_res = 0
  }
  # if (pp4_res >= 0.8) {
  #   fwrite(coloc_save, paste0('/project/xuanyao/jinghui/pqtl/08_gwas_coloc/ukb_pco_coloc/', 
  #                             gwas_trait, '/coloc_', loci_file[i]), sep = '\t')
  # }
  pp4_all[i] = pp4_res
  print(i)
}

pp4_table = data.frame(file = loci_file, pp4 = pp4_all)
fwrite(pp4_table, paste0('pqtl/08_gwas_coloc/pco_mcl_coloc/', 
                      gwas_trait, '/pp4_all.txt'), sep = '\t')




