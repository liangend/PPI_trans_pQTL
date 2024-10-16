library(coloc)
library(data.table)

setwd('/project/xuanyao/jinghui')
args = commandArgs(trailingOnly = T)
trait_index = as.numeric(args[1])
chr_index = as.numeric(args[2])

gwas_loci_files = list.files('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/gwas_loci_500kb/')
gwas_loci = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/cc_trait/gwas_loci_500kb/', 
                         gwas_loci_files[trait_index]))
gwas_loci$CHR = as.numeric(gwas_loci$CHR)
gwas_loci = gwas_loci[gwas_loci$CHR == chr_index, ]

gwas_trait = sub('.txt', '', gwas_loci_files[trait_index])
print(gwas_trait)

if (!dir.exists(paste0('pqtl/08_gwas_coloc/gwas_coloc_mcl/', gwas_trait))) {
  dir.create(paste0('pqtl/08_gwas_coloc/gwas_coloc_mcl/', gwas_trait))
}
mcl_mod = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')
n_mod = max(mcl_mod$mod)

pqtl_maf_file = fread(paste0('/project2/xuanyao/jinghui/UKB_PPP_combined/A1BG_P04217_OID30771_v1_Inflammation_II/combined_chr', 
                             chr_index, '_A1BG:P04217:OID30771:v1:Inflammation_II.gz'))
pqtl_sample_size = pqtl_maf_file$N[1]

pp4_all = c()
mod_all = c()
gwas_loci_all = c()
minp_all = c()
for (i in 1:n_mod) {
  pco_loci = readRDS(paste0('pqtl/03_p/ukb_mcl/p.mod', i, '.chr', chr_index, '.rds'))
  if (length(pco_loci) < 1) {
    next
  }
  pco_loci = data.frame(P = unname(pco_loci), ID = names(pco_loci))
  pco_loci$A1FREQ = pqtl_maf_file$A1FREQ[match(pco_loci$ID, pqtl_maf_file$ID)]
  pco_loci$hg38_bp = pqtl_maf_file$GENPOS[match(pco_loci$ID, pqtl_maf_file$ID)]
  pco_loci$hg37_bp = as.numeric(sapply(strsplit(pco_loci$ID, ':', fixed = T), '[', 2))
  pco_loci = na.omit(pco_loci)
  if (length(pco_loci) < 1) {
    next
  }
  
  for (j in 1:nrow(gwas_loci)) {
    gwas_loci_j = fread(paste0('gtex/10_GWAS_coloc/gwas_loci/cc_trait/', gwas_trait, '/chr', 
                               chr_index, '_', gwas_loci$loci[j], '.txt.gz'))
    gwas_loci_j = na.omit(gwas_loci_j)
    gwas_loci_j$Freq_Tested_Allele = as.numeric(gwas_loci_j$Freq_Tested_Allele)
    gwas_loci_j = gwas_loci_j[gwas_loci_j$Freq_Tested_Allele > 0 & 
                                gwas_loci_j$Freq_Tested_Allele < 1, ]
    if (gwas_trait %in% c('CD', 'IBD', 'MS', 'stroke', 'T1D', 'T2D')) {  ## gwas with hg38 bp
      pco_loci$POS = pco_loci$hg38_bp
    } else {
      pco_loci$POS = pco_loci$hg37_bp
    }
    common_rsid = intersect(gwas_loci_j$POS, pco_loci$POS)
    if (length(common_rsid) > 0) {
      pco_loci_j = pco_loci[match(common_rsid, pco_loci$POS), ]
      gwas_loci_j = gwas_loci_j[match(common_rsid, gwas_loci_j$POS), ]
      
      gwas_frq = gwas_loci_j$Freq_Tested_Allele
      gwas_maf = pmin(gwas_frq, 1 - gwas_frq)
      
      pqtl_frq = pco_loci_j$A1FREQ
      pqtl_maf = pmin(pqtl_frq, 1 - pqtl_frq)
      
      gwas_loci_j$P[gwas_loci_j$P == 0] = 1e-310
      pco_loci_j$P[pco_loci_j$P == 0] = 1e-310
      
      gwas_region = list('type' = 'cc', 
                         'snp' = gwas_loci_j$POS, 
                         'MAF' = gwas_maf,
                         'N' = gwas_loci_j$N,
                         "pvalues" = gwas_loci_j$P,
                         "s" = gwas_loci_j$N_case / gwas_loci_j$N)
      
      pqtl_region = list('pvalues' = pco_loci_j$P,
                         'N' = pqtl_sample_size,
                         'MAF' = pqtl_maf,
                         'type' = 'quant',
                         'snp' = pco_loci_j$POS)
      
      coloc_res = coloc.abf(pqtl_region, gwas_region)
      pp4_res = coloc_res$summary["PP.H4.abf"]
      minp_j = min(pco_loci_j$P)
    } else {
      pp4_res = 0
      minp_j = NA
    }
    
    pp4_all = c(pp4_all, pp4_res)
    mod_all = c(mod_all, i)
    gwas_loci_all = c(gwas_loci_all, gwas_loci$loci[j])
    minp_all = c(minp_all, minp_j)
  }
}
pp4_table = data.frame(mcl_mod = mod_all,
                       gwas_loci = gwas_loci_all, 
                       pp4 = pp4_all,
                       min_pco_p = minp_all)
fwrite(pp4_table, paste0('pqtl/08_gwas_coloc/gwas_coloc_mcl/', 
                         gwas_trait, '/pp4_chr', chr_index, '.txt'), sep = '\t')




