library(coloc)
library(data.table)
args = commandArgs(trailingOnly = T)
chr_i = as.numeric(args[1])
setwd('/project/xuanyao/jinghui')

ukb_loci = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
ukb_loci = ukb_loci[ukb_loci$cis_trans == 'trans', ]
ukb_loci = ukb_loci[ukb_loci$CHROM == chr_i, ]

eqtl_summ = fread(paste0('gtex/00_ref/eQTLGen/cis_chr', chr_i, '.txt.gz'))

# eqtl_summ = fread(paste0('pqtl/16_pco_cis_coloc/INTERVAL/eQTL_summ_stats/INTERVAL_eQTL_nominal_chr',
#                          chr_i, '.tsv'))

print(paste0('total loci: ', nrow(ukb_loci)))
for (i in 1:nrow(ukb_loci)) {
  loci_i = fread(paste0('pqtl/08_gwas_coloc/ukb_qtl_loci/', 
                        ukb_loci$file[i], '_', 'chr', chr_i, '_trans_',
                        ukb_loci$loci[i], '.txt.gz'))
  loci_i$hg37_bp = as.numeric(sapply(strsplit(loci_i$ID, ':', fixed = T), '[', 2))
  loci_i$P = 10^-(loci_i$LOG10P)
  
  eqtl_sub = eqtl_summ[which(eqtl_summ$SNPPos >= min(loci_i$hg37_bp) &
                               eqtl_summ$SNPPos <= max(loci_i$hg37_bp)), ]
  # eqtl_sub = eqtl_summ[which(eqtl_summ$pos_b37 >= min(loci_i$hg37_bp) &
  #                              eqtl_summ$pos_b37 <= max(loci_i$hg37_bp)), ]
  if (nrow(eqtl_sub) < 1) {
    next
  }
  
  uniq_gene = unique(eqtl_sub$GeneSymbol)
  pp4_all = c()
  minp = c()
  for (j in uniq_gene) {
    eqtl_j = eqtl_sub[eqtl_sub$GeneSymbol == j, ]
    common_rsid = intersect(loci_i$hg37_bp, eqtl_j$SNPPos)
    if (length(common_rsid) > 0) {
      loci_j = loci_i[match(common_rsid, loci_i$hg37_bp), ]
      pqtl_frq = loci_j$A1FREQ
      pqtl_maf = pmin(pqtl_frq, 1 - pqtl_frq)
      
      eqtl_j = eqtl_j[match(common_rsid, eqtl_j$SNPPos), ]
      # eqtl_frq = eqtl_j$af
      # eqtl_maf = pmin(eqtl_frq, 1 - eqtl_frq)
      
      ### coloc has an error when all the p values are 0
      loci_j$P[loci_j$P == 0] = 1e-310
      eqtl_j$Pvalue[eqtl_j$Pvalue == 0] = 1e-310
      
      eqtl_region = list('pvalues' = eqtl_j$Pvalue,
                         'type' = 'quant', 
                         'snp' = eqtl_j$SNPPos, 
                         'MAF' = pqtl_maf,
                         'N' = eqtl_j$NrSamples)
                         # 'N' = 4372)
      pqtl_region = list('pvalues' = loci_j$P,
                         'type' = 'quant',
                         'snp' = loci_j$hg37_bp,
                         'MAF' = pqtl_maf,
                         'N' = loci_j$N)
      
      coloc_res = coloc.abf(pqtl_region, eqtl_region)
      pp4_res = coloc_res$summary["PP.H4.abf"]
      minp_j = min(eqtl_j$Pvalue)
    } else {
      pp4_res = 0
      minp_j = NA
    }
    pp4_all = c(pp4_all, pp4_res)
    minp = c(minp, minp_j)
  }
  
  pp4_table = data.frame(eqtl_gene = uniq_gene, ukb_loci = ukb_loci$loci[i], 
                         pp4 = pp4_all, eqtl_minp = minp)
  fwrite(pp4_table, paste0('pqtl/14_cis_trans_coloc/ukb_trans_eqtlgen_eqtl_coloc/', 
                           ukb_loci$loci[i], '.txt'), 
         sep = '\t', col.names = F)
  print(i)
}



