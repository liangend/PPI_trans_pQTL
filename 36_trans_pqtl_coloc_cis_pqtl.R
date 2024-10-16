library(data.table)
library(coloc)
setwd('/project/xuanyao/jinghui')

args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])

ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_qtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
ukb_trans = ukb_qtl[ukb_qtl$cis_trans == 'trans', ]

ukb_trans = ukb_trans[ukb_trans$CHROM == i, ]
ukb_prot = ukb_prot[ukb_prot$chr == i, ]
chr_i = paste0('chr', i, '_')
cis_size = 1000000
for (j in 1:nrow(ukb_trans)) {
  bp_j = ukb_trans$bp_hg38[j]
  cis_prot = ukb_prot[abs(ukb_prot$gene_start - bp_j) < cis_size, ]
  if (nrow(cis_prot) < 1) {
    next
  } else {
    trans_file = list.files(paste0('pqtl/UKB_PPP_combined/', ukb_trans$file[j]))
    trans_file = trans_file[grep(chr_i, trans_file)]
    trans_loci = fread(paste0('pqtl/UKB_PPP_combined/', ukb_trans$file[j], '/', trans_file), 
                       select = c(1:3,6,8,10:11,13))
    trans_loci = trans_loci[trans_loci$GENPOS >= ukb_trans$loci_start[j] &
                              trans_loci$GENPOS <= ukb_trans$loci_end[j], ]
    pp4_all = c()
    lead_snp_all = c()
    trans_beta = c()
    trans_se = c()
    cis_beta = c()
    cis_se = c()
    for (k in 1:nrow(cis_prot)) {
      cis_file = list.files(paste0('pqtl/UKB_PPP_combined/', cis_prot$file[k]))
      cis_file = cis_file[grep(chr_i, cis_file)]
      cis_loci = fread(paste0('pqtl/UKB_PPP_combined/', cis_prot$file[k], '/', cis_file), 
                         select = c(1:3, 8, 10:11, 13))
      cis_loci = cis_loci[cis_loci$GENPOS >= ukb_trans$loci_start[j] &
                            cis_loci$GENPOS <= ukb_trans$loci_end[j], ]
      common_snp = intersect(trans_loci$ID, cis_loci$ID)
      
      trans_loci_k = trans_loci[match(common_snp, trans_loci$ID), ]
      pqtl_frq = trans_loci_k$A1FREQ
      pqtl_maf = pmin(pqtl_frq, 1 - pqtl_frq)
      
      cis_loci = cis_loci[match(common_snp, cis_loci$ID), ]
      
      trans_region = list('beta' = trans_loci_k$BETA,
                         'varbeta' = (trans_loci_k$SE)^2, 
                         'type' = 'quant', 
                         'snp' = trans_loci_k$ID, 
                         'MAF' = pqtl_maf,
                         'N' = trans_loci_k$N)
      cis_region = list('beta' = cis_loci$BETA,
                        'varbeta' = (cis_loci$SE)^2, 
                        'type' = 'quant', 
                        'snp' = cis_loci$ID, 
                        'MAF' = pqtl_maf,
                        'N' = cis_loci$N)
      
      coloc_res = coloc.abf(trans_region, cis_region)
      coloc_tab = coloc_res$results
      
      pp4_res = coloc_res$summary["PP.H4.abf"]
      lead_snp = coloc_tab$snp[which.max(coloc_tab$SNP.PP.H4)]
      lead_trans_beta = trans_loci$BETA[which(trans_loci$ID == lead_snp)]
      lead_trans_se = trans_loci$SE[which(trans_loci$ID == lead_snp)]
      lead_cis_beta = cis_loci$BETA[which(cis_loci$ID == lead_snp)]
      lead_cis_se = cis_loci$SE[which(cis_loci$ID == lead_snp)]
      
      pp4_all = c(pp4_all, pp4_res)
      lead_snp_all = c(lead_snp_all, lead_snp)
      trans_beta = c(trans_beta, lead_trans_beta)
      trans_se = c(trans_se, lead_trans_se)
      cis_beta = c(cis_beta, lead_cis_beta)
      cis_se = c(cis_se, lead_cis_se)
    }
    trans_j = ukb_trans[rep(j, nrow(cis_prot)), ]
    trans_j$cis_prot = cis_prot$file
    trans_j$pp4 = pp4_all
    trans_j$lead_snp = lead_snp_all
    trans_j$trans_beta = trans_beta
    trans_j$trans_se = trans_se
    trans_j$cis_beta = cis_beta
    trans_j$cis_se = cis_se
    fwrite(trans_j, paste0('pqtl/14_cis_trans_coloc/each_trans_coloc/', ukb_trans$loci[j], '.txt'),
           sep = '\t')
  }
}



