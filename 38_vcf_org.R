##### make vcf inputs for snpEff
setwd('/project/xuanyao/jinghui')
library(openxlsx)
library(data.table)
library(dplyr)
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
ukb_qtl = fread('pqtl/04_fdr/ukb/ukb_500kb_loci.txt')
ukb_trans = ukb_qtl[ukb_qtl$cis_trans == 'trans', ]

pco_all_vcf = data.frame(CHROM = as.numeric(sub('chr', '', pco_all$chr)), 
                         POS = pco_all$bp_hg37, 
                         ID = '.', REF = 'A', ALT = 'C', QUAL = '.', 
                         FILTER = '.', INFO = '.')
ukb_trans_vcf = data.frame(CHROM = ukb_trans$CHROM, 
                           POS = ukb_trans$bp_hg19, 
                           ID = '.', REF = 'A', ALT = 'C', QUAL = '.', 
                           FILTER = '.', INFO = '.')

#### DGN trans SNPs
## DGN snps from the paper
dgn_qtl_intra = read.xlsx('pqtl/00_ref/dgn_qtl.xlsx', sheet = 4)
dgn_qtl_inter = read.xlsx('pqtl/00_ref/dgn_qtl.xlsx', sheet = 6)
dgn_qtl_inter$LOG_PVAL = as.numeric(dgn_qtl_inter$LOG_PVAL)
dgn_qtl_inter_sub = dgn_qtl_inter %>%
  group_by(GENE_NAME) %>%
  slice(which.max(LOG_PVAL))
dgn_trans_gene = rbind(dgn_qtl_intra, dgn_qtl_inter_sub)

dgn_trans_gene_vcf = data.frame(CHROM = dgn_trans_gene$SNP_CHROM, 
                                POS = dgn_trans_gene$SNP_POS, 
                                ID = '.', REF = 'A', ALT = 'C', QUAL = '.', 
                                FILTER = '.', INFO = '.')

## DGN snps from my analysis
dgn_trans = fread('gtex/06_qtl_z/dgn_trans/fdr005.trans_qtl_pairs.txt.gz')
## make dgn trans loci
# significant loci prot pairs
loci_window = 250000
dgn_trans_loci = c()
sig_gene = unique(dgn_trans$phenotype_id)
loci_i = 1
for (i in sig_gene) {
  dgn_i = dgn_trans[dgn_trans$phenotype_id == i, ]
  while (nrow(dgn_i > 0)) {
    lead_sig = which.min(dgn_i$pval)
    lead_snp = dgn_i[lead_sig, ]
    lead_snp$loci = loci_i
    lead_snp$loci_start = lead_snp$var_bp - loci_window
    lead_snp$loci_end = lead_snp$var_bp + loci_window
    row_rm = which(dgn_i$var_chr == lead_snp$var_chr & dgn_i$var_bp >= lead_snp$loci_start &
                     dgn_i$var_bp <= lead_snp$loci_end)
    dgn_trans_loci = rbind(dgn_trans_loci, lead_snp)
    loci_i = loci_i + 1
    dgn_i = dgn_i[-row_rm, ]
  }
  print(i)
}
dgn_trans_loci = dgn_trans_loci[order(dgn_trans_loci$phenotype_id,
                                              dgn_trans_loci$var_bp,
                                              dgn_trans_loci$var_bp), ]
# merge overlapping windows
dgn_trans_new = dgn_trans_loci[1, ]
for (i in 2:nrow(dgn_trans_loci)) {
  chr_i = dgn_trans_loci$var_chr[i]
  start_i = dgn_trans_loci$loci_start[i]
  gene_i = dgn_trans_loci$phenotype_id[i]
  bp19_i = dgn_trans_loci$var_bp[i]
  chr_prev = dgn_trans_new$var_chr[nrow(dgn_trans_new)]
  end_prev = dgn_trans_new$loci_end[nrow(dgn_trans_new)]
  gene_prev = dgn_trans_new$phenotype_id[nrow(dgn_trans_new)]
  bp19_prev = dgn_trans_new$var_bp[nrow(dgn_trans_new)]
  # merge MHC regions
  if (gene_prev == gene_i & chr_i == 6 & chr_prev == 6 & bp19_i > 28477797 & bp19_i < 33448354 &
      bp19_prev > 28477797 & bp19_prev < 33448354) {
    # merge two overlapping loci
    add_pool_i = rbind(dgn_trans_new[nrow(dgn_trans_new), ], dgn_trans_loci[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$pval), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    dgn_trans_new[nrow(dgn_trans_new), ] = add_loci_i
    next
  }
  # keep non-overlapping windows
  if (gene_prev != gene_i | chr_i != chr_prev |
      (gene_prev == gene_i & chr_i == chr_prev & start_i > end_prev)) {
    # add the loci directly if two loci do not overlap
    dgn_trans_new = rbind(dgn_trans_new, dgn_trans_loci[i, ])
  } else {
    # merge two overlapping loci
    add_pool_i = rbind(dgn_trans_new[nrow(dgn_trans_new), ], dgn_trans_loci[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$pval), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    dgn_trans_new[nrow(dgn_trans_new), ] = add_loci_i
  }
  if (i %% 100 == 0) {
    print(i)
  }
}
dgn_trans_snp_vcf = data.frame(CHROM = as.numeric(sub('chr', '', dgn_trans_new$var_chr)), 
                                POS = dgn_trans_new$var_bp, 
                                ID = '.', REF = dgn_trans_new$A0, ALT = dgn_trans_new$A1, QUAL = '.', 
                                FILTER = '.', INFO = '.')

## DGN trans-pco snps
dgn_pco = read.xlsx('gtex/00_ref/Wang_transPCO_supp.xlsx', sheet = 3)
dgn_pco$var_chr = as.numeric(sapply(strsplit(dgn_pco$SNP, ':', fixed = T), '[', 1))
dgn_pco$var_bp = as.numeric(sapply(strsplit(dgn_pco$SNP, ':', fixed = T), '[', 2))
# significant loci prot pairs
loci_window = 250000
dgn_pco_loci = c()
sig_gene = unique(dgn_pco$gene_module)
loci_i = 1
for (i in sig_gene) {
  dgn_i = dgn_pco[dgn_pco$gene_module == i, ]
  while (nrow(dgn_i > 0)) {
    lead_sig = which.min(dgn_i$P)
    lead_snp = dgn_i[lead_sig, ]
    lead_snp$loci = loci_i
    lead_snp$loci_start = lead_snp$var_bp - loci_window
    lead_snp$loci_end = lead_snp$var_bp + loci_window
    row_rm = which(dgn_i$var_chr == lead_snp$var_chr & dgn_i$var_bp >= lead_snp$loci_start &
                     dgn_i$var_bp <= lead_snp$loci_end)
    dgn_pco_loci = rbind(dgn_pco_loci, lead_snp)
    loci_i = loci_i + 1
    dgn_i = dgn_i[-row_rm, ]
  }
  print(i)
}


dgn_pco_loci = dgn_pco_loci[order(dgn_pco_loci$gene_module,
                                      dgn_pco_loci$var_bp,
                                      dgn_pco_loci$var_bp), ]
# merge overlapping windows
dgn_pco_new = dgn_pco_loci[1, ]
for (i in 2:nrow(dgn_pco_loci)) {
  chr_i = dgn_pco_loci$var_chr[i]
  start_i = dgn_pco_loci$loci_start[i]
  gene_i = dgn_pco_loci$gene_module[i]
  bp19_i = dgn_pco_loci$var_bp[i]
  chr_prev = dgn_pco_new$var_chr[nrow(dgn_pco_new)]
  end_prev = dgn_pco_new$loci_end[nrow(dgn_pco_new)]
  gene_prev = dgn_pco_new$gene_module[nrow(dgn_pco_new)]
  bp19_prev = dgn_pco_new$var_bp[nrow(dgn_pco_new)]
  # merge MHC regions
  if (gene_prev == gene_i & chr_i == 6 & chr_prev == 6 & bp19_i > 28477797 & bp19_i < 33448354 &
      bp19_prev > 28477797 & bp19_prev < 33448354) {
    # merge two overlapping loci
    add_pool_i = rbind(dgn_pco_new[nrow(dgn_pco_new), ], dgn_pco_loci[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$P), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    dgn_pco_new[nrow(dgn_pco_new), ] = add_loci_i
    next
  }
  # keep non-overlapping windows
  if (gene_prev != gene_i | chr_i != chr_prev |
      (gene_prev == gene_i & chr_i == chr_prev & start_i > end_prev)) {
    # add the loci directly if two loci do not overlap
    dgn_pco_new = rbind(dgn_pco_new, dgn_pco_loci[i, ])
  } else {
    # merge two overlapping loci
    add_pool_i = rbind(dgn_pco_new[nrow(dgn_pco_new), ], dgn_pco_loci[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$P), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    dgn_pco_new[nrow(dgn_pco_new), ] = add_loci_i
  }
  if (i %% 100 == 0) {
    print(i)
  }
}

dgn_pco_vcf = data.frame(CHROM = dgn_pco_new$var_chr, 
                               POS = dgn_pco_new$var_bp, 
                               ID = '.', REF = 'A', ALT = 'C', QUAL = '.', 
                               FILTER = '.', INFO = '.')
colnames(pco_all_vcf)[1] = '#CHROM'
colnames(ukb_trans_vcf)[1] = '#CHROM'
colnames(dgn_trans_gene_vcf)[1] = '#CHROM'
colnames(dgn_trans_snp_vcf)[1] = '#CHROM'
colnames(dgn_pco_vcf)[1] = '#CHROM'

fwrite(pco_all_vcf, 'pqtl/17_snpEff/ukb_pco.vcf', sep = '\t')
fwrite(ukb_trans_vcf, 'pqtl/17_snpEff/ukb_trans.vcf', sep = '\t')
fwrite(dgn_trans_gene_vcf, 'pqtl/17_snpEff/dgn_trans_gene.vcf', sep = '\t')
fwrite(dgn_trans_snp_vcf, 'pqtl/17_snpEff/dgn_trans_snp.vcf', sep = '\t')
fwrite(dgn_pco_vcf, 'pqtl/17_snpEff/dgn_pco_snp.vcf', sep = '\t')

### all of ukb qtls (cis and trans)
ukb_qtl_vcf = data.frame(CHROM = ukb_qtl$CHROM, 
                         POS = ukb_qtl$bp_hg19, 
                         ID = '.', REF = ukb_qtl$ALLELE0, ALT = ukb_qtl$ALLELE1, 
                         QUAL = '.', FILTER = ukb_qtl$cis_trans, INFO = '.')
colnames(ukb_qtl_vcf)[1] = '#CHROM'
fwrite(ukb_qtl_vcf, 'pqtl/17_snpEff/ukb_all.vcf', sep = '\t')

## dgn cis
dgn_cis = fread('gtex/06_qtl_z/dgn_cis/all.cis_gene.txt.gz')
dgn_cis_vcf = data.frame(CHROM = sapply(strsplit(dgn_cis$variant_id, ':', fixed = T), '[', 1), 
                         POS = sapply(strsplit(dgn_cis$variant_id, ':', fixed = T), '[', 2), 
                         ID = '.', REF = dgn_cis$A0, ALT = dgn_cis$A1, 
                         QUAL = '.', FILTER = '.', INFO = '.')
colnames(dgn_cis_vcf)[1] = '#CHROM'
fwrite(dgn_cis_vcf, 'pqtl/17_snpEff/dgn_cis.vcf', sep = '\t')

## ukb all snp
ukb_all = fread('pqtl/00_ref/ukb_geno/all.txt')
set.seed(101)
ukb_sub = ukb_all[sample(1:nrow(ukb_all), 2300000), ]
ukb_sub = ukb_sub[ukb_sub$CHROM != 'CHROM', ]
ukb_sub$bp_hg37 = sapply(strsplit(ukb_sub$ID, ':', fixed = T), '[', 2)
ukb_sub = ukb_sub[, c(1:2,7,4:6)]
colnames(ukb_sub)[2] = 'bp_hg38'
fwrite(ukb_sub, 'pqtl/00_ref/ukb_geno/rand_sub.txt', sep = '\t')

ukb_rand_vcf = data.frame(CHROM = ukb_sub$CHROM, 
                          POS = ukb_sub$bp_hg37, 
                          ID = '.', REF = ukb_sub$ALLELE0, ALT = ukb_sub$ALLELE1, 
                          QUAL = '.', FILTER = '.', INFO = '.')
colnames(ukb_rand_vcf)[1] = '#CHROM'
fwrite(ukb_rand_vcf, 'pqtl/17_snpEff/ukb_rand.vcf', sep = '\t')


