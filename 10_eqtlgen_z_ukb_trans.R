library(data.table)
library(openxlsx)
setwd('/project/xuanyao/jinghui')
gene_meta = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)

ukb_pqtl = read.xlsx('pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                     startRow = 5, cols = c(1:11, 20), sheet = 10)
colnames(ukb_pqtl)[1] = 'var_id_hg19'
colnames(ukb_pqtl)[12] = 'cis_trans'
ukb_pqtl$gene_id = gene_meta$gene_id[match(ukb_pqtl$Assay.Target, gene_meta$gene_name)]

eqtl_gen_genes = fread('/project2/xuanyao/data/eQTLGen/eQTLGen.gene.txt', header = F)

ukb_trans = ukb_pqtl[ukb_pqtl$gene_id %in% eqtl_gen_genes$V1 & ukb_pqtl$cis_trans == 'trans', ]

eqtlgen_snp = fread('/project2/xuanyao/data/eQTLGen/eQTLGen.snp.txt', header = F)
ukb_trans = ukb_trans[ukb_trans$rsID %in% eqtlgen_snp$V1, ]

eqtl_gen_trans = fread('/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz')
pval = c()
bon_pval = c()
z = c()
A0 = c()
A1 = c()
fdr = c()
for (i in 1:nrow(ukb_trans)) {
  gene_i = ukb_trans$gene_id[i]
  snp_i = ukb_trans$rsID[i]
  row_i = eqtl_gen_trans[which(eqtl_gen_trans$Gene == gene_i & 
                                 eqtl_gen_trans$SNP == snp_i), ]
  if (nrow(row_i) > 0) {
    pval[i] = row_i$Pvalue[1]
    bon_pval[i] = row_i$BonferroniP[1]
    z[i] = row_i$Zscore[1]
    A0[i] = row_i$AssessedAllele[1]
    A1[i] = row_i$OtherAllele[1]
    fdr[i] = row_i$FDR[1]
  } else {
    pval[i] = 1
    bon_pval[i] = 1
    z[i] = 0
    A0[i] = NA
    A1[i] = NA
    fdr[i] = 1
  }
  print(i)
}

ukb_trans$eqtlgen_p = pval
ukb_trans$eqtlgen_bon_pval = bon_pval
ukb_trans$eqtlgen_z = z
ukb_trans$eqtlgen_A0 = A0
ukb_trans$eqtlgen_A1 = A1
ukb_trans$eqtlgen_fdr = fdr

fwrite(ukb_trans, 'pqtl/12_beta_across_two_data/eqtlgen_z_ukb_trans_all.txt', sep = '\t')



