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

ukb_pqtl_sub = ukb_pqtl[ukb_pqtl$gene_id %in% eqtl_gen_genes$V1 & 
                          ukb_pqtl$cis_trans == 'cis', ]

eqtl_gen_cis = fread('/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
pval = c()
bon_pval = c()
z = c()
for (i in 1:nrow(ukb_pqtl_sub)) {
  gene_i = ukb_pqtl_sub$gene_id[i]
  chr_i = ukb_pqtl_sub$CHROM[i]
  bp_i = as.numeric(strsplit(ukb_pqtl_sub$var_id_hg19[i], ':', fixed = T)[[1]][2])
  row_i = eqtl_gen_cis[which(eqtl_gen_cis$Gene == gene_i & 
                               eqtl_gen_cis$SNPChr == chr_i &
                               eqtl_gen_cis$SNPPos == bp_i), ]
  if (nrow(row_i) > 0) {
    pval[i] = row_i$Pvalue
    bon_pval[i] = row_i$BonferroniP
    z[i] = row_i$Zscore
  } else {
    pval[i] = 1
    bon_pval[i] = 1
    z[i] = 0
  }
  print(i)
}

ukb_pqtl_sub$eqtlgen_p = pval
ukb_pqtl_sub$eqtlgen_bon_pval = bon_pval
ukb_pqtl_sub$eqtlgen_z = z

fwrite(ukb_pqtl_sub, 'pqtl/12_beta_across_two_data/eqtlgen_beta_ukb_pqtl_all.txt', sep = '\t')





