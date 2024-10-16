library(data.table)
setwd('/project/xuanyao/jinghui')
eqtlgen_trans = fread('/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz')
dgn_meta = fread('gtex/06_qtl_z/dgn_gene_meta.txt')
dgn_meta$GeneNameConv = sapply(strsplit(dgn_meta$GeneNameConv, '.', fixed = T), '[', 1)
eqtlgen_trans = eqtlgen_trans[eqtlgen_trans$Gene %in% dgn_meta$GeneNameConv, ]
uniq_gene = unique(eqtlgen_trans$Gene)

args = commandArgs(trailingOnly = T)
k = as.numeric(args[1])
i = uniq_gene[k]

eqtlgen_trans_sub = eqtlgen_trans[eqtlgen_trans$Gene == i, ]
eqtlgen_trans_sub = eqtlgen_trans_sub[order(eqtlgen_trans_sub$SNPChr), ]

gene_name_i = dgn_meta$gene_id[which(dgn_meta$GeneNameConv == i)]
set_i = dgn_meta$gene_set[which(dgn_meta$gene_id == gene_name_i)]

uniq_chr = unique(eqtlgen_trans_sub$SNPChr)
eqtlgen_trans_i = c()
for (j in uniq_chr) {
  dgn_summ = fread(paste0('gtex/06_qtl_z/dgn_z/geneSet', set_i, 
                          '.chr', j, '.trans_qtl_pairs.txt.gz'))
  dgn_summ = dgn_summ[dgn_summ$phenotype_id == gene_name_i, ]
  dgn_summ$bp = as.numeric(sapply(strsplit(dgn_summ$variant_id, ':', fixed = T), '[', 2))
  eqtlgen_trans_j = eqtlgen_trans_sub[eqtlgen_trans_sub$SNPChr == j, ]
  eqtlgen_trans_j = cbind(eqtlgen_trans_j, 
                          dgn_summ[match(eqtlgen_trans_j$SNPPos, dgn_summ$bp), 
                                   c('b', 'b_se', 'pval')])
  eqtlgen_trans_i = rbind(eqtlgen_trans_i, eqtlgen_trans_j)
}

fwrite(eqtlgen_trans_i, paste0('pqtl/12_beta_across_two_data/dgn_beta_eqtlgen_trans/', i, '.txt'),
       sep = '\t', col.names = F)


