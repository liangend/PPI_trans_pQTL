library(data.table)
setwd('/project/xuanyao/jinghui')
eqtlgen_cis = fread('/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
dgn_meta = fread('gtex/06_qtl_z/dgn_gene_meta.txt')
dgn_meta$GeneNameConv = sapply(strsplit(dgn_meta$GeneNameConv, '.', fixed = T), '[', 1)
eqtlgen_cis = eqtlgen_cis[eqtlgen_cis$Gene %in% dgn_meta$GeneNameConv, ]
uniq_gene = unique(eqtlgen_cis$Gene)

args = commandArgs(trailingOnly = T)
k = as.numeric(args[1])
i = uniq_gene[k]

eqtlgen_cis_sub = eqtlgen_cis[eqtlgen_cis$Gene == i, ]
eqtlgen_cis_sub = eqtlgen_cis_sub[order(eqtlgen_cis_sub$SNPChr), ]

gene_name_i = dgn_meta$gene_id[which(dgn_meta$GeneNameConv == i)]
chr_i = dgn_meta$chr[which(dgn_meta$gene_id == gene_name_i)]
set_i = dgn_meta$gene_set[which(dgn_meta$gene_id == gene_name_i)]

dgn_summ = fread(paste0('gtex/06_qtl_z/dgn_z/geneSet', set_i, 
                        '.', chr_i, '.trans_qtl_pairs.txt.gz'))
dgn_summ = dgn_summ[dgn_summ$phenotype_id == gene_name_i, ]
dgn_summ$bp = as.numeric(sapply(strsplit(dgn_summ$variant_id, ':', fixed = T), '[', 2))

eqtlgen_cis_sub = cbind(eqtlgen_cis_sub, 
                        dgn_summ[match(eqtlgen_cis_sub$SNPPos, dgn_summ$bp), 
                                       c('b', 'b_se', 'pval')])

fwrite(eqtlgen_cis_sub, paste0('pqtl/12_beta_across_two_data/dgn_beta_eqtlgen_cis/', i, '.txt'),
       sep = '\t', col.names = F)


