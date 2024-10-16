library(data.table)
setwd('/project/xuanyao/jinghui')
dgn_cis_gene = fread('gtex/06_qtl_z/dgn_cis/all.cis_gene.txt.gz')
eqtlgen_gene = fread('/project2/xuanyao/data/eQTLGen/eQTLGen.gene.txt', header = F)

dgn_meta = fread('gtex/06_qtl_z/dgn_gene_meta.txt')
dgn_meta$GeneNameConv = sapply(strsplit(dgn_meta$GeneNameConv, '.', fixed = T), '[', 1)
dgn_cis_gene$gene_id = dgn_meta$GeneNameConv[match(dgn_cis_gene$phenotype_id, dgn_meta$gene_id)]

dgn_cis_gene = dgn_cis_gene[dgn_cis_gene$gene_id %in% eqtlgen_gene$V1, ]
uniq_gene = unique(dgn_cis_gene$phenotype_id)
cis_window = 1000000

args = commandArgs(trailingOnly = T)
k = as.numeric(args[1])
i = uniq_gene[k]

eqtlgen_cis = fread('/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')

gene_id_i = dgn_meta$GeneNameConv[which(dgn_meta$gene_id == i)]
chr_i = dgn_meta$chr[which(dgn_meta$gene_id == i)]
set_i = dgn_meta$gene_set[which(dgn_meta$gene_id == i)]
tss_i = dgn_meta$start[which(dgn_meta$gene_id == i)]
pval_thre = dgn_cis_gene$pval_nominal_threshold[dgn_cis_gene$phenotype_id == i]

dgn_summ = fread(paste0('gtex/06_qtl_z/dgn_z/geneSet', set_i, 
                        '.', chr_i, '.trans_qtl_pairs.txt.gz'))
dgn_summ = dgn_summ[dgn_summ$phenotype_id == i, ]
dgn_sig = dgn_summ[dgn_summ$pval < pval_thre, ]
dgn_sig$chr = sapply(strsplit(dgn_sig$variant_id, ':', fixed = T), '[', 1)
dgn_sig$bp = as.numeric(sapply(strsplit(dgn_sig$variant_id, ':', fixed = T), '[', 2))
dgn_cis = dgn_sig[dgn_sig$bp > tss_i - cis_window & 
                    dgn_sig$bp < tss_i + cis_window, ]

dgn_snp = fread(paste0('/project2/xuanyao/data/DGN/genotype/', chr_i, '_QCed.bim'))
dgn_cis = cbind(dgn_cis, dgn_snp[match(dgn_cis$variant_id, dgn_snp$V2), 5:6])
colnames(dgn_cis)[9:10] = c('A0', 'A1')

eqtlgen_sub = eqtlgen_cis[eqtlgen_cis$Gene == gene_id_i, ]
dgn_cis = cbind(dgn_cis, 
                eqtlgen_sub[match(dgn_cis$bp, eqtlgen_sub$SNPPos), 
                            c('AssessedAllele', 'OtherAllele', 'Zscore', 'Pvalue', 'FDR', 'BonferroniP')])

fwrite(dgn_cis, paste0('pqtl/12_beta_across_two_data/eqtlgen_z_dgn_cis/', i, '.txt'),
       sep = '\t', col.names = F)



