library(data.table)
setwd('/project/xuanyao/jinghui')
pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
ppi_list = readRDS('pqtl/11_prot_complex/ppi_list.rds')
univar_prot = strsplit(pco_all$univar_trans_prot, ',', fixed = T)
# background null enrichment in PPI
n_ite = 1000
n_ppi = c()
for (i in 1:n_ite) {
  sample_i = sample(pco_all$nearest_gene)
  prot_pair_i = sapply(1:length(univar_prot), function(x){paste0(univar_prot[[x]], ',', 
                                                                 sample_i[x])})
  ppi_i = sapply(prot_pair_i, function(x){sum(x %in% ppi_list)})
  n_ppi[i] = sum(ppi_i > 0)
  print(i)
}
fwrite(as.data.frame(n_ppi), 'pqtl/11_prot_complex/backgroud_n_ppi.txt')
