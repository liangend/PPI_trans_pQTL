### this script generates the z matrix based on summary stats, which is the input for trans PCO

prot.meta = snakemake@input[['prot_meta']]
#prot.meta = '/project/xuanyao/jinghui/pqtl/00_ref/prot_meta.txt'

file.zmat = snakemake@output[['file_zmat']]
#file.zmat = paste0('/project/xuanyao/jinghui/pqtl/01_zmat/mod1.chr', 1:22, '.txt.gz')
file.sigma = snakemake@output[['file_sigma']]
#file.sigma = '/project/xuanyao/jinghui/pqtl/02_sigma/mod1.rds'

module = as.numeric(snakemake@params[['module']])
#module = 1
dir_sum_stats = snakemake@params[['dir_sum_stats']]
#dir_sum_stats = '/project/xuanyao/jinghui/pqtl/AGES_sum_stats'
indep.snp = snakemake@params[['indep_snp']]
#indep.snp = '/project/xuanyao/jinghui/pqtl/00_ref/dgn_geno/indep_snp.txt'
p_thre = snakemake@params[['p_thre']]
#p_thre = 0.0001

## read files
library(data.table)
prot_meta = as.data.frame(fread(prot.meta))
meta_mod_i = prot_meta[which(prot_meta$module == module), ]

# remove files not downloaded
all_file = list.files(dir_sum_stats)
meta_mod_i = meta_mod_i[meta_mod_i$file %in% all_file, ]

# read the whole z matrix for module i
mod_i_file = meta_mod_i$file
mod_i_file_path = paste0(dir_sum_stats, '/', mod_i_file, '/harmonised/')
z_mat = c()
for (i in mod_i_file_path) {
  sum_stats_i = list.files(i)
  prot_i = fread(paste0(i, sum_stats_i), select = c('chromosome', 'base_pair_location', 'beta', 'standard_error'))
  prot_i$beta = as.numeric(prot_i$beta)
  prot_i$standard_error = as.numeric(prot_i$standard_error)
  z_i = prot_i$beta / prot_i$standard_error
  z_mat = cbind(z_mat, z_i)
  gc()
}
colnames(z_mat) = meta_mod_i$gene
prot_snp = paste0(prot_i$chromosome, ":", prot_i$base_pair_location)
rownames(z_mat) = prot_snp
z_mat = as.data.frame(z_mat)

## save zmat by chromosome
for (i in 1:22) {
  chr_i = which(prot_i$chromosome == i)
  z_mat_i = z_mat[chr_i, ]
  fwrite(z_mat_i, file.zmat[i], sep = '\t', row.names = T)
}

## calculate sigma based on independent null snps
# independent SNPs based on --indep-pairwise 50 5 0.2 in plink
indep_snp = fread(indep.snp, header = F)
indep_prot_snp = which(prot_snp %in% indep_snp$V1)

# null snps are not significant associated with any proteins in the module
indep_z_mat = z_mat[indep_prot_snp, ]
rm(z_mat)
n_sig = apply(indep_z_mat, 1, function(x){sum(abs(x) > qnorm(p_thre/2, lower.tail = F))})
null_z_mat = indep_z_mat[which(n_sig == 0), ]
snp_prot_ratio = nrow(null_z_mat) / ncol(null_z_mat)
fwrite(as.data.frame(snp_prot_ratio), 
       paste0('/project/xuanyao/jinghui/pqtl/02_sigma/ratio_mod', module,'.txt'))

# write a sigma for each module
sigma_mat = cor(null_z_mat)
saveRDS(sigma_mat, file.sigma)



