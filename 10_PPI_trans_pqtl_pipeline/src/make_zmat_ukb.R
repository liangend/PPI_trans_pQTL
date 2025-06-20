### this script generates the z matrix based on summary stats, which is the input for trans PCO

prot.meta = snakemake@input[['prot_meta']]
#prot.meta = '/project/xuanyao/jinghui/pqtl/11_prot_complex/tf_mod.txt'

file.zmat = snakemake@output[['file_zmat']]
#file.zmat = paste0('/project/xuanyao/jinghui/pqtl/01_zmat/ukb/mod1.chr', 1:22, '.txt.gz')
file.sigma = snakemake@output[['file_sigma']]
#file.sigma = '/project/xuanyao/jinghui/pqtl/02_sigma/ukb/mod1.rds'

module = as.numeric(snakemake@params[['module']])
#module = 1
dir_sum_stats = snakemake@params[['dir_sum_stats']]
#dir_sum_stats = '/project/xuanyao/jinghui/pqtl/UKB_PPP/'
indep.snp = snakemake@params[['indep_snp']]
#indep.snp = '/project/xuanyao/jinghui/pqtl/00_ref/dgn_geno/indep_snp.txt'
p_thre = snakemake@params[['p_thre']]
#p_thre = 0.0001

## read files
library(data.table)
prot_meta = as.data.frame(fread(prot.meta))
meta_mod_i = prot_meta[which(prot_meta$mod == module), ]

# read the whole z matrix for module i by chromosome
mod_i_file = meta_mod_i$file
mod_i_file_path = paste0(dir_sum_stats, '/', mod_i_file)

z_mat = c()
for (i in 1:22) {
  chr_i = paste0('chr', i, '_')
  z_mat_i = c()
  for (j in mod_i_file_path) {
    sum_stats_j = list.files(j)
    sum_stats_j = sum_stats_j[grep(chr_i, sum_stats_j)]
    prot_j = fread(paste0(j, '/', sum_stats_j), select = c('ID', 'BETA', 'SE'))
    prot_j$BETA = as.numeric(prot_j$BETA)
    prot_j$SE = as.numeric(prot_j$SE)
    z_j = prot_j$BETA / prot_j$SE
    if (length(z_mat_i) == 0) {
      z_mat_i = data.frame(snp = prot_j$ID, z_j)
    } else {
      #z_mat_j = data.frame(snp_j, z_j)
      z_mat_i = data.frame(z_mat_i, z_j[match(z_mat_i$snp, prot_j$ID)])
    }
  }
  colnames(z_mat_i) = c('snp', meta_mod_i$prot)
  rownames(z_mat_i) = z_mat_i$snp
  
  z_mat_i = z_mat_i[, -1]
  z_mat_i = na.omit(z_mat_i)
  z_mat = rbind(z_mat, z_mat_i)
  fwrite(z_mat_i, file.zmat[i], sep = '\t', row.names = T)
  gc()
}

## calculate sigma based on independent null snps
# independent SNPs based on --indep-pairwise 50 5 0.2 in plink
indep_snp = fread(indep.snp, header = F)
prot_snp_chr = sapply(strsplit(rownames(z_mat), ':', fixed = T), '[', 1)
prot_snp_bp = sapply(strsplit(rownames(z_mat), ':', fixed = T), '[', 2)
prot_snp = paste0(prot_snp_chr, ':', prot_snp_bp)
indep_prot_snp = which(prot_snp %in% indep_snp$V1)

# null snps are not significant associated with any proteins in the module
indep_z_mat = z_mat[indep_prot_snp, ]
rm(z_mat)
n_sig = apply(indep_z_mat, 1, function(x){sum(abs(x) > qnorm(p_thre/2, lower.tail = F))})
null_z_mat = indep_z_mat[which(n_sig == 0), ]
snp_prot_ratio = nrow(null_z_mat) / ncol(null_z_mat)
fwrite(as.data.frame(snp_prot_ratio), 
       paste0('/project/xuanyao/jinghui/pqtl/02_sigma/ukb/ratio_mod', module,'.txt'))

# write a sigma for each module
sigma_mat = cor(null_z_mat)
saveRDS(sigma_mat, file.sigma)



