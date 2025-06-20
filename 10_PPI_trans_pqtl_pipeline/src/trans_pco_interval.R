### this script calculate the p value using the trans PCO method

prot.meta = snakemake@input[['prot_meta']]
file.zmat = snakemake@input[['file_zmat']]
#file.zmat = paste0('/project/xuanyao/jinghui/pqtl/01_zmat/interval/mod1.chr1.txt.gz')
file.sigma = snakemake@input[['file_sigma']]
#file.sigma = '/project/xuanyao/jinghui/pqtl/02_sigma/interval/mod1.rds'

file.p = snakemake@output[['file_p']]
#file.p = '/project/xuanyao/jinghui/pqtl/03_p/interval/p.mod1.chr1.rds'

params1 = snakemake@params[['dir_script']]
#params1 = '/project/xuanyao/jinghui/pqtl/99_scripts/src/'
chr = as.numeric(snakemake@params[['chr']])
#chr = 1
module = as.numeric(snakemake@params[['module']])
#module = 1

library(data.table)
library(ACAT)
## read files
# protein meta files
prot_meta = as.data.frame(fread(prot.meta))
prot_meta = na.omit(prot_meta)

# read z-scores of snps on proteins in the module
z.mat = fread(file.zmat)
z.mat = as.matrix(z.mat, rownames = TRUE)

# read pco gadgets
source(paste0(params1, "/pco_acat.R"))
source(paste0(params1, "/liu.R"))
source(paste0(params1, "/liumod.R"))
source(paste0(params1, "/davies.R"))
dyn.load(paste0(params1, "/qfc.so"))

## extract proteins in the module
gene_w_pos = data.frame(gene = colnames(z.mat))
gene_w_pos$chr = prot_meta$chr[match(gene_w_pos$gene, prot_meta$prot)]
# remove proteins/protein complex on the same chromosome as the snps
gene_trans = gene_w_pos$gene[gene_w_pos$chr != chr]

## run PCO to calculate p-values
if(length(gene_trans) > 1){
  Sigma = readRDS(file.sigma)
  Sigma = Sigma[gene_trans, gene_trans]
  z.mat_trans = z.mat[, gene_trans]; rm(z.mat)
  p.all = PCO_acat(Z.mat = z.mat_trans, Sigma = Sigma)
}else{
  p.all <- NULL
}

# save results -----
saveRDS(p.all, file = file.p)


