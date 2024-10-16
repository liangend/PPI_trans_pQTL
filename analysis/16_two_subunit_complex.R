library(data.table)
library(openxlsx)
library(ggplot2)
library(coloc)
setwd('/project/xuanyao/jinghui/')
gene_meta = fread('/project/xuanyao/jinghui/gtex/00_ref/genecode.GRCh38.gene.meta.gtf')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)

prot_module = fread('pqtl/11_prot_complex/prot_module_corum_2022.txt')
corum = read.xlsx('pqtl/11_prot_complex/CORUM_2022_09_12.xlsx')
corum = corum[corum$Organism == 'Human', ]

## remove modules having exactly the same proteins
repeat_mod = names(which(table(corum$ComplexName) > 1))
repeat_mod = prot_module[which(prot_module$complex %in% repeat_mod), ]
repeat_mod = repeat_mod[order(repeat_mod$complex), ]
uniq_mod = unique(repeat_mod$mod)
uniq_mod_prot = list()
for (i in 1:length(uniq_mod)) {
  mod_i = repeat_mod[repeat_mod$mod == uniq_mod[i], ]
  uniq_mod_prot[[i]] = mod_i$prot
}

rm_mod = uniq_mod[c(2,3,5,7,9,11,14,16:18,23,25,27,29,32)]
prot_module = prot_module[!(prot_module$mod %in% rm_mod), ]
prot_module$gene_tss = gene_meta$start[match(prot_module$prot, gene_meta$gene_name)]

## direct comparison of p values
two_unit_mod = names(table(prot_module$mod)[(table(prot_module$mod) == 2)])
two_unit_mod = as.numeric(two_unit_mod)
window_size = 1000000
cis_snp = c()
cis_p = c()
cis_beta = c()
cis_p_other = c()
cis_beta_other = c()
for (i in two_unit_mod) {
  prot_module_i = prot_module[prot_module$mod == i, ]
  file_i = prot_module_i$file
  chr_i = paste0('chr', prot_module_i$chr, '_')
  tss_i = prot_module_i$gene_tss
  
  file1 = list.files(paste0('pqtl/UKB_PPP/', file_i[1]))
  file1 = file1[grep(chr_i[1], file1)]
  file2 = list.files(paste0('pqtl/UKB_PPP/', file_i[2]))
  file2 = file2[grep(chr_i[1], file2)]
  summ_stat1 = fread(paste0('pqtl/UKB_PPP/', file_i[1], '/', file1))
  summ_stat2 = fread(paste0('pqtl/UKB_PPP/', file_i[2], '/', file2))
  ## cis p value and beta of one subunit
  cis1 = summ_stat1[which(summ_stat1$GENPOS > tss_i[1] - window_size &
                            summ_stat1$GENPOS < tss_i[1] + window_size), ]
  cis_p1 = max(cis1$LOG10P)[1]
  cis_snp1 = cis1$ID[which.max(cis1$LOG10P)[1]]
  cis_beta1 = cis1$BETA[which.max(cis1$LOG10P)[1]]
  ## corresponding beta and p value of the other subunit
  cis_snp1_beta2 = summ_stat2$BETA[which(summ_stat2$ID == cis_snp1)]
  cis_snp1_p2 = summ_stat2$LOG10P[which(summ_stat2$ID == cis_snp1)]
  
  file2 = list.files(paste0('pqtl/UKB_PPP/', file_i[2]))
  file2 = file2[grep(chr_i[2], file2)]
  file1 = list.files(paste0('pqtl/UKB_PPP/', file_i[1]))
  file1 = file1[grep(chr_i[2], file1)]
  summ_stat2 = fread(paste0('pqtl/UKB_PPP/', file_i[2], '/', file2))
  summ_stat1 = fread(paste0('pqtl/UKB_PPP/', file_i[1], '/', file1))
  cis2 = summ_stat2[which(summ_stat2$GENPOS > tss_i[2] - window_size &
                            summ_stat2$GENPOS < tss_i[2] + window_size), ]
  cis_p2 = max(cis2$LOG10P)[1]
  cis_snp2 = cis2$ID[which.max(cis2$LOG10P)[1]]
  cis_beta2 = cis2$BETA[which.max(cis2$LOG10P)[1]]
  cis_snp2_beta1 = summ_stat1$BETA[which(summ_stat1$ID == cis_snp2)]
  cis_snp2_p1 = summ_stat1$LOG10P[which(summ_stat1$ID == cis_snp2)]
  
  cis_snp = c(cis_snp, cis_snp1, cis_snp2)
  cis_p = c(cis_p, cis_p1, cis_p2)
  cis_beta = c(cis_beta, cis_beta1, cis_beta2)
  cis_p_other = c(cis_p_other, cis_snp1_p2, cis_snp2_p1)
  cis_beta_other = c(cis_beta_other, cis_snp1_beta2, cis_snp2_beta1)
  print(i)
}

p_plot = data.frame(p = c(cis_p, cis_p_other), 
                    group = c(rep('extreme cis p value of one subunit', length(cis_p)), 
                              rep('p of the other subunit', length(cis_p_other))))
ggplot(p_plot, aes(x=group, y=p, fill=group)) +
  geom_boxplot() + ylim(c(0,50)) +
  labs(x = "", y = '-log10(p)', fill = '') +
  geom_hline(yintercept = -log10(0.05/866), linetype = "dashed") + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')





