library(data.table)
library(ggplot2)
library(gridExtra)
## filter out cis signals of univariate analysis to get trans only
prot_meta = fread('/project/xuanyao/jinghui/pqtl/00_ref/prot_meta.txt')
univar_sig = fread('/project/xuanyao/jinghui/pqtl/04_fdr/sig_univar_all.txt')
univar_sig$snp = paste0(univar_sig$V1, ':', univar_sig$V2)
univar_trans = univar_sig[univar_sig$V1 != 'X', ]
univar_trans = univar_trans[univar_trans$V4 %in% prot_meta$file, ]
univar_trans$prot_chr = prot_meta$chr[match(univar_trans$V4, prot_meta$file)]
univar_trans$V1 = paste0('chr', univar_trans$V1)
univar_trans = univar_trans[univar_trans$V1 != univar_trans$prot_chr, ]
univar_trans = univar_trans[-which(univar_trans$V1 =='chr5' & grepl('chr5', univar_trans$prot_chr)), ]
univar_trans = univar_trans[-which(univar_trans$V1 =='chr1' & grepl('chr1;', univar_trans$prot_chr)), ]
univar_trans = univar_trans[-which(univar_trans$V1 =='chr19' & grepl('chr19', univar_trans$prot_chr)), ]
univar_trans$V3 = as.numeric(univar_trans$V3)
summary(univar_trans)

## unique loci across all proteins in univariate trans analysis
window_size = 100000
loci_index_univar = 1   # record the region
#univar_uniq_loci = rep(0, nrow(univar_trans))
univar_uniq_loci = c()
univar_chr = c()
univar_start = c()
univar_end = c()
univar_trans_copy = univar_trans

while (nrow(univar_trans_copy) > 0) {
  ## define leading SNP and the surrounding window
  lead_snp = which.min(univar_trans_copy$V3)[1]
  lead_chr = univar_trans_copy$V1[lead_snp]
  lead_pos = univar_trans_copy$V2[lead_snp]
  start_pos = lead_pos - window_size
  end_pos = lead_pos + window_size
  
  ## remove signals already included in the region
  rm_index = which(univar_trans_copy$V1 == lead_chr & univar_trans_copy$V2 >= start_pos &
                     univar_trans_copy$V2 <= end_pos)
  univar_trans_copy = univar_trans_copy[-rm_index, ]
  
  #loci_index = which(univar_trans$V1 == lead_chr & univar_trans$V2 >= start_pos &
  #                     univar_trans$V2 <= end_pos)
  #univar_uniq_loci[loci_index] = loci_index_univar
  univar_uniq_loci = c(univar_uniq_loci, loci_index_univar)
  univar_chr = c(univar_chr, lead_chr)
  univar_start = c(univar_start, start_pos)
  univar_end = c(univar_end, end_pos)
  loci_index_univar = loci_index_univar + 1
}
#univar_trans$uniq_loci = univar_uniq_loci
univar_uniq_loci = data.frame(loci = univar_uniq_loci, chr = univar_chr, 
                              start = univar_start, end = univar_end)

## define loci-protein pairs in univariate trans analysis
prot_uniq = unique(univar_trans$V4)
loci_index_univar = 1   # record the region
univar_trans_loci = rep(0, nrow(univar_trans))
univar_trans_target = univar_trans$V4
summary(univar_trans)
for (m in prot_uniq) {
  trans_i_prot_m = univar_trans[univar_trans$V4 == m, ]
  while (nrow(trans_i_prot_m) > 0) {
    ## define leading SNP and the surrounding window
    lead_snp = which.min(trans_i_prot_m$V3)[1]
    lead_chr = trans_i_prot_m$V1[lead_snp]
    lead_pos = trans_i_prot_m$V2[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(trans_i_prot_m$V1 == lead_chr & trans_i_prot_m$V2 >= start_pos &
                       trans_i_prot_m$V2 <= end_pos)
    trans_i_prot_m = trans_i_prot_m[-rm_index, ]
    
    loci_index = which(univar_trans$V4 == m & univar_trans$V1 == lead_chr & univar_trans$V2 >= start_pos &
                         univar_trans$V2 <= end_pos)
    univar_trans$V4[loci_index] = NA
    univar_trans_loci[loci_index] = loci_index_univar
    loci_index_univar = loci_index_univar + 1
  }
  print(m)
}
univar_trans$V4 = univar_trans_target
univar_trans$loci = univar_trans_loci

## lead SNP in each locus
univar_trans_minp = c()
for (i in 1:max(univar_trans$loci)) {
  trans_loci_i = univar_trans[univar_trans$loci == i, ]
  minp_row = trans_loci_i[which.min(trans_loci_i$V3), ]
  univar_trans_minp = rbind(univar_trans_minp, minp_row)
}

## min p value of PCO for the corresponding protein-loci pair
prot_meta = fread('/project/xuanyao/jinghui/pqtl/00_ref/prot_meta.txt')
pco_mod = c()
pco_minp = c()
pco_snp = c()
window_size = 100000

for (i in 1:nrow(univar_trans_minp)) {
  targe_i = univar_trans_minp$V4[i]
  mod_i = prot_meta$module[which(prot_meta$file == targe_i)]
  chr_i = univar_trans_minp$V1[i]
  bp_center = univar_trans_minp$V2[i]
  bp_start = bp_center - window_size
  bp_end = bp_center + window_size
  min_p = c()
  for (j in mod_i) {
    pco_p_j = readRDS(paste0('/project/xuanyao/jinghui/pqtl/03_p/', 
                           'p.mod', j, '.', chr_i, '.rds'))
    pco_p_bp = as.numeric(sapply(strsplit(names(pco_p_j), ':', fixed = T), '[', 2))
    pco_p_loci = pco_p_j[which(pco_p_bp >= bp_start & pco_p_bp <= bp_end)]
    min_p_j = which.min(pco_p_loci)[1]
    min_p = c(min_p, pco_p_loci[min_p_j])
  }
  minp_uniq = which.min(min_p)[1]
  pco_mod[i] = mod_i[minp_uniq]
  pco_minp[i] = min_p[minp_uniq]
  pco_snp[i] = names(min_p)[minp_uniq]
  print(i)
}
univar_trans_minp$pco_minp = pco_minp
univar_trans_minp$pco_mod = pco_mod
univar_trans_minp$pco_snp = pco_snp

## assign a unique loci to each loci-protein pair
univar_minp_loci = c()
for (i in 1:nrow(univar_trans_minp)) {
  chr_i = univar_trans_minp$V1[i]
  bp_i = univar_trans_minp$V2[i]
  univar_minp_loci[i] = min(univar_uniq_loci$loci[which(univar_uniq_loci$chr == chr_i &
                                                          univar_uniq_loci$start < bp_i &
                                                          univar_uniq_loci$end > bp_i)])
}
univar_trans_minp$uniq_loci = univar_minp_loci
colnames(univar_trans_minp)[c(1:4, 7, 11)] = c('chr', 'bp', 'lead_p', 'target', 'loci_prot_pair', 'genome_loci')
univar_trans_minp = univar_trans_minp[, c(5,1,2,4,3,6,7,11,8:10)]
fwrite(univar_trans_minp, '/project/xuanyao/jinghui/pqtl/04_fdr/univar_pco_comp.txt', sep = '\t')

## find protein-loci pairs with window size of 200 kb around the leading SNP for pco
pco_sig = fread('/project/xuanyao/jinghui/pqtl/04_fdr/ukb/pco_small_p.txt')
pco_sig = pco_sig[pco_sig$V2 < 5e-8/598, ]
prot_meta = fread('/project/xuanyao/jinghui/pqtl/00_ref/prot_meta.txt')

mod_uniq = unique(pco_sig$mod)
window_size = 100000
loci_index_pco = 1   # record the region
pco_trans_loci = rep(0, nrow(pco_sig))
pco_mod = pco_sig$mod
for (m in mod_uniq) {
  trans_i_prot_m = pco_sig[pco_sig$mod == m, ]
  while (nrow(trans_i_prot_m) > 0) {
    ## define leading SNP and the surrounding window
    lead_snp = which.min(trans_i_prot_m$p)[1]
    lead_chr = trans_i_prot_m$chr[lead_snp]
    lead_pos = trans_i_prot_m$bp[lead_snp]
    start_pos = lead_pos - window_size
    end_pos = lead_pos + window_size
    
    ## remove signals already included in the region
    rm_index = which(trans_i_prot_m$chr == lead_chr & 
                       trans_i_prot_m$bp >= start_pos &
                       trans_i_prot_m$bp <= end_pos)
    trans_i_prot_m = trans_i_prot_m[-rm_index, ]
    
    loci_index = which(pco_sig$mod == m & pco_sig$chr == lead_chr & 
                         pco_sig$bp >= start_pos & pco_sig$bp <= end_pos)
    pco_sig$mod[loci_index] = NA
    pco_trans_loci[loci_index] = loci_index_pco
    loci_index_pco = loci_index_pco + 1
  }
  print(m)
}
pco_sig$mod = pco_mod
pco_sig$loci = pco_trans_loci

## min p value of univariate analysis for the corresponding loci-module pair
pco_trans_minp = c()
for (i in 1:max(pco_sig$loci)) {
  trans_loci_i = pco_sig[pco_sig$loci == i, ]
  minp_row = trans_loci_i[which.min(trans_loci_i$p), ]
  pco_trans_minp = rbind(pco_trans_minp, minp_row)
}
fwrite(pco_trans_minp, '/project/xuanyao/jinghui/pqtl/04_fdr/pco_univar_comp.txt', sep = '\t')

pco_trans_minp = fread('/project/xuanyao/jinghui/pqtl/04_fdr/pco_univar_comp.txt')
trans_gene = c()
univar_minp = c()
univar_snp = c()
window_size = 100000
for (i in 1:nrow(pco_trans_minp)) {
  mod_i = pco_trans_minp$mod[i]
  chr_i = pco_trans_minp$chr[i]
  bp_center = pco_trans_minp$bp[i]
  bp_start = bp_center - window_size
  bp_end = bp_center + window_size
  
  z_mod_i = fread(paste0('/project/xuanyao/jinghui/pqtl/01_zmat/', mod_i, 
                         '.chr', chr_i, '.txt.gz'))
  z_mod_i_sub = z_mod_i[z_mod_i$V1 %in% paste0(chr_i, ':', bp_start:bp_end), ]
  z_mod_i_sub = as.data.frame(z_mod_i_sub)
  z_mod_i_sub[, -1] = abs(z_mod_i_sub[, -1])
  minp_index = which(z_mod_i_sub == max(z_mod_i_sub[,-1]), arr.ind = T)[1,]

  trans_gene[i] = colnames(z_mod_i_sub)[minp_index[2]]
  univar_snp[i] = z_mod_i_sub$V1[minp_index[1]]
  univar_minp[i] = pnorm(max(z_mod_i_sub[,-1]), lower.tail = F) * 2
}
pco_trans_minp$trans_gene = trans_gene
pco_trans_minp$univar_minp = univar_minp
pco_trans_minp$univar_snp = univar_snp
fwrite(pco_trans_minp, '/project/xuanyao/jinghui/pqtl/04_fdr/pco_univar_comp.txt', sep = '\t')


## compare pco minP and univar P
# summary of number of signals
n_signal = data.frame(n = c(84291, 172608, 1595, 5077, 657, 1029),
                      cat = c('SNP-target pair', 'SNP-target pair', 'loci-target pair', 
                              'loci-target pair', 'Genome loci', 'Genome loci'),
                      group = rep(c('PCO', 'Univar'), 3))
ggplot(n_signal, aes(x = group, y = n, fill = reorder(cat, -n))) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) + 
  geom_text(aes(label = n), vjust=-0.3, size=4, position = position_dodge(0.5)) +
  labs(x = "", y = "N") + 
  guides(fill=guide_legend(title="")) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


univar_trans_minp = fread('/project/xuanyao/jinghui/pqtl/04_fdr/univar_pco_comp.txt')
sum(univar_trans_minp$pco_minp < univar_trans_minp$V3)
sum(univar_trans_minp$pco_minp < 5e-8 / 27) / nrow(univar_trans_minp)

# Summary of univariate signal replication
comp_dat = data.frame(n_loci = c(5077, 4581, 4239, 3366), 
                      cat = c('All', 'Replicated', 'Significant in PCO', 'Smaller p in PCO'))
ggplot(comp_dat, aes(x = cat, y= n_loci)) + 
  geom_bar(stat="identity", fill = 'steelblue', width = 0.5) + 
  geom_text(aes(label = n_loci), vjust=-0.3, size=4, position = position_dodge(0.5)) +
  labs(x = "", y = "# loci") + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# minp comparison
univar_trans_minp_sub = univar_trans_minp[univar_trans_minp$pco_minp > 0, ]
univar_trans_minp_sub$neg_log_uni_p = -log10(univar_trans_minp_sub$V3)
univar_trans_minp_sub$neg_log_pco_p = -log10(univar_trans_minp_sub$pco_minp)
ggplot(univar_trans_minp_sub, aes(x = neg_log_uni_p, y= neg_log_pco_p)) + 
  geom_point(size = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(x = "-logP, univariate ", y = "-logP, PCO") + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# signal distribution
univar_pco_sig = univar_trans_minp[univar_trans_minp$pco_minp < 5e-8/27, ]
univar_pco_nonsig = univar_trans_minp[univar_trans_minp$pco_minp > 5e-8/27, ]
length(unique(univar_pco_sig$genome_loci))
length(unique(univar_pco_nonsig$genome_loci))

accu_sig = as.data.frame(table(univar_pco_sig$genome_loci))
accu_nonsig = as.data.frame(table(univar_pco_nonsig$genome_loci))
max(accu_sig$Freq)

p1 = ggplot(accu_sig, aes(Freq)) + stat_ecdf(geom = "point") + 
  xlim(0, 985) + 
  labs(x = "# sig proteins per locus", y = "Accumulated frequency", title = 'PCO significant') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2 = ggplot(accu_nonsig, aes(Freq)) + stat_ecdf(geom = "point") + 
  xlim(0, 10) + 
  labs(x = "# sig proteins per locus", y = "Accumulated frequency", title = 'PCO non-significant') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, nrow = 1)
quantile(accu_sig$Freq)
quantile(accu_nonsig$Freq)

# replication of PCO signals using univariate method
pco_trans_minp = fread('/project/xuanyao/jinghui/pqtl/04_fdr/pco_univar_comp.txt')
length(unique(pco_trans_minp$univar_snp))
sum(pco_trans_minp$univar_minp < 5e-8/4782) / nrow(pco_trans_minp)
kkk = pco_trans_minp[pco_trans_minp$univar_minp > 0 & pco_trans_minp$p > 0, ]
plot(-log10(pco_trans_minp$univar_minp), -log10(pco_trans_minp$p))
abline(0,1)
sum(pco_trans_minp$univar_minp > pco_trans_minp$p)




