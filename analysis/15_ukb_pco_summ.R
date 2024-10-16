library(data.table)
library(openxlsx)
library(ggplot2)
setwd('/project/xuanyao/jinghui/')
gene_meta = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)

prot_module = fread('pqtl/11_prot_complex/prot_module_corum_2022.txt')
corum = read.xlsx('pqtl/11_prot_complex/CORUM_2022_09_12.xlsx')
corum = corum[corum$Organism == 'Human', ]

pco_sig = fread('pqtl/04_fdr/ukb/pco_univar_comp.txt')
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
pco_sig = pco_sig[!(pco_sig$mod %in% paste0('mod', rm_mod)), ]

# find if any proteins of the target protein complex are nearby genes (on the same chr) of trans-pQTL
near_gene = c()
gene_tss = c()
for (i in 1:nrow(pco_sig)) {
  mod_i = pco_sig$mod[i]
  chr_i = paste0('chr', pco_sig$chr[i])
  pos_i = pco_sig$bp[i]
  prot_complex = unique(prot_module$complex[which(paste0('mod', prot_module$mod) == mod_i)])
  prot_list = unlist(corum$`subunits(Gene.name)`[which(corum$ComplexName == prot_complex)])
  prot_list = unique(unlist(strsplit(prot_list, ';', fixed = T)))
  prot_meta_i = gene_meta[match(prot_list, gene_meta$gene_name), ]
  prot_meta_i = prot_meta_i[prot_meta_i$chr == chr_i, ]
  if (nrow(prot_meta_i) > 0) {
    gene_index = which.min(abs(prot_meta_i$start - pos_i))
    gene_tss[i] = prot_meta_i$start[gene_index]
    near_gene[i] = prot_meta_i$gene_name[gene_index]
  } else {
    gene_tss[i] = NA
    near_gene[i] = NA
  }
  print(i)
}
pco_sig$prot_complex = prot_module$complex[match(pco_sig$mod, 
                                                 paste0('mod', prot_module$mod))]
pco_sig$subunit = corum$`subunits(Gene.name)`[match(pco_sig$prot_complex, 
                                                    corum$ComplexName)]

uniq_mod = unique(prot_module$mod)
gene_list = c()
for (i in uniq_mod) {
  mod_i = prot_module[prot_module$mod == i, ]
  gene_i = mod_i$prot
  gene_list = c(gene_list, paste(gene_i, collapse = ';'))
}
uniq_mod = paste0('mod', uniq_mod)

pco_sig$subunit_in_pco = gene_list[match(pco_sig$mod, uniq_mod)]
pco_sig$near_gene = near_gene
pco_sig$gene_tss = gene_tss
fwrite(pco_sig, 'pqtl/11_prot_complex/pco_sig_loci.txt', sep = '\t')

pco_sig = fread('pqtl/11_prot_complex/pco_sig_loci.txt')

# distance of nearby gene to the leading SNP of trans-pQTL
distance = gene_tss - pco_sig$bp
log_distance = log10(abs(distance))
log_distance[which(distance < 0)] = -log_distance[which(distance < 0)]
distance = na.omit(distance)
sum(abs(distance) < 1000000)

dist_plot = data.frame(distance = distance, log_dist = log_distance)
p1 = ggplot(dist_plot, aes(x=distance)) + 
  geom_density() + 
  labs(x = "distance, bp", y = 'density') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

dist_count = as.data.frame(table(cut(log_distance, breaks = c(-9,-8,-7,-6,6,7,8,9))))
names(dist_count) = c('log_dist', 'N')
p2 = ggplot(data=dist_count, aes(x=log_dist, y=N)) +
  geom_bar(stat="identity", fill = 'steelblue') + 
  labs(x = "distance, +/-log(bp)", y = '# loci') +
  geom_text(aes(label = N), vjust=-0.3, size=4, position = position_dodge(0.5)) +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)

# near genes also included in UKB dataset
ukb_prot_list = fread('pqtl/05_h2/04_h2_summ/prot_h2_ukb_1mb_5mb.txt')
ukb_prot_list$gene_name = gene_meta$gene_name[match(ukb_prot_list$gene_id, gene_meta$gene_id)]
ukb_prot_list = na.omit(ukb_prot_list)
ukb_overlap = pco_sig[which(pco_sig$near_gene %in% ukb_prot_list$gene_name), ]
ukb_overlap = ukb_overlap[which(abs(ukb_overlap$gene_tss - ukb_overlap$bp) < 1000000), ]

# if trans-pQTL is the cis-pQTL of the nearby gene
cis_p = c() # cis p value of the leading SNP in trans-pQTL
cis_minp = c() # min cis p value of trans-pQTL region
window_size = 100000
for (i in 1:nrow(ukb_overlap)) {
  prot_i = ukb_overlap$near_gene[i]
  chr_i = ukb_overlap$chr[i]
  bp_i = ukb_overlap$bp[i]
  #start_i = ukb_overlap$loci_start[i]
  #end_i = ukb_overlap$loci_end[i]
  file_i = prot_module$file[which(prot_module$prot == prot_i)][1]
  summ_file_list = list.files(paste0('/project/xuanyao/jinghui/pqtl/UKB_PPP/',
                                     file_i))
  summ_file_list = summ_file_list[grep(paste0('chr', chr_i, '_'), summ_file_list)]
  summ_i = fread(paste0('/project/xuanyao/jinghui/pqtl/UKB_PPP/',
                                     file_i, '/', summ_file_list))
  
  cis_pos = which(grepl(paste0(chr_i, ':', bp_i), summ_i$ID))
  cis_start = summ_i$GENPOS[cis_pos] - window_size
  cis_end = summ_i$GENPOS[cis_pos] + window_size
  cis_region_p = summ_i$LOG10P[which(summ_i$GENPOS >= cis_start & summ_i$GENPOS <= cis_end)]
  cis_minp[i] = max(cis_region_p)
  cis_p[i] = summ_i$LOG10P[cis_pos]
  print(i)
}

ukb_overlap$cis_p = 10^(-cis_p)
ukb_overlap$cis_minp = 10^(-cis_minp)

# effect direction on protein complexes for SNPs that are both significant cis- and trans-pQTL
sig_ukb_overlap = ukb_overlap[ukb_overlap$cis_p < 0.05/157, ]
beta_list = list()
p_list = list()
for (i in 1:nrow(sig_ukb_overlap)) {
  chr_i = sig_ukb_overlap$chr[i]
  bp_i = sig_ukb_overlap$bp[i]
  prot_list_i = unlist(strsplit(sig_ukb_overlap$subunit_in_pco[i], ';', fixed = T))
  file_list_i = prot_module$file[match(prot_list_i, prot_module$prot)]
  beta_i = c()
  p_i = c()
  for (j in file_list_i) {
    summ_file_list = list.files(paste0('/project/xuanyao/jinghui/pqtl/UKB_PPP/',
                                       j))
    summ_file_list = summ_file_list[grep(paste0('chr', chr_i, '_'), summ_file_list)]
    summ_j = fread(paste0('/project/xuanyao/jinghui/pqtl/UKB_PPP/',
                          j, '/', summ_file_list))
    cis_pos = which(grepl(paste0(chr_i, ':', bp_i), summ_j$ID))
    beta_i = c(beta_i, summ_j$BETA[cis_pos])
    p_i = c(p_i, summ_j$LOG10P[cis_pos])
  }
  beta_list[[i]] = beta_i
  p_list[[i]] = p_i
  print(i)
}

sig_p_list = sapply(p_list, function(x){which(x > -log10(0.05/sum(sapply(p_list, length))))})
same_dir = c()
for (i in 1:length(beta_list)) {
  same_dir[i] = ifelse(length(unique(sign(beta_list[[i]][sig_p_list[[i]]]))) > 1, 0, 1)
}
sig_ukb_overlap$same_dir = same_dir
fwrite(sig_ukb_overlap, 'pqtl/11_prot_complex/shared_cis_trans.txt', sep = '\t')



