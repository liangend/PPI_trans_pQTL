library(data.table)
setwd('/project/xuanyao/jinghui/')
pco_sig = fread('pqtl/11_prot_complex/pco_sig_loci.txt')
gene_meta = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)

prot_module = fread('pqtl/11_prot_complex/prot_module_corum_2022.txt')
prot_module = prot_module[paste0('mod', prot_module$mod) %in% pco_sig$mod, ]
prot_module$tss = gene_meta$start[match(prot_module$prot, gene_meta$gene_name)]

#### enrichment analysis of trans-pQTL of PCO are within 1 Mb of one subunit of the targeting complex
## distance of trans-pQTL to subunits on the same chromosome
sig_mod = unique(prot_module$mod)
loci_dist = c()
for (i in 1:nrow(pco_sig)) {
  chr_i = pco_sig$chr[i]
  pos_i = pco_sig$bp[i]
  dist_i = abs(pos_i - prot_module$tss)
  dist_i[which(prot_module$chr != chr_i)] = NA
  min_dist_i = aggregate(dist_i ~ prot_module$mod, FUN = min)
  all_min_i = min_dist_i$dist_i[match(sig_mod, min_dist_i$`prot_module$mod`)]
  loci_dist = rbind(loci_dist, all_min_i)
  print(i)
}

loci_dist_1mb = loci_dist
loci_dist_1mb[which(loci_dist_1mb > 1000000)] = NA

## random background
n_1mb = c()
for (i in 1:100) {
  perm_i = sample(1:ncol(loci_dist_1mb), nrow(loci_dist_1mb), replace = T)
  dist_i = loci_dist_1mb[(row(loci_dist_1mb) == 1:nrow(loci_dist_1mb) & 
                            col(loci_dist_1mb) == perm_i)]
  n_1mb[i] = sum(!(is.na(dist_i)))
}
## enrichment
sum(abs(pco_sig$bp - pco_sig$gene_tss) < 1000000, na.rm = T) / mean(n_1mb)

#### enrichment of trans-pQTL is also a cis-pQTL of the nearby subunit
pco_sig_near = pco_sig[which(abs(pco_sig$bp - pco_sig$gene_tss) < 1000000), ]
ukb_prot_list = fread('05_h2/04_h2_summ/prot_h2_ukb_1mb_5mb.txt')
ukb_prot_list$gene_name = gene_meta$gene_name[match(ukb_prot_list$gene_id, gene_meta$gene_id)]
ukb_prot_list = na.omit(ukb_prot_list)
ukb_overlap = pco_sig_near[which(pco_sig_near$near_gene %in% ukb_prot_list$gene_name), ]

## all cis genes for each pQTL
cis_p = list()
for (i in 134:nrow(ukb_overlap)) {
  chr_i = ukb_overlap$chr[i]
  bp_i = ukb_overlap$bp[i]
  # cis genes
  cis_prot = ukb_prot_list[ukb_prot_list$chr == chr_i & 
                             abs(ukb_prot_list$gene_start - bp_i) < 1000000, ]
  if (nrow(cis_prot) < 1) {
    cis_p[[i]] = 0
    next
  }
  cis_p_i = c()
  print(paste0('n cis prot: ', nrow(cis_prot)))
  for (j in 1:nrow(cis_prot)) {
    file_j = cis_prot$file[j]
    chr_all = list.files(paste0('UKB_PPP/', file_j))
    chr_j = chr_all[grep(paste0('chr', chr_i, '_'), chr_all)]
    summ_stat_j = fread(paste0('UKB_PPP/', file_j, '/', chr_j))
    bp_pos = grep(bp_i, summ_stat_j$ID)[1]
    if (length(bp_pos) > 0) {
      cis_p_i = c(cis_p_i, summ_stat_j$LOG10P[bp_pos])
    }
  }
  cis_p[[i]] = cis_p_i
  print(paste0('loci', i, ' finished'))
}

## random background
n_sig = c()
for (i in 1:100) {
  sample_p = sapply(cis_p, function(x){
    return(sample(x, 1))})
  n_sig[i] = sum(sample_p > -log10(0.05/157))
}
76 / mean(n_sig)




