library(data.table)
library(ggplot2)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')
sun_prot = fread('pqtl/05_h2/00_ref/Sun_2018_prot_w_coor.txt')
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
gtex_gene = fread('gtex/06_qtl_z/gtex_blood_meta.txt')
dgn_gene = fread('gtex/06_qtl_z/dgn_gene_meta.txt')

shared_gene = Reduce(intersect, list(sun_prot$gene, ukb_prot$gene_name, 
                                     gtex_gene$gene_name, dgn_gene$gene_id))
n_gene = length(shared_gene)

sun_small_p = fread('pqtl/04_fdr/sun_2018/small_p.txt')
sun_small_p$gene = sun_prot$gene[match(sun_small_p$file, sun_prot$target)]
sun_small_p = sun_small_p[sun_small_p$`log(P)` < log10(5e-8/n_gene) & 
                            sun_small_p$gene %in% shared_gene, ]
sun_small_p$LOG10P = -sun_small_p$`log(P)`
colnames(sun_small_p)[2:3] = c('chr', 'bp_hg19')
sun_small_p$file = sun_small_p$gene

ukb_small_p = fread('pqtl/04_fdr/ukb/ukb_univar_small_p.txt')
ukb_small_p$gene = ukb_prot$gene_name[match(ukb_small_p$file, ukb_prot$file)]
ukb_small_p = ukb_small_p[ukb_small_p$LOG10P > -log10(5e-8/n_gene) &
                            ukb_small_p$gene %in% shared_gene, ]
colnames(ukb_small_p)[1] = 'chr'

ukb_sub_small_p = fread('pqtl/04_fdr/ukb/ukb_discovery/all_prot_univar.txt')
ukb_sub_small_p$gene = ukb_prot$gene_name[match(ukb_sub_small_p$file, ukb_prot$file)]
ukb_sub_small_p = ukb_sub_small_p[ukb_sub_small_p$LOG10P > -log10(5e-8/n_gene) &
                                    ukb_sub_small_p$gene %in% shared_gene, ]
ukb_sub_small_p = na.omit(ukb_sub_small_p)
colnames(ukb_sub_small_p)[1] = 'chr'

gtex_small_p = fread('gtex/06_qtl_z/gtex_blood_small_p.txt')
colnames(gtex_small_p) = c('snp', 'gene_id', 'p', 'beta', 'beta_se', 'af')
gtex_small_p$gene = gtex_gene$gene_name[match(gtex_small_p$gene_id, gtex_gene$gene_id)]
gtex_small_p = gtex_small_p[gtex_small_p$p < 5e-8/n_gene &
                              gtex_small_p$gene %in% shared_gene, ]
gtex_small_p$LOG10P = -log10(gtex_small_p$p)
gtex_small_p$chr = sapply(strsplit(gtex_small_p$snp, '_', fixed = T), '[', 1)
gtex_small_p$chr = as.numeric(sub('chr', '', gtex_small_p$chr))
gtex_small_p$bp_hg38 = as.numeric(sapply(strsplit(gtex_small_p$snp, '_', fixed = T), '[', 2))
gtex_small_p$bp_hg19 = gtex_small_p$bp_hg38
gtex_small_p$file = gtex_small_p$gene

dgn_small_p = fread('gtex/06_qtl_z/dgn_small_p.txt')
colnames(dgn_small_p) = c('snp', 'gene', 'p', 'beta', 'beta_se', 'af')
dgn_small_p = dgn_small_p[dgn_small_p$p < 5e-8/n_gene &
                            dgn_small_p$gene %in% shared_gene, ]
dgn_small_p$LOG10P = -log10(dgn_small_p$p)
dgn_small_p$chr = as.numeric(sapply(strsplit(dgn_small_p$snp, ':', fixed = T), '[', 1))
dgn_small_p$bp_hg19 = sapply(strsplit(dgn_small_p$snp, ':', fixed = T), '[', 2)
dgn_small_p$bp_hg19 = as.numeric(sub('_2', '', dgn_small_p$bp_hg19))
dgn_small_p$file = dgn_small_p$gene

make_loci = function(data, loci_window = 250000, gtex = F) {
  sig_loci = c()
  sig_file = unique(data$file)
  loci_i = 1
  for (i in sig_file) {
    sig_i = data[data$file == i, ]
    while (nrow(sig_i > 0)) {
      lead_sig = which.max(sig_i$LOG10P)
      lead_snp = sig_i[lead_sig, ]
      lead_snp$loci = loci_i
      lead_snp$loci_start = lead_snp$bp_hg19 - loci_window
      lead_snp$loci_end = lead_snp$bp_hg19 + loci_window
      row_rm = which(sig_i$chr == lead_snp$chr & sig_i$bp_hg19 >= lead_snp$loci_start &
                       sig_i$bp_hg19 <= lead_snp$loci_end)
      sig_loci = rbind(sig_loci, lead_snp)
      loci_i = loci_i + 1
      sig_i = sig_i[-row_rm, ]
    }
    print(i)
  }
  sig_loci = sig_loci[order(sig_loci$file, sig_loci$chr, 
                            sig_loci$bp_hg19), ]
  ## merge overlapping windows
  sig_loci_new = sig_loci[1, ]
  for (i in 2:nrow(sig_loci)) {
    chr_i = sig_loci$chr[i]
    start_i = sig_loci$loci_start[i]
    file_i = sig_loci$file[i]
    bp19_i = sig_loci$bp_hg19[i]
    chr_prev = sig_loci_new$chr[nrow(sig_loci_new)]
    end_prev = sig_loci_new$loci_end[nrow(sig_loci_new)]
    file_prev = sig_loci_new$file[nrow(sig_loci_new)]
    bp19_prev = sig_loci_new$bp_hg19[nrow(sig_loci_new)]
    # merge MHC regions
    mhc_lower = ifelse(gtex, 29602228, 28477797)
    mhc_upper = ifelse(gtex, 33410226, 33448354)
    if (file_prev == file_i & chr_i == 6 & chr_prev == 6 & bp19_i > mhc_lower & bp19_i < mhc_upper &
        bp19_prev > mhc_lower & bp19_prev < mhc_upper) {
      # merge two overlapping loci
      add_pool_i = rbind(sig_loci_new[nrow(sig_loci_new), ], sig_loci[i, ])
      add_loci_i = add_pool_i[which.max(add_pool_i$LOG10P), ]
      add_loci_i$loci_start = min(add_pool_i$loci_start)
      add_loci_i$loci_end = max(add_pool_i$loci_end)
      sig_loci_new[nrow(sig_loci_new), ] = add_loci_i
      next
    }
    # keep non-overlapping windows
    if (file_prev != file_i | chr_i != chr_prev | 
        (file_prev == file_i & chr_i == chr_prev & start_i > end_prev)) {
      # add the loci directly if two loci do not overlap
      sig_loci_new = rbind(sig_loci_new, sig_loci[i, ])
    } else {
      # merge two overlapping loci
      add_pool_i = rbind(sig_loci_new[nrow(sig_loci_new), ], sig_loci[i, ])
      add_loci_i = add_pool_i[which.max(add_pool_i$LOG10P), ]
      add_loci_i$loci_start = min(add_pool_i$loci_start)
      add_loci_i$loci_end = max(add_pool_i$loci_end)
      sig_loci_new[nrow(sig_loci_new), ] = add_loci_i
    }
  }
  return(sig_loci_new)
}

sun_loci = make_loci(sun_small_p)
ukb_loci = make_loci(ukb_small_p)
ukb_sub_loci = make_loci(ukb_sub_small_p)
gtex_loci = make_loci(gtex_small_p, gtex = T)
dgn_loci = make_loci(dgn_small_p)

gene_meta_hg19 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta_hg38 = fread('gtex/00_ref/genecode.GRCh38.gene.meta.gtf')

sun_loci$gene_chr = gene_meta_hg19$chr[match(sun_loci$gene, gene_meta_hg19$gene_name)]
sun_loci$chr = paste0('chr', sun_loci$chr)
sun_loci$gene_tss = gene_meta_hg19$start[match(sun_loci$gene, gene_meta_hg19$gene_name)]

ukb_loci$gene_chr = gene_meta_hg19$chr[match(ukb_loci$gene, gene_meta_hg19$gene_name)]
ukb_loci$chr = paste0('chr', ukb_loci$chr)
ukb_loci$gene_tss = gene_meta_hg19$start[match(ukb_loci$gene, gene_meta_hg19$gene_name)]

ukb_sub_loci$gene_chr = gene_meta_hg19$chr[match(ukb_sub_loci$gene, gene_meta_hg19$gene_name)]
ukb_sub_loci$chr = paste0('chr', ukb_sub_loci$chr)
ukb_sub_loci$gene_tss = gene_meta_hg19$start[match(ukb_sub_loci$gene, gene_meta_hg19$gene_name)]

gtex_loci$gene_chr = gene_meta_hg38$chr[match(gtex_loci$gene, gene_meta_hg38$gene_name)]
gtex_loci$chr = paste0('chr', gtex_loci$chr)
gtex_loci$gene_tss = gene_meta_hg38$start[match(gtex_loci$gene, gene_meta_hg38$gene_name)]

dgn_loci$gene_chr = gene_meta_hg38$chr[match(dgn_loci$gene, gene_meta_hg38$gene_name)]
dgn_loci$chr = paste0('chr', dgn_loci$chr)
dgn_loci$gene_tss = gene_meta_hg38$start[match(dgn_loci$gene, gene_meta_hg38$gene_name)]

sun_loci$cis_trans = 'trans'
sun_loci$cis_trans[which(sun_loci$chr == sun_loci$gene_chr & 
                           abs(sun_loci$bp_hg19 - sun_loci$gene_tss) < 1000000)] = 'cis'
ukb_loci$cis_trans = 'trans'
ukb_loci$cis_trans[which(ukb_loci$chr == ukb_loci$gene_chr & 
                           abs(ukb_loci$bp_hg19 - ukb_loci$gene_tss) < 1000000)] = 'cis'
ukb_sub_loci$cis_trans = 'trans'
ukb_sub_loci$cis_trans[which(ukb_sub_loci$chr == ukb_sub_loci$gene_chr & 
                           abs(ukb_sub_loci$bp_hg19 - ukb_sub_loci$gene_tss) < 1000000)] = 'cis'
gtex_loci$cis_trans = 'trans'
gtex_loci$cis_trans[which(gtex_loci$chr == gtex_loci$gene_chr & 
                            abs(gtex_loci$bp_hg19 - gtex_loci$gene_tss) < 1000000)] = 'cis'
dgn_loci$cis_trans = 'trans'
dgn_loci$cis_trans[which(dgn_loci$chr == dgn_loci$gene_chr & 
                           abs(dgn_loci$bp_hg19 - dgn_loci$gene_tss) < 1000000)] = 'cis'

## cis trans snp proportion
qtl_table = data.frame(n_qtl = c(545, 9499, 4943, 185, 401),
                       n_cis = c(186, 576, 548, 182, 308), 
                       n_trans = c(349, 8923, 4395, 3, 93),
                       data = c('INTERVAL', 'UKB-PPP', 'UKB-PPP discovery', 'GTEx', 'DGN'),
                       type = c('protein', 'protein', 'protein', 'mRNA', 'mRNA'))
qtl_table$cis_prop = qtl_table$n_cis / qtl_table$n_qtl
qtl_table$trans_prop = qtl_table$n_trans / qtl_table$n_qtl
qtl_prop_plot = reshape(qtl_table[,-(1:3)], direction = "long", idvar = c('data', 'type'), 
                   varying = c('cis_prop', 'trans_prop'), v.names = 'proportion',
                   timevar = 'cis_trans', times = c('cis', 'trans'))

ggplot(qtl_prop_plot, aes(x = factor(data, levels = c('INTERVAL', 'UKB-PPP', 'UKB-PPP discovery', 'GTEx', 'DGN')), 
                     y = proportion, 
                     fill = cis_trans)) + 
  geom_bar(stat="identity", color="black", width = 0.6) +
  labs(x = "", fill = '', y = 'Proportion') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


## cis trans gene proportion
eqtlgen_cis = fread('gtex/00_ref/eQTLGen/cis_small_p.txt.gz', select = c(1,8,9))
colnames(eqtlgen_cis) = c('p',	'Gene', 'GeneSymbol')
eqtlgen_cis = eqtlgen_cis[eqtlgen_cis$p < 5e-8/n_gene, ]
eqtlgen_cis_gene = unique(eqtlgen_cis$GeneSymbol)

eqtlgen_trans = fread('gtex/00_ref/eQTLGen/trans_500k_loci.txt')
eqtlgen_trans = eqtlgen_trans[eqtlgen_trans$Pvalue < 5e-8/n_gene, ]
eqtlgen_trans_gene = unique(eqtlgen_trans$GeneSymbol)

gene_table = data.frame(gene = shared_gene)
sun_group = c()
ukb_group = c()
ukb_sub_group = c()
gtex_group = c()
dgn_group = c()
eqtlgen_group = c()
for (i in 1:nrow(gene_table)) {
  gene_i = gene_table$gene[i]
  sun_cis_trans = unique(sun_loci$cis_trans[which(sun_loci$gene == gene_i)])
  if (length(sun_cis_trans) == 0) {
    sun_group[i] = 'no signal'
  } else if (length(sun_cis_trans) == 1){
    sun_group[i] = paste0(sun_cis_trans, ' only')
  } else {
    sun_group[i] = 'cis & trans'
  }
  
  ukb_cis_trans = unique(ukb_loci$cis_trans[which(ukb_loci$gene == gene_i)])
  if (length(ukb_cis_trans) == 0) {
    ukb_group[i] = 'no signal'
  } else if (length(ukb_cis_trans) == 1){
    ukb_group[i] = paste0(ukb_cis_trans, ' only')
  } else {
    ukb_group[i] = 'cis & trans'
  }
  
  ukb_sub_cis_trans = unique(ukb_sub_loci$cis_trans[which(ukb_sub_loci$gene == gene_i)])
  if (length(ukb_sub_cis_trans) == 0) {
    ukb_sub_group[i] = 'no signal'
  } else if (length(ukb_sub_cis_trans) == 1){
    ukb_sub_group[i] = paste0(ukb_sub_cis_trans, ' only')
  } else {
    ukb_sub_group[i] = 'cis & trans'
  }
  
  gtex_cis_trans = unique(gtex_loci$cis_trans[which(gtex_loci$gene == gene_i)])
  if (length(gtex_cis_trans) == 0) {
    gtex_group[i] = 'no signal'
  } else if (length(gtex_cis_trans) == 1){
    gtex_group[i] = paste0(gtex_cis_trans, ' only')
  } else {
    gtex_group[i] = 'cis & trans'
  }
  
  dgn_cis_trans = unique(dgn_loci$cis_trans[which(dgn_loci$gene == gene_i)])
  if (length(dgn_cis_trans) == 0) {
    dgn_group[i] = 'no signal'
  } else if (length(dgn_cis_trans) == 1){
    dgn_group[i] = paste0(dgn_cis_trans, ' only')
  } else {
    dgn_group[i] = 'cis & trans'
  }
  
  if (gene_i %in% eqtlgen_cis_gene & gene_i %in% eqtlgen_trans_gene) {
    eqtlgen_group[i] = 'cis & trans'
  } else if (gene_i %in% eqtlgen_cis_gene & !(gene_i %in% eqtlgen_trans_gene)){
    eqtlgen_group[i] = 'cis only'
  } else if (gene_i %in% eqtlgen_trans_gene & !(gene_i %in% eqtlgen_cis_gene)){
    eqtlgen_group[i] = 'trans only'
  } else {
    eqtlgen_group[i] = 'no signal'
  }
}
gene_table$sun_signal = sun_group
gene_table$ukb_signal = ukb_group
gene_table$ukb_sub_signal = ukb_sub_group
gene_table$gtex_signal = gtex_group
gene_table$dgn_signal = dgn_group
gene_table$eqtlgen_signal = eqtlgen_group
fwrite(gene_table, 'pqtl/04_fdr/46_gene_cis_trans.txt', sep = '\t')

gene_table = fread('pqtl/04_fdr/46_gene_cis_trans.txt')
sun_signal = as.data.frame(table(gene_table$sun_signal))
ukb_signal = as.data.frame(table(gene_table$ukb_signal))
ukb_sub_signal = as.data.frame(table(gene_table$ukb_sub_signal))
gtex_signal = as.data.frame(table(gene_table$gtex_signal))
dgn_signal = as.data.frame(table(gene_table$dgn_signal))
eqtlgen_signal = as.data.frame(table(gene_table$eqtlgen_signal))

colnames(sun_signal) = c('signal', 'proportion')
colnames(ukb_signal) = c('signal', 'proportion')
colnames(ukb_sub_signal) = c('signal', 'proportion')
colnames(gtex_signal) = c('signal', 'proportion')
colnames(dgn_signal) = c('signal', 'proportion')
colnames(eqtlgen_signal) = c('signal', 'proportion')

gene_cat = rbind(sun_signal, ukb_signal, ukb_sub_signal, 
                 gtex_signal, dgn_signal, eqtlgen_signal)
gene_cat$data = rep(c('INTERVAL', 'UKB-PPP', 'UKB-PPP discovery', 
                      'GTEx', 'DGN', 'eQTLGen'), each = 4)
gene_cat$group = c(rep('Protein', 12), rep('mRNA', 12))

ggplot(gene_cat, aes(x = factor(data, levels = c('INTERVAL', 'UKB-PPP', 'UKB-PPP discovery', 
                                                 'GTEx', 'DGN', 'eQTLGen')), 
                     y = proportion, 
                     fill = factor(signal, levels = c('cis only', 'trans only', 
                                                      'cis & trans', 'no signal')))) + 
  geom_bar(stat="identity", color="black", width = 0.6) +
  labs(x = "", fill = '', y = 'Number of genes') +
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




