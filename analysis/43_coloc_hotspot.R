library(data.table)
library(ggplot2)
library(RColorBrewer)
setwd('/project/xuanyao/jinghui')

all_trait = list.files('pqtl/08_gwas_coloc/gwas_coloc_mcl/')
all_trait = all_trait[!grepl('.txt', all_trait)]
immune_trait = c('CD', 'IBD', 'MS', 'UC', 'lupus', 'RA', 'asthma')
blood_trait = c("BASO","BASO-P","EO","EO-P","HCT",
                "HGB","HLSR","HLSR-P","IRF","LYMPH",
                "LYMPH-P","MCH","MCHC","MCV","MONO","MONO-P","MPV","MRV",
                "MSCV","NEUT","NEUT-P","PDW","PLT","PLTC","RBC","RBC-W",
                "RET","RET-P","WBC")
other_trait = c("AD","BMI","HDL","LDL","PD","T1D","T2D","WHR",
                "bipolar","height","hypertension",
                "schizo","stroke","weight")

gwas_i = fread(paste0('gtex/10_GWAS_coloc/GWAS_sum_stats/gwas_loci_500kb/', 
                      trait_i,'.txt'))
pp4_trans_pqtl_i = pp4_pqtl_trans[pp4_pqtl_trans$trait == trait_i, ]
pp4_mcl_i = pp4_mcl_sig[pp4_mcl_sig$trait == trait_i, ]
pp4_mcl_i = cbind(pp4_mcl_i, gwas_i[match(pp4_mcl_i$gwas_loci, gwas_i$loci), 
                                    c('CHR', 'loci_start', 'loci_end', 'P')])
table(pp4_mcl_i$CHR)
immune_full_i = fread(paste0('pqtl/08_gwas_coloc/pco_mcl_coloc/pp4_immune/', 
                             trait_i, '_full.txt'))

pp4_mcl_i_nov = pp4_mcl_i[!pp4_mcl_i$gwas_loci %in% pp4_trans_pqtl_i$gwas_loci, ]

immune_full_all = c()
for (i in immune_trait) {
  immune_full_i = fread(paste0('pqtl/08_gwas_coloc/pco_mcl_coloc/pp4_immune/', 
                               i, '_full.txt'))
  immune_full_i$trait = i
  immune_full_all = rbind(immune_full_all, immune_full_i)
}

pco_all = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')
pp4_all = fread('pqtl/08_gwas_coloc/pco_mcl_coloc/mcl_pp4.txt')
pp4_all = as.data.frame(pp4_all)
pp4_sig = c()
for (t in 3:52) {
  trait_i = pp4_all[, c(2, t)]
  colnames(trait_i)[2] = 'pp4'
  trait_i = trait_i[trait_i$pp4 > 0.75, ]
  trait_i$nearest_gene = pco_all$nearest_gene[match(trait_i$loci, pco_all$loci)]
  trait_i$trait = colnames(pp4_all)[t]
  pp4_sig = rbind(pp4_sig, trait_i)
  print(t)
}

top_gene = names(sort(table(pp4_sig$nearest_gene), 
                      decreasing = T)[1:200])
all_trait = unique(pp4_sig$trait)
top_n_all = c()
for (t in all_trait) {
  pp4_sig_i = pp4_sig[pp4_sig$trait == t, ]
  pp4_sig_i$mod = pco_all$mod[match(pp4_sig_i$loci, pco_all$loci)]
  mod_count = aggregate(mod ~ nearest_gene, data = pp4_sig_i, function(x){length(unique(x))})
  n_mod = mod_count$mod
  names(n_mod) = mod_count$nearest_gene
  top_n_i = n_mod[top_gene]
  top_n_i[is.na(top_n_i)] = 0
  top_n_i = unname(top_n_i)
  top_n_all = c(top_n_all, top_n_i)
}

top_tab = data.frame(n_coloc = top_n_all, 
                       gene = rep(top_gene, length(all_trait)),
                       trait = rep(all_trait, each = length(top_gene)))

gene_meta_hg37 = fread('pqtl/05_h2/00_ref/gene.meta.hg37.gft')
top_tab$chr = pco_all$chr[match(top_tab$gene, pco_all$nearest_gene)]
top_tab$tss = gene_meta_hg37$start[match(top_tab$gene, gene_meta_hg37$gene_name)]

chr_pos = readRDS('pqtl/00_ref/chromosome_location.rds')
top_tab$tot = chr_pos$tot[match(top_tab$chr, chr_pos$CHR)]
top_tab$def_pos = top_tab$tss + top_tab$tot

top_sub = top_tab[top_tab$trait %in% c(immune_trait), ]

mark_tab = data.frame(mark_gene = top_gene[1:20])
mark_tab$chr = pco_all$chr[match(mark_tab$mark_gene, pco_all$nearest_gene)]
mark_tab$tss = gene_meta_hg37$start[match(mark_tab$mark_gene, gene_meta_hg37$gene_name)]
mark_tab = mark_tab[order(mark_tab$chr, mark_tab$tss), ]
rownames(mark_tab) = 1:nrow(mark_tab)
mark_tab = mark_tab[c(1,3,4,10,13,15,19), ]

ggplot(data = top_sub) + geom_rect(data = chr_pos, aes(xmin = tot, xmax = xmax, ymin = -Inf,
                                                       ymax = Inf, fill = factor(CHR))) + 
  geom_bar(aes(x = def_pos, y = n_coloc), color = 'red', 
           stat="identity") +
  scale_fill_manual(values = rep(c("#e5e5e5", "#ffffff"), 11), guide = "none") +
  facet_grid(trait~., scales="free") + 
  labs(x = 'Chromosome', y = 'Number of PPI clusters') +
  geom_vline(xintercept = top_sub$def_pos[top_sub$gene %in% mark_tab$mark_gene],
             linetype="dotted", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
                     label = chr_pos$CHR, breaks = chr_pos$center, expand = c(0, 0)) + 
  theme(text = element_text(size=13, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.text.y = element_text(colour = 'black', size = 10),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')






