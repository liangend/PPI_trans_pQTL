library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
### enhancer-gene pairs
# liver
# state6 = fread('/project/xuanyao/jinghui/pqtl/09_enhancer/RoadmapLinks/links_E066_6_2.5.txt')
# state7 = fread('/project/xuanyao/jinghui/pqtl/09_enhancer/RoadmapLinks/links_E066_7_2.5.txt')
# state12 = fread('/project/xuanyao/jinghui/pqtl/09_enhancer/RoadmapLinks/links_E066_12_2.5.txt')
# 
# enh_region = rbind(state6, state7, state12)
# enh_region = enh_region[order(enh_region$V1, enh_region$V2), ]
# colnames(enh_region) = c('chr', 'start', 'end', 'gene_id', 'pred_score', 'distance_to_TSS')
# 
# gene_meta = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/gene.meta.hg37.gft')
# gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
# enh_region$gene = gene_meta$gene_name[match(enh_region$gene_id, gene_meta$gene_id)]

# across tissues
enh_list = fread('/project/xuanyao/jinghui/gtex/99_scripts/analysis/06_sig_enrich/roadmap_tissue_ID.txt',
                 header = F)
enh_len = list()
for (i in 1:nrow(enh_list)) {
  tissue_i = enh_list$V1[i]
  state6 = fread(paste0('/project/xuanyao/jinghui/pqtl/09_enhancer/RoadmapLinks/links_', 
                        tissue_i, '_6_2.5.txt'))
  state7 = fread(paste0('/project/xuanyao/jinghui/pqtl/09_enhancer/RoadmapLinks/links_', 
                        tissue_i, '_7_2.5.txt'))
  state12 = fread(paste0('/project/xuanyao/jinghui/pqtl/09_enhancer/RoadmapLinks/links_', 
                         tissue_i, '_12_2.5.txt'))
  enh_region = rbind(state6, state7, state12)
  enh_region = enh_region[order(enh_region$V1, enh_region$V2), ]
  colnames(enh_region) = c('chr', 'start', 'end', 'gene_id', 'pred_score', 'distance_to_TSS')
  enh_len_i = as.data.frame(table(enh_region$gene_id))
  enh_len[[i]] = enh_len_i
}

all_gene_id = unlist(sapply(enh_len, function(x){return(x$Var1)}))
all_gene_id = unique(as.character(all_gene_id))
enh_len_uniq = sapply(enh_len, function(x){return(x$Freq[match(all_gene_id, x$Var1)])})
enh_len_uniq = as.data.frame(enh_len_uniq)
enh_len_uniq[is.na(enh_len_uniq)] = 0
enh_len_uniq = enh_len_uniq * 200    # each enh annotation is 200 bp length
ave_enh_len = apply(enh_len_uniq, 1, mean)
n_enh_tissue = apply(enh_len_uniq, 1, function(x){sum(x > 0)})

enh_table = data.frame(gene_id = all_gene_id, ave_enh_len = ave_enh_len,
                       n_enh_tissue = n_enh_tissue)

gene_meta = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/gene.meta.hg37.gft')
gene_meta$gene_id = sapply(strsplit(gene_meta$gene_id, '.', fixed = T), '[', 1)
enh_table$gene = gene_meta$gene_name[match(enh_table$gene_id, gene_meta$gene_id)]

## 2018 Sun
sun_pqtl = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2018_Sun_supplement.xlsx', 
                     sheet = 4, startRow = 5, cols = 1:16)
sun_prot_list = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/Sun_prot_w_coor.txt')
sun_prot_list$is_sig = sun_prot_list$target %in% sun_pqtl$SOMAmer.ID

uniq_targe = unique(sun_pqtl$SOMAmer.ID)
signal = c()
for (i in 1:length(uniq_targe)) {
  signal_i = unique(sun_pqtl$`cis/.trans`[which(sun_pqtl$SOMAmer.ID == uniq_targe[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
sun_prot_list$signal = signal[match(sun_prot_list$target, uniq_targe)]
sun_prot_list$signal[which(is.na(sun_prot_list$signal))] = 'no'

sun_prot_list$enh_length = enh_table$ave_enh_len[match(sun_prot_list$gene, enh_table$gene)]
sun_prot_list$n_enh_tissue = enh_table$n_enh_tissue[match(sun_prot_list$gene, enh_table$gene)]
sun_prot_list$is_cis = 1
sun_prot_list$is_cis[which(sun_prot_list$signal %in% c('no', 'trans'))] = 0

## 2022 Gudjonsson
gud_pqtl = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2022_Gudjonsson_pqtl.xlsx', 
                     startRow = 4, cols = 1:15)
gud_prot_list = fread('/project/xuanyao/jinghui/pqtl/05_h2/04_h2_summ/prot_h2_Gudjonsson.txt')
gud_prot_info = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2022_Gudjonsson_prot_info.xlsx', 
                          startRow = 3)
gud_prot_list$study = gud_prot_info$Study.tag[match(gud_prot_list$target, gud_prot_info$Study.Accession)]
gud_prot_list = gud_prot_list[,c(1,3,4,10)]

gud_pqtl_merge = merge(gud_pqtl, gud_prot_list, by.x = 'SOMAmer', by.y = 'study', all.x = T)
gud_pqtl_trans = gud_pqtl_merge[gud_pqtl_merge$`cis/trans` == 'trans', ]
gud_pqtl_merge$Chr = paste0('chr', gud_pqtl_merge$Chr)
cis_index = which(gud_pqtl_merge$`cis/trans` == 'trans' &
                    gud_pqtl_merge$Chr == gud_pqtl_merge$chr &
                    abs(gud_pqtl_merge$tss - gud_pqtl_merge$`Pos.(GRCh37)`) < 1000000)
gud_pqtl_merge$`cis/trans`[cis_index] = 'cis'

gud_prot_list$is_sig = gud_prot_list$study %in% gud_pqtl_merge$SOMAmer
uniq_targe = unique(gud_pqtl_merge$`Protein.(Entrez.symbol)`)
signal = c()
for (i in 1:length(uniq_targe)) {
  signal_i = unique(gud_pqtl_merge$`cis/trans`[which(gud_pqtl_merge$`Protein.(Entrez.symbol)` == uniq_targe[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
gud_prot_list$signal = signal[match(gud_prot_list$prot, uniq_targe)]
gud_prot_list$signal[which(is.na(gud_prot_list$signal))] = 'no'

gud_prot_list$enh_length = enh_table$ave_enh_len[match(gud_prot_list$prot, enh_table$gene)]
gud_prot_list$n_enh_tissue = enh_table$n_enh_tissue[match(gud_prot_list$prot, enh_table$gene)]
gud_prot_list$is_cis = 1
gud_prot_list$is_cis[which(gud_prot_list$signal %in% c('no', 'trans'))] = 0

# gud_enh_length = c()
# for (i in 1:nrow(gud_prot_list)) {
#   gene_i = gud_prot_list$prot[i]
#   enh_i = enh_region[enh_region$gene == gene_i, ]
#   gud_enh_length[i] = 200 * nrow(enh_i)
# }
# gud_prot_list$enh_length = gud_enh_length


## 2022 ukb
ukb_pqtl = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                     startRow = 3, cols = 1:11, sheet = 7)
ukb_prot_list = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2022_ukbiobank_prot.xlsx', 
                          startRow = 3, sheet = 4)
ukb_prot_list$is_sig = ukb_prot_list$UniProt %in% ukb_pqtl$Target.UniProt
ukb_pqtl$gene_chr = ukb_prot_list$Gene.CHROM[match(ukb_pqtl$Assay.Target, ukb_prot_list$Assay.Target)]
ukb_pqtl$gene_tss = ukb_prot_list$Gene.start[match(ukb_pqtl$Assay.Target, ukb_prot_list$Assay.Target)]
ukb_pqtl$gene_tss = sapply(strsplit(ukb_pqtl$gene_tss, ';', fixed = T), '[', 1)
ukb_pqtl$gene_tss = as.numeric(ukb_pqtl$gene_tss)

ukb_pqtl$cis_trans = 'trans'
ukb_pqtl$cis_trans[which(ukb_pqtl$gene_chr == ukb_pqtl$CHROM &
                           ukb_pqtl$`GENPOS.(hg38)` < ukb_pqtl$gene_tss + 1000000 &
                           ukb_pqtl$`GENPOS.(hg38)` > ukb_pqtl$gene_tss - 1000000)] = 'cis'

uniq_target = unique(ukb_pqtl$Assay.Target)
signal = c()
for (i in 1:length(uniq_target)) {
  signal_i = unique(ukb_pqtl$cis_trans[which(ukb_pqtl$Assay.Target == uniq_target[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
ukb_prot_list$signal = signal[match(ukb_prot_list$Assay.Target, uniq_target)]
ukb_prot_list$signal[which(is.na(ukb_prot_list$signal))] = 'no'

ukb_prot_list$enh_length = enh_table$ave_enh_len[match(ukb_prot_list$Gene.symbol, enh_table$gene)]
ukb_prot_list$n_enh_tissue = enh_table$n_enh_tissue[match(ukb_prot_list$Gene.symbol, enh_table$gene)]
ukb_prot_list$is_cis = 1
ukb_prot_list$is_cis[which(ukb_prot_list$signal %in% c('no', 'trans'))] = 0

# ukb_enh_length = c()
# for (i in 1:nrow(ukb_prot_list)) {
#   gene_i = ukb_prot_list$Gene.symbol[i]
#   enh_i = enh_region[enh_region$gene == gene_i, ]
#   ukb_enh_length[i] = 200 * nrow(enh_i)
# }
# ukb_prot_list$enh_length = ukb_enh_length

## generate se through bootstrapping
library(boot)
boot_mean = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, mean)
  return(aggre_mean$enh_length)
}
sun_boot = boot(sun_prot_list, boot_mean, 1000, formula = enh_length~is_cis)
gud_boot = boot(gud_prot_list, boot_mean, 1000, formula = enh_length~is_cis)
ukb_boot = boot(ukb_prot_list, boot_mean, 1000, formula = enh_length~is_cis)

enh_len_plot = data.frame(enh_len = c(apply(sun_boot$t, 2, mean), apply(gud_boot$t, 2, mean),
                                      apply(ukb_boot$t, 2, mean)),
                          se = c(apply(sun_boot$t, 2, sd), apply(gud_boot$t, 2, sd), 
                                 apply(ukb_boot$t, 2, sd)), 
                          is_cis = rep(c('no', 'yes'), 3), 
                          data = rep(c('Sun', 'Gudjonsson', 'UKB'), each = 2))
library(scales)
p1 = ggplot(enh_len_plot, aes(x=is_cis, y=enh_len, colour=data)) + 
  geom_errorbar(aes(ymin=enh_len-se, ymax=enh_len+se), width=.1, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  geom_hline(yintercept = mean(gud_prot_list$enh_length, na.rm = T), lty = 2, color = hue_pal()(3)[1]) +
  geom_hline(yintercept = mean(sun_prot_list$enh_length, na.rm = T), lty = 2, color = hue_pal()(3)[2]) +
  geom_hline(yintercept = mean(ukb_prot_list$enh_length, na.rm = T), lty = 2, color = hue_pal()(3)[3]) +
  labs(x = "If containing cis signal", y = 'Enhancer length, bp', title = 'Ave enh length across all tissues') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = 'none')

## construct logistic regression using enhancer length and # active tissues as covariates
# normalize enhancer length and # active tissues
pLI = fread('/project/xuanyao/jinghui/pqtl/06_LoF/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz')

sun_prot_list$gene_length = pLI$gene_length[match(sun_prot_list$gene, pLI$gene)]
sun_prot_list = na.omit(sun_prot_list)
sun_prot_list$norm_enh_len = (sun_prot_list$enh_length - mean(sun_prot_list$enh_length)) / 
  sd(sun_prot_list$enh_length)
sun_prot_list$norm_n_tissue = (sun_prot_list$n_enh_tissue - mean(sun_prot_list$n_enh_tissue)) / 
  sd(sun_prot_list$n_enh_tissue)
sun_prot_list$norm_gene_len = (sun_prot_list$gene_length - mean(sun_prot_list$gene_length)) / 
  sd(sun_prot_list$gene_length)

gud_prot_list$gene_length = pLI$gene_length[match(gud_prot_list$prot, pLI$gene)]
gud_prot_list = na.omit(gud_prot_list)
gud_prot_list$norm_enh_len = (gud_prot_list$enh_length - mean(gud_prot_list$enh_length)) / 
  sd(sun_prot_list$enh_length)
gud_prot_list$norm_n_tissue = (gud_prot_list$n_enh_tissue - mean(gud_prot_list$n_enh_tissue)) / 
  sd(gud_prot_list$n_enh_tissue)
gud_prot_list$norm_gene_len = (gud_prot_list$gene_length - mean(gud_prot_list$gene_length)) / 
  sd(gud_prot_list$gene_length)

ukb_prot_list = na.omit(ukb_prot_list)
ukb_prot_list$gene_length = as.numeric(ukb_prot_list$Gene.end) - 
  as.numeric(ukb_prot_list$Gene.start) + 1
ukb_prot_list$norm_enh_len = (ukb_prot_list$enh_length - mean(ukb_prot_list$enh_length)) / 
  sd(sun_prot_list$enh_length)
ukb_prot_list$norm_n_tissue = (ukb_prot_list$n_enh_tissue - mean(ukb_prot_list$n_enh_tissue)) / 
  sd(sun_prot_list$n_enh_tissue)
ukb_prot_list$norm_gene_len = (ukb_prot_list$gene_length - mean(ukb_prot_list$gene_length)) / 
  sd(ukb_prot_list$gene_length)

fit_sun = glm(is_cis ~ norm_enh_len + norm_n_tissue + norm_gene_len, 
              data = sun_prot_list, family = "binomial")
fit_gud = glm(is_cis ~ norm_enh_len + norm_n_tissue + norm_gene_len, 
              data = gud_prot_list, family = "binomial")
fit_ukb = glm(is_cis ~ norm_enh_len + norm_n_tissue + norm_gene_len, 
              data = ukb_prot_list, family = "binomial")

logit_reg_coef = as.data.frame(rbind(summary(fit_sun)$coefficients[, 1:2], 
                                     summary(fit_gud)$coefficients[, 1:2],
                                     summary(fit_ukb)$coefficients[, 1:2]))
logit_reg_coef$data = rep(c('Sun', 'Gudjonsson', 'UKB'), each = 4)
logit_reg_coef$coef = rep(c('Intercept', 'Enhancer length', 
                            'Count of active tissues', 'Gene length'), 3)

p2 = ggplot(logit_reg_coef[logit_reg_coef$coef != 'Intercept', ], aes(x=coef, y=Estimate, colour=data)) + 
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), width=.1, 
                position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = 'coeff', title = 'If a gene contains cis signal as response variable') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)

### enhancer length, gene length and pLI
# gene_count = as.data.frame(table(enh_table$gene))
# colnames(gene_count) = c('gene', 'enh_length')
# gene_count$enh_length = gene_count$enh_length * 200
enh_table$gene_length = pLI$gene_length[match(enh_table$gene, pLI$gene)]
enh_table$pLI = pLI$pLI[match(enh_table$gene, pLI$gene)]
enh_table = na.omit(enh_table)

enh_table$pLI_bin = cut(enh_table$pLI, breaks=seq(0, 1, 0.1), right = T, labels = F)
enh_table$gene_len_bin = cut(enh_table$gene_length, 
                             breaks=quantile(enh_table$gene_length, seq(0, 1, 0.1)), 
                             right = T, labels = F)

gene_len_summ = as.data.frame(aggregate(gene_length ~ pLI_bin, data = enh_table, mean))
gene_len_summ$se = aggregate(gene_length ~ pLI_bin, data = enh_table, sd)[,2] / table(enh_table$pLI_bin)^0.5
gene_len_summ$pLI_bin = gene_len_summ$pLI_bin/10
p1 = ggplot(gene_len_summ, aes(x=pLI_bin, y=gene_length)) + 
  geom_point(color = 'steelblue') +
  geom_errorbar(aes(ymin=gene_length-se, ymax=gene_length+se), width=.01, color = 'steelblue') + 
  labs(x = "pLI", y = 'Gene length, bp') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

enh_len_summ = as.data.frame(aggregate(ave_enh_len ~ pLI_bin, data = enh_table, mean))
enh_len_summ$se = aggregate(ave_enh_len ~ pLI_bin, data = enh_table, sd)[,2] / table(enh_table$pLI_bin)^0.5
enh_len_summ$pLI_bin = enh_len_summ$pLI_bin/10
p2 = ggplot(enh_len_summ, aes(x=pLI_bin, y=ave_enh_len)) + 
  geom_point(color = 'steelblue') +
  geom_errorbar(aes(ymin=ave_enh_len-se, ymax=ave_enh_len+se), width=.01, color = 'steelblue') + 
  labs(x = "pLI", y = 'Enhancer length, bp') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

enh_len_summ2 = as.data.frame(aggregate(ave_enh_len ~ gene_len_bin, data = enh_table, mean))
enh_len_summ2$se = aggregate(ave_enh_len ~ gene_len_bin, data = enh_table, sd)[,2] / table(enh_table$gene_len_bin)^0.5
enh_len_summ2$gene_len_bin = enh_len_summ2$gene_len_bin/10
p3 = ggplot(enh_len_summ2, aes(x=gene_len_bin, y=ave_enh_len)) + 
  geom_point(color = 'steelblue') +
  geom_errorbar(aes(ymin=ave_enh_len-se, ymax=ave_enh_len+se), width=0.01, color = 'steelblue') +
  labs(x = "Gene length quantile", y = 'Enhancer length, bp') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, p3, nrow = 2)
cor.test(enh_table$pLI, enh_table$gene_length)
cor.test(enh_table$ave_enh_len, enh_table$pLI)
cor.test(enh_table$ave_enh_len, enh_table$gene_length)


