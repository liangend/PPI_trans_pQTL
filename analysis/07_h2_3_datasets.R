library(data.table)
library(ggplot2)
library(openxlsx)
library(gridExtra)
library(ggsignif)
h2_summ_Gud = fread('/project/xuanyao/jinghui/pqtl/05_h2/04_h2_summ/prot_h2_gud_1mb_5mb.txt')
h2_summ_sun = fread('/project/xuanyao/jinghui/pqtl/05_h2/04_h2_summ/prot_h2_sun_1mb_5mb.txt')
h2_summ_ukb = fread('/project/xuanyao/jinghui/pqtl/05_h2/04_h2_summ/prot_h2_ukb_1mb_5mb.txt')

trans_prop_gud = mean(h2_summ_Gud$trans_5mb_h2_ind) / 
  mean(h2_summ_Gud$trans_5mb_h2_ind + h2_summ_Gud$cis_1mb_h2)
trans_prop_sun = mean(h2_summ_sun$trans_5mb_h2_ind) / 
  mean(h2_summ_sun$trans_5mb_h2_ind + h2_summ_sun$cis_1mb_h2)
trans_prop_ukb = mean(h2_summ_ukb$trans_5mb_h2_ind) / 
  mean(h2_summ_ukb$trans_5mb_h2_ind + h2_summ_ukb$cis_1mb_h2)

trans_se_gud = sd(h2_summ_Gud$trans_5mb_h2_ind / (h2_summ_Gud$trans_5mb_h2_ind + h2_summ_Gud$cis_1mb_h2)) / 
  sqrt(nrow(h2_summ_Gud))
trans_se_sun = sd(h2_summ_sun$trans_5mb_h2_ind / (h2_summ_sun$trans_5mb_h2_ind + h2_summ_sun$cis_1mb_h2)) /
  sqrt(nrow(h2_summ_sun))
trans_se_ukb = sd(h2_summ_ukb$trans_5mb_h2_ind / (h2_summ_ukb$trans_5mb_h2_ind + h2_summ_ukb$cis_1mb_h2)) /
  sqrt(nrow(h2_summ_ukb))

trans_prop_all = data.frame(trans_prop = c(trans_prop_sun, trans_prop_gud, trans_prop_ukb),
                            se = c(trans_se_sun, trans_se_gud, trans_se_ukb), 
                            data = c('Sun', 'Gudjonsson', 'UKB'))
ggplot(trans_prop_all, aes(x=factor(data, levels = c('Sun', 'Gudjonsson', 'UKB')), 
                           y=trans_prop, color=data)) +
  geom_errorbar(aes(ymin=trans_prop-se, ymax=trans_prop+se), width=.05) +
  geom_point(size = 3) + 
  labs(x = "", y = "trans/(cis + trans) h2") + 
  ylim(c(0,1.3)) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = 'none')

## Sun et al., 2018
sun_pqtl = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2018_Sun_supplement.xlsx', 
                     sheet = 4, startRow = 5, cols = 1:16)
h2_summ_sun$is_sig = h2_summ_sun$target %in% sun_pqtl$SOMAmer.ID

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
h2_summ_sun$signal = signal[match(h2_summ_sun$target, uniq_targe)]
h2_summ_sun$signal[which(is.na(h2_summ_sun$signal))] = 'no'

## Gudjonsson et al., 2022
gud_pqtl = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2022_Gudjonsson_pqtl.xlsx', 
                     startRow = 4, cols = 1:12)
gud_prot_info = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2022_Gudjonsson_prot_info.xlsx', 
                          startRow = 3)
gud_pqtl$target = gud_prot_info$Study.Accession[match(gud_pqtl$SOMAmer, gud_prot_info$Study.tag)]
h2_summ_Gud$is_sig = h2_summ_Gud$target %in% gud_pqtl$target

uniq_targe = unique(gud_pqtl$`Protein.(Entrez.symbol)`)
signal = c()
for (i in 1:length(uniq_targe)) {
  signal_i = unique(gud_pqtl$`cis/trans`[which(gud_pqtl$`Protein.(Entrez.symbol)` == uniq_targe[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
h2_summ_Gud$signal = signal[match(h2_summ_Gud$prot, uniq_targe)]
h2_summ_Gud$signal[which(is.na(h2_summ_Gud$signal))] = 'no'

  
## ukb, 2022
ukb_pqtl = fread('/project/xuanyao/jinghui/pqtl/04_fdr/ukb/all_prot.txt')
ukb_prot_list = fread('/project/xuanyao/jinghui/pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_pqtl$gene_chr = ukb_prot_list$chr[match(ukb_pqtl$file, ukb_prot_list$file)]
ukb_pqtl$tss = ukb_prot_list$gene_start[match(ukb_pqtl$file, ukb_prot_list$file)]
ukb_pqtl$cis_trans = 'trans'
ukb_pqtl$cis_trans[which(ukb_pqtl$CHROM == ukb_pqtl$gene_chr & 
                           ukb_pqtl$GENPOS > ukb_pqtl$tss - 1000000 & 
                           ukb_pqtl$GENPOS < ukb_pqtl$tss + 1000000)] = 'cis'

ukb_pqtl = ukb_pqtl[ukb_pqtl$LOG10P > -log10(1e-8 / 2940), ]


h2_summ_ukb$is_sig = h2_summ_ukb$file %in% ukb_pqtl$file

uniq_targe = unique(ukb_pqtl$file)
signal = c()
for (i in 1:length(uniq_targe)) {
  signal_i = unique(ukb_pqtl$cis_trans[which(ukb_pqtl$file == uniq_targe[i])])
  if (length(signal_i) > 1) {
    signal[i] = 'both'
  } else {
    signal[i] = signal_i
  }
}
h2_summ_ukb$signal = signal[match(h2_summ_ukb$file, uniq_targe)]
h2_summ_ukb$signal[which(is.na(h2_summ_ukb$signal))] = 'no'

## generate se through bootstrapping
library(boot)
boot_mean = function(formula, data, indices){
  d <- data[indices,]
  aggre_mean = aggregate(formula, data = d, mean)
  return(aggre_mean[,2])
}

h2_summ_sun$trans_prop = h2_summ_sun$trans_5mb_h2_ind / 
  (h2_summ_sun$trans_5mb_h2_ind + h2_summ_sun$cis_1mb_h2)
h2_summ_Gud$trans_prop = h2_summ_Gud$trans_5mb_h2_ind / 
  (h2_summ_Gud$trans_5mb_h2_ind + h2_summ_Gud$cis_1mb_h2)
h2_summ_ukb$trans_prop = h2_summ_ukb$trans_5mb_h2_ind / 
  (h2_summ_ukb$trans_5mb_h2_ind + h2_summ_ukb$cis_1mb_h2)

sun_cis_boot = boot(h2_summ_sun, boot_mean, 1000, formula = cis_1mb_h2~signal)
gud_cis_boot = boot(h2_summ_Gud, boot_mean, 1000, formula = cis_1mb_h2~signal)
ukb_cis_boot = boot(h2_summ_ukb, boot_mean, 1000, formula = cis_1mb_h2~signal)

sun_trans_boot = boot(h2_summ_sun, boot_mean, 1000, formula = trans_5mb_h2_ind~signal)
gud_trans_boot = boot(h2_summ_Gud, boot_mean, 1000, formula = trans_5mb_h2_ind~signal)
ukb_trans_boot = boot(h2_summ_ukb, boot_mean, 1000, formula = trans_5mb_h2_ind~signal)

sun_trans_prop_boot = boot(h2_summ_sun, boot_mean, 1000, formula = trans_prop~signal)
gud_trans_prop_boot = boot(h2_summ_Gud, boot_mean, 1000, formula = trans_prop~signal)
ukb_trans_prop_boot = boot(h2_summ_ukb, boot_mean, 1000, formula = trans_prop~signal)

h2_plot = data.frame(cis_h2 = c(apply(sun_cis_boot$t, 2, mean), apply(gud_cis_boot$t, 2, mean),
                                apply(ukb_cis_boot$t, 2, mean)),
                     cis_se = c(apply(sun_cis_boot$t, 2, sd), apply(gud_cis_boot$t, 2, sd), 
                                apply(ukb_cis_boot$t, 2, sd)),
                     trans_h2 = c(apply(sun_trans_boot$t, 2, mean), apply(gud_trans_boot$t, 2, mean),
                                  apply(ukb_trans_boot$t, 2, mean)),
                     trans_se = c(apply(sun_trans_boot$t, 2, sd), apply(gud_trans_boot$t, 2, sd),
                                  apply(ukb_trans_boot$t, 2, sd)),
                     trans_prop = c(apply(sun_trans_prop_boot$t, 2, mean), 
                                    apply(gud_trans_prop_boot$t, 2, mean),
                                    apply(ukb_trans_prop_boot$t, 2, mean)),
                     trans_prop_se = c(apply(sun_trans_prop_boot$t, 2, sd), 
                                       apply(gud_trans_prop_boot$t, 2, sd),
                                       apply(ukb_trans_prop_boot$t, 2, sd)),
                     group = rep(c('both', 'cis', 'no', 'trans'), 3), 
                     data = rep(c('Sun', 'Gudjonsson', 'UKB'), each = 4))

p1 = ggplot(h2_plot, aes(x=data, y=cis_h2, color=group)) +
  geom_errorbar(aes(ymin=cis_h2-cis_se, ymax=cis_h2+cis_se), width=.1, 
                position=position_dodge(0.5)) +
  geom_point(position=position_dodge(0.5), size = 3) + 
  labs(x = "", y = "cis h2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2 = ggplot(h2_plot, aes(x=data, y=trans_h2, color=group)) +
  geom_errorbar(aes(ymin=trans_h2-trans_se, ymax=trans_h2+trans_se), width=.1, 
                position=position_dodge(0.5)) +
  geom_point(position=position_dodge(0.5), size = 3) + 
  labs(x = "", y = "trans h2") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p3 = ggplot(h2_plot, aes(x=data, y=trans_prop, color=group)) +
  geom_errorbar(aes(ymin=trans_prop-trans_prop_se, ymax=trans_prop+trans_prop_se), width=.1, 
                position=position_dodge(0.5)) +
  geom_point(position=position_dodge(0.5), size = 3) + 
  labs(x = "", y = "trans prop") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, p3, nrow = 2)




