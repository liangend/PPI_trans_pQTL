library(data.table)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui/')

## Roadmap annotated tissues
chrom_state = fread("pqtl/07_chr_state/E062_25_imputed12marks_hg38lift_stateno.bed") # blood chromatin states
chrom_state$V3 = as.numeric(chrom_state$V3)
chrom_state$V2 = as.numeric(chrom_state$V2)
chrom_state = na.omit(chrom_state)
chrom_state$n_snp = chrom_state$V3 - chrom_state$V2 + 1
n_snp_by_state = aggregate(n_snp ~ V4, data = chrom_state, sum)
n_snp_by_state$prop = n_snp_by_state$n_snp / sum(n_snp_by_state$n_snp)
names(n_snp_by_state)[1] = "state"

## UKB cis and trans pqtl
ukb_w_state = fread("pqtl/07_chr_state/bed_file/ukb_pqtl_blood_w_state.bed")
ukb_all_prop = as.data.frame(table(ukb_w_state$V9) / nrow(ukb_w_state))
ukb_cis_prop = as.data.frame(table(ukb_w_state$V9[which(ukb_w_state$V4 == 'cis')]) / 
                               sum(ukb_w_state$V4 == 'cis'))
ukb_trans_prop = as.data.frame(table(ukb_w_state$V9[which(ukb_w_state$V4 == 'trans')]) / 
                                 sum(ukb_w_state$V4 == 'trans'))
ukb_state_prop = data.frame(state = 1:25, 
                            all = ukb_all_prop$Freq[match(1:25, ukb_all_prop$Var1)],
                            cis = ukb_cis_prop$Freq[match(1:25, ukb_cis_prop$Var1)],
                            trans = ukb_trans_prop$Freq[match(1:25, ukb_trans_prop$Var1)])

ukb_state_comp = merge(n_snp_by_state, ukb_state_prop, by = "state", all = T)
ukb_state_comp[is.na(ukb_state_comp)] = 0


## DGN trans eqtl
dgn_trans = fread('/project/xuanyao/jinghui/pqtl/07_chr_state/bed_file/dgn_trans_w_state.bed')
dgn_trans_prop = as.data.frame(table(dgn_trans$V7) / nrow(dgn_trans))
colnames(dgn_trans_prop) = c('state', 'dgn_trans')
dgn_trans_prop = rbind(dgn_trans_prop, data.frame(state = as.factor(13), dgn_trans = 0))

state_comp_all = merge(ukb_state_comp, dgn_trans_prop, by = 'state')

chromHMM = c("active TSS", "prom upstream TSS", "prom downstrem TSS 1", "prom downstrem TSS 2", 
             "transcribed 5' prefer", "strong transcribed", "transcribed 3' prefer", "weak transcribed",
             "transcribed & reg", "transcribed 5' prefer & enh", "transcribed 3' prefer & enh", 
             "transcribed & weak enh", "active enh 1", "active enh 2", "active enh flank", "weak enh 1", 
             "weak enh 2", "H3K27ac possible enh", "DNase", "ZNF genes & repeats", "hetero", "poised prom",
             "bivalent prom", "repressed polycomb", "quiescent")

all_plot = data.frame(state = rep(c(paste0(0, 1:9), 10:25), 4), 
                      enrich = c(state_comp_all$all/state_comp_all$prop, 
                                 state_comp_all$cis/state_comp_all$prop,
                                 state_comp_all$trans/state_comp_all$prop,
                                 state_comp_all$dgn_trans/state_comp_all$prop),
                      group = rep(c("pqtl", 'cis-pqtl', 'trans-pqtl', 'trans-eqtl'), each = 25), 
                      cat = rep(c(rep("prom", 4), rep("transcription", 3), 
                              rep("else", 5), rep("enh", 3), rep("else", 9), "quiescent"), 4))


ggplot(all_plot, 
       aes(x=factor(state, levels = rev(levels(factor(state)))), y=enrich, fill = group)) + 
  coord_flip() +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(title = 'UKB pQTL chromHMM (blood)', y = 'enrichment') +
  #scale_fill_manual(values=c(hue_pal()(5)[1:4], “grey”)) +
  scale_x_discrete('', labels = rev(chromHMM)) +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## trans pqtl vs trans eqtl
trans_plot = data.frame(state = chromHMM, 
                        enrich = state_comp_all$trans/state_comp_all$dgn_trans,
                        cat = c(rep("prom", 4), rep("transcription", 3), 
                                rep("else", 5), rep("enh", 3), rep("else", 9), "quiescent"))
ggplot(trans_plot, 
       aes(x= reorder(state, enrich), y = enrich)) + 
  coord_flip() +
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue') +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(x = '', title = 'trans-pQTL vs trans-eQTL', y = 'Enrichment') +
  #scale_fill_manual(values=c(hue_pal()(5)[1:4], “grey”)) +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## trans pqtl less eqtl vs trans pqtl more eqtl
ukb_less_eqtl_prop = as.data.frame(table(ukb_w_state$V9[which(ukb_w_state$V4 == 'trans' & 
                                                                ukb_w_state$V5 >= 0.05)]) / 
                               sum(ukb_w_state$V4 == 'trans' & ukb_w_state$V5 >= 0.05))
ukb_more_eqtl_prop = as.data.frame(table(ukb_w_state$V9[which(ukb_w_state$V4 == 'trans' & 
                                                                ukb_w_state$V5 < 0.05)]) / 
                                     sum(ukb_w_state$V4 == 'trans' & ukb_w_state$V5 < 0.05))
ukb_trans_prop = data.frame(state = 1:25, 
                            less_eqtl = ukb_less_eqtl_prop$Freq[match(1:25, ukb_less_eqtl_prop$Var1)],
                            more_eqtl = ukb_more_eqtl_prop$Freq[match(1:25, ukb_more_eqtl_prop$Var1)])

ukb_trans_state_comp = merge(n_snp_by_state, ukb_trans_prop, by = "state", all = T)
ukb_trans_state_comp[is.na(ukb_trans_state_comp)] = 0

ukb_trans_plot = data.frame(state = rep(c(paste0(0, 1:9), 10:25), 2), 
                      enrich = c(ukb_trans_state_comp$less_eqtl/ukb_state_comp$prop,
                                 ukb_trans_state_comp$more_eqtl/ukb_state_comp$prop),
                      group = rep(c('dgn p >= 0.05', 'dgn p < 0.05'), each = 25), 
                      cat = rep(c(rep("prom", 4), rep("transcription", 3), 
                                  rep("else", 5), rep("enh", 3), rep("else", 9), "quiescent"), 2))

ggplot(ukb_trans_plot, 
       aes(x=factor(state, levels = rev(levels(factor(state)))), y=enrich, fill = group)) + 
  coord_flip() +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(title = 'UKB trans-pQTL chromHMM (blood)', y = 'enrichment') +
  #scale_fill_manual(values=c(hue_pal()(5)[1:4], “grey”)) +
  scale_x_discrete('', labels = rev(chromHMM)) +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)



