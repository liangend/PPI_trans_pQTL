library(data.table)
library(ggplot2)
setwd('/project/xuanyao/jinghui/pqtl/')

## Roadmap annotated tissues
chrom_state = fread("07_chr_state/E066_25_imputed12marks_stateno.bed") # liver chromatin states
chrom_state$V3 = as.numeric(chrom_state$V3)
chrom_state$V2 = as.numeric(chrom_state$V2)
chrom_state = na.omit(chrom_state)
chrom_state$n_snp = chrom_state$V3 - chrom_state$V2 + 1
n_snp_by_state = aggregate(n_snp ~ V4, data = chrom_state, sum)
n_snp_by_state$prop = n_snp_by_state$n_snp / sum(n_snp_by_state$n_snp)
names(n_snp_by_state)[1] = "state"

## Sun et al., 2018
sun_w_state = fread("07_chr_state/sun_w_liver_state.bed")
sun_all_prop = as.data.frame(table(sun_w_state$V8) / nrow(sun_w_state))
sun_cis_prop = as.data.frame(table(sun_w_state$V8[which(sun_w_state$V4 == 'cis')]) / 
                               sum(sun_w_state$V4 == 'cis'))
sun_trans_prop = as.data.frame(table(sun_w_state$V8[which(sun_w_state$V4 == 'trans')]) / 
                                 sum(sun_w_state$V4 == 'trans'))
sun_state_prop = data.frame(state = 1:25, 
                            all = sun_all_prop$Freq[match(1:25, sun_all_prop$Var1)],
                            cis = sun_cis_prop$Freq[match(1:25, sun_cis_prop$Var1)],
                            trans = sun_trans_prop$Freq[match(1:25, sun_trans_prop$Var1)])

sun_state_comp = merge(n_snp_by_state, sun_state_prop, by = "state", all = T)
sun_state_comp[is.na(sun_state_comp)] = 0

chromHMM = c("active TSS", "prom upstream TSS", "prom downstrem TSS 1", "prom downstrem TSS 2", 
             "transcribed 5' prefer", "strong transcribed", "transcribed 3' prefer", "weak transcribed",
             "transcribed & reg", "transcribed 5' prefer & enh", "transcribed 3' prefer & enh", 
             "transcribed & weak enh", "active enh 1", "active enh 2", "active enh flank", "weak enh 1", 
             "weak enh 2", "H3K27ac possible enh", "DNase", "ZNF genes & repeats", "hetero", "poised prom",
             "bivalent prom", "repressed polycomb", "quiescent")
sun_plot = data.frame(state = rep(c(paste0(0, 1:9), 10:25), 3), 
                      enrich = c(sun_state_comp$all/sun_state_comp$prop, 
                                 sun_state_comp$cis/sun_state_comp$prop,
                                 sun_state_comp$trans/sun_state_comp$prop),
                      group = rep(c("all", 'cis', 'trans'), each = 25), 
                      cat = rep(c(rep("prom", 4), rep("transcription", 3), 
                              rep("else", 5), rep("enh", 3), rep("else", 9), "quiescent"), 3))


# dat_plot_sub = dat_plot[dat_plot$cat != "else" & dat_plot$cat != "quiescent", ]
# library(scales)

p1 = ggplot(sun_plot, aes(x=factor(state, levels = rev(levels(factor(state)))), y=enrich, fill = group)) + 
  coord_flip() +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(title = 'Sun et al., 2018', y = 'enrichment') +
  #scale_fill_manual(values=c(hue_pal()(5)[1:4], “grey”)) +
  scale_x_discrete('', labels = rev(chromHMM)) +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

## Gudjonsson et al., 2022
gud_w_state = fread("07_chr_state/gud_w_liver_state.bed")
gud_all_prop = as.data.frame(table(gud_w_state$V8) / nrow(gud_w_state))
gud_cis_prop = as.data.frame(table(gud_w_state$V8[which(gud_w_state$V4 == 'cis')]) / 
                               sum(gud_w_state$V4 == 'cis'))
gud_trans_prop = as.data.frame(table(gud_w_state$V8[which(gud_w_state$V4 == 'trans')]) / 
                                 sum(gud_w_state$V4 == 'trans'))
gud_state_prop = data.frame(state = 1:25, 
                            all = gud_all_prop$Freq[match(1:25, gud_all_prop$Var1)],
                            cis = gud_cis_prop$Freq[match(1:25, gud_cis_prop$Var1)],
                            trans = gud_trans_prop$Freq[match(1:25, gud_trans_prop$Var1)])

gud_state_comp = merge(n_snp_by_state, gud_state_prop, by = "state", all = T)
gud_state_comp[is.na(gud_state_comp)] = 0

gud_plot = data.frame(state = rep(c(paste0(0, 1:9), 10:25), 3), 
                      enrich = c(gud_state_comp$all/gud_state_comp$prop, 
                                 gud_state_comp$cis/gud_state_comp$prop,
                                 gud_state_comp$trans/gud_state_comp$prop),
                      group = rep(c("all", 'cis', 'trans'), each = 25), 
                      cat = rep(c(rep("prom", 4), rep("transcription", 3), 
                                  rep("else", 5), rep("enh", 3), rep("else", 9), "quiescent"), 3))
p2 = ggplot(gud_plot, aes(x=factor(state, levels = rev(levels(factor(state)))), y=enrich, fill = group)) + 
  coord_flip() +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs(title = 'Gudjonsson et al., 2022', y = 'enrichment') +
  #scale_fill_manual(values=c(hue_pal()(5)[1:4], “grey”)) +
  scale_x_discrete('', labels = rep(' ', 25)) +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, nrow = 1)




