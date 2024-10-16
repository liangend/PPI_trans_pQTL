library(data.table)
library(openxlsx)
library(ggplot2)
library(gridExtra)
## 2018 Sun
sun_pqtl = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2018_Sun_supplement.xlsx', 
                     sheet = 4, startRow = 5, cols = 1:16)
sun_pqtl$gene = sun_prot_list$gene[match(sun_pqtl$SOMAmer.ID, sun_prot_list$target)]
sun_pqtl$MAF = pmin(sun_pqtl$EAF, 1 - sun_pqtl$EAF)

sun_pqtl_eff = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2018_Sun_supplement.xlsx', 
                         sheet = 4, startRow = 6, cols = 23:25)
sun_pqtl = cbind(sun_pqtl, sun_pqtl_eff)

# make a pqtl bed
sun_pqtl_bed = data.frame(chr = paste0('chr', sun_pqtl$Chr), 
                          start = sun_pqtl$Pos-1, end = sun_pqtl$Pos,
                          type = sun_pqtl$`cis/.trans`)
fwrite(sun_pqtl_bed, '/project/xuanyao/jinghui/pqtl/07_chr_state/sun_pqtl.bed', col.names = F,
       sep = '\t')

# maf difference between cis and trans
p1 = ggplot(sun_pqtl, aes(x = `cis/.trans`, y = MAF, fill = `cis/.trans`)) + 
  geom_boxplot(width=0.5) +
  labs(x = "", y = "MAF", title = 'Sun et al., 2018') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
# maf vs |effect size|
p2 = ggplot(sun_pqtl, aes(x = MAF, y = abs(beta), color = `cis/.trans`)) + 
  geom_point(size = 1.5) +
  labs(x = "MAF", y = "|effect size|", title = 'Sun et al., 2018') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# if trans signal overlap with cis signal
sun_pqtl_cis = sun_pqtl[sun_pqtl$`cis/.trans` == 'cis', ]
sun_pqtl_trans = sun_pqtl[sun_pqtl$`cis/.trans` == 'trans', ]
is_cis = c()
for (i in 1:nrow(sun_pqtl_trans)) {
  chr_i = sun_pqtl_trans$Chr[i]
  start_i = sun_pqtl_trans$Region.start[i]
  end_i = sun_pqtl_trans$Region.end[i]
  cis_i = which(sun_pqtl_cis$Chr == chr_i & (
    (sun_pqtl_cis$Region.start <= start_i &  sun_pqtl_cis$Region.end >= start_i) |
      (sun_pqtl_cis$Region.start <= end_i &  sun_pqtl_cis$Region.end >= end_i) |
      (sun_pqtl_cis$Region.start >= start_i &  sun_pqtl_cis$Region.end <= end_i) |
      (sun_pqtl_cis$Region.start <= start_i &  sun_pqtl_cis$Region.end >= end_i)
  ))
  is_cis[i] = length(cis_i) > 0
}
sun_pqtl_trans$is_cis = is_cis

ggplot(sun_pqtl_trans, aes(x = is_cis, y = MAF, fill = is_cis)) + 
  geom_boxplot(width=0.5) +
  labs(x = "trans overlap with cis", y = "MAF", title = 'Sun et al., 2018') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

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

# make a pqtl bed
gud_pqtl_bed = data.frame(chr = gud_pqtl_merge$Chr, 
                          start = gud_pqtl_merge$`Pos.(GRCh37)`-1, 
                          end = gud_pqtl_merge$`Pos.(GRCh37)`,
                          type = gud_pqtl_merge$`cis/trans`)
fwrite(gud_pqtl_bed, '/project/xuanyao/jinghui/pqtl/07_chr_state/gud_pqtl.bed', col.names = F,
       sep = '\t')

# maf difference between cis and trans
p3 = ggplot(gud_pqtl_merge, aes(x = `cis/trans`, y = EAF, fill = `cis/trans`)) + 
  geom_boxplot(width=0.5) +
  labs(x = "", y = "MAF", title = 'Gudjonsson et al., 2022') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')

# maf vs |effect size|
p4 = ggplot(gud_pqtl_merge, aes(x = EAF, y = abs(beta), color = `cis/trans`)) + 
  geom_point(size = 1.5) +
  labs(x = "MAF", y = "|effect size|", title = 'Gudjonsson et al., 2022') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p3, p2, p4, nrow = 2)
table(gud_pqtl_merge$`cis/trans`)
















