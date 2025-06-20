library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
## strongest cis-pQTL vs trans-pQTL, using UKB as discovery set and Sun et al., 2018 as replication set
setwd('/project/xuanyao/jinghui/')
sun_beta_ukb_extr = fread('pqtl/12_beta_across_two_data/sun_beta_ukb_extreme.txt')
sun_beta_ukb_extr = na.omit(sun_beta_ukb_extr)
sun_beta_ukb_extr$cis_z = abs(sun_beta_ukb_extr$cis_beta/sun_beta_ukb_extr$cis_se)
sun_beta_ukb_extr$trans_z = abs(sun_beta_ukb_extr$trans_beta/sun_beta_ukb_extr$trans_se)

sun_z_ukb_extr = data.frame(cis_trans = rep(c('cis', 'trans'), each = nrow(sun_beta_ukb_extr)),
                            abs_z = c(sort(sun_beta_ukb_extr$cis_z, decreasing = T),
                                      sort(sun_beta_ukb_extr$trans_z, decreasing = T)))
sun_z_ukb_extr$rank = rep(1:nrow(sun_beta_ukb_extr), 2)
sun_z_ukb_extr$rank_percent = sun_z_ukb_extr$rank / 
  max(sun_z_ukb_extr$rank) * 100

## strongest cis-eQTL vs trans-eQTL, using Wright et al., 2014 as discovery set and DGN as replication set
cis_all=NULL
trans_all=NULL
for(i in c(1:22)){
  dat1=read.table(paste("pqtl/00_ref/xuanyao_ominigenic/chr",i,"_cis_topallNTR_sign_p.txt",sep=""))
  dat2=read.table(paste("pqtl/00_ref/xuanyao_ominigenic/chr",i,"_trans_topallNTR_sign_p.txt",sep=""))
  cis_all=c(cis_all,abs(as.numeric(dat1[,1])))
  trans_all=c(trans_all,abs(as.numeric(dat2[,1])))
}

logp_cis=-log(cis_all)
logp_trans=-log(trans_all)

temp_cis=abs(qnorm(-logp_cis-log(2),log=T))
temp_trans=abs(qnorm(-logp_trans-log(2),log=T))

mrna_z_extr = data.frame(cis_trans = rep(c('cis', 'trans'), each = length(temp_cis)),
                            abs_z = c(sort(temp_cis, decreasing = T), 
                                      sort(temp_trans, decreasing = T)))

mrna_z_extr$rank = rep(1:length(temp_cis), 2)
mrna_z_extr$rank_percent = mrna_z_extr$rank / 
  max(mrna_z_extr$rank) * 100

## plot protein and mRNA together
mrna_z_extr$mole = 'mRNA'
sun_z_ukb_extr$mole = 'protein'
z_all = rbind(mrna_z_extr[c(seq(1, 50, 12), seq(51, nrow(mrna_z_extr)/2, 200),
                            seq(nrow(mrna_z_extr)/2+1, nrow(mrna_z_extr)/2+50, 12), 
                            seq(nrow(mrna_z_extr)/2+51, nrow(mrna_z_extr), 200)),], 
              sun_z_ukb_extr[c(seq(1, 50, 5), seq(51, nrow(sun_z_ukb_extr)/2, 20),
                               seq(nrow(sun_z_ukb_extr)/2+1, nrow(sun_z_ukb_extr)/2+50, 5), 
                               seq(nrow(sun_z_ukb_extr)/2+51, nrow(sun_z_ukb_extr), 20)),])
ggplot() +
  geom_point(data = z_all, aes(x = rank_percent, y = abs_z, shape = cis_trans, color = mole), 
             size = 2) +
  labs(x = "QTL rank percentile", y = '|Z|', 
       title = 'Strongest cis- and trans-QTLs effect size', 
       color = "", shape = '') +
  scale_color_manual(values = brewer.pal(5,"Set1")[c(5)]) +
  scale_shape_manual(values = c(16,4)) +
  ylim(c(0,100)) + 
  theme(text = element_text(size=14, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### plot trans/cis ratio
ratio_rank_percent = data.frame(ratio = c(sort(sun_beta_ukb_extr$trans_z, decreasing = T) / 
                                             sort(sun_beta_ukb_extr$cis_z, decreasing = T),
                                           sort(temp_trans, decreasing = T) / 
                                             sort(temp_cis, decreasing = T)),
                                rank_percent = c(sun_z_ukb_extr$rank_percent[1:nrow(sun_beta_ukb_extr)],
                                                 mrna_z_extr$rank_percent[1:length(temp_cis)]),
                                type = c(rep('protein', nrow(sun_beta_ukb_extr)),
                                         rep('mRNA', length(temp_cis))))
ggplot() +
  geom_point(data = ratio_rank_percent, aes(x = rank_percent, y = ratio, color = type), size=1) +
  labs(x = "QTL rank percentile", y = '|trans Z / cis Z|', title = '', color = "") +
  scale_color_manual(values = brewer.pal(5,"Set1")[c(2,5)]) +
  theme(text = element_text(size=15, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 13, angle = 0, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')


