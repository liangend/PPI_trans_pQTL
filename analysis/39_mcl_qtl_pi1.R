library(data.table)
library(ggplot2)
library(qvalue)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
mcl_univar_comp = fread('pqtl/04_fdr/ukb/pco_mcl_meta.txt')

univar_z = sapply(strsplit(mcl_univar_comp$univar_z, ',', fixed = T), as.numeric)
univar_p = sapply(univar_z, function(x){
  pnorm(abs(x), lower.tail = F) * 2
})
pi1_all = sapply(univar_p, function(x){
  qval_i = try(qvalue(p = x), silent = T)  # pi0 estimate may not be available for small
  if (length(qval_i) == 8) {
    return(1 - qval_i$pi0)
  } else {
    return(NA)
  }
})
mcl_univar_comp$pi1 = pi1_all
mcl_univar_comp$n_prot = sapply(univar_z, length)
mcl_univar_comp$n_univar_p_5e2 = sapply(univar_p, function(x){sum(x < 0.05)})
mcl_univar_comp$n_univar_p_1e2 = sapply(univar_p, function(x){sum(x < 0.01)})
mcl_univar_comp$n_univar_p_1e3 = sapply(univar_p, function(x){sum(x < 1e-3)})
mcl_univar_comp$n_univar_p_1e5 = sapply(univar_p, function(x){sum(x < 1e-5)})
mcl_univar_comp$n_univar_p_1e7 = sapply(univar_p, function(x){sum(x < 1e-7)})

mcl_univar_comp$prop_univar_p_5e2 = mcl_univar_comp$n_univar_p_5e2/mcl_univar_comp$n_prot
mcl_univar_comp$prop_univar_p_1e2 = mcl_univar_comp$n_univar_p_1e2/mcl_univar_comp$n_prot
mcl_univar_comp$prop_univar_p_1e3 = mcl_univar_comp$n_univar_p_1e3/mcl_univar_comp$n_prot
mcl_univar_comp$prop_univar_p_1e5 = mcl_univar_comp$n_univar_p_1e5/mcl_univar_comp$n_prot

mcl_univar_w_pi = na.omit(mcl_univar_comp)
p1 = ggplot(mcl_univar_w_pi, aes(x=pi1)) + 
  geom_histogram(fill = 'steelblue') +
  labs(x = "pi1", title = 'pi1 distribution (all)') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p2 = ggplot(mcl_univar_w_pi[mcl_univar_w_pi$n_prot > 3, ], aes(x=pi1)) + 
  geom_histogram(fill = 'steelblue') +
  labs(x = "pi1", title = 'pi1 distribution (# prot > 3)') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)

p3=ggplot(mcl_univar_comp, aes(x=n_univar_p_1e3/n_prot)) + 
  geom_histogram(fill = 'steelblue') +
  labs(x = "# univar p < 1e-3", title = 'prop univar p < 1e-3 distribution (all)') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p4=ggplot(mcl_univar_comp[mcl_univar_comp$n_prot > 3, ], aes(x=n_univar_p_1e3/n_prot)) + 
  geom_histogram(fill = 'steelblue') +
  labs(x = "# univar p < 1e-3", title = 'prop univar p < 1e-3 distribution (# prot > 3)') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p3, p4, nrow = 1)


mcl_univar_w_pi_nov = mcl_univar_w_pi[mcl_univar_w_pi$is_nov, ]
mcl_univar_w_pi_share = mcl_univar_w_pi[!mcl_univar_w_pi$is_nov, ]

mcl_univar_nov = mcl_univar_comp[mcl_univar_comp$is_nov, ]
mcl_univar_share = mcl_univar_comp[!mcl_univar_comp$is_nov, ]

pi_tab = data.frame(prot = c(sum(mcl_univar_w_pi$pi1 > 1/mcl_univar_w_pi$n_prot)/nrow(mcl_univar_w_pi),
                             sum(mcl_univar_w_pi_nov$pi1 > 1/mcl_univar_w_pi_nov$n_prot)/nrow(mcl_univar_w_pi_nov),
                             sum(mcl_univar_w_pi_share$pi1 > 1/mcl_univar_w_pi_share$n_prot)/nrow(mcl_univar_w_pi_share),
                             sum(mcl_univar_comp$n_univar_p_1e5 > 1) / nrow(mcl_univar_comp),
                             sum(mcl_univar_nov$n_univar_p_1e5 > 1) / nrow(mcl_univar_nov),
                             sum(mcl_univar_share$n_univar_p_1e5 > 1) / nrow(mcl_univar_share),
                             sum(mcl_univar_comp$n_univar_p_1e7 > 1) / nrow(mcl_univar_comp),
                             sum(mcl_univar_nov$n_univar_p_1e7 > 1) / nrow(mcl_univar_nov),
                             sum(mcl_univar_share$n_univar_p_1e7 > 1) / nrow(mcl_univar_share),
                             sum(mcl_univar_comp$n_univar_p_5e2 > 1) / nrow(mcl_univar_comp),
                             sum(mcl_univar_nov$n_univar_p_5e2 > 1) / nrow(mcl_univar_nov),
                             sum(mcl_univar_share$n_univar_p_5e2 > 1) / nrow(mcl_univar_share),
                             sum(mcl_univar_comp$n_univar_p_1e2 > 1) / nrow(mcl_univar_comp),
                             sum(mcl_univar_nov$n_univar_p_1e2 > 1) / nrow(mcl_univar_nov),
                             sum(mcl_univar_share$n_univar_p_1e2 > 1) / nrow(mcl_univar_share),
                             sum(mcl_univar_comp$n_univar_p_1e3 > 1) / nrow(mcl_univar_comp),
                             sum(mcl_univar_nov$n_univar_p_1e3 > 1) / nrow(mcl_univar_nov),
                             sum(mcl_univar_share$n_univar_p_1e3 > 1) / nrow(mcl_univar_share)),
                    group = rep(c('pi1', 'p < 1e-5', 'p < 1e-7', 'p < 5e-2',
                                  'p < 1e-2', 'p < 1e-3'), each = 3),
                    cat = rep(c('all', 'novel', 'shared'), 6))
ggplot(pi_tab, aes(x = factor(group, levels = c('p < 5e-2', 'p < 1e-2', 'p < 1e-3', 'p < 1e-5', 'p < 1e-7',
                                                'pi1')), y = prot, fill = cat, label = round(prot, 2))) + 
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.6) +
  geom_text(position=position_dodge(0.6), vjust=-0.2, size = 4.5) + 
  labs(x = "", fill = '', title = 'PCO pQTL pleiotropy') +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())







