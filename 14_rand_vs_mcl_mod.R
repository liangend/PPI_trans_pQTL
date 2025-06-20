library(data.table)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
### random background
pco_rand = fread('pqtl/04_fdr/ukb/ukb_discovery/pco_rand_univar_comp_all.txt')
rand_mod = fread('pqtl/11_prot_complex/random_mod.txt')
rand_mod$mod = paste0('mod', rand_mod$mod)
pco_rand$is_nov = (pco_rand$univar_max_logP < -log10(5e-8/2923))

### pco signals
pco_ppi = fread('pqtl/04_fdr/ukb/ukb_discovery/pco_mcl_meta.txt')
ppi_mod = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')
ppi_mod$mod = paste0('mod', ppi_mod$mod)
ppi_mod$mod_ori = sapply(strsplit(ppi_mod$mod_ori, '_', fixed = T), '[', 1)

total_nov_ppi = list()
total_nov_rand = list()
for (i in 2:10) {
  mod_size_i = names(which(table(ppi_mod$mod) == i))
  
  #### number of novel loci in total
  pco_ppi_i = pco_ppi[pco_ppi$mod %in% mod_size_i, ]
  n_nov_ppi_size_i = aggregate(is_nov ~ mod, data = pco_ppi_i, sum)
  
  pco_rand_i = pco_rand[pco_rand$size == i, ]
  n_nov_rand_size_i = aggregate(is_nov ~ mod, data = pco_rand_i, sum)
  
  total_nov_ppi[[i-1]] = c(n_nov_ppi_size_i$is_nov,
                           rep(0, length(mod_size_i) - nrow(n_nov_ppi_size_i))) # 0 novel loci for modules with no pQTL
  total_nov_rand[[i-1]] = c(n_nov_rand_size_i$is_nov, 
                             rep(0, 100 - nrow(n_nov_rand_size_i)))
  
  print(i)
}
n_nov_comp = data.frame(n_nov = c(unlist(total_nov_ppi), unlist(total_nov_rand)),
                        mod_size = c(rep(2:10, times = sapply(total_nov_ppi, length)),
                                     rep(2:10, each = 100)),
                        group = c(rep('PPI', length(unlist(total_nov_ppi))),
                                  rep('random', length(unlist(total_nov_rand)))))
n_nov_p = c()
for (i in 2:10) {
  n_nov_sub = n_nov_comp[n_nov_comp$mod_size == i, ]
  fit_i = wilcox.test(n_nov ~ group, data = n_nov_sub)
  n_nov_p = c(n_nov_p, fit_i$p.value)
}


n_nov_boot_plot = data.frame(n_nov = c(sapply(total_nov_ppi, mean), 
                                       sapply(total_nov_rand, mean)),
                             se = c(sapply(total_nov_ppi, function(x){sd(x)/sqrt(length(x))}), 
                                    sapply(total_nov_rand, function(x){sd(x)/sqrt(length(x))})),
                             mod_size = rep(as.factor(c(2:10)), 2),
                             group = rep(c('PPI', 'random'), each = 9))

ggplot(n_nov_boot_plot, aes(x = mod_size, y = n_nov, color = group)) + 
  geom_point(stat="identity", position=position_dodge(0.6), size = 1.5) +
  geom_errorbar(position=position_dodge(0.6), aes(ymin=n_nov-se, ymax=n_nov+se), width=.05) +
  labs(x = "Cluster size", y = 'Number novel loci per cluster', title = '', color = 'Protein cluster') + 
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


### correlation comparison
pco_cor = list()
rand_cor = list()
for (i in 2:10) {
  ppi_mod_i = names(which(table(ppi_mod$mod) == i))
  rand_mod_i = names(which(table(rand_mod$mod) == i))
  
  pco_cor_i = c()
  for (j in 1:length(ppi_mod_i)) {
    cor_j = readRDS(paste0('pqtl/02_sigma/ukb/', ppi_mod_i[j], '.rds'))
    pco_cor_i[j] = mean(abs(cor_j[lower.tri(cor_j)]))
  }
  
  rand_cor_i = c()
  for (j in 1:length(rand_mod_i)) {
    cor_j = readRDS(paste0('pqtl/02_sigma/ukb_random/', rand_mod_i[j], '.rds'))
    rand_cor_i[j] = mean(abs(cor_j[lower.tri(cor_j)]))
  }
  
  pco_cor[[i-1]] = pco_cor_i
  rand_cor[[i-1]] = rand_cor_i
  
  print(i)
}

cor_plot = data.frame(cor = c(sapply(pco_cor, mean), 
                              sapply(rand_cor, mean)),
                      se = c(sapply(pco_cor, function(x){sd(x)/sqrt(length(x))}), 
                             sapply(rand_cor, function(x){sd(x)/sqrt(length(x))})),
                      mod_size = rep(as.factor(c(2:10)), 2),
                      group = rep(c('PPI', 'random'), each = 9))

cor_comp = data.frame(cor = c(unlist(pco_cor), unlist(rand_cor)),
                      mod_size = c(rep(2:10, times = sapply(pco_cor, length)),
                                   rep(2:10, times = sapply(rand_cor, length))),
                      group = c(rep('PPI', length(unlist(pco_cor))),
                                rep('random', length(unlist(rand_cor)))))
wilcox.test(cor ~ group, data = cor_comp)
aggregate(cor ~ group, data = cor_comp, mean)

ggplot(cor_plot, aes(x = mod_size, y = cor, color = group)) + 
  geom_point(stat="identity", position=position_dodge(0.6), size = 1.5) +
  geom_errorbar(position=position_dodge(0.6), aes(ymin=cor-se, ymax=cor+se), width=.05) +
  labs(x = "Cluster size", y = 'Cluster correlation', title = '', color = 'Protein cluster') + 
  theme(text = element_text(size=13, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')










