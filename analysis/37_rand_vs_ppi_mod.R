library(data.table)
library(ggplot2)
library(gridExtra)
setwd('/project/xuanyao/jinghui')
### random background
pco_rand1 = fread('pqtl/04_fdr/ukb/pco_rand_univar_comp1.txt')
colnames(pco_rand1) = c('chr', 'bp_hg37', 'mod', 'pco_p', 'loci', 
                        'bp_start', 'bp_end', 'univar_max_logP')

pco_rand2 = fread('pqtl/04_fdr/ukb/pco_rand_univar_comp2.txt')
colnames(pco_rand2) = c('chr', 'bp_hg37', 'mod', 'pco_p', 'loci', 
                        'bp_start', 'bp_end', 'univar_max_logP')

pco_rand2$mod = as.numeric(sub('mod', '', pco_rand2$mod))
pco_rand2$mod = pco_rand2$mod + 150
pco_rand2$mod = paste0('mod', pco_rand2$mod)

pco_rand = rbind(pco_rand1, pco_rand2)
pco_rand$is_nov = (pco_rand$univar_max_logP < -log10(5e-8/2923))

### pco signals
pco_ppi = fread('pqtl/04_fdr/ukb/pco_nonoverlap_meta_share.txt')
pco_ppi$mod_ori = sapply(strsplit(pco_ppi$mod_ori, '_', fixed = T), '[', 1)

ppi_mod = fread('pqtl/11_prot_complex/4_database_mod.txt')
ppi_mod$mod = paste0('mod', ppi_mod$mod)
mod_w_p = list.files('pqtl/03_p/ukb_4_database')
mod_w_p = sapply(strsplit(mod_w_p, '.', fixed = T), '[', 2)
mod_w_p = unique(mod_w_p)
ppi_mod = ppi_mod[ppi_mod$mod %in% mod_w_p, ]

### modules with 3 proteins
mod_size3 = names(which(table(ppi_mod$mod) == 3))
ppi_mod_size3 = ppi_mod[ppi_mod$mod %in% mod_size3, ]
# non-overlapping modules with 3 proteins
non_overlap_ppi_mod_size3 = ppi_mod_size3[ppi_mod_size3$mod == mod_size3[1], ]
for (i in 2:length(mod_size3)) {
  prot_i = ppi_mod_size3$prot[ppi_mod_size3$mod == mod_size3[i]]
  if (sum(prot_i %in% non_overlap_ppi_mod_size3$prot) == 0) {
    non_overlap_ppi_mod_size3 = rbind(non_overlap_ppi_mod_size3, 
                                      ppi_mod_size3[ppi_mod_size3$mod == mod_size3[i], ])
  }
}
non_overlap_mod_size3 = unique(non_overlap_ppi_mod_size3$mod)

# number of novel loci for overlapping vs random modules
pco_ppi_size3 = pco_ppi[pco_ppi$mod %in% mod_size3, ]
n_nov_ppi_size3 = aggregate(is_nov ~ mod, data = pco_ppi_size3, sum)

pco_rand_size3 = pco_rand[pco_rand$mod %in% paste0('mod', c(1:50, 151:200)), ]
n_nov_rand_size3 = aggregate(is_nov ~ mod, data = pco_rand_size3, sum)

n_nov_size3_tab = data.frame(mod = mod_size3, 
                             mod_ori = ppi_mod_size3$mod_ori[match(mod_size3, ppi_mod_size3$mod)],
                             n_nov = n_nov_ppi_size3$is_nov[match(mod_size3, n_nov_ppi_size3$mod)])
n_nov_rand_size3_tab = data.frame(mod = paste0('mod', c(1:50, 151:200)),
                                  mod_ori = rep('random', 100),
                                  n_nov = n_nov_rand_size3$is_nov[match(paste0('mod', c(1:50, 151:200)),
                                                                        n_nov_rand_size3$mod)])
n_nov_size3_tab = rbind(n_nov_size3_tab, n_nov_rand_size3_tab)
n_nov_size3_tab[is.na(n_nov_size3_tab)] = 0
n_nov_size3_tab$group = 'PPI'
n_nov_size3_tab$group[which(n_nov_size3_tab$mod_ori == 'random')] = 'random'

wilcox.test(n_nov_size3_tab$n_nov[n_nov_size3_tab$mod_ori == 'random'], 
            n_nov_size3_tab$n_nov[n_nov_size3_tab$mod_ori != 'random'])

p1=ggplot(n_nov_size3_tab, aes(x=n_nov, fill=group)) +
  geom_density(alpha=0.4) + 
  labs(title = "Module size = 3, Mann-Whitney p = 0.21", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p2=ggplot(n_nov_size3_tab, aes(x=n_nov, color=mod_ori)) +
  geom_density() + 
  labs(title = "Module size = 3", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# number of novel loci for non-overlapping vs random modules
pco_ppi_size3 = pco_ppi[pco_ppi$mod %in% non_overlap_mod_size3, ]
n_nov_ppi_size3 = aggregate(is_nov ~ mod, data = pco_ppi_size3, sum)

pco_rand_size3 = pco_rand[pco_rand$mod %in% paste0('mod', c(1:50, 151:200)), ]
n_nov_rand_size3 = aggregate(is_nov ~ mod, data = pco_rand_size3, sum)

n_nov_size3_tab = data.frame(mod = non_overlap_mod_size3, 
                             mod_ori = ppi_mod_size3$mod_ori[match(non_overlap_mod_size3, ppi_mod_size3$mod)],
                             n_nov = n_nov_ppi_size3$is_nov[match(non_overlap_mod_size3, n_nov_ppi_size3$mod)])

n_nov_size3_tab = rbind(n_nov_size3_tab, n_nov_rand_size3_tab)
n_nov_size3_tab[is.na(n_nov_size3_tab)] = 0
n_nov_size3_tab$group = 'PPI'
n_nov_size3_tab$group[which(n_nov_size3_tab$mod_ori == 'random')] = 'random'

wilcox.test(n_nov_size3_tab$n_nov[n_nov_size3_tab$mod_ori == 'random'], 
            n_nov_size3_tab$n_nov[n_nov_size3_tab$mod_ori != 'random'])

p3=ggplot(n_nov_size3_tab, aes(x=n_nov, fill=group)) +
  geom_density(alpha=0.4) + 
  labs(title = "Module size = 3, Mann-Whitney p = 0.32, non-overlapping", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p4=ggplot(n_nov_size3_tab, aes(x=n_nov, color=mod_ori)) +
  geom_density() + 
  labs(title = "Module size = 3, non-overlapping", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

### modules with 10 proteins
mod_size10 = names(which(table(ppi_mod$mod) == 22))
ppi_mod_size10 = ppi_mod[ppi_mod$mod %in% mod_size10, ]
# non-overlapping modules with 10 proteins
non_overlap_ppi_mod_size10 = ppi_mod_size10[ppi_mod_size10$mod == mod_size10[1], ]
for (i in 2:length(mod_size10)) {
  prot_i = ppi_mod_size10$prot[ppi_mod_size10$mod == mod_size10[i]]
  if (sum(prot_i %in% non_overlap_ppi_mod_size10$prot) == 0) {
    non_overlap_ppi_mod_size10 = rbind(non_overlap_ppi_mod_size10, 
                                      ppi_mod_size10[ppi_mod_size10$mod == mod_size10[i], ])
  }
}
non_overlap_mod_size10 = unique(non_overlap_ppi_mod_size10$mod)

# number of novel loci for overlapping vs random modules
pco_ppi_size10 = pco_ppi[pco_ppi$mod %in% mod_size10, ]
n_nov_ppi_size10 = aggregate(is_nov ~ mod, data = pco_ppi_size10, sum)

pco_rand_size10 = pco_rand[pco_rand$mod %in% paste0('mod', c(51:100, 201:250)), ]
n_nov_rand_size10 = aggregate(is_nov ~ mod, data = pco_rand_size10, sum)

n_nov_size10_tab = data.frame(mod = mod_size10, 
                             mod_ori = ppi_mod_size10$mod_ori[match(mod_size10, ppi_mod_size10$mod)],
                             n_nov = n_nov_ppi_size10$is_nov[match(mod_size10, n_nov_ppi_size10$mod)])
n_nov_rand_size10_tab = data.frame(mod = paste0('mod', c(51:100, 201:250)),
                                  mod_ori = rep('random', 100),
                                  n_nov = n_nov_rand_size10$is_nov[match(paste0('mod', c(51:100, 201:250)),
                                                                        n_nov_rand_size10$mod)])
n_nov_size10_tab = rbind(n_nov_size10_tab, n_nov_rand_size10_tab)
n_nov_size10_tab[is.na(n_nov_size10_tab)] = 0
n_nov_size10_tab$group = 'PPI'
n_nov_size10_tab$group[which(n_nov_size10_tab$mod_ori == 'random')] = 'random'

wilcox.test(n_nov_size10_tab$n_nov[n_nov_size10_tab$mod_ori == 'random'], 
            n_nov_size10_tab$n_nov[n_nov_size10_tab$mod_ori != 'random'])

p5=ggplot(n_nov_size10_tab, aes(x=n_nov, fill=group)) +
  geom_density(alpha=0.4) + 
  labs(title = "Module size = 10, Mann-Whitney p = 0.02", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p6=ggplot(n_nov_size10_tab, aes(x=n_nov, color=mod_ori)) +
  geom_density() + 
  labs(title = "Module size = 10", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# number of novel loci for non-overlapping vs random modules
pco_ppi_size10 = pco_ppi[pco_ppi$mod %in% non_overlap_mod_size10, ]
n_nov_ppi_size10 = aggregate(is_nov ~ mod, data = pco_ppi_size10, sum)

pco_rand_size10 = pco_rand[pco_rand$mod %in% paste0('mod', c(51:100, 201:250)), ]
n_nov_rand_size10 = aggregate(is_nov ~ mod, data = pco_rand_size10, sum)

n_nov_size10_tab = data.frame(mod = non_overlap_mod_size10, 
                             mod_ori = ppi_mod_size10$mod_ori[match(non_overlap_mod_size10, ppi_mod_size10$mod)],
                             n_nov = n_nov_ppi_size10$is_nov[match(non_overlap_mod_size10, n_nov_ppi_size10$mod)])

n_nov_size10_tab = rbind(n_nov_size10_tab, n_nov_rand_size10_tab)
n_nov_size10_tab[is.na(n_nov_size10_tab)] = 0
n_nov_size10_tab$group = 'PPI'
n_nov_size10_tab$group[which(n_nov_size10_tab$mod_ori == 'random')] = 'random'

wilcox.test(n_nov_size10_tab$n_nov[n_nov_size10_tab$mod_ori == 'random'], 
            n_nov_size10_tab$n_nov[n_nov_size10_tab$mod_ori != 'random'])

p7=ggplot(n_nov_size10_tab, aes(x=n_nov, fill=group)) +
  geom_density(alpha=0.4) + 
  labs(title = "Module size = 10, Mann-Whitney p = 0.08, non-overlapping", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p8=ggplot(n_nov_size10_tab, aes(x=n_nov, color=mod_ori)) +
  geom_density() + 
  labs(title = "Module size = 10, non-overlapping", x = '# novel loci', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1,p2,p5,p6, nrow=2)
grid.arrange(p3,p4,p7,p8, nrow=2)


#### number of novel loci in total
library(boot)
boot_mean = function(data, i) {
  df = data[i]
  sum(df)
}

rand_boot = function(data, n_sample, rep) {
  rand_total_nov = c()
  for (i in 1:rep) {
    df = sample(x = data, size = n_sample, replace = T)
    rand_total_nov[i] = sum(df)
  }
  return(rand_total_nov)
}

# modules with 3 proteins
pco_ppi_size3 = pco_ppi[pco_ppi$mod %in% mod_size3, ]
n_nov_ppi_size3 = aggregate(is_nov ~ mod, data = pco_ppi_size3, sum)

pco_rand_size3 = pco_rand[pco_rand$mod %in% paste0('mod', c(1:50, 151:200)), ]
n_nov_rand_size3 = aggregate(is_nov ~ mod, data = pco_rand_size3, sum)

n_nov_size3_tab = data.frame(mod = mod_size3, 
                             mod_ori = ppi_mod_size3$mod_ori[match(mod_size3, ppi_mod_size3$mod)],
                             n_nov = n_nov_ppi_size3$is_nov[match(mod_size3, n_nov_ppi_size3$mod)])
n_nov_rand_size3_tab = data.frame(mod = paste0('mod', c(1:50, 151:200)),
                                  mod_ori = rep('random', 100),
                                  n_nov = n_nov_rand_size3$is_nov[match(paste0('mod', c(1:50, 151:200)),
                                                                        n_nov_rand_size3$mod)])
n_nov_size3_tab = rbind(n_nov_size3_tab, n_nov_rand_size3_tab)
n_nov_size3_tab[is.na(n_nov_size3_tab)] = 0
n_nov_size3_tab$group = 'PPI'
n_nov_size3_tab$group[which(n_nov_size3_tab$mod_ori == 'random')] = 'random'

total_nov_ppi_size3 = boot(n_nov_size3_tab$n_nov[which(n_nov_size3_tab$group == 'PPI')], 
                           boot_mean, 1000)
total_nov_rand_size3 = rand_boot(n_nov_size3_tab$n_nov[which(n_nov_size3_tab$group == 'random')], 
                                 length(mod_size3), 1000)
wilcox.test(total_nov_ppi_size3$t[,1], total_nov_rand_size3)

# modules with 10 proteins
pco_ppi_size10 = pco_ppi[pco_ppi$mod %in% mod_size10, ]
n_nov_ppi_size10 = aggregate(is_nov ~ mod, data = pco_ppi_size10, sum)

pco_rand_size10 = pco_rand[pco_rand$mod %in% paste0('mod', c(51:100, 201:250)), ]
n_nov_rand_size10 = aggregate(is_nov ~ mod, data = pco_rand_size10, sum)

n_nov_size10_tab = data.frame(mod = mod_size10, 
                             mod_ori = ppi_mod_size10$mod_ori[match(mod_size10, ppi_mod_size10$mod)],
                             n_nov = n_nov_ppi_size10$is_nov[match(mod_size10, n_nov_ppi_size10$mod)])
n_nov_rand_size10_tab = data.frame(mod = paste0('mod', c(51:100, 201:250)),
                                  mod_ori = rep('random', 100),
                                  n_nov = n_nov_rand_size10$is_nov[match(paste0('mod', c(51:100, 201:250)),
                                                                        n_nov_rand_size10$mod)])
n_nov_size10_tab = rbind(n_nov_size10_tab, n_nov_rand_size10_tab)
n_nov_size10_tab[is.na(n_nov_size10_tab)] = 0
n_nov_size10_tab$group = 'PPI'
n_nov_size10_tab$group[which(n_nov_size10_tab$mod_ori == 'random')] = 'random'

total_nov_ppi_size10 = boot(n_nov_size10_tab$n_nov[which(n_nov_size10_tab$group == 'PPI')], 
                           boot_mean, 1000)
total_nov_rand_size10 = rand_boot(n_nov_size10_tab$n_nov[which(n_nov_size10_tab$group == 'random')], 
                                 length(mod_size10), 1000)
wilcox.test(total_nov_ppi_size10$t[,1], total_nov_rand_size10)

n_nov_boot_plot = data.frame(n_nov = c(mean(total_nov_ppi_size3$t), 
                                       mean(total_nov_rand_size3),
                                       mean(total_nov_ppi_size10$t),
                                       mean(total_nov_rand_size10)),
                     lower = c(quantile(total_nov_ppi_size3$t, 0.025), 
                               quantile(total_nov_rand_size3, 0.025),
                               quantile(total_nov_ppi_size10$t, 0.025),
                               quantile(total_nov_rand_size10, 0.025)), 
                     upper = c(quantile(total_nov_ppi_size3$t, 0.975),
                               quantile(total_nov_rand_size3, 0.975),
                               quantile(total_nov_ppi_size10$t, 0.975),
                               quantile(total_nov_rand_size10, 0.975)),
                     mod_size = as.factor(c(3, 3, 10, 10)),
                     group = c('PPI', 'random', 'PPI', 'random'))

ggplot(n_nov_boot_plot, aes(x = mod_size, y = n_nov, fill = group)) + 
  geom_bar(stat="identity", position=position_dodge(0.6), width = 0.5) +
  geom_errorbar(position=position_dodge(0.6), aes(ymin=lower, ymax=upper), width=.05) +
  labs(x = "module size", y = 'total # novel loci', title = '') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


#### number of unique novel loci
# modules with 3 proteins
pco_ppi_size3 = pco_ppi[pco_ppi$mod %in% mod_size3, ]
pco_rand_size3 = pco_rand[pco_rand$mod %in% paste0('mod', c(1:50, 151:200)), ]

ppi_uniq_loci_size3 = c()
pco_ppi_size3$snp = paste0(pco_ppi_size3$chr, ':', pco_ppi_size3$bp_hg37)
uniq_snp = unique(pco_ppi_size3$snp)
for (i in 1:length(uniq_snp)) {
  pco_loci_i = pco_ppi_size3[pco_ppi_size3$snp == uniq_snp[i], ]
  lead_sig = which.min(pco_loci_i$pco_p)
  lead_snp = pco_loci_i[lead_sig, ]
  lead_snp$loci_start = min(pco_loci_i$loci_start)
  lead_snp$loci_end = max(pco_loci_i$loci_end)
  lead_snp$is_nov_all = (sum(lead_snp$is_nov) == nrow(lead_snp))
  lead_snp$is_nov_any = (sum(lead_snp$is_nov) > 0)
  ppi_uniq_loci_size3 = rbind(ppi_uniq_loci_size3, lead_snp)
}
ppi_uniq_loci_size3 = ppi_uniq_loci_size3[order(ppi_uniq_loci_size3$chr, 
                                                ppi_uniq_loci_size3$loci_start), ]
ppi_uniq_final_size3 = ppi_uniq_loci_size3[1, ]
for (i in 2:nrow(ppi_uniq_loci_size3)) {
  chr_i = ppi_uniq_loci_size3$chr[i]
  start_i = ppi_uniq_loci_size3$loci_start[i]
  chr_prev = ppi_uniq_final_size3$chr[nrow(ppi_uniq_final_size3)]
  end_prev = ppi_uniq_final_size3$loci_end[nrow(ppi_uniq_final_size3)]
  if (chr_i != chr_prev | (chr_i == chr_prev & start_i > end_prev)) {
    # add the loci directly if two loci do not overlap
    ppi_uniq_final_size3 = rbind(ppi_uniq_final_size3, ppi_uniq_loci_size3[i, ])
  } else {
    # merge two overlapping loci
    add_pool_i = rbind(ppi_uniq_final_size3[nrow(ppi_uniq_final_size3), ], ppi_uniq_loci_size3[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$pco_p), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    add_loci_i$is_nov_all = (sum(add_pool_i$is_nov_all) == nrow(add_pool_i))
    add_loci_i$is_nov_any = (sum(add_pool_i$is_nov_any) > 0)
    ppi_uniq_final_size3[nrow(ppi_uniq_final_size3), ] = add_loci_i
  }
}

rand_uniq_loci_size3 = c()
pco_rand_size3$snp = paste0(pco_rand_size3$chr, ':', pco_rand_size3$bp_hg37)
uniq_snp = unique(pco_rand_size3$snp)
for (i in 1:length(uniq_snp)) {
  pco_loci_i = pco_rand_size3[pco_rand_size3$snp == uniq_snp[i], ]
  lead_sig = which.min(pco_loci_i$pco_p)
  lead_snp = pco_loci_i[lead_sig, ]
  lead_snp$bp_start = min(pco_loci_i$bp_start)
  lead_snp$bp_end = max(pco_loci_i$bp_end)
  lead_snp$is_nov_all = (sum(lead_snp$is_nov) == nrow(lead_snp))
  lead_snp$is_nov_any = (sum(lead_snp$is_nov) > 0)
  rand_uniq_loci_size3 = rbind(rand_uniq_loci_size3, lead_snp)
}
rand_uniq_loci_size3 = rand_uniq_loci_size3[order(rand_uniq_loci_size3$chr,
                                                  rand_uniq_loci_size3$bp_start), ]
rand_uniq_final_size3 = rand_uniq_loci_size3[1, ]
for (i in 2:nrow(rand_uniq_loci_size3)) {
  chr_i = rand_uniq_loci_size3$chr[i]
  start_i = rand_uniq_loci_size3$bp_start[i]
  chr_prev = rand_uniq_final_size3$chr[nrow(rand_uniq_final_size3)]
  end_prev = rand_uniq_final_size3$bp_end[nrow(rand_uniq_final_size3)]
  if (chr_i != chr_prev | (chr_i == chr_prev & start_i > end_prev)) {
    # add the loci directly if two loci do not overlap
    rand_uniq_final_size3 = rbind(rand_uniq_final_size3, rand_uniq_loci_size3[i, ])
  } else {
    # merge two overlapping loci
    add_pool_i = rbind(rand_uniq_final_size3[nrow(rand_uniq_final_size3), ], rand_uniq_loci_size3[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$pco_p), ]
    add_loci_i$bp_start = min(add_pool_i$bp_start)
    add_loci_i$bp_end = max(add_pool_i$bp_end)
    add_loci_i$is_nov_all = (sum(add_pool_i$is_nov_all) == nrow(add_pool_i))
    add_loci_i$is_nov_any = (sum(add_pool_i$is_nov_any) > 0)
    rand_uniq_final_size3[nrow(rand_uniq_final_size3), ] = add_loci_i
  }
}
sum(rand_uniq_final_size3$is_nov_all) / nrow(rand_uniq_final_size3)
sum(ppi_uniq_final_size3$is_nov_all) / nrow(ppi_uniq_final_size3)


# modules with 10 proteins
pco_ppi_size10 = pco_ppi[pco_ppi$mod %in% mod_size10, ]
pco_rand_size10 = pco_rand[pco_rand$mod %in% paste0('mod', c(51:100, 201:250)), ]

ppi_uniq_loci_size10 = c()
pco_ppi_size10$snp = paste0(pco_ppi_size10$chr, ':', pco_ppi_size10$bp_hg37)
uniq_snp = unique(pco_ppi_size10$snp)
for (i in 1:length(uniq_snp)) {
  pco_loci_i = pco_ppi_size10[pco_ppi_size10$snp == uniq_snp[i], ]
  lead_sig = which.min(pco_loci_i$pco_p)
  lead_snp = pco_loci_i[lead_sig, ]
  lead_snp$loci_start = min(pco_loci_i$loci_start)
  lead_snp$loci_end = max(pco_loci_i$loci_end)
  lead_snp$is_nov_all = (sum(lead_snp$is_nov) == nrow(lead_snp))
  lead_snp$is_nov_any = (sum(lead_snp$is_nov) > 0)
  ppi_uniq_loci_size10 = rbind(ppi_uniq_loci_size10, lead_snp)
}
ppi_uniq_loci_size10 = ppi_uniq_loci_size10[order(ppi_uniq_loci_size10$chr, 
                                                ppi_uniq_loci_size10$loci_start), ]
ppi_uniq_final_size10 = ppi_uniq_loci_size10[1, ]
for (i in 2:nrow(ppi_uniq_loci_size10)) {
  chr_i = ppi_uniq_loci_size10$chr[i]
  start_i = ppi_uniq_loci_size10$loci_start[i]
  chr_prev = ppi_uniq_final_size10$chr[nrow(ppi_uniq_final_size10)]
  end_prev = ppi_uniq_final_size10$loci_end[nrow(ppi_uniq_final_size10)]
  if (chr_i != chr_prev | (chr_i == chr_prev & start_i > end_prev)) {
    # add the loci directly if two loci do not overlap
    ppi_uniq_final_size10 = rbind(ppi_uniq_final_size10, ppi_uniq_loci_size10[i, ])
  } else {
    # merge two overlapping loci
    add_pool_i = rbind(ppi_uniq_final_size10[nrow(ppi_uniq_final_size10), ], ppi_uniq_loci_size10[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$pco_p), ]
    add_loci_i$loci_start = min(add_pool_i$loci_start)
    add_loci_i$loci_end = max(add_pool_i$loci_end)
    add_loci_i$is_nov_all = (sum(add_pool_i$is_nov_all) == nrow(add_pool_i))
    add_loci_i$is_nov_any = (sum(add_pool_i$is_nov_any) > 0)
    ppi_uniq_final_size10[nrow(ppi_uniq_final_size10), ] = add_loci_i
  }
}

rand_uniq_loci_size10 = c()
pco_rand_size10$snp = paste0(pco_rand_size10$chr, ':', pco_rand_size10$bp_hg37)
uniq_snp = unique(pco_rand_size10$snp)
for (i in 1:length(uniq_snp)) {
  pco_loci_i = pco_rand_size10[pco_rand_size10$snp == uniq_snp[i], ]
  lead_sig = which.min(pco_loci_i$pco_p)
  lead_snp = pco_loci_i[lead_sig, ]
  lead_snp$bp_start = min(pco_loci_i$bp_start)
  lead_snp$bp_end = max(pco_loci_i$bp_end)
  lead_snp$is_nov_all = (sum(lead_snp$is_nov) == nrow(lead_snp))
  lead_snp$is_nov_any = (sum(lead_snp$is_nov) > 0)
  rand_uniq_loci_size10 = rbind(rand_uniq_loci_size10, lead_snp)
}
rand_uniq_loci_size10 = rand_uniq_loci_size10[order(rand_uniq_loci_size10$chr, 
                                                  rand_uniq_loci_size10$bp_start), ]
rand_uniq_final_size10 = rand_uniq_loci_size10[1, ]
for (i in 2:nrow(rand_uniq_loci_size10)) {
  chr_i = rand_uniq_loci_size10$chr[i]
  start_i = rand_uniq_loci_size10$bp_start[i]
  chr_prev = rand_uniq_final_size10$chr[nrow(rand_uniq_final_size10)]
  end_prev = rand_uniq_final_size10$bp_end[nrow(rand_uniq_final_size10)]
  if (chr_i != chr_prev | (chr_i == chr_prev & start_i > end_prev)) {
    # add the loci directly if two loci do not overlap
    rand_uniq_final_size10 = rbind(rand_uniq_final_size10, rand_uniq_loci_size10[i, ])
  } else {
    # merge two overlapping loci
    add_pool_i = rbind(rand_uniq_final_size10[nrow(rand_uniq_final_size10), ], rand_uniq_loci_size10[i, ])
    add_loci_i = add_pool_i[which.min(add_pool_i$pco_p), ]
    add_loci_i$bp_start = min(add_pool_i$bp_start)
    add_loci_i$bp_end = max(add_pool_i$bp_end)
    add_loci_i$is_nov_all = (sum(add_pool_i$is_nov_all) == nrow(add_pool_i))
    add_loci_i$is_nov_any = (sum(add_pool_i$is_nov_any) > 0)
    rand_uniq_final_size10[nrow(rand_uniq_final_size10), ] = add_loci_i
  }
}
sum(rand_uniq_final_size10$is_nov_any) / nrow(rand_uniq_final_size10)
sum(ppi_uniq_final_size10$is_nov_any) / nrow(ppi_uniq_final_size10)

uniq_nov_tab = data.frame(prop_any = c(sum(rand_uniq_final_size3$is_nov_any) / nrow(rand_uniq_final_size3),
                                   sum(ppi_uniq_final_size3$is_nov_any) / nrow(ppi_uniq_final_size3),
                                   sum(rand_uniq_final_size10$is_nov_any) / nrow(rand_uniq_final_size10),
                                   sum(ppi_uniq_final_size10$is_nov_any) / nrow(ppi_uniq_final_size10)),
                          prop_all = c(sum(rand_uniq_final_size3$is_nov_all) / nrow(rand_uniq_final_size3),
                                       sum(ppi_uniq_final_size3$is_nov_all) / nrow(ppi_uniq_final_size3),
                                       sum(rand_uniq_final_size10$is_nov_all) / nrow(rand_uniq_final_size10),
                                       sum(ppi_uniq_final_size10$is_nov_all) / nrow(ppi_uniq_final_size10)),
                          count_any = c(sum(rand_uniq_final_size3$is_nov_any),
                                        sum(ppi_uniq_final_size3$is_nov_any),
                                        sum(rand_uniq_final_size10$is_nov_any),
                                        sum(ppi_uniq_final_size10$is_nov_any)),
                          count_all = c(sum(rand_uniq_final_size3$is_nov_all),
                                        sum(ppi_uniq_final_size3$is_nov_all),
                                        sum(rand_uniq_final_size10$is_nov_all),
                                        sum(ppi_uniq_final_size10$is_nov_all)),
                          mod_size = as.factor(c(3,3,10,10)),
                          group = c('random', 'PPI', 'random', 'PPI'))
p1 = ggplot(uniq_nov_tab, aes(x = mod_size, y = count_any, fill = group, label = count_any)) + 
  geom_bar(stat="identity", position=position_dodge(0.6), width = 0.5) +
  geom_text(vjust=-0.2, size = 4.5, position = position_dodge(0.6)) +
  labs(x = "module size", y = '# of novel loci', title = 'Novel loci: any locus-mod pairs included are novel') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p2 = ggplot(uniq_nov_tab, aes(x = mod_size, y = count_all, fill = group, label = count_all)) + 
  geom_bar(stat="identity", position=position_dodge(0.6), width = 0.5) +
  geom_text(vjust=-0.2, size = 4.5, position = position_dodge(0.6)) + 
  labs(x = "module size", y = '# of novel loci', title = 'Novel loci: all locus-mod pairs included are novel') + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)


## minP distribution
ppi_minp = aggregate(pco_p ~ mod, data = pco_ppi, min)
rand_minp = aggregate(pco_p ~ mod, data = pco_rand, min)
ppi_minp$group = 'PPI'
rand_minp$group = 'random'

minp_size3 = rbind(ppi_minp[ppi_minp$mod %in% mod_size3, ],
                   rand_minp[rand_minp$mod %in% paste0('mod', c(1:50, 151:200)), ])
minp_size3$pco_p[minp_size3$pco_p == 0] = 1e-320
minp_size3$pco_p = minp_size3$pco_p * 1e6
minp_size3$log_pco_p = -log10(minp_size3$pco_p)
wilcox.test(minp_size3$log_pco_p[minp_size3$group == 'random'],
            minp_size3$log_pco_p[minp_size3$group == 'PPI'])
p1=ggplot(minp_size3, aes(x=log_pco_p, fill=group)) +
  geom_density(alpha=0.4) + 
  labs(title = "Module size = 3, Mann-Whitney p = 0.19", x = '-log10(minP * 1e6)', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


minp_size10 = rbind(ppi_minp[ppi_minp$mod %in% mod_size10, ],
                   rand_minp[rand_minp$mod %in% paste0('mod', c(51:100, 201:250)), ])
minp_size10$pco_p[minp_size10$pco_p == 0] = 1e-320
minp_size10$pco_p = minp_size10$pco_p * 1e6
minp_size10$log_pco_p = -log10(minp_size10$pco_p)
wilcox.test(minp_size10$log_pco_p[minp_size10$group == 'random'],
            minp_size10$log_pco_p[minp_size10$group == 'PPI'])
p2=ggplot(minp_size10, aes(x=log_pco_p, fill=group)) +
  geom_density(alpha=0.4) + 
  labs(title = "Module size = 10, Mann-Whitney p = 0.22", x = '-log10(minP * 1e6)', y = 'Density') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 13),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1,p2, nrow = 1)

