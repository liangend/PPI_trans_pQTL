library(data.table)
library(ggplot2)
library(gridExtra)
trans_pqtl = fread('/project/xuanyao/jinghui/pqtl/04_fdr/sig_pair_all.txt')
mean(trans_pqtl$pi1, na.rm = T)
prot_meta = fread('/project/xuanyao/jinghui/pqtl/00_ref/prot_meta.txt')
mod_prot_count = as.data.frame(table(prot_meta$module))
names(mod_prot_count) = c('mod', 'n_prot')
mod_prot_count$mod = paste0('mod', mod_prot_count$mod)

trans_pqtl$n_prot = mod_prot_count$n_prot[match(trans_pqtl$mod, mod_prot_count$mod)]
trans_pqtl$n_prot = as.factor(trans_pqtl$n_prot)
p1 = ggplot(trans_pqtl, aes(x = n_prot, y = pi1, color = n_prot)) + 
  geom_boxplot(width = 0.4) +
  geom_hline(yintercept = mean(trans_pqtl$pi1, na.rm = T), linetype = 2) + 
  geom_hline(yintercept = median(trans_pqtl$pi1, na.rm = T), linetype = 2, color = 'red') + 
  labs(x = "", y = "pi1") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")


na_prop = aggregate(pi1 ~ n_prot, data = trans_pqtl, function(x){sum(is.na(x))}, na.action = na.pass)
na_prop$pi1 = na_prop$pi1/table(trans_pqtl$n_prot)
na_prop$pi1 = round(na_prop$pi1, 2)
p2 = ggplot(na_prop, aes(x = n_prot, y = pi1)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = 'steelblue') +
  geom_text(aes(label = pi1), vjust=-0.3, size=4, position = position_dodge(0.5)) + 
  labs(x = "# proteins in module", y = "Proportion of NA pi1") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, nrow = 2)


