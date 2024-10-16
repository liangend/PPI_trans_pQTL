library(data.table)
library(ggplot2)
library(gridExtra)

small_p = fread('/project/xuanyao/jinghui/pqtl/04_fdr/small_p.txt')
n_mod = 27
# Bonferroni correction
sig_p = small_p[small_p$V2 <= 5e-8/n_mod, ]

sig_mod = sapply(strsplit(sig_p$V1, ':', fixed = T), '[', 1)
sig_chr = sapply(strsplit(sig_p$V1, ':', fixed = T), '[', 2)
sig_bp = sapply(strsplit(sig_p$V1, ':', fixed = T), '[', 3)
sig_all = data.frame(mod = sig_mod, snp = paste0(sig_chr, ':', sig_bp), 
                     chr = sig_chr, bp = as.numeric(sig_bp), p = sig_p$V2)
fwrite(sig_all, '/project/xuanyao/jinghui/pqtl/04_fdr/sig_pair_all.txt', sep = '\t')

sig_all = fread('/project/xuanyao/jinghui/pqtl/04_fdr/sig_pair_all.txt')
prot_meta = fread('/project/xuanyao/jinghui/pqtl/00_ref/prot_meta.txt')
# signal distribution
mod_meta_table = as.data.frame(table(prot_meta$module))
names(mod_meta_table) = c("module", "n_protein")
ggplot(mod_meta_table, aes(x=reorder(module, -n_protein), y=n_protein)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title = 'Module division', x = 'module', y = '# protein') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

mod_table = as.data.frame(table(sig_all$mod))
names(mod_table) = c("module", "n_SNP")
f1 = ggplot(mod_table, aes(x=reorder(module, n_SNP), y=n_SNP)) + 
  coord_flip() +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title = 'Module distribution', x = '', y = '# SNP') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

chr_table = as.data.frame(table(sig_all$chr))
names(chr_table) = c("chr", "n_SNP")
f2 = ggplot(chr_table, aes(x=reorder(chr, n_SNP), y=n_SNP)) + 
  coord_flip() +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title = 'Chr distribution', x = '', y = '# SNP') +
  theme(text = element_text(size=15, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15),
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(f1, f2, nrow = 1)


## relationship between module size, module cor and # sig associations
mod_table = as.data.frame(table(sig_all$mod))
colnames(mod_table) = c('mod', 'n_sig')
mod_table$mod = as.numeric(sub('mod', '', mod_table$mod))
mod_table = mod_table[order(mod_table$mod), ]

mod_meta = as.data.frame(table(prot_meta$module))
colnames(mod_meta) = c('mod', 'n_gene')
mod_table$n_gene = mod_meta$n_gene

mod_cor = matrix(rep(0, 27*22), nrow = 27, ncol = 22)
for (i in 1:27) {
  for (j in 1:22) {
    cor_ij = readRDS(paste0('/project/xuanyao/jinghui/pqtl/02_sigma/mod', i, '.chr', j, '.rds'))
    # average correlation of module i
    mod_cor[i,j] = (sum(abs(cor_ij)) - nrow(cor_ij)) / (nrow(cor_ij) * (nrow(cor_ij) - 1))
  }
}
mod_table$cor = apply(mod_cor, 1, mean)
mod_table$name = paste0('mod', mod_table$mod)
mod_table$name[-c(3,5,7,9,26,27)] = ''
cor(mod_table$cor, mod_table$n_sig)
p1 = ggplot(mod_table, aes(x = n_gene, y = n_sig, label = name)) + geom_point(size = 2) +
  labs(x = "# genes", y = "# sig asso", title = 'Module size') + 
  geom_text(hjust=0.8, vjust=1.5, size = 4.5) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2 = ggplot(mod_table, aes(x = cor, y = n_sig, label = name)) + geom_point(size = 2) +
  labs(x = "# genes", y = "# sig asso", title = 'Module cor') + 
  geom_text(hjust=0.5, vjust=1.5, size = 4.5) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)


