library(igraph)
library(data.table)
library(ggplot2)
library(gridExtra)
library(openxlsx)
setwd('/project/xuanyao/jinghui')

### using igraph to do clustering
# bioplex = fread('pqtl/11_prot_complex/BioPlex_293T_Network_10K_Dec_2019.tsv')
# bioplex_ppi = data.frame(prot1 = bioplex$SymbolA, prot2 = bioplex$SymbolB)
# 
# g = graph.edgelist(as.matrix(test_ppi), directed = F)
# prot_clust1 = fastgreedy.community(g)
# prot_clust2 = walktrap.community(g)
# 
# mod1 = c()
# for (i in 1:length(prot_clust1)) {
#   mod_i = data.frame(prot = prot_clust1[[i]], mod = i)
#   mod1 = rbind(mod1, mod_i)
# }
# mod2 = c()
# for (i in 1:length(prot_clust2)) {
#   mod_i = data.frame(prot = prot_clust2[[i]], mod = i)
#   mod2 = rbind(mod2, mod_i)
# }
# 
# n_prot_bioplex1 = as.data.frame(table(sizes(prot_clust1)))
# colnames(n_prot_bioplex1) = c('n_prot', 'Freq')
# n_prot_bioplex1$n_prot = as.numeric(as.character(n_prot_bioplex1$n_prot))
# 
# n_prot_bioplex2 = as.data.frame(table(sizes(prot_clust2)))
# colnames(n_prot_bioplex2) = c('n_prot', 'Freq')
# n_prot_bioplex2$n_prot = as.numeric(as.character(n_prot_bioplex2$n_prot))

#### make PPI modules based on Ahn
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')
ukb_prot$gene_name = sapply(strsplit(ukb_prot$file, '_', fixed = T), '[', 1)

## BioPlex
ahn_bioplex = fread('pqtl/11_prot_complex/ppi_input/Ahn_result/bioplex_ppi_maxS0.299065_maxD0.054792.edge2comm.txt', 
                    header = F)
uniq_mod = unique(ahn_bioplex$V3)
bioplex_mod = c()
for (i in 1:length(uniq_mod)) {
  mod_i = ahn_bioplex[ahn_bioplex$V3 == uniq_mod[i], 1:2]
  prot_i = unique(unlist(mod_i))
  bioplex_mod = rbind(bioplex_mod, data.frame(prot = prot_i, 
                                              mod = rep(uniq_mod[i], length(prot_i))))
}
bioplex_mod$mod = paste0('bioplex_', bioplex_mod$mod)
# keep prot in UKB
bioplex_mod_sub = bioplex_mod[bioplex_mod$prot %in% ukb_prot$gene_name, ]
# remove modules with only 1 prot
bioplex_mod_sub = bioplex_mod_sub[bioplex_mod_sub$mod %in% 
                                    names(which(table(bioplex_mod_sub$mod) > 1)), ]

## STRING
ahn_string = fread('pqtl/11_prot_complex/ppi_input_pc2p/string_ppi_maxS0.590062_maxD0.379506.edge2comm.txt', 
                    header = F)
uniq_mod = unique(ahn_string$V3)
string_mod = c()
for (i in 1:length(uniq_mod)) {
  mod_i = ahn_string[ahn_string$V3 == uniq_mod[i], 1:2]
  prot_i = unique(unlist(mod_i))
  string_mod = rbind(string_mod, data.frame(prot = prot_i, 
                                            mod = rep(uniq_mod[i], length(prot_i))))
}
string_mod$mod = paste0('string_', string_mod$mod)
string_mod_sub = string_mod[string_mod$prot %in% ukb_prot$gene_name, ]
string_mod_sub = string_mod_sub[string_mod_sub$mod %in% 
                                    names(which(table(string_mod_sub$mod) > 1)), ]

## HIPPIE
ahn_hippie = fread('pqtl/11_prot_complex/ppi_input_pc2p/hippie_ppi_maxS0.251462_maxD0.067695.edge2comm.txt', 
                   header = F)
uniq_mod = unique(ahn_hippie$V3)
hippie_mod = c()
for (i in 1:length(uniq_mod)) {
  mod_i = ahn_hippie[ahn_hippie$V3 == uniq_mod[i], 1:2]
  prot_i = unique(unlist(mod_i))
  hippie_mod = rbind(hippie_mod, data.frame(prot = prot_i, 
                                            mod = rep(uniq_mod[i], length(prot_i))))
}
hippie_mod$mod = paste0('hippie_', hippie_mod$mod)
hippie_mod_sub = hippie_mod[hippie_mod$prot %in% ukb_prot$gene_name, ]
hippie_mod_sub = hippie_mod_sub[hippie_mod_sub$mod %in% 
                                  names(which(table(hippie_mod_sub$mod) > 1)), ]

## CORUM (not clustered by Ahn)
corum = read.xlsx('pqtl/11_prot_complex/CORUM_2022_09_12.xlsx')
corum = corum[corum$Organism == 'Human', ]
prot_list = strsplit(corum$`subunits(Gene.name)`, ';', fixed = T)
corum_size = sapply(prot_list, length)
corum_mod = data.frame(prot = unlist(prot_list), 
                       mod = rep(1:length(prot_list), corum_size))
corum_mod$mod = paste0('corum_', corum_mod$mod)
corum_mod_sub = corum_mod[corum_mod$prot %in% ukb_prot$gene_name, ]
corum_mod_sub = corum_mod_sub[corum_mod_sub$mod %in% 
                                names(which(table(corum_mod_sub$mod) > 1)), ]

### merge 4 databases together
ahn_ppi_all = rbind(bioplex_mod_sub, string_mod_sub, 
                    hippie_mod_sub, corum_mod_sub)
ahn_ppi_all = ahn_ppi_all[order(ahn_ppi_all$mod, ahn_ppi_all$prot), ]

# remove repeated modules
prot_each_mod = aggregate(prot ~ mod, data = ahn_ppi_all, unique)
prot_vec = sapply(prot_each_mod$prot, function(x){paste(x, collapse = ',')})
ppi_tab = data.frame(mod = prot_each_mod$mod, prot = prot_vec)
prot_uniq = unique(ppi_tab$prot)
mod_uniq = ppi_tab$mod[match(prot_uniq, ppi_tab$prot)]
ppi_tab_uniq = c()
for (i in 1:length(mod_uniq)) {
  prot_i = unlist(strsplit(prot_uniq[i], ',', fixed = T))
  mod_i = data.frame(mod = rep(mod_uniq[i], length(prot_i)), prot = prot_i)
  ppi_tab_uniq = rbind(ppi_tab_uniq, mod_i)
}

# remove modules that are sub-modules of larger ones
is_sub = c()
uniq_mod = unique(ppi_tab_uniq$mod)
for (i in 1:length(uniq_mod)) {
  prot_i = ppi_tab_uniq$prot[which(ppi_tab_uniq$mod == uniq_mod[i])]
  n_prot_i = length(prot_i)
  n_overlap_i = table(ppi_tab_uniq$mod[which(ppi_tab_uniq$prot %in% prot_i)])
  is_sub[i] = sum(n_overlap_i == n_prot_i)
  if (i %% 100 == 0) {
    print(i)
  }
}
ahn_ppi_mod = ppi_tab_uniq[ppi_tab_uniq$mod %in% 
                             (uniq_mod[is_sub == 1]), ]
colnames(ahn_ppi_mod)[1] = 'mod_ori'
n_prot = table(ahn_ppi_mod$mod_ori)
ahn_ppi_mod$mod = rep(1:length(n_prot), n_prot)

ahn_ppi_mod$gene_id = ukb_prot$gene_id[match(ahn_ppi_mod$prot, ukb_prot$gene_name)]
ahn_ppi_mod$chr = ukb_prot$chr[match(ahn_ppi_mod$prot, ukb_prot$gene_name)]
ahn_ppi_mod$file = ukb_prot$file[match(ahn_ppi_mod$prot, ukb_prot$gene_name)]

uniq_mod_ori = unique(ahn_ppi_mod$mod_ori)
all_ppi = rbind(bioplex_mod, string_mod, hippie_mod, corum_mod)
prot_list = c()
for (i in uniq_mod_ori) {
  prot_i = all_ppi$prot[which(all_ppi$mod == i)]
  prot_list = c(prot_list, paste(prot_i, collapse = ','))
}
ahn_mod_all_prot = data.frame(mod = uniq_mod_ori, prot = prot_list)

fwrite(ahn_ppi_mod, 'pqtl/11_prot_complex/ppi_input/Ahn_result/ahn_ppi_mod.txt', sep = '\t')
fwrite(ahn_mod_all_prot, 'pqtl/11_prot_complex/ppi_input/Ahn_result/ahn_mod_all_prot.txt', sep = '\t')

n_prot_ahn = as.data.frame(table(table(ahn_ppi_mod$mod)))
colnames(n_prot_ahn) = c('n_prot', 'Freq')
n_prot_ahn$n_prot = as.numeric(as.character(n_prot_ahn$n_prot))
n_prot_tab = rbind(n_prot_ahn[n_prot_ahn$n_prot < 10, ],
                   data.frame(n_prot = '> 10', 
                              Freq = sum(n_prot_ahn$Freq[which(n_prot_ahn$n_prot >= 10)])))

p1 = ggplot(n_prot_tab, aes(x = factor(n_prot, 
                                  levels = c(as.character(2:9),'> 10')), 
                               y = Freq, label = Freq)) + 
  geom_bar(stat="identity", fill = 'steelblue', width = 0.7) +
  geom_text(vjust=-0.2, size = 4.5) + 
  labs(x = "# protein", y = "# module",
       title = paste0('Ahn cluster, # modules = ', max(ahn_ppi_mod$mod))) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

overlap_prot_ahn = as.data.frame(table(table(ahn_ppi_mod$prot)))
colnames(overlap_prot_ahn) = c('n_prot', 'Freq')
overlap_prot_ahn$n_prot = as.numeric(as.character(overlap_prot_ahn$n_prot))
overlap_prot_tab = rbind(overlap_prot_ahn[overlap_prot_ahn$n_prot < 10, ],
                         data.frame(n_prot = '> 10', 
                              Freq = sum(overlap_prot_ahn$Freq[which(overlap_prot_ahn$n_prot >= 10)])))

p2 = ggplot(overlap_prot_tab, aes(x = factor(n_prot, 
                                  levels = c(as.character(1:9),'> 10')), 
                       y = Freq, label = Freq)) + 
  geom_bar(stat="identity", fill = 'steelblue', width = 0.7) +
  geom_text(vjust=-0.2, size = 4.5) + 
  labs(x = "# module", y = "# protein",
       title = paste0('Ahn cluster, # uniq protein = ', length(unique(ahn_ppi_mod$prot)))) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)

### my method
ppi_mod = fread('pqtl/11_prot_complex/4_database_mod.txt')
n_prot_old = as.data.frame(table(table(ppi_mod$mod)))
colnames(n_prot_old) = c('n_prot', 'Freq')
n_prot_old$n_prot = as.numeric(as.character(n_prot_old$n_prot))
n_prot_tab = rbind(n_prot_old[n_prot_old$n_prot < 10, ],
                   data.frame(n_prot = '> 10', 
                              Freq = sum(n_prot_old$Freq[which(n_prot_old$n_prot >= 10)])))

p1 = ggplot(n_prot_tab, aes(x = factor(n_prot, 
                                       levels = c(as.character(2:9),'> 10')), 
                            y = Freq, label = Freq)) + 
  geom_bar(stat="identity", fill = 'steelblue', width = 0.7) +
  geom_text(vjust=-0.2, size = 4.5) + 
  labs(x = "# protein", y = "# module",
       title = paste0('Old cluster, # modules = ', max(ppi_mod$mod))) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

overlap_prot_old = as.data.frame(table(table(ppi_mod$prot)))
colnames(overlap_prot_old) = c('n_prot', 'Freq')
overlap_prot_old$n_prot = as.numeric(as.character(overlap_prot_old$n_prot))
overlap_prot_old = rbind(overlap_prot_old[overlap_prot_old$n_prot < 10, ],
                         data.frame(n_prot = '> 10', 
                                    Freq = sum(overlap_prot_old$Freq[which(overlap_prot_old$n_prot >= 10)])))

p2 = ggplot(overlap_prot_old, aes(x = factor(n_prot, 
                                             levels = c(as.character(1:9),'> 10')), 
                                  y = Freq, label = Freq)) + 
  geom_bar(stat="identity", fill = 'steelblue', width = 0.7) +
  geom_text(vjust=-0.2, size = 4.5) + 
  labs(x = "# module", y = "# protein",
       title = paste0('Old cluster, # uniq protein = ', length(unique(ppi_mod$prot)))) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)

##### make PPI modules based on MCL
## BioPlex
mcl_bioplex = fread('pqtl/11_prot_complex/ppi_input/mcl_result/bioplex_I20', 
                    header = F, sep = ',')
bioplex_mod = c()
for (i in 1:nrow(mcl_bioplex)) {
  prot_i = unlist(strsplit(mcl_bioplex$V1[i], '\t', fixed = T))
  bioplex_mod = rbind(bioplex_mod, data.frame(prot = prot_i, 
                                              mod = rep(i, length(prot_i))))
}
bioplex_mod$mod = paste0('bioplex_', bioplex_mod$mod)
# keep prot in UKB
bioplex_mod_sub = bioplex_mod[bioplex_mod$prot %in% ukb_prot$gene_name, ]
# remove modules with only 1 prot
bioplex_mod_sub = bioplex_mod_sub[bioplex_mod_sub$mod %in% 
                                    names(which(table(bioplex_mod_sub$mod) > 1)), ]


## STRING
mcl_string = fread('pqtl/11_prot_complex/ppi_input/mcl_result/string_I20', 
                    header = F, sep = ',')
string_mod = c()
for (i in 1:nrow(mcl_string)) {
  prot_i = unlist(strsplit(mcl_string$V1[i], '\t', fixed = T))
  string_mod = rbind(string_mod, data.frame(prot = prot_i, 
                                              mod = rep(i, length(prot_i))))
}
string_mod$mod = paste0('string_', string_mod$mod)
# keep prot in UKB
string_mod_sub = string_mod[string_mod$prot %in% ukb_prot$gene_name, ]
# remove modules with only 1 prot
string_mod_sub = string_mod_sub[string_mod_sub$mod %in% 
                                    names(which(table(string_mod_sub$mod) > 1)), ]

## HIPPIE
mcl_hippie = fread('pqtl/11_prot_complex/ppi_input/mcl_result/hippie_I20', 
                   header = F, sep = ',')
hippie_mod = c()
for (i in 1:nrow(mcl_hippie)) {
  prot_i = unlist(strsplit(mcl_hippie$V1[i], '\t', fixed = T))
  hippie_mod = rbind(hippie_mod, data.frame(prot = prot_i, 
                                            mod = rep(i, length(prot_i))))
}
hippie_mod$mod = paste0('hippie_', hippie_mod$mod)
# keep prot in UKB
hippie_mod_sub = hippie_mod[hippie_mod$prot %in% ukb_prot$gene_name, ]
# remove modules with only 1 prot
hippie_mod_sub = hippie_mod_sub[hippie_mod_sub$mod %in% 
                                  names(which(table(hippie_mod_sub$mod) > 1)), ]

### merge 4 databases together
mcl_ppi_all = rbind(bioplex_mod_sub, string_mod_sub, 
                    hippie_mod_sub, corum_mod_sub)
mcl_ppi_all = mcl_ppi_all[order(mcl_ppi_all$mod, mcl_ppi_all$prot), ]

# remove repeated modules
prot_each_mod = aggregate(prot ~ mod, data = mcl_ppi_all, unique)
prot_vec = sapply(prot_each_mod$prot, function(x){paste(x, collapse = ',')})
ppi_tab = data.frame(mod = prot_each_mod$mod, prot = prot_vec)
prot_uniq = unique(ppi_tab$prot)
mod_uniq = ppi_tab$mod[match(prot_uniq, ppi_tab$prot)]
ppi_tab_uniq = c()
for (i in 1:length(mod_uniq)) {
  prot_i = unlist(strsplit(prot_uniq[i], ',', fixed = T))
  mod_i = data.frame(mod = rep(mod_uniq[i], length(prot_i)), prot = prot_i)
  ppi_tab_uniq = rbind(ppi_tab_uniq, mod_i)
}

# remove modules that are sub-modules of larger ones
is_sub = c()
uniq_mod = unique(ppi_tab_uniq$mod)
for (i in 1:length(uniq_mod)) {
  prot_i = ppi_tab_uniq$prot[which(ppi_tab_uniq$mod == uniq_mod[i])]
  n_prot_i = length(prot_i)
  n_overlap_i = table(ppi_tab_uniq$mod[which(ppi_tab_uniq$prot %in% prot_i)])
  is_sub[i] = sum(n_overlap_i == n_prot_i)
  if (i %% 100 == 0) {
    print(i)
  }
}
mcl_ppi_mod = ppi_tab_uniq[ppi_tab_uniq$mod %in% 
                             (uniq_mod[is_sub == 1]), ]
colnames(mcl_ppi_mod)[1] = 'mod_ori'
n_prot = table(mcl_ppi_mod$mod_ori)
mcl_ppi_mod$mod = rep(1:length(n_prot), n_prot)

mcl_ppi_mod$gene_id = ukb_prot$gene_id[match(mcl_ppi_mod$prot, ukb_prot$gene_name)]
mcl_ppi_mod$chr = ukb_prot$chr[match(mcl_ppi_mod$prot, ukb_prot$gene_name)]
mcl_ppi_mod$file = ukb_prot$file[match(mcl_ppi_mod$prot, ukb_prot$gene_name)]

uniq_mod_ori = unique(mcl_ppi_mod$mod_ori)
all_ppi = rbind(bioplex_mod, string_mod, hippie_mod, corum_mod)
prot_list = c()
for (i in uniq_mod_ori) {
  prot_i = all_ppi$prot[which(all_ppi$mod == i)]
  prot_list = c(prot_list, paste(prot_i, collapse = ','))
}
mcl_mod_all_prot = data.frame(mod = uniq_mod_ori, prot = prot_list)

fwrite(mcl_ppi_mod, 'pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt', sep = '\t')
fwrite(mcl_mod_all_prot, 'pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_all_prot.txt', sep = '\t')

mcl_ppi_mod = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_ppi_mod.txt')
n_prot_mcl = as.data.frame(table(table(mcl_ppi_mod$mod)))
colnames(n_prot_mcl) = c('n_prot', 'Freq')
n_prot_mcl$n_prot = as.numeric(as.character(n_prot_mcl$n_prot))
n_prot_tab = rbind(n_prot_mcl[n_prot_mcl$n_prot < 10, ],
                   data.frame(n_prot = '> 10', 
                              Freq = sum(n_prot_mcl$Freq[which(n_prot_mcl$n_prot >= 10)])))

p1 = ggplot(n_prot_tab, aes(x = factor(n_prot, 
                                       levels = c(as.character(2:9),'> 10')), 
                            y = Freq, label = Freq)) + 
  geom_bar(stat="identity", fill = 'steelblue', width = 0.7) +
  geom_text(vjust=-0.2, size = 4.5) + 
  labs(x = "Number of proteins", y = "Number of clusters",
       title = paste0('Total number of PPI clusters = ', max(mcl_ppi_mod$mod))) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

overlap_prot_mcl = as.data.frame(table(table(mcl_ppi_mod$prot)))
colnames(overlap_prot_mcl) = c('n_prot', 'Freq')
overlap_prot_mcl$n_prot = as.numeric(as.character(overlap_prot_mcl$n_prot))
overlap_prot_tab = rbind(overlap_prot_mcl[overlap_prot_mcl$n_prot < 10, ],
                         data.frame(n_prot = '> 10', 
                                    Freq = sum(overlap_prot_mcl$Freq[which(overlap_prot_mcl$n_prot >= 10)])))

p2 = ggplot(overlap_prot_tab, aes(x = factor(n_prot, 
                                             levels = c(as.character(1:9),'> 10')), 
                                  y = Freq, label = Freq)) + 
  geom_bar(stat="identity", fill = 'steelblue', width = 0.7) +
  geom_text(vjust=-0.2, size = 4.5) + 
  labs(x = "Number of clusters", y = "Number of proteins",
       title = paste0('Total number of uniq protein = ', length(unique(mcl_ppi_mod$prot)))) + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
grid.arrange(p1, p2, nrow = 1)


