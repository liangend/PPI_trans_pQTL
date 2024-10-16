library(igraph)
library(data.table)
setwd('/project/xuanyao/jinghui')
ukb_prot = fread('pqtl/05_h2/00_ref/ukb_prot_w_coor.txt')

prot_mod = fread('pqtl/11_prot_complex/ppi_input/mcl_result/mcl_mod_all_prot.txt')
prot_plot = prot_mod$prot[which(prot_mod$mod == 'corum_1219')]
complex_topology = data.frame(from = 'Nop56p-associated \n pre-rRNA \n protein complex',
                              to = unlist(strsplit(prot_plot, ',', fixed = T)))
complex_node = data.frame(node = c(complex_plot$from[1], complex_plot$to))
complex_node$is_ukb = (complex_node$node %in% ukb_prot$gene_name)
complex_node$color = 'lightblue'
complex_node$color[which(complex_node$is_ukb)] = 'orange'
complex_node$color[1] = 'pink'
complex_node$color[which(complex_node$node %in% c('RPL6', 'RPL7A'))] = 'red'

node.size = setNames(c(40, rep(15, nrow(complex_node)-1)), complex_node$node)

complex_g = graph.data.frame(complex_plot, vertices = complex_node, directed = F)
g = simplify(complex_g)
V(g)$color = complex_node$color

plot(g, vertex.label.cex = 0.5, vertex.size=as.matrix(node.size))

