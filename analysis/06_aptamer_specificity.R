library(ggplot2)
library(openxlsx)
library(data.table)
## aptamer specificity assessed by mass spectrometry (Emilsson et al., 2018)
ms_val = data.frame(n_aptamer = c(734, 102, 773, 734-5, 102-2,  773-5), 
                    method = rep(c('DDA', 'MMR', 'Either'), 2), 
                    validation = rep(c('total', 'true'), each = 3))

ggplot(ms_val, aes(x = reorder(method, -n_aptamer), y = n_aptamer, fill = validation)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  geom_text(aes(label = n_aptamer), vjust=-0.3, size=4, position = position_dodge(0.5)) + 
  labs(x = "method", y = "# aptamers") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## aptamer specificity assessed by Sun et al., 2018
react = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2018_Sun_supplement.xlsx', sheet = 3, startRow = 3)

react_sum = as.data.frame(table(react$Binding.result, react$Product.of.same.gene))
colnames(react_sum) = c('result', 'same_gene', 'n_pair')
ggplot(react_sum, aes(x = reorder(result, -n_pair), y = n_pair, fill = same_gene)) + 
  geom_bar(stat="identity", width = 0.3) +
  scale_fill_brewer(palette="Paired") + 
  #geom_text(aes(label = n_pair), vjust=-1, size=3.5) + 
  labs(x = "", y = "# pairs") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

react_comp_bind = react[react$Binding.result == 'Comparable binding observed', ]
length(unique(react_comp_bind$SOMAmer.ID))
table(react_comp_bind$Product.of.same.gene)

problem_apatamer = react_comp_bind[react_comp_bind$Product.of.same.gene == 'No', ]
problem_somamer = unique(problem_apatamer$SOMAmer.ID)
problem_somamer = sub('-', '_', problem_somamer)

gud_prot_info = read.xlsx('/project/xuanyao/jinghui/pqtl/00_ref/2022_Gudjonsson_prot_info.xlsx', 
                          startRow = 3)
sum(gud_prot_info$Study.tag %in% problem_somamer)



