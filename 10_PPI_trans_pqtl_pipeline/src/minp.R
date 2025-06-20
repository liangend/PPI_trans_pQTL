library(data.table)
#setwd('/project/xuanyao/jinghui/pqtl/03_p/')
input_p = snakemake@input[['file_p']]
output_small_p = snakemake@output[['small_p']]

for (i in 1:length(input_p)) {
  tmp_y=readRDS(input_p[i])
  tmp_y = tmp_y[tmp_y < 1e-6]
  if(!is.null(tmp_y) & length(tmp_y) > 0){
    tmp_p = as.data.table(setNames(tmp_y, paste0(strsplit(input_p[i], '.', fixed = T)[[1]][2], ":", 
                                                 names(tmp_y))), keep.rownames=T)
    fwrite(tmp_p, output_small_p, sep = '\t', append = T)
  }
  if (i %% 500 == 0) {
    print(i)
    gc()
  }
}


