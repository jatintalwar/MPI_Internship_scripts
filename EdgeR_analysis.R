## edge R analysis script:

library(edgeR)

######################     the count matrix     ################

raw.data <- read.table("count_table_complete_HTSeq_tophat.txt",header=T,row.names=1)
View(raw.data)
colnames(raw.data)
count_matrix <- raw.data
View(count_matrix)
