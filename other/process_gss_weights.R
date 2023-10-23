library(dplyr)
library(reshape2)

breast <- read.table('~/phd/data/slant_cancer/weights/gss_weights_breast.csv', sep=',', header=T) %>%
  select(-X)
kidney <- read.table('~/phd/data/slant_cancer/weights/gss_weights_kidney.csv', sep=',', header=T) %>%
  select(-X)

weights <- breast %>%
  left_join(kidney)

head(weights)

weights = melt(weights, id.vars = c('gene1', 'gene2')) 
head(weights)
weights <- weights %>%
  select(gene1, gene2, cl=variable, weight=value)

head(weights)


write.csv(weights, 'gosem_weights.csv')
