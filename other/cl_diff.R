
dependencies <- read.table(sprintf('%s/dependency/processed/dependencies.csv', data_dir), stringsAsFactors = F)


dependencies <- dependencies %>% filter(dependency_p > 0.65)
cl1 <- dependencies %>% filter(cell_line == 'HCC1395_BREAST') 
cl2 <- dependencies %>% filter(cell_line == 'CAKI2_KIDNEY') 

total_mean <- mean(nrow(cl1), nrow(cl2))

shared <- length(intersect(cl1$gene_id, cl2$gene_id))

shared/total_mean
