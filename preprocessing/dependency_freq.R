# DRB check dependencies
head(dependencies)

# add disease column based on cl column (everything after the underscore
dependencies$disease <-  gsub('^[^_]*_', '\\1', dependencies$cl)

# Subset dependencies, omit NA's and filter only rows with dependency_p > 0.65, and where cl name is in train_cell_lines vector
# Then count how many times gene is dependent in the 39 training cell lines
# Switched from 'omit.na()' to 'filter(rowSums(is.na(select(., -gene_id))) == 0)' 
# This drops rows if anything other gene_id is na, rather than droping rows where any col is na
# NOTE: removed part of pipreline that maps gene_id to dependecies df, so omit.na() has been restored
d <- dependencies %>%
  #filter(rowSums(is.na(select(., -gene_id))) == 0) %>% 
  na.omit() %>%
  filter(dependency_p > 0.65) %>%
  filter(cl %in% train_cell_lines) %>%
  #group_by(cl, gene_id) %>%
  group_by(cl, gene_name) %>%
  dplyr::summarise(count=n()) %>%
  arrange(desc(count))
d
  
# Group d by gene_name. Save the count (n) as occurrences_in_cell_lines. Assign to g and sort on count column decending.
# Returns list of gene_ids with a count of how many of the training cell lines they are dependent in
g <- d %>%
  group_by(gene_name) %>%
  dplyr::summarise(occurrences_in_cell_lines=n()) %>%
  arrange(desc(occurrences_in_cell_lines))

g

# get gene1 and gene2 columns from processed base PPI data (interactions) and assign to other_genes
# filter out duplicates
# keep only the genes where the gene_name is not already in g$gene_name
# add occurances_in_cell_lines column and set to zero
# combine this with the g dataframe
other_genes <- data.frame(gene_name= c(interactions$gene1, interactions$gene2)) %>%
  distinct() %>%
  filter(!gene_name %in% g$gene_name) %>%
  mutate(occurrences_in_cell_lines = 0)
g <- rbind(g, other_genes)

g <- g %>% rename(gene_id = gene_name)

g 

# #########
# # Test: see which genes are in the dep_freq list but not in the ppi genes and vice versa
# # Note: all genes in the PPI are in the dep_freq list, but not viec vera
# # Checke din the Pyhton weights pipeline and this is because they are in exp and dep, but
# # never in the PPI to ebgin with 
# ########
# dep_freq_genes <- g %>%
#   select(gene_id) %>%
#   distinct() %>%
#   unlist() %>%
#   as.vector()
# 
# ppi_genes <- data.frame(gene_name= c(interactions$gene1, interactions$gene2)) %>%
#   distinct() %>%
#   unlist() %>%
#   as.vector()
# 
# not_in_dep_freq <- ppi_genes[!ppi_genes %in% dep_freq_genes]
# not_in_ppi <- dep_freq_genes[!dep_freq_genes %in% ppi_genes]
# 
# not_in_dep_freq
# not_in_ppi

# # Test: How many dependent genes are there for 'HCC1428_BREAST' in the depmap data?
# # Check how many end up in the dependency frequency file
# dep_test <- dependencies %>%
#   filter(cl == 'HCC1428_BREAST') %>%
#   filter(dependency_p > 0.65)
# 
# dep_test

write.table(g, sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir), sep=',', row.names = F)

# ggplot(g) +
#   aes(occurrences_in_cell_lines) +
#   geom_histogram(bins=42, color='cornflowerblue', fill=NA)
# 
# g %>% filter(occurrences_in_cell_lines <= 1 )
# 
# table(g$occurrences_in_cell_lines)
# 
# # cl <- dependencies %>%
# #   na.omit() %>%
# #   filter(dependency_p > 0.75) %>%
# #   group_by(cell_line) %>%
# #   summarise(count=n()) %>%
# #   arrange(desc(count))
# # cl
# 
# g
# 
#
# sum(g$gene_id %in% other_genes$gene)
# #raw_dependencies <- read.csv(sprintf('%s/dependency/gene_dependency.csv', data_dir), stringsAsFactors = F) ### https://portals.broadinstitute.org/achilles/datasets/18/download/gene_dependency.csv 

rm(list=c('d', 'g', 'other_genes'))
