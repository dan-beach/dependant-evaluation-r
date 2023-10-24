# DRB check dependencies
head(dependencies)

# add disease column based on cl column (everything after the underscore
dependencies$disease <-  gsub('^[^_]*_', '\\1', dependencies$cl)

# subset dependencies
# omit NAs
# filter only rows with dependency_p > 0.5, and where cl name is in train_cell_lines vector
# count how many times gene is dependent in the 39 training cell lines
d <- dependencies %>%
  na.omit() %>%
  filter(dependency_p > 0.65) %>%
  filter(cl %in% train_cell_lines) %>%
  group_by(cl, gene_id) %>%
  dplyr::summarise(count=n()) %>%
  arrange(desc(count))
d

# Group d by gene_id. Save the count (n) as occurrences_in_cell_lines. Assign to g and sort on count column decending.
# Returns list of gene_ids with a count of how many of the training cell lines they are dependent in
g <- d %>%
  group_by(gene_id) %>%
  dplyr::summarise(occurrences_in_cell_lines=n()) %>%
  arrange(desc(occurrences_in_cell_lines))

g

# get gene1 and gene2 columns from processed base PPI data (interactions) and assign to other_genes
# filter out duplicates
# keep only the genes where the gene_id is not already in g$gene_id
# add occurances_in_cell_lines column and set to zero
# combine this with the g dataframe
other_genes <- data.frame(gene_id= c(interactions$gene1, interactions$gene2)) %>%
  distinct() %>%
  filter(!gene_id %in% g$gene_id) %>%
  mutate(occurrences_in_cell_lines = 0)
g <- rbind(g, other_genes)

g 

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
# 
# 
# sum(g$gene_id %in% other_genes$gene)
# #raw_dependencies <- read.csv(sprintf('%s/dependency/gene_dependency.csv', data_dir), stringsAsFactors = F) ### https://portals.broadinstitute.org/achilles/datasets/18/download/gene_dependency.csv 

rm(list=c('d', 'g'))
