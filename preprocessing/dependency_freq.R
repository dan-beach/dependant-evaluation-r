
dependencies$disease <-  gsub('^[^_]*_', '\\1', dependencies$cl)

d <- dependencies %>%
  na.omit() %>%
  filter(dependency_p > 0.65) %>%
  filter(cl %in% train_cell_lines) %>%
  group_by(cl, gene_id) %>%
  dplyr::summarise(count=n()) %>%
  arrange(desc(count))
d

g <- d %>%
  group_by(gene_id) %>%
  dplyr::summarise(occurrences_in_cell_lines=n()) %>%
  arrange(desc(occurrences_in_cell_lines))

g

other_genes <- data.frame(gene_id= c(interactions$gene1, interactions$gene2)) %>%
  distinct() %>%
  filter(!gene_id %in% g$gene_id) %>%
  mutate(occurrences_in_cell_lines = 0)
g <- rbind(g, other_genes)

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
