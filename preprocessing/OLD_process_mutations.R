
mutations <- read.csv(sprintf('%s/mutations/%s_mutations.csv', data_dir, cell_line))

# Add gene IDs based on protein ID
gene_symbol_map <- id_map %>% select(c(Gene.ID, Associated.Gene.Name))

mutations <- mutations %>%
  left_join(gene_symbol_map, by=c('Gene'='Associated.Gene.Name')) %>%
  select(gene = Gene.ID) %>%
  na.omit() %>%
  distinct()

mutations <- mutations %>%
  mutate(loss_gain=0)

write.table(mutations, sprintf('%s/mutations/processed/%s_mutations.csv', data_dir, cell_line), sep=',', row.names = F)

