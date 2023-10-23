
# Download dependency data https://ndownloader.figshare.com/files/10416306 
# and add data to [data_dir]/dependency/gene_dependency.csv

raw_dependencies <- read.csv(sprintf('%s/dependency/gene_dependency.csv', data_dir), stringsAsFactors = F, sep=',') ### https://portals.broadinstitute.org/achilles/datasets/18/download/gene_dependency.csv 
nrow(raw_dependencies)
dependencies_melt <- melt(raw_dependencies, 'line') %>%
  select(cell_line=line, gene=variable, dependency_p=value) %>% 
  distinct()

length(unique(dependencies_melt$cell_line))
nrow(dependencies_melt)

dependencies_melt$disease <-  gsub('^[^_]*_', '\\1', dependencies_melt$cell_line)
nrow(dependencies_melt)

diseases <- c('BREAST', 'KIDNEY', 'PANCREAS')
d <- dependencies_melt %>%
  select(-gene, -dependency_p) %>%
  distinct() %>%
  filter(disease %in% diseases)

str(dependencies_melt)
dependencies_melt <- dependencies_melt %>%
  mutate(gene=gsub("\\..*","",gene))

gene_map <- id_map %>%
  select(Gene.ID, Associated.Gene.Name)

dependencies_melt <- dependencies_melt %>%
  left_join(gene_map, c('gene' = 'Associated.Gene.Name')) %>%
  select(cl=cell_line, gene_id=Gene.ID, gene_name=gene, dependency_p) %>%
  distinct()

write.table(dependencies_melt, sprintf('%s/dependency/processed/dependencies.csv', data_dir))

rm(list = c('dependencies_melt', 'raw_dependencies', 'd', 'gene_map'))







