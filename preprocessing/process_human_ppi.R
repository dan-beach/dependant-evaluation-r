#Read ID map

## Read interaction data
print('Reading and preprocess PPI data')
raw_interactions <- read.table(sprintf('%s/supporting_files/9606.protein.links.detailed.v10.txt', data_dir), sep = ' ', header=T, stringsAsFactors = F)

## preprocess interaction data
# Only keep interactions with experimental evidence > 800
# Remove '9606.' from protein names in the two protein columns
# Select just the protein1 and protein2 cols
raw_interactions <- raw_interactions %>% 
  filter(experimental > 800)%>% 
  mutate(protein1 = gsub('9606.', "", protein1)) %>%
  mutate(protein2 = gsub('9606.', "", protein2)) %>%
  select(c(protein1, protein2))

# DRB see what's in raw_interaction
#head(raw_interactions, n=10)
nrow(raw_interactions)

# Add reverse interactions - string data already bidirectional?
print('Add reverse direction interactions')
vv <- data.frame('protein1'=raw_interactions$protein2, 'protein2'=raw_interactions$protein1)
raw_interactions <- rbind(raw_interactions, vv)

# DRB check content
head(raw_interactions)

# Add gene IDs based on protein ID
str(id_map)
gene_protein_map <- id_map %>% select(c(Gene.ID, Protein.ID))
raw_interactions <- raw_interactions %>%
                left_join(gene_protein_map, by=c('protein1'='Protein.ID')) %>%
                left_join(gene_protein_map, by=c('protein2'='Protein.ID')) %>%
                select(gene1 = Gene.ID.x, gene2 = Gene.ID.y) %>%
                na.omit()

raw_interactions
write.table(raw_interactions, sprintf('%s/interactions/processed/base_ppi.csv', data_dir), sep=',', row.names = F)

rm(list=c('raw_interactions', 'gene_protein_map', 'vv'))
