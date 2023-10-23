library(dplyr)

data_dir <- '/home/g/gb/gb293/phd/data/slorth'


human_predictions <- read.csv('results/recent/human_predictions.csv', sep=' ')
human_ids <- read.csv(sprintf('%s/id_maps/human_id_map.txt', data_dir), sep='\t') %>%
  select(Gene.ID, Associated.Gene.Name)

named_human_predictions <- human_predictions %>%
  left_join(human_ids, c('gene1' = 'Gene.ID')) %>%
  left_join(human_ids, c('gene2' = 'Gene.ID')) %>%
  mutate(gene1 = Associated.Gene.Name.x) %>%
  mutate(gene2 = Associated.Gene.Name.y) %>%
  select(gene1, gene2, consensus) %>%
  filter(consensus>0.75) %>%
  distinct()


genes <- c('PBRM1')

human_interest <- named_human_predictions %>%
  filter(gene1 %in% genes | gene2 %in% genes)

  

nodes <- unique(as.vector(as.matrix(human_interest[,c("gene1", "gene2")])))
interest_graph <- graph_from_data_frame(d=human_interest[,c("gene1", "gene2")], vertices=nodes, directed=F)
E(interest_graph)$weight <- human_interest$consensus
write_graph(interest_graph, 'results/recent/human_interest_smc.gml', format = c("gml"))

human_interest <- human_interest %>%
  arrange(desc(consensus))

human_interest

write.table(human_interest, 'results/recent/humantop_predictions.csv')


####   get drug targets

drug_targets <- read.csv(sprintf('%s/drug_data/druggables.csv',data_dir), stringsAsFactors = F, header = F, col.names = c('uniprot'))
head(drug_targets)
uniprot_key <- read.csv(sprintf('%s/drug_data/hugo_uniprot.csv', data_dir), stringsAsFactors = F, col.names=c('hugo', 'uniprot'))
head(uniprot_key)
drug_targets <- drug_targets %>%
                left_join(uniprot_key)
targets <- drug_targets$hugo

human_druggable <- named_human_predictions %>%
  filter(gene1 %in% targets | gene2 %in% targets)

genes <- c('PBRM1')

human_druggable_interest <- human_druggable %>%
  filter(gene1 %in% genes | gene2 %in% genes) %>%
  arrange(desc(consensus))

human_druggable_interest




