print(sprintf('Processing %s ppi pertubation', cell_line))
print('read PPI')
cl_interactions <- interactions ### Base PPI loaded in pipeline
head(cl_interactions)

if (tag != 'base') {
  weights_raw <- read.csv(sprintf('%s/weights/%s_weights.csv', data_dir, tag))
  print(head(weights_raw))
  
  print('melt data')
  print(cell_line)
  
  weights <- melt(weights_raw, id.vars = 'Name') %>%
    select(gene=Name, cl=variable, weight=value) %>%
    filter(cl==cell_line) 
  print(head(weights))
  
  print('Create reverse interactions')
  cl_interactions_rev <- data.frame('gene1' = cl_interactions$gene2, 'gene2' = cl_interactions$gene1)
  cl_interactions <- rbind(cl_interactions, cl_interactions_rev) %>%
                  distinct()
  
  cl_interactions <- cl_interactions %>%
    left_join(weights, by=c('gene1'='gene')) %>%
    distinct()
  
  if(tag == 'tt_tanh' | tag == 'raw_all_tanh'){ 
    cl_interactions <- cl_interactions %>%
      filter(weight < 1)
  }

} else {
  cl_interactions <- cl_interactions %>%
    mutate(weight = 1, cl=cell_line)
}


print(sprintf('Saving ppi with weights for %s', cell_line))   
print(head(cl_interactions))
write.table(cl_interactions, sprintf('%s/interactions/processed/%s_ppi_%s.csv', data_dir, cell_line, tag), sep=',', row.names = F)

