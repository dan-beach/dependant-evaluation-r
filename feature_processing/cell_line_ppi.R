## vars below used to run file outside of pipline.R loop
#cell_line <- 'HCC1428_BREAST'
#tag = "raw_0"
  
print(sprintf('Processing %s ppi pertubation', cell_line))
print('read PPI')
cl_interactions <- interactions ### Base PPI loaded in pipeline
head(cl_interactions)
  
# if not using the base ppi (tag == base), load the correct weights file
if (tag != 'base' && tag != 'reactome_base') {
  weights_raw <- read.csv(sprintf('%s/weights/%s_weights.csv', data_dir, tag))
  print(head(weights_raw))
  
  print('melt data')
  print(cell_line)

# Convert raw weights matrix to long format
# rename columns gene, cl, weight
  weights <- melt(weights_raw, id.vars = 'Name') %>%
    select(gene=Name, cl=variable, weight=value) %>%
    filter(cl==cell_line) 
  print(head(weights))
  
# Only add the reverse interactions if not reactome
  if (!grepl("reactome", tag)) {
    
    # add the reverse interactions and assign to cl_interactions_rev
    print('Create reverse interactions')
    cl_interactions_rev <- data.frame('gene1' = cl_interactions$gene2, 'gene2' = cl_interactions$gene1)
    cl_interactions <- rbind(cl_interactions, cl_interactions_rev) %>%
      distinct()
  }

# join the weights df onto the cl df using the gene1 & genes columns
# remove duplicates
  cl_interactions <- cl_interactions %>%
    left_join(weights, by=c('gene1'='gene')) %>%
    distinct()

# Filter out interactions where weight = 1
# This will remove any nodes with LOF or where expression data was missing
# Need to adjust conditional statement here (or use one of the tags below) to ensure these nodes are removed
# EDIT - adjusted to always remove nodes unless tag = raw_1 
#  if(tag == 'tt_tanh' | tag == 'raw_all_tanh'){ 
  if(tag != 'raw_1'){     
    cl_interactions <- cl_interactions %>%
      filter(weight < 1)
  }
  
} else {
  cl_interactions <- cl_interactions %>%
    mutate(weight = 1, cl=cell_line)
}
  
# save to file
print(sprintf('Saving ppi with weights for %s', cell_line))   
print(head(cl_interactions))
write.table(cl_interactions, sprintf('%s/interactions/processed/%s_ppi_%s.csv', data_dir, cell_line, tag), sep=',', row.names = F)
