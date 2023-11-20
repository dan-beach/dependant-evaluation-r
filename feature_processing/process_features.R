# "Process features" function
# 1) Reads the PPI interactions for this cell line
# 2) Builds the PPI graph network
# 3) Extracts topological features for each node
#cell_line <- 'HCC1428_BREAST'
#tag = 'raw_2'

process_features <- function(cell_line, data_dir) {
  
  print(sprintf('Processing features for %s', cell_line))

  # read the personalised ppi for this cell line
  print('read PPI')
  cl_interactions <- read.csv(sprintf('%s/interactions/processed/%s_ppi_%s.csv', data_dir, cell_line, tag), stringsAsFactors = F)
  head(cl_interactions)

if(tag == 'raw_all_tanh_gss') {
  #### Just for go sem sim score weights
  print('Load GSS')
  gosemsim_weights <- read.csv(sprintf('%s/weights/raw_all_tanh_gss_pair_weights.csv', data_dir), stringsAsFactors=F) %>%
  		  filter(cl == cell_line) %>%
		  select(-X)
  print(head(gosemsim_weights))
  print('Left join')
  cl_interactions <- cl_interactions %>%
  		     select(-weight) %>%
		     left_join(gosemsim_weights)
  cl_interactions <- cl_interactions %>%
  		     mutate(weight= signif(weight*1000000, 2)) %>%
		     na.omit()
  str(cl_interactions)
  print(head(cl_interactions))
  rm(gosemsim_weights)
}

  # build graph
  print('Building graph')
  
  # get list of genes and assign to 'nodes' 
  nodes <- unique(as.vector(as.matrix(cl_interactions[,c("gene1", "gene2")])))
  
  # build a directed graph from the PPI interacion data for this cell line
  g <- graph_from_data_frame(d=cl_interactions, vertices=nodes, directed=T)
  E(g)$weight = cl_interactions$weight
  rm(cl_interactions)
  
  print(sprintf('nodes remaining %s', length(nodes)))
  
  features <- data.frame(gene=nodes)
  
  # Try inverting some of the weights
  if (tag == 'raw_2') {
  
    # Create a modified graph for calculting some values with inverse weights
    g_modified <- g
    E(g_modified)$weight <- 1 - E(g)$weight
    
    # Calculations features with a mix of original and modifeid (inverted) weights
    print('Extract node-wise features')
    features$betweenness <- betweenness(g, weights = E(g)$weight)
    features$constraint <- constraint(g_modified)
    features$closeness <- closeness(g, weights = E(g)$weight)
    features$coreness <- coreness(g)
    features$degree <- degree(g)
    features$eccentricity <- eccentricity(g)
    features$eigen_centrality <- eigen_centrality(g_modified)$vector
    features$hub_score <- hub_score(g_modified)$vector
    features$neighborhood1.size <- neighborhood.size(g, 1)
    features$neighborhood2.size <- neighborhood.size(g, 2)
    features$neighborhood6.size <- neighborhood.size(g, 6)
    features$d6_to_d2_neighbours <- features$neighborhood6.size / features$neighborhood2.size 
    features$page_rank <- page_rank(g_modified)$vector
    
  } else {
    
    print('Extract node-wise features')
    features$betweenness <- betweenness(g)
    features$constraint <- constraint(g)
    features$closeness <- closeness(g)
    features$coreness <- coreness(g)
    features$degree <- degree(g)
    features$eccentricity <- eccentricity(g)
    features$eigen_centrality <- eigen_centrality(g)$vector
    features$hub_score <- hub_score(g)$vector
    features$neighborhood1.size <- neighborhood.size(g, 1)
    features$neighborhood2.size <- neighborhood.size(g, 2)
    #features$neighborhood5.size <- neighborhood.size(g, 5)
    features$neighborhood6.size <- neighborhood.size(g, 6)
    #features$d6neighbours <- features$neighborhood6.size - features$neighborhood5.size
    features$d6_to_d2_neighbours <- features$neighborhood6.size / features$neighborhood2.size
    # Other neighbourhood stuff
    features$page_rank <- page_rank(g)$vector

  }

  print('Label genes')
  
  # If the current cell line is in the list of training cell lines, 
  # get the dependencies df and: filter for the current cell line, and dependency score > 0.65, and remove NAs
  # Then convert to a list of the dependent gene_ids for this cell line
  if (cell_line %in% train_cell_lines) {
    print('In training cell lines!')
  c_dep <- dependencies %>%
    filter(cl == cell_line) %>%
    filter(dependency_p > 0.65) %>%
    na.omit()
  c_dep <- as.vector(unlist(c_dep$gene_id))
  print(head(c_dep))
  print(head(features))
  print(str(features))
  
  # add a 'dependent' column to the features df
  # if the gene column value is in the list of dependenrt genes for this cell line, dependent col = dependent, otherwise = non_dependent
  features <- features %>%
    mutate(dependent = case_when(gene %in% c_dep ~ 'dependent', TRUE ~ 'non_dependent'))
  
  print(table(features$dependent))
  } else {
    # current cell line isn't in the training set so create a dependent column and set value to zero
    print('Not in training cell lines')
    features <- features %>%
      mutate(dependent=0)
  }
  
  print('Done')
  
  
  #print('Add some data to graph and save for visualisation')
  #V(g)$pagerank <- features$page_rank
  #V(g)$dependent <- features$dependent
  #V(g)$Label <- V(g)$name
  #write_graph(g, file=sprintf('%s/results/graph_plots/%s_graph_%s.gml', data_dir, cell_line, tag), format = 'gml')

  return(features)
  
}

# Call the 'process_features' function
# Define the output file and write features to file
features <- process_features(cell_line, data_dir)
output_features_file <- sprintf("%s/training/%s_features_%s.csv", data_dir, cell_line, tag)
write.table(features, file=output_features_file, row.names = F, quote = F, sep='\t')


