

process_features <- function(cell_line, data_dir) {
  
  print(sprintf('Processing features for %s', cell_line))

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
  nodes <- unique(as.vector(as.matrix(cl_interactions[,c("gene1", "gene2")])))
  g <- graph_from_data_frame(d=cl_interactions, vertices=nodes, directed=T)
  E(g)$weight = cl_interactions$weight
  rm(cl_interactions)
  
  print(sprintf('nodes remaining %s', length(nodes)))
  
  print('Extract node-wise features')
  features <- data.frame(gene=nodes)
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

  print('Label genes')
  
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
  
  features <- features %>%
    mutate(dependent = case_when(gene %in% c_dep ~ 'dependent', TRUE ~ 'non_dependent'))
  
  print(table(features$dependent))
  } else {
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

#Read ID map and patient mutations

features <- process_features(cell_line, data_dir)
output_features_file <- sprintf("%s/training/%s_features_%s.csv", data_dir, cell_line, tag)
write.table(features, file=output_features_file, row.names = F, quote = F, sep='\t')


