


node_features <- c('betweenness', 'constraint', 'closeness', 'coreness', 'degree', 'eccentricity', 'eigen_centrality', 'hub_score',  'neighborhood2.size', 'page_rank')
feature_list <- c('gene', node_features)

process_training_data <- function(features) {
  
  print('preprocess features')
  print(str(features))   
  
  print('Normalise feature data')
  genes <- features %>% select(gene)
  features <- features %>%
    select(-gene, -dependent)
  #preprocessParams <- preProcess(features, method=c("range"))
  #features <- predict(preprocessParams, features)

  print('create training data')
  training <- cbind(genes, features)
  
  training <- training %>%
    na.omit()
  print(table(training$dependent))
  
  print(str(training))
  
  return(list('training'=training))   
}

feature_processing <- function(training) {
  ### feature selection
  x <- training %>%
    select(-dependent)
  y <- training %>%
    select(dependent)
  control <- rfeControl(functions=rfFuncs, method="cv", number=5)
  results <- rfe(x, y$dependent, rfeControl=control)
  print(results)
  predictors(results)
  plot(results)
  important <- as.vector(unlist(results$optVariables[1:5]))
  return(important)
}


features <- read.table(sprintf('%s/training/%s_features_%s.csv', data_dir, cell_line, tag), sep='\t', header=TRUE, fill=T, stringsAsFactors = F)
print(table(features$dependent))
feature_data <- process_training_data(features)
training <- feature_data[['training']]
head(training)
testing <- feature_data[['testing']]

write.table(training, sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)









