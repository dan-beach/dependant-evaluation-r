# Processes testing data (gene features) for cell lines in test_cell_lines var
# Basically the same as the process_training_data.R file, but doesn't split data into 80/20 training/testing subsets
# As with process_training_data.R the data normalisation functions are commented out currently 
# Also includes feature_processing function for feature selection which is currently unused

#cell_line <- 'CAL51_BREAST' # needs to be in 'testy_cell_line'
#tag <- c('base')

node_features <- c('betweenness', 'constraint', 'closeness', 'coreness', 'degree', 'eccentricity', 'eigen_centrality', 'hub_score',  'neighborhood2.size', 'page_rank')
feature_list <- c('gene', node_features)

process_training_data <- function(features) {
  
  print(head(features))
  
  print('preprocess features')
  
  print('Normalise feature data')
  genes <- features %>% select(gene)
  features <- features %>%
    select(-gene, -dependent)
  # preprocessParams <- preProcess(features, method=c("range"))
  # features <- predict(preprocessParams, features)
  
  print('create training data')
  training <- cbind(genes, features)
  
  training <- training %>%
    na.omit()
  
  return(training)   
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
training <- process_training_data(features)
write.table(training, sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)
rm(features)
