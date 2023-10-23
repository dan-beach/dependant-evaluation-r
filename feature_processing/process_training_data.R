
# set.seed(0)

node_features <- c('betweenness', 'constraint', 'closeness', 'coreness', 'degree', 'eccentricity', 'eigen_centrality', 'hub_score',  'neighborhood2.size', 'page_rank')
#pair_features <- c('cohesion')
feature_list <- c('gene', node_features, 'dependent')

process_training_data <- function(features) {
  
  print('preprocess features')
  # print('select features')
  # # features <- features %>%
  # #   select(one_of(feature_list))
  # print('completed preprocessing features')
  print(str(features))   
  print(table(features$dependent)) 
  
  print('Normalise feature data')
  genes <- features %>% select(gene)
  labels <- features %>% select(dependent)
  features <- features %>%
    select(-gene, -dependent)
  # preprocessParams <- preProcess(features, method=c("range"))
  # features <- predict(preprocessParams, features)

  print('create training data')
  training <- cbind(genes, features, labels)
  
  training <- training %>%
    na.omit()
  print(table(training$dependent))
  
  print(str(training))
  
   print('Undersample to create balanced sets')
   print('make non_sl')
   ratio <- 1 ###############  :1 --- To change ratio between dependent and non-dependent
   non_dependent <- filter(training, dependent=='non_dependent') 
   print('make non_dependent_pairs')
   non_dependent_genes <- sample_n(non_dependent, nrow(filter(training, dependent=='dependent'))*ratio)
   print('make dependent_pairs')
   dependent_genes <- filter(training, dependent=='dependent') 
   print('Make training')
   training <- rbind(non_dependent_genes, dependent_genes) 
   print('randomise order of training data')
   training <- training[sample(nrow(training)),]
  
  
  print('Create holdout test set')
  print(table(training$dependent))
  
  # Make test subset
  trainIndex <- createDataPartition(training$dependent, p = .2, 
                                    list = FALSE, 
                                    times = 1)
  testing  <- training[-trainIndex,]
  training <- training[ trainIndex,]
  
  #training_genes <- training %>% select(gene)
  #testing_genes <- testing %>% select(gene)
  # training <- training %>% select(-gene)
  # testing <- testing %>% select(-gene)
  
  training$dependent <- factor(make.names(training$dependent))
  testing$dependent <- factor(make.names(testing$dependent))

  
  

  return(list('training'=training, 'testing'=testing))   
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

#print('Create cls test set')
#cls_dep <- read.table(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir), sep=',', header = T)
#print(head(cls_dep))

#cls_dep <- cls_dep %>%
#  filter(occurrences_in_cell_lines <= cell_line_count) %>%
#  select(gene_id) %>%
#  unlist() %>%
#  as.vector()

#cls_testing <- testing %>%
#  filter(gene %in% cls_dep)

#print(table(training$dependent))
#print(table(testing$dependent))

write.table(training, sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)
write.table(testing, sprintf('%s/training/%s_testing_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)
#write.table(testing, sprintf('%s/training/%s_testing_cls_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)
rm(features)








