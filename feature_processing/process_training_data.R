# Assign cell_line and tag to run script outside of pipline loop

#cell_line <- 'HCC1428_BREAST' # needs to be in 'train_cell_line'
#tag <- c('base')

# Set.seed(0) makes any randomised functions reproducible (but was initially commented out here). 
# Later on random subsets of the non_dependent genes are selected. Without this, will a different random subset be selected each time?

set.seed(0)

# Generate a list of features (to be used as column names)
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
  
  # seperate the genes, features and labels columns
  # apply preprocessing to the features columns (commented out here)
  # recombine the genes, features and labels columns
  # (this step doing nothing currently due to commented out code, I think)
  print('Normalise feature data')
  genes <- features %>% select(gene)
  labels <- features %>% select(dependent)
  features <- features %>%
    select(-gene, -dependent)
  # preprocessParams <- preProcess(features, method=c("range"))
  # features <- predict(preprocessParams, features)

  print('create training data')
  training <- cbind(genes, features, labels)

  # remove any rows with missing values
  training <- training %>%
    na.omit()
  
  # print table showing number of dependent/non-dependent rows (to quickly see class distribution) 
  print(table(training$dependent))
  
  print(str(training))
  
  # Create balanced sets of dependent and non-dependent genes - currently set to do this in 1:1 ratio
  # Filter out all genes labelled as 'non_dependent' and assign to non_dependent
  # Randomly select a subset of the non_dependent dataframe, which contains the same number of rows as 
  # there are 'dependent' labelled genes, and satisying the specified ratio of dependent: non-dependent genes. Assign to 'non_dependent_genes' var.
  # Filter genes labelled as 'dependent' and assign to dependent_genes variable
  # Combine dependent_genes and non_dependent_genes in one dataframe called training and randomise the order 
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
  
  # Use the caret createDataPartition function to split the data into training (20%) and testing (80%) subsets
  # See https://www.rdocumentation.org/packages/caret/versions/6.0-94/topics/createDataPartition
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
  
  # make the dependent columns of the training and testing dataframes factors
  training$dependent <- factor(make.names(training$dependent))
  testing$dependent <- factor(make.names(testing$dependent))

  # return the testing and training data in two variables (testing/training)
  return(list('training'=training, 'testing'=testing))   
}

# Feature_processing function doesn't appear to be used in this file
# Looks like it's performing the leave one out analysis (recursive feature elimination/RFE)
# to determine the importnace of each feature.
# Was this previously used for feature selection?
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

# Open the training features csv file for this cell line
# pass the features to the 'process_training_data' function
# Assign output to 'feature_data'
# Assign "training" and "testing" dataframes (returned by function) to training/testing variables

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

# write training and testing dataframes to files
write.table(training, sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)
write.table(testing, sprintf('%s/training/%s_testing_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)
#write.table(testing, sprintf('%s/training/%s_testing_cls_%s.csv', data_dir, cell_line, tag), row.names = F, sep='\t', quote = F)
rm(features)