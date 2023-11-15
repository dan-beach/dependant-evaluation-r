# DRB assign cell_line and tag to run script outside of pipeline loop
#cell_line <- 'HCC1428_BREAST' # needs to be in 'train_cell_line'
#tag <- c('raw_0')

# Load training data for the current cell line
training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)

#### Randomise labels!
#dep_temp <- training['dependent']
#rand_dep  <- training$dependent
#rand_dep <- sample(as.vector(rand_dep))
#training['dependent'] <- rand_dep
#print(dep_temp==rand_dep)

#print(str(training))
print('Validation')

# get list of genes in the training data set
training_genes <- training %>% select(gene) %>% unlist() %>% as.vector()

# set up control parameters for training  model
# use cross validation to evaluate the model
# summaryFunction=twoClassSummary: summary function for evaluating the model performance  - typically used for binary classification problems. It calculates performance metrics such as sensitivity, specificity, AUC, etc.
# classProbs=T: set to TRUE to calculate class probabilities - required for ROC curve analysis (AUC metric in particular).
# savePredictions = T: This is set to TRUE to save the predictions made during the training process, which can be useful for later analysis.

print('Model being trained...')
control <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                        savePredictions = T)

#str(training)

# This line trains an AdaBoost model - not RF?
# Train a model on training data (minus gene column), using 'dependent' as the target variable
# Use controls define above. TuneLength specifies that 5 different values of the tuning parameter should be evaluated.
# Model selection should be based on the Receiver Operating Characteristic (ROC) curve area (common performance metric for binary classification models).
mdl <- train(dependent ~ ., data=select(training, -gene), method='ada', trControl=control, tuneLength=5, metric="ROC")

print('CV results')
print(mdl)


# preProc=c("center", "scale"))

auc <- c()
acc <- c()
spec <- c()
sens <- c()

# load the list of dependency frequencies for each gene
print('Create cls test set')
cls_dep <- read.table(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir), sep=',', header = T)

# Do this for each number of cell lines specified in the 'count' variable (assigned in multi_pipe.R)
# Currently set to 1, 10, 20, 30, 39
# This tests model performance when predicting genes that are dependent in X number of cell lines

for(cell_line_count in counts) {
  print(sprintf('%s test set', cell_line_count))
  #print('cell line specific validation')

# Initialise variables for ROC metrics  
  auc <- c()
  acc <- c()
  spec <- c()
  sens <- c()

  #print(train_cell_lines )

  # For each cell line in the train_cell_lines list
  # Get testing data file contents and assign to 'testing'
  for(test_cell_line in train_cell_lines) {
    print(sprintf('validating against  %s', test_cell_line ))
    testing <- read.table(sprintf('%s/training/%s_testing_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
  
    # Filter for genes that are dependent in a certain number of cell lines
    # DRB - get list of genes from dependency frequencies file where X number of cell lines are dependent on the gene 
    # Assign list to cls_genes 
    cls_genes <- cls_dep %>%
    filter(occurrences_in_cell_lines <= cell_line_count) %>%
    select(gene_id) %>%
    unlist() %>%
    as.vector()

    # Filter testing data to include only genes that are in cls_genes
    cls_testing <- testing %>%
    		filter(gene %in% cls_genes)

  #print(str(testing))
  ### Remove genes in test set that are present in training set
  cls_testing <- cls_testing %>% 
    filter(!(gene %in% training_genes)) %>% 
    select(-gene) 
  #print(str(cls_testing))
  #print('performing predictions for validation...')
  
  
  
  # Make predictions on the new data (cls_testing) using the AdaBoost model trained above
  # type="prob": specifies the type of prediction. When set to "prob", the function will return 
  # the predicted class probabilities instead of the class labels. For a binary classification problem, 
  # this will be a two-column data frame with the probability of each observation being in the first and second class, respectively.
    validation <- predict(mdl, newdata=cls_testing, type="prob")
  #print('done validation predictions!')
  
  # add 'pred' column to validation dataframe
  # if predicted dependent prob > 0.5 then set the value to 'dependent', otherwise, set to non-dependent
  validation <- validation %>% 
    mutate(pred = case_when(dependent > 0.5 ~ 'dependent', 
                            TRUE ~ 'non_dependent'))
  
  
  # If there are two unique predicted classes and two unique actual classes:
  # Calculate the confusion matrix, comparing the actual outcomes (cls_testing$dependent) with the predicted outcomes (validation$pred). The results are stored in cm.
  # Calculate the Receiver Operating Characteristic (ROC) curve and the Area Under the Curve (AUC) for the model’s predictions. The AUC value is then appended to an existing vector auc.
  # Extract the accuracy, specificity, and sensitivity from the confusion matrix and append them to existing vectors acc, spec, and sens.
  
  #print(validation$pred)
  #print(testing$dependent)
  if (length(unique(as.factor(validation$pred))) == 2 & length(unique(as.factor(cls_testing$dependent))) == 2) {

     cm <- confusionMatrix(as.factor(cls_testing$dependent), as.factor(validation$pred))
     roc <- roc(cls_testing$dependent, validation$dependent)
  
	  auc <- c(roc$auc, auc)
  	acc <- c(cm$overall['Accuracy'], acc)
  	spec <- c(cm$byClass['Specificity'], spec)
  	sens <- c(cm$byClass['Sensitivity'], sens)
  
  	# If there are not two unique predicted classes or two unique actual classes:
  	# Assign default values to the vectors auc, acc, spec, and sens. 
  	# These default values seem to represent a scenario where the model is not able to distinguish between the two classes at all, 
  	# as indicated by an AUC of 0.5 and an accuracy of 0.5. Specificity is set to 1 and sensitivity to 0, indicating that all predictions are of a single class.
	} else {

  	auc <- c(0.5, auc)
  	acc <- c(0.5, acc)
  	spec <- c(1, spec)
  	sens <- c(0, sens)
  	}

}

# write the results to file 
# location: /results/cell_line_results/{cell_line}_rocs_cls_{tag}_{cell_line_count}.csv
### 
print(sprintf('Saving results for %s - %s - %s', tag, cell_line, cell_line_count))
write.table(data.frame('cl'=train_cell_lines , 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec), sprintf('%s/results/cell_line_results/%s_rocs_cls_%s_%s.csv', data_dir,  cell_line, tag, cell_line_count), row.names=F, sep=',')
print('Done!')

}