
training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)

#### Randomise labels!
#dep_temp <- training['dependent']
#rand_dep  <- training$dependent
#rand_dep <- sample(as.vector(rand_dep))
#training['dependent'] <- rand_dep
#print(dep_temp==rand_dep)

#print(str(training))
print('Validation')
training_genes <- training %>% select(gene) %>% unlist() %>% as.vector()


print('Model being trained...')
control <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                        savePredictions = T)


#str(training)

mdl <- train(dependent ~ ., data=select(training, -gene), method='ada', trControl=control, tuneLength=5, metric="ROC")

print('CV results')
print(mdl)


# preProc=c("center", "scale"))

auc <- c()
acc <- c()
spec <- c()
sens <- c()


print('Create cls test set')
cls_dep <- read.table(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir), sep=',', header = T)

for(cell_line_count in counts) {
  print(sprintf('%s test set', cell_line_count))
  #print('cell line specific validation')

  auc <- c()
  acc <- c()
  spec <- c()
  sens <- c()

  #print(train_cell_lines )

  for(test_cell_line in train_cell_lines) {
    print(sprintf('validating against  %s', test_cell_line ))
    testing <- read.table(sprintf('%s/training/%s_testing_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
  
    # Filter for genes in a certain number of cell lines

    cls_genes <- cls_dep %>%
    filter(occurrences_in_cell_lines <= cell_line_count) %>%
    select(gene_id) %>%
    unlist() %>%
    as.vector()

    cls_testing <- testing %>%
    		filter(gene %in% cls_genes)

  #print(str(testing))
  ### Remove genes in test set that are present in training set
  cls_testing <- cls_testing %>% 
    filter(!(gene %in% training_genes)) %>% 
    select(-gene) 
  #print(str(cls_testing))
  #print('performing predictions for validation...')
  
  
  
  
  validation <- predict(mdl, newdata=cls_testing, type="prob")
  #print('done validation predictions!')
  
  validation <- validation %>% 
    mutate(pred = case_when(dependent > 0.5 ~ 'dependent', 
                            TRUE ~ 'non_dependent'))
  

  #print(validation$pred)
  #print(testing$dependent)
  if (length(unique(as.factor(validation$pred))) == 2 & length(unique(as.factor(cls_testing$dependent))) == 2) {

     cm <- confusionMatrix(as.factor(cls_testing$dependent), as.factor(validation$pred))
     roc <- roc(cls_testing$dependent, validation$dependent)
  
	  auc <- c(roc$auc, auc)
  	acc <- c(cm$overall['Accuracy'], acc)
  	spec <- c(cm$byClass['Specificity'], spec)
  	sens <- c(cm$byClass['Sensitivity'], sens)
  

	} else {

  	auc <- c(0.5, auc)
  	acc <- c(0.5, acc)
  	spec <- c(1, spec)
  	sens <- c(0, sens)
  	}

}

### 
print(sprintf('Saving results for %s - %s - %s', tag, cell_line, cell_line_count))
write.table(data.frame('cl'=train_cell_lines , 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec), sprintf('%s/results/cell_line_results/%s_rocs_cls_%s_%s.csv', data_dir,  cell_line, tag, cell_line_count), row.names=F, sep=',')
print('Done!')

}