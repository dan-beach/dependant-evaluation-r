#tags
#tag='raw_0'

# Initialize a dataframe with cols for feature names 
training = data.frame("gene"=character(), "betweenness"=numeric(), "constraint"=numeric(), "closeness"=numeric(), "coreness"=numeric(), "degree"=numeric(), "eccentricity"=numeric(), "eigen_centrality"=numeric(),  "hub_score"=numeric(), "neighborhood1.size"=numeric(), "neighborhood2.size"=numeric(), "neighborhood6.size"=numeric(), "d6_to_d2_neighbours"=numeric(), "page_rank"=numeric(), "dependent"=factor())

# For each cell line in the training set, open it's training data csv and add it onto the 'training' dataframe
# This will produce training super-set
for(training_cell_line in train_cell_lines) {
  training_temp <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, training_cell_line, tag), sep='\t', header = T)
  training = rbind(training, training_temp)
}

str(training)

# select a subset (not used here)
n = as.integer(nrow(training)*1)
training = sample_n(training, size =n)

print('Validation')
print('Model being trained...')

# Set controls and train model
# Settings are same as those used in analysis.R and use AdaBoost classifier
control <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                        savePredictions = T)
rf <- train(dependent ~ ., data=select(training, -gene), method='ada', trControl=control, tuneLength=5, metric="ROC", preProc=c("center", "scale"))

auc <- c()
acc <- c()
spec <- c()
sens <- c()

# For each cell line in train_cell_lines get a list of genes from the testing data and assign to 'testing'
# This is a validation of the model trained on the super traning set - using testing data form the train_cell_lines lines 
for(test_cell_line in train_cell_lines) {
  print('validation for:')
  print(test_cell_line)
  testing <- read.table(sprintf('%s/training/%s_testing_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
  testing <- testing %>% 
    select(-gene)
  
  # Predict whether the genes in the 'testing' set are dependent or not using the model trained above
  # add a 'pred' column to the results and set to 'dependent' if prob > 0.5 | non_dependent
  # Does not filter genes out of test set if they are also in the training set
  print('performing predictions for validation...')
  validation <- predict(rf, newdata=testing, type="prob")
  print('completed predictions')
  validation <- validation %>% 
    mutate(pred = case_when(dependent > 0.5 ~ 'dependent', 
                            TRUE ~ 'non_dependent'))
  
  # Generate confusion matrix and roc curve
  cm <- confusionMatrix(as.factor(testing$dependent), as.factor(validation$pred))
  roc <- roc(testing$dependent, validation$dependent)
  
  # save metrics 
  auc <- c(roc$auc, auc)
  acc <- c(cm$overall['Accuracy'], acc)
  spec <- c(cm$byClass['Specificity'], spec)
  sens <- c(cm$byClass['Sensitivity'], sens)
  
  # generate ROC plot
  # save to results/{cell_line}_super_roc_{tag}.pdf
  interval = 0.2
  breaks = seq(0, 1, interval)
  roc_plot <- ggplot() +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) +
    geom_step() +
    geom_path(aes(x=roc$specificities, y=roc$sensitivities)) +
    scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    ggtitle(sprintf('Dependency validation')) +
    annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC \t",sprintf("%.3f",roc$auc)))
  
  plot(roc_plot)
  print(sprintf('saving %s', test_cell_line ))
  ggsave(sprintf('%s/results/%s_super_roc_%s.pdf', data_dir, test_cell_line, tag), roc_plot)
}

# Save results to file 
# location: results/super_rocs_{tag}.csv
print(mean(auc))
print(data.frame(train_cell_lines, auc))
auc
train_cell_lines
write.table(data.frame('cl'=train_cell_lines, 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec), sprintf('%s/results/super_rocs%s.csv', data_dir, tag), row.names=F, sep=',')



#####   Grades of cell line specific genes  

# Repeat the analysis to predict subset of genes that are dependent in 1, 10 20, 30 39 cell lines
# Does not filter genes out of test set if they are also in the training set

print('Create cls test set')
cls_dep <- read.table(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir), sep=',', header = T)

for(cell_line_count in c(counts)) {
  print(sprintf('%s test set', cell_line_count))
  #print('cell line specific validation')
  
  auc <- c()
  acc <- c()
  spec <- c()
  sens <- c()
  
  for(test_cell_line in train_cell_lines) {
    
    # Filter for genes in a certain number of cell lines
    
    print('validation for:')
    print(test_cell_line)
    testing <- read.table(sprintf('%s/training/%s_testing_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
    
    if (cell_line_count != 0) {
    cls_genes <- cls_dep %>%
      filter(occurrences_in_cell_lines <= cell_line_count) %>%
      select(gene_id) %>%
      unlist() %>%
      as.vector()
    
    testing <- testing %>%
      filter(gene %in% cls_genes)
    }
    
    testing <- testing %>% 
      select(-gene)
    
    print('performing predictions for validation...')
    validation <- predict(rf, newdata=testing, type="prob")
    print('completed predictions')
    validation <- validation %>% 
      mutate(pred = case_when(dependent > 0.5 ~ 'dependent', 
                              TRUE ~ 'non_dependent'))
    
    # If there are two unique predicted classes and two unique actual classes:
    # Calculate the ROC curve and the Area Under the Curve (AUC)
    # Extract the accuracy, specificity, and sensitivity from the confusion matrix
    # plot the ROC curve and save to /results/{cell_line}_super_roc_{tag}.pdf
    
    if (length(unique(as.factor(validation$pred))) == 2 & length(unique(as.factor(testing$dependent))) == 2) {
    
      cm <- confusionMatrix(as.factor(testing$dependent), as.factor(validation$pred))
      roc <- roc(testing$dependent, validation$dependent)
      
      auc <- c(roc$auc, auc)
      acc <- c(cm$overall['Accuracy'], acc)
      spec <- c(cm$byClass['Specificity'], spec)
      sens <- c(cm$byClass['Sensitivity'], sens)
      
      interval = 0.2
      breaks = seq(0, 1, interval)
      roc_plot <- ggplot() +
        geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) +
        geom_step() +
        geom_path(aes(x=roc$specificities, y=roc$sensitivities)) +
        scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) +
        scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
        ggtitle(sprintf('Dependency validation')) +
        annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC \t",sprintf("%.3f",roc$auc)))
      
      plot(roc_plot)
      print(sprintf('saving %s', test_cell_line ))
      ggsave(sprintf('%s/results/%s_super_roc_%s.pdf', data_dir, test_cell_line, tag), roc_plot)
    } else {
      
      auc <- c(0.5, auc)
      acc <- c(0.5, acc)
      spec <- c(1, spec)
      sens <- c(0, sens)
    }
  } 
  
  ### 
  print(cell_line_count)
  super_rocs <- data.frame('cl'=train_cell_lines , 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec)
  print(super_rocs)
  write.table(super_rocs, sprintf('%s/results/cell_line_results/super_rocs_cls_%s_%s.csv', data_dir, tag, cell_line_count), row.names=F, sep=',')
  print('Done!')

}
