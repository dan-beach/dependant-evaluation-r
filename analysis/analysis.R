# Assign cell_line and tag to run script outside of pipeline loop
# cell_line <- 'HCC1428_BREAST' # needs to be in 'train_cell_line'
# tag <- c('base')
set.seed(0)

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
# use cross validation to evaluate the model (default is 10 folds as number not specified here)
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

# New metrics added (not available in useful form for CV)
prec <- c()
f1 <- c()
npv <- c()
pr_auc <- c()

head(mdl)
mdl$results

# ##########
# # Get precision and recall metrics for CV
# # NOTE: it is not possible to get this. You can extract Precision and Recall for ALL of 
# # the cross-validation folds and parameter combinations, but not just for the best model. 
# # This isn't very useful so I won't save it. For now I will just collect the ROC metrics for 
# # the best model (i.e. parameter combination that produced the highest ROC AUC) 
# ##########
# 
# # Access the saved predictions
# predictions <- mdl$pred
# 
# # Calculate precision and recall for each fold using mdl$pred
# new_metrics <- mdl$pred %>%
#   group_by(Resample) %>%
#   summarise(Precision = posPredValue(as.factor(obs), as.factor(pred)),
#             Recall = sensitivity(as.factor(obs), as.factor(pred)))
# 
# # Viewing the metrics
# print(new_metrics)
# 
# # This is what is added to csv file in pipe
# max(mdl$results$ROC)
# 
# ##########

# load the list of dependency frequencies for each gene
print('Create cls test set')
cls_dep <- read.table(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir), sep=',', header = T)

# Do this for each number of cell lines specified in the 'count' variable (assigned in multi_pipe.R)
# E.g. 1, 10, 20, 30, 39
# This tests model performance when predicting genes that are dependent in X number of cell lines

for(cell_line_count in counts) {
  print(sprintf('%s test set', cell_line_count))
  #print('cell line specific validation')

# Initialise variables for ROC metrics  
  auc <- c()
  acc <- c()
  spec <- c()
  sens <- c()
  
# Other metrics 
  prec <- c()
  f1 <- c()
  npv <- c()
  pr_auc <- c()
  
  #print(train_cell_lines )

  # For each cell line in the train_cell_lines list
  # Get testing data file contents and assign to 'testing'
  
  for(test_cell_line in train_cell_lines) {
    
    print(sprintf('validating against  %s', test_cell_line ))
    testing <- read.table(sprintf('%s/training/%s_testing_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
  
    # Filter for genes that are dependent in a certain number of cell lines
    # Get list of genes from dependency frequencies file where X number of cell lines are dependent on the gene 
    # Assign list to cls_genes 
    cls_genes <- cls_dep %>%
    filter(occurrences_in_cell_lines <= cell_line_count) %>%
    select(gene_id) %>%
    unlist() %>%
    as.vector()

    # Filter testing data to include only genes that are in cls_genes
    cls_testing <- testing %>%
    		filter(gene %in% cls_genes)

    # Remove genes in test set that are present in training set
    # Beacuse of how the test/train data are split this will not remove anything when testing the cell line against itself 
    cls_testing <- cls_testing %>% 
      filter(!(gene %in% training_genes)) %>% 
      select(-gene)
  
    
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
      
      # Create cm, roc and prroc
      cm <- confusionMatrix(as.factor(cls_testing$dependent), as.factor(validation$pred))
      roc <- roc(cls_testing$dependent, validation$dependent)
      prc <- pr.curve(scores.class0 = validation$dependent, weights.class0=cls_testing$dependent == "dependent", curve=T)
      
      # Plot the ROC and PR curves
      #plot(roc, main="ROC Curve")
      #plot(prc, main="PR Curve")
      
      
      ################################
      
      # Get values of acc, sens, spec, prec, f1, npv, pr_auc. 
      # Acc, sens and spec are are named vectors and won't behave corretcly in the calcs below (esp f1)
      auc_value <- roc$auc
      acc_value <- as.numeric(cm$overall['Accuracy'])
      spec_value <- as.numeric(cm$byClass['Specificity'])
      sens_value <- as.numeric(cm$byClass['Sensitivity'])
      prec_value <- posPredValue(as.factor(validation$pred), as.factor(cls_testing$dependent))
      f1_value <- 2 * (prec_value * sens_value) / (prec_value + sens_value)
      npv_value <- negPredValue(as.factor(validation$pred), as.factor(cls_testing$dependent))
      pr_auc_value <- prc$auc.integral
      
      # Save metrics
      auc <- c(auc, auc_value)
      acc <- c(acc, acc_value)
      spec <- c(spec, spec_value)
      sens <- c(sens, sens_value)
      prec <- c(prec, prec_value)
      f1 <- c(f1, f1_value)
      npv <- c(npv, npv_value)
      pr_auc <- c(pr_auc, pr_auc_value)
      
      ################################
      
      
  #     # ROC derived metrics
  # 	  auc <- c(roc$auc, auc)
  #   	acc <- c(cm$overall['Accuracy'], acc)
  #   	spec <- c(cm$byClass['Specificity'], spec)
  #   	sens <- c(cm$byClass['Sensitivity'], sens)
  #   	
  #   	# Get numeric value of sens. It is named vector and won't behave corretdly in the f1 calcs below
  #   	sens_value <- as.numeric(cm$byClass['Sensitivity'])
  #   	
  #   	# New PR derived metrics for evaluation
  #   	prec <- c(posPredValue(as.factor(validation$pred), as.factor(cls_testing$dependent)), prec)
  #   	f1 <- c((2 * ((prec * sens_value) / (prec + sens_value))), f1)
  #   	npv <- c(negPredValue(as.factor(validation$pred), as.factor(cls_testing$dependent)), npv)
  #   	pr_auc <- c(prc$auc.integral, pr_auc)
    	

    	#########
    	# Generate plots for ROC and PR and save to file with data
    	#########
    	
    	# Define the interval for breaks
    	# Used for both plots
    	interval = 0.2
    	breaks = seq(0, 1, interval)
    	    	
    	#####
    	# ROC
    	#####
    	
    	# Extract sens and spec values
    	specificity_vals <- roc$specificities
    	sensitivity_vals <- roc$sensitivities
    	
    	# Creating a dataframe for ggplot
    	roc_data <- data.frame(Specificity = specificity_vals, Sensitivity = sensitivity_vals)
    	
    	roc_plot <- ggplot() +
    	  geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) +
    	  geom_step() +
    	  geom_path(aes(x=specificity_vals, y=sensitivity_vals)) +
    	  scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) +
    	  scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    	  ggtitle(sprintf('Dependency validation')) +
    	  annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC: ",sprintf("%.3f",auc_value)))
    	
    	print(roc_plot)
    	
    	#####
    	# PR
    	#####
    	
    	# Extracting PR curve data
    	pr_curve_data <- prc$curve
    	
    	# Assuming the first column is Recall, the second column is Precision
    	recall_vals <- pr_curve_data[, 1]
    	precision_vals <- pr_curve_data[, 2]
    	
    	# Creating a dataframe for ggplot
    	pr_data <- data.frame(Recall = recall_vals, Precision = precision_vals)
    	
    	# Creating the PR plot using ggplot2
    	pr_plot <- ggplot(pr_data, aes(x = Recall, y = Precision)) +
    	  geom_path() +
    	  scale_x_continuous(name = "Recall", limits = c(0, 1), breaks = breaks, expand = c(0.001, 0.001)) +
    	  scale_y_continuous(name = "Precision", limits = c(0, 1), breaks = breaks, expand = c(0.001, 0.001)) +
    	  ggtitle('Precision-Recall Curve') +
    	  annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC: ",sprintf("%.3f", pr_auc_value))) +
    	  theme_minimal()
    	
    	# Display the plot
    	print(pr_plot)
    	
    	#####
    	# Save the plots and data
    	#####
    	
    	print(sprintf('saving %s plots', test_cell_line ))
    	
    	# Save ROC plot and data
    	ggsave(sprintf('%s/results/cell_line_results/plots/%s_%s_roc_%s_%s.pdf', data_dir, cell_line, test_cell_line, tag, cell_line_count), roc_plot)
    	write.csv(roc_data, sprintf("%s/results/cell_line_results/plot_data/%s_%s_roc_%s_%s_data.csv", data_dir, cell_line, test_cell_line, tag, cell_line_count), row.names = FALSE)
    	
    	# Save PR plot and data
    	ggsave(sprintf('%s/results/cell_line_results/plots/%s_%s_pr_%s_%s.pdf', data_dir, cell_line, test_cell_line, tag, cell_line_count), pr_plot)
    	write.csv(pr_data, sprintf("%s/results/cell_line_results/plot_data/%s_%s_pr_%s_%s_data.csv", data_dir, cell_line, test_cell_line, tag, cell_line_count), row.names = FALSE)
    	
	} else {
	  
  	  # If there are not two unique predicted classes or two unique actual classes:
  	  # Assign default values to the vectors auc, acc, spec, and sens. 
  	  # These default values seem to represent a scenario where the model is not able to distinguish between the two classes at all, 
  	  # as indicated by an AUC of 0.5 and an accuracy of 0.5. Specificity is set to 1 and sensitivity to 0, indicating that all predictions are of a single class.
    	
    	auc <- c(auc, 0.5)
    	acc <- c(acc, 0.5)
    	spec <- c(spec, 1)
    	sens <- c(sens, 0)
    	prec <- c(prec, 0.5)
    	f1 <- c(f1, 0.5)
    	npv <- c(npv, 0.5)
    	pr_auc <- c(pr_auc, 0)
  }

    # Debugging
    print(train_cell_lines)
    print(auc)
    print(acc)
    print(spec)
    print(sens)
    print(prec)
    print(f1)
    print(npv)
    print(pr_auc)
}

# Debugging  
print(train_cell_lines)
print(auc)
print(acc)
print(spec)
print(sens)
print(prec)
print(f1)
print(npv)
print(pr_auc)
  
  
# Write the results to file 
# Location: /results/cell_line_results/{cell_line}_rocs_cls_{tag}_{cell_line_count}.csv
print(sprintf('Saving results for %s - %s - %s', tag, cell_line, cell_line_count))
write.table(data.frame('cl'=train_cell_lines , 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec, 'prec'=prec, 'f1'=f1, 'npv'=npv, 'pr_auc'=pr_auc), sprintf('%s/results/cell_line_results/%s_rocs_cls_%s_%s.csv', data_dir,  cell_line, tag, cell_line_count), row.names=F, sep=',')
print('Done!')

}

# Remove some variables that are used again later
rm(list=c('auc', 'acc', 'spec', 'sens', 'prec', 'f1', 'npv', 'pr_auc', 'auc_value', 'acc_value', 'spec_value', 'sens_value', 'prec_value', 'f1_value', 'npv_value', 'pr_auc_value'))
rm(list=c('interval', 'breaks', 'specificity_vals', 'sensitivity_vals', 'roc_data', 'roc_plot', 'pr_curve_data', 'recall_vals', 'precision_vals', 'pr_data', 'pr_plot'))