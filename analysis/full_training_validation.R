#tag='base'

set.seed(0)

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

# Initialise ROC metrics 
auc <- c()
acc <- c()
spec <- c()
sens <- c()

# New metrics added
prec <- c()
f1 <- c()
npv <- c()
pr_auc <- c()


# For each cell line in train_cell_lines get a list of genes from the testing data and assign to 'testing'
# This is a validation of the model trained on the super traning set - using testing data from the train_cell_lines lines 
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
  
  # some debugging
  print('prc first round')
  print(validation$dependent)
  print(testing$dependent == "dependent")
  print(test_cell_line)
  print('end of debug printout')
  
  # Generate confusion matrix and roc curve
  cm <- confusionMatrix(as.factor(testing$dependent), as.factor(validation$pred))
  roc <- roc(testing$dependent, validation$dependent)
  prc <- pr.curve(scores.class0 = validation$dependent, weights.class0=testing$dependent == "dependent", curve=T)
  
  # Get values of acc, sens, spec, prec, f1, npv, pr_auc. 
  # Acc, sens and spec are are named vectors and won't behave corretcly in the calcs below (esp f1)
  auc_value <- roc$auc
  acc_value <- as.numeric(cm$overall['Accuracy'])
  spec_value <- as.numeric(cm$byClass['Specificity'])
  sens_value <- as.numeric(cm$byClass['Sensitivity'])
  prec_value <- posPredValue(as.factor(validation$pred), as.factor(testing$dependent))
  f1_value <- 2 * (prec_value * sens_value) / (prec_value + sens_value)
  npv_value <- negPredValue(as.factor(validation$pred), as.factor(testing$dependent))
  pr_auc_value <- prc$auc.integral
  
  # Save metrics
  # auc <- c(auc_value, auc)
  # acc <- c(acc_value, acc)
  # spec <- c(spec_value, spec)
  # sens <- c(sens_value, sens)
  # prec <- c(prec_value, prec)
  # f1 <- c(f1_value, f1)
  # npv <- c(npv_value, npv)
  # pr_auc <- c(pr_auc_value, pr_auc)
  
  auc <- c(auc, auc_value)
  acc <- c(acc, acc_value)
  spec <- c(spec, spec_value)
  sens <- c(sens, sens_value)
  prec <- c(prec, prec_value)
  f1 <- c(f1, f1_value)
  npv <- c(npv, npv_value)
  pr_auc <- c(pr_auc, pr_auc_value)
  
  # auc
  # acc
  # spec
  # sens
  # prec
  # f1
  # npv
  # pr_auc
  
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
  ggsave(sprintf('%s/results/plots/%s_super_roc_%s.pdf', data_dir, test_cell_line, tag), roc_plot)
  write.csv(roc_data, sprintf("%s/results/plot_data/%s_super_roc_%s_data.csv", data_dir, test_cell_line, tag), row.names = FALSE)
  
  # Save PR plot and data
  ggsave(sprintf('%s/results/plots/%s_super_pr_%s.pdf', data_dir, test_cell_line, tag), pr_plot)
  write.csv(pr_data, sprintf("%s/results/plot_data/%s_super_pr_%s_data.csv", data_dir, test_cell_line, tag), row.names = FALSE)

  
  # # Original ROC plot code here
  # # generate ROC plot
  # # save to results/{cell_line}_super_roc_{tag}.pdf
  # interval = 0.2
  # breaks = seq(0, 1, interval)
  # roc_plot <- ggplot() +
  #   geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) +
  #   geom_step() +
  #   geom_path(aes(x=roc$specificities, y=roc$sensitivities)) +
  #   scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) +
  #   scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
  #   ggtitle(sprintf('Dependency validation')) +
  #   annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC \t",sprintf("%.3f",roc$auc)))
  # 
  # plot(roc_plot)
  # print(sprintf('saving %s', test_cell_line ))
  # ggsave(sprintf('%s/results/%s_super_roc_%s.pdf', data_dir, test_cell_line, tag), roc_plot)

}

# Save results to file 
# location: results/super_rocs_{tag}.csv
print(mean(auc))
print(data.frame(train_cell_lines, auc))
auc
train_cell_lines
write.table(data.frame('cl'=train_cell_lines, 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec, 'prec'=prec, 'f1'=f1, 'npv'=npv, 'pr_auc'=pr_auc), sprintf('%s/results/super_rocs_%s.csv', data_dir, tag), row.names=F, sep=',')


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

  # New metrics added
  prec <- c()
  f1 <- c()
  npv <- c()
  pr_auc <- c()
    
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
      
      # some debugging
      print('prc second round')
      print(validation$dependent)
      print(testing$dependent == "dependent")
      print(test_cell_line)
      print('end of debug printout')
      
      cm <- confusionMatrix(as.factor(testing$dependent), as.factor(validation$pred))
      roc <- roc(testing$dependent, validation$dependent)
      prc <- pr.curve(scores.class0 = validation$dependent, weights.class0=testing$dependent == "dependent", curve=T)
      
      # Get values of acc, sens, spec, prec, f1, npv, pr_auc. 
      # Acc, sens and spec are are named vectors and won't behave corretcly in the calcs below (esp f1)
      auc_value <- roc$auc
      acc_value <- as.numeric(cm$overall['Accuracy'])
      spec_value <- as.numeric(cm$byClass['Specificity'])
      sens_value <- as.numeric(cm$byClass['Sensitivity'])
      prec_value <- posPredValue(as.factor(validation$pred), as.factor(testing$dependent))
      f1_value <- 2 * (prec_value * sens_value) / (prec_value + sens_value)
      npv_value <- negPredValue(as.factor(validation$pred), as.factor(testing$dependent))
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
      
      # auc
      # acc
      # spec
      # sens
      # prec
      # f1
      # npv
      # pr_auc

      # sens
      # sens_value
      # prec
      # f1
      # 
      # sens <- c()
      # sens_value <- c()
      # prec <- c()
      # f1 <- c()
      
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
      ggsave(sprintf('%s/results/plots/%s_super_roc_%s_%s.pdf', data_dir, test_cell_line, tag, cell_line_count), roc_plot)
      write.csv(roc_data, sprintf("%s/results/plot_data/%s_super_roc_%s_%s_data.csv", data_dir, test_cell_line, tag, cell_line_count), row.names = FALSE)
      
      # Save PR plot and data
      ggsave(sprintf('%s/results/plots/%s_super_pr_%s_%s.pdf', data_dir, test_cell_line, tag, cell_line_count), pr_plot)
      write.csv(pr_data, sprintf("%s/results/plot_data/%s_super_pr_%s_%s_data.csv", data_dir, test_cell_line, tag, cell_line_count), row.names = FALSE)
  
            
      # # Original plot code below
      # 
      # interval = 0.2
      # breaks = seq(0, 1, interval)
      # roc_plot <- ggplot() +
      #   geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) +
      #   geom_step() +
      #   geom_path(aes(x=roc$specificities, y=roc$sensitivities)) +
      #   scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) +
      #   scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
      #   ggtitle(sprintf('Dependency validation')) +
      #   annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC \t",sprintf("%.3f",roc$auc)))
      # 
      # plot(roc_plot)
      # print(sprintf('saving %s', test_cell_line ))
      # ggsave(sprintf('%s/results/%s_super_roc_%s_%s.pdf', data_dir, test_cell_line, tag, cell_line_count), roc_plot)
    
    } else {
      
      auc <- c(auc, 0.5)
      acc <- c(acc, 0.5)
      spec <- c(spec, 1)
      sens <- c(sens, 0)
      prec <- c(prec, 0.5)
      f1 <- c(f1, 0.5)
      npv <- c(npv, 0.5)
      pr_auc <- c(pr_auc, 0)
      
    }
    
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
  
  ### 
  #print(cell_line_count)
  
  print(train_cell_lines)
  print(auc)
  print(acc)
  print(spec)
  print(sens)
  print(prec)
  print(f1)
  print(npv)
  print(pr_auc)
  
  super_rocs <- data.frame('cl'=train_cell_lines , 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec, 'prec'=prec, 'f1'=f1, 'npv'=npv, 'pr_auc'=pr_auc)
  print(super_rocs)
  write.table(super_rocs, sprintf('%s/results/cell_line_results/super_rocs_cls_%s_%s.csv', data_dir, tag, cell_line_count), row.names=F, sep=',')
  print('Done!')

}