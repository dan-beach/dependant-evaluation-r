# Set this here so I can run the file directly instead of inside the loop in multi_pipe.R 
tag <- c('raw_0')

# Original loop generated personalised PPI and processed features all together
# I have seperated these into two different loops below to evaluate more easily

#for(cell_line in all_cell_lines) {  
#  print('Generate features for ')
#  print(cell_line)
#  print('Generate weighted PPI')
#  source('feature_processing/cell_line_ppi.R')
#  print('Processing features')
#  source('feature_processing/process_features.R')
#}

# Generate weighted PPI networks for all the cell lines
for(cell_line in all_cell_lines) {  
  print('Generate weighted PPI for ')#
  print(cell_line)
  source('feature_processing/cell_line_ppi.R')
}

# Build PPI network graphs for each cell line and calculate features
# also labels data when cell line is in training set
for(cell_line in all_cell_lines) {  
  print('Processing features for ')
  print(cell_line)
  source('feature_processing/process_features.R')
}

print('Completed feature generation')

print('performing analysis')

for(cell_line in train_cell_lines) {
  print('######## Training data processing for:')
  print(cell_line)
  source('feature_processing/process_training_data.R')
}

for(cell_line in test_cell_lines) {
  print('######## Training data processing for:')
  print(cell_line)
  source('feature_processing/process_testing_data.R')
}
# Train RF model on each cell line individually and estimate feature importance for each cell line
# Then test train will one cell line and test against all other cell lines one by one, and for each also 
# test perfommace when prediciting rare and common essential genes (1, 10, 20, 30, 39)
# The model trained in analysis.R tests 5 iterations and delcted the best one based on the ROC
# Select the highest ROC here and add to cv_roc for each cell line 
cv_roc <- c()
for(cell_line in train_cell_lines) {
  print('######## Analysis data for:')
  print(cell_line)
  source('analysis/feature_importance_analysis.R')
  source('analysis/analysis.R')
  cv_roc <- c(cv_roc, cell_line = max(mdl$results$ROC))
}

# Write ROC results to file (1 line for each cell line) to file
cv_roc <- data.frame(train_cell_lines, cv_roc)
write.csv(cv_roc, sprintf('%s/results/cv_roc_%s.csv', data_dir, tag))
cv_roc

# Concatenate all training data and retrain the model
# Use new model to predict dependent genes in each cell line
# Repeat to see how well it predicts genes that are essential in 1, 10, 20, 30, 39 cell lines
source('analysis/full_training_validation.R')

# Predict cell essential genes in the 35 unlabeled cell lines (test_cell_lines)
source('analysis/prediction.R')