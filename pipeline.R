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

cv_roc <- c()
for(cell_line in train_cell_lines) {
  print('######## Analysis data for:')
  print(cell_line)
  source('analysis/feature_importance_analysis.R')
  source('analysis/analysis.R')
  cv_roc <- c(cv_roc, cell_line = max(mdl$results$ROC))
}

cv_roc <- data.frame(train_cell_lines, cv_roc)
write.csv(cv_roc, sprintf('%s/results/cv_roc_%s.csv', data_dir, tag))
cv_roc

source('analysis/full_training_validation.R')


source('analysis/prediction.R')