# DRB assign cell_line and tag to run script outside of pipeline loop
#cell_line <- 'HCC1428_BREAST' # needs to be in 'train_cell_line'
#tag <- c('raw_0')

# Read the training data file for this cell line and assign to 'training'
training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)

print('Calculate important features')

# set up control parameters for training rf model: repeated cross validation, 10 folds, 3 repeats
# number = number of folds (k) data are split into for cross-validation. Here k = 10. Model will be trained on 9 and then validated on the one remaining set
# This is repeated 10 times with each fold being used as the validation set exactly once.
# repeats = 3 means the whole cross-validation process is repeated 3 times. This is different to boot-strapping in which a random subset is selected each time
control <- trainControl(method="repeatedcv", number=10, repeats=3)

# Train a random forest model using the 'train' function of the caret package
# defines 'dependent' as the target variable -everything other than 'gene' id a training feature
# 'training' is the data source, method is random forest (rf), feature values are scaled first (between 0-1)
# training contols (trControl) = the 'controls' object defined above
model <- train(dependent~. - gene, data=training, method="rf", preProcess="scale", trControl=control)

# Estimate variable importance
# uses varImp function from caret package to estimate
# scale =  FALSE means do not scale importance scores
importance <- varImp(model, scale=FALSE)

# summarize importance
print(importance)

# plot importance
plot(importance)

# save results to file
fn <- sprintf('%s/results/%s_importance_%s.csv', data_dir, cell_line, tag)
write.table(importance$importance, fn, sep=',', row.names=F)