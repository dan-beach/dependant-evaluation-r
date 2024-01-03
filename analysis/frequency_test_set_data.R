# This runs through the test set construction post-hoc and saves meta deta 
# relating to the number of genes in the test set and their class balance.
# This is to assess whether the results obtained from test sets that have been 
# filtered by frequency can be relied upon

cell_line <- 'HCC1428_BREAST' # needs to be in 'train_cell_lines'
tag <- c('base')

# Set TRUE to get meta data for super analysis test data, FALSE for cell-line specific
super_analysis = TRUE
  
# Initialise some vectors to store data about the test sets
cell_line_vec <- c()
train_cell_line_vec <- c()
cell_line_count_vec <- c()
num_test_genes_vec <- c()
num_dependent_vec <- c()
num_non_dependent_vec <- c()
percent_minority_vec <- c()

for (cell_line in train_cell_lines) {
  
  print(sprintf('Analysing test set for %s', cell_line ))
  
  # Load training data for the current cell line
  training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)
  
  # get list of genes in the training data set
  training_genes <- training %>% select(gene) %>% unlist() %>% as.vector()
  
  # load the list of dependency frequencies for each gene
  cls_dep <- read.table(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir), sep=',', header = T)
  
  # Do this for each number of cell lines specified in the 'count' variable (assigned in multi_pipe.R)
  # Currently set to 1, 10, 20, 30, 39
  # This tests model performance when predicting genes that are dependent in X number of cell lines
  
  for(cell_line_count in counts) {
  
    # For each cell line in the train_cell_lines list
    # Get testing data file contents and assign to 'testing'
    for(test_cell_line in train_cell_lines) {
      
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
      # Don't do this if getting data for the super analysis training sets
      if (super_analysis != TRUE) {
        
        cls_testing <- cls_testing %>% 
          filter(!(gene %in% training_genes)) %>% 
          select(-gene)
      }
      
      # Save some data about the test set
      num_test_genes <- nrow(cls_testing)
      num_dependent <- nrow(subset(cls_testing, dependent == 'dependent'))
      num_non_dependent <- nrow(subset(cls_testing, dependent == 'non_dependent'))
      percent_minority <- (num_dependent / num_test_genes) * 100
      
      # Append the calculated values to their respective vectors
      cell_line_vec <- c(cell_line_vec, cell_line)
      train_cell_line_vec <- c(train_cell_line_vec, test_cell_line)
      cell_line_count_vec <- c(cell_line_count_vec, cell_line_count)
      num_test_genes_vec <- c(num_test_genes_vec, num_test_genes)
      num_dependent_vec <- c(num_dependent_vec, num_dependent)
      num_non_dependent_vec <- c(num_non_dependent_vec, num_non_dependent)
      percent_minority_vec <- c(percent_minority_vec, percent_minority)
      
      print(sprintf('Saving test set meta data for %s - %s - %s', tag, cell_line, cell_line_count))
    }
  }
}

# Create a dataframe from the vectors
test_set_info <- data.frame(cell_line=cell_line_vec,
                        train_cell_line=train_cell_line_vec,
                        cell_line_count=cell_line_count_vec,
                        num_test_genes=num_test_genes_vec,
                        num_dependent=num_dependent_vec,
                        num_non_dependent=num_non_dependent_vec,
                        percent_minority=percent_minority_vec)

# write the results to file 
# location: /results/cell_line_results/cl_frequency_test_set_data_{tag}.csv
###
print('Saving to file')

if (super_analysis != TRUE) {
  write.table(test_set_info, sprintf('%s/results/cell_line_results/test_set_meta_data/cl_frequency_test_set_data_%s.csv', data_dir, tag), row.names=F, sep=',')
} else {
  write.table(test_set_info, sprintf('%s/results/cell_line_results/test_set_meta_data/super_frequency_test_set_data_%s.csv', data_dir, tag), row.names=F, sep=',')
}

print('Done!')