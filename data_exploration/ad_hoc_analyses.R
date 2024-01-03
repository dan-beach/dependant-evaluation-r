#####
# Tests
#####

#####
# 1
#####

# Testing the effect of leaving "occurances_in_cell_lines = 0" in test set when filtering for rarity of gene dependency 

# Get data and set params
cell_line <- 'HCC1428_BREAST' # needs to be in 'train_cell_lines'
tag <- c('base')
cell_line_count = 1

# Get training and testing data for this cell line
HCC1428_training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
HCC1428_testing <- read.table(sprintf('%s/training/%s_testing_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)

# Get list of genes in the dependency data
# Filter out genes that are not dependent in any cell lines (occurrences_in_cell_lines != 0)
# Retain only those dependencies with frequency less than or equal to the cutoff
cls_genes <- cls_dep %>%
  filter(occurrences_in_cell_lines != 0) %>%
  filter(occurrences_in_cell_lines <= cell_line_count) %>%
  select(gene_id) %>%
  unlist() %>%
  as.vector()

# How many genes in traingin and testing sets? 
nrow(HCC1428_training)
nrow(HCC1428_testing)

head(HCC1428_testing)

# Filter testing data to include only genes that are in cls_genes (dependency genes)
cls_testing <- HCC1428_testing  %>%
  filter(gene %in% cls_genes)

# Get list of testing genes
cls_testing_genes <- HCC1428_testing  %>%
  select(gene) %>%
  unlist() %>%
  as.vector()

# Get list of training genes
cls_training_genes <- HCC1428_training  %>%
  select(gene) %>%
  unlist() %>%
  as.vector()

# How many genes in the test set?
nrow(cls_testing)

# Remove genes in test set that are present in training set
cls_testing <- cls_testing %>%
  filter(!(gene %in% cls_training_genes)) %>% 
  select(-gene)

# How many left in the test set?
nrow(cls_testing)



#####
# 2
#####

# Checking effect of "na.omit()" in dependencies_freq.R
# If dependencies not set yet, uncomment below to get data
dependencies <- read.table(sprintf('%s/dependency/processed/dependencies.csv', data_dir), stringsAsFactors = F, sep=' ')

d <- dependencies

# How many dependencies?
nrow(d)

# Filter for one cell line only 
d <- dependencies %>%
  filter(dependency_p > 0.65) %>%
  filter(cl== 'HCC1428_BREAST')

# How may now?
nrow(d)

# Filter rows have na and and see which cols they are in
# NOTE: it was all missing gene_ids, so no need to remove these in dependency_freq.R
rows_with_na <- d[apply(d, 1, function(x) any(is.na(x))), ]
rows_with_na

# How many rows with na?
nrow(rows_with_na)

# Omit the na rows, or rows where anoything other than gene_id is na (uncomment as appropriate)
d <- d %>%
  #na.omit()
  filter(rowSums(is.na(select(., -gene_id))) == 0)

# How many rows without the na rows?
nrow(d)


