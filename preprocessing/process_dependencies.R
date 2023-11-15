
# Download dependency data https://ndownloader.figshare.com/files/10416306 
# and add data to [data_dir]/dependency/gene_dependency.csv

# read raw depmap data from csv file
raw_dependencies <- read.csv(sprintf('%s/dependency/gene_dependency.csv', data_dir), stringsAsFactors = F, sep=',') ### https://portals.broadinstitute.org/achilles/datasets/18/download/gene_dependency.csv 

# print to screen num rows
nrow(raw_dependencies)
ncol(raw_dependencies)

# DRB view first 10 rows or forst 3 cols (very large dataframe)
raw_dependencies[1:10, 1:3]

# convert data frame from wide to long format
# rename columns to cell_line, gene, dependency_p
# filter duplciate rows
dependencies_melt <- melt(raw_dependencies, 'line') %>%
  select(cell_line=line, gene=variable, dependency_p=value) %>% 
  distinct()

# DRB view first 10 rows
dependencies_melt[1:10,]

# check number of unique cell_lines and number of rows
length(unique(dependencies_melt$cell_line))
nrow(dependencies_melt)

# add 'disease' column containing the cell line name minus the first part up to the underscore using gsub "global substitution" function and regex 
# then check number of rows again
dependencies_melt$disease <-  gsub('^[^_]*_', '\\1', dependencies_melt$cell_line)
nrow(dependencies_melt)

# DRB view first 10 rows
head(dependencies_melt, n=10)

# Filter dependencies_melt so it includes only cell_lines and disease columns,
# containing names of BREAST, KIDNEY and PANCREAS cell lines only
# store result as "d"
diseases <- c('BREAST', 'KIDNEY', 'PANCREAS')
d <- dependencies_melt %>%
  select(-gene, -dependency_p) %>%
  distinct() %>%
  filter(disease %in% diseases)

# DRB show first 10 rows
head(d, n=10)
View(d)

# add new column 'gene', containing existing gene column value, up to the forst period
# (gene names currently in the format ZZZ3...26009 -  this will save just ZZZ3 in the new column)
str(dependencies_melt)
dependencies_melt <- dependencies_melt %>%
  mutate(gene=gsub("\\..*","",gene))

# DRB check the first 10 rows
head(dependencies_melt, n=10)

# Select the Gene.ID and Assoiated.Gene.Name cols from id_map dataframe
# and save as gene_map dataframe 
gene_map <- id_map %>%
  select(Gene.ID, Associated.Gene.Name)

# DRB check the first 10 rows
head(gene_map, n=10)

# Left join the gene_map df to dependencies_melt on dependencies_melt$genes = gene_map$Associated.Gene.Name
# rename cell_line = cl, Gene.ID = gene_id, gene = gene_name
# select these and the depdendency_p column and remove duplicates
dependencies_melt <- dependencies_melt %>%
  left_join(gene_map, c('gene' = 'Associated.Gene.Name')) %>%
  select(cl=cell_line, gene_id=Gene.ID, gene_name=gene, dependency_p) %>%
  distinct()

# save this as the processed dependancies csv file
write.table(dependencies_melt, sprintf('%s/dependency/processed/dependencies.csv', data_dir))

# remove these objects from workspace to free up memory
rm(list = c('dependencies_melt', 'raw_dependencies', 'd', 'gene_map'))
