# Download dependency data from DepMap: https://depmap.org/portal/download/all/
# and add data to [data_dir]/dependency/CRISPRGeneDependency.csv

library(stringr)

# read raw depmap data from csv file
raw_dependencies <- read.csv(sprintf('%s/dependency/CRISPRGeneDependency.csv', data_dir), stringsAsFactors = F, sep=',') ### DepMap 23Q4 release
model_data <- read.csv(sprintf('%s/dependency/Model.csv', data_dir), stringsAsFactors = F, sep=',')

# print to screen num rows
nrow(raw_dependencies)
ncol(raw_dependencies)
View(head(model_data))

# DRB view first 10 rows of first 3 cols (very large dataframe)
raw_dependencies[1:10, 1:3]

# For ease, get data into same format as the original dataset, with cellLineName_tissue as line ID instead of BROAD ID
# Join the dataframes by ModelID
dep_model_joined <- left_join(raw_dependencies, model_data, by = "ModelID")

# Replace ModelID in raw_dependencies with the CCLEName 
raw_dependencies$ModelID <- dep_model_joined$CCLEName

# Optionally, remove rows where CombinedName is NA (if there are any)
#df1 <- df1[!is.na(df1$ModelID), ]

# View the modified dataframe
View(raw_dependencies[1:10, 1:3])

# convert data frame from wide to long format
# rename columns to cell_line, gene, dependency_p
# filter duplciate rows
dependencies_melt <- melt(raw_dependencies, 'ModelID') %>%
  select(cell_line=ModelID, gene=variable, dependency_p=value) %>% 
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

# add new column 'gene', containing existing gene column value, up to the first period
# (gene names currently in the format ZZZ3...26009 -  this will save just ZZZ3 in the new column)
str(dependencies_melt)
dependencies_melt <- dependencies_melt %>%
  mutate(gene=gsub("\\..*","",gene))

# DRB check the first 10 rows
head(dependencies_melt, n=10)

#####
# Left join to gene_map section below commented out
# We are not using ENSG gene_id for Reactome and this join was adding nearly 2 million more rows of data!
#####

# # Check num rows before left join
# nrow(dependencies_melt)

# # This step not needed for reactome!
# # Select the Gene.ID and Associated.Gene.Name cols from id_map dataframe
# # and save as gene_map dataframe 
# gene_map <- id_map %>%
#   select(Gene.ID, Associated.Gene.Name)

# # DRB check the first 10 rows
# head(gene_map, n=10)

# # Left join the gene_map df to dependencies_melt on dependencies_melt$genes = gene_map$Associated.Gene.Name
# # rename cell_line = cl, Gene.ID = gene_id, gene = gene_name
# # select these and the depdendency_p column and remove duplicates
# dependencies_melt <- dependencies_melt %>%
#   left_join(gene_map, c('gene' = 'Associated.Gene.Name')) %>%
#   select(cl=cell_line, gene_id=Gene.ID, gene_name=gene, dependency_p) %>%
#   distinct()
# 
# head(dependencies_melt)
# 
# # Check num rows after left join
# nrow(dependencies_melt)


# If left join above commented out, we need to select and rename the columns
 dependencies_melt <- dependencies_melt %>%
   select(cl=cell_line, gene_name=gene, dependency_p) %>%
   distinct()

# save this as the processed dependencies csv file
write.table(dependencies_melt, sprintf('%s/dependency/processed/dependencies.csv', data_dir))

# remove these objects from workspace to free up memory
#rm(list = c('d', 'gene_map', 'dep_model_joined'))
rm(list = c('dependencies_melt', 'raw_dependencies', 'd', 'dep_model_joined'))
