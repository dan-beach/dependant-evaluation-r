# Read the base_weights.csv file and store as weights_raw
# Where does this come from?
weights_raw <- read.csv(sprintf('%s/weights/base_weights_v1.csv', data_dir))

# Check col names
colnames(weights_raw)

# Need to rename some cols becasue the new data has different headings
# Using the original col heading here so there is less code to update
# The old col name are also more useful becasue they identify the tissue of origin 

# Load the model data
model_data <- read.csv(sprintf('%s/dependency/Model.csv', data_dir), stringsAsFactors = F, sep=',')

head(model_data)

# Map the CCLE name in the model data to the model id in the weights_raw col names
# R converts hyphens in col names to dots, so this swicthes them back to hypens before mapping to CCLE names 
names(weights_raw) <- model_data$CCLEName[match(gsub("\\.", "-", names(weights_raw)), model_data$ModelID)]

# Rename the first col to "Name"
names(weights_raw)[1] <- 'Name'

colnames(weights_raw)

# DRB check head 
head(weights_raw, n=10)

# DRB see what this melt function is doing before it is used below
weights_cls <- melt(weights_raw, id.vars = 'Name')
head(weights_cls, n=10)

# Gets a list of distinct cell lines from the weights_raw dataframe
# First melts from wide to long format, keeping name as ID value
# then selects just the cell_line column, gets distinct values and converts to a vector
weights_cls <- melt(weights_raw, id.vars = 'Name') %>%
  select(gene=Name, cl=variable, weight=value) %>%
  select(cell_line = cl) %>%
  distinct() %>%
  unlist() %>%
  as.vector()

# print to screen
head(weights_cls)

# selects list of cell lines from dependencies dataframe
# gets distinct values and saves as a vector
dep_cls <- dependencies %>%
  select(cl) %>%
  distinct() %>%
  unlist() %>%
  as.vector()

# print to screen
head(dep_cls)

# Create dataframes for training cell lines and testing cell lines
# training_cls = weights_cls cell-lines that are also in the dep_cls list
# testing_cls = weights_cls cell lines that are not in the dep_cls list
# set the cell line column name to cl and add a train column, set to 1 for training cell lines and 0 for testing cell lines
training_cls <- data.frame(cl=weights_cls[weights_cls %in% dep_cls], train=1)
testing_cls <- data.frame(cl=weights_cls[!weights_cls %in% dep_cls], train=0)

# DRB see what these look like
training_cls
testing_cls

# append testing dataframe to training dataframe
# add a disease column based on cl column (everything after the underscore)
# sort (arrange) by disease column in ascending order
cl_list <- rbind(training_cls, testing_cls) %>%
  mutate(disease = gsub('^[^_]*_', '\\1', cl)) %>%
  arrange(disease)

cl_list

### a few manual tweeks!

cl_list <- cl_list %>%
  filter(cl != 'NCIH460_LUNG')

# DRB see contents of cl_list
cl_list

# training_cl <-  c("X769P_KIDNEY", "X786O_KIDNEY", "A498_KIDNEY", "A704_KIDNEY", "ACHN_KIDNEY", "ASPC1_PANCREAS", "AU565_BREAST", "BFTC909_KIDNEY", "BT20_BREAST",  "BT474_BREAST","BXPC3_PANCREAS", "CAKI2_KIDNEY", "CAL120_BREAST","CAL51_BREAST", "CAL54_KIDNEY","CFPAC1_PANCREAS", "DANG_PANCREAS","EFM19_BREAST", "HCC1143_BREAST",  "HCC1187_BREAST",    "HCC1395_BREAST", "HCC1428_BREAST",    "HCC1500_BREAST", "HCC1569_BREAST","HCC1599_BREAST", "HCC1806_BREAST", "HCC1937_BREAST", "HCC1954_BREAST", "HCC202_BREAST","HCC2218_BREAST", "HCC38_BREAST", "HCC70_BREAST", "HEKTE_KIDNEY", "HPAC_PANCREAS","HPAFII_PANCREAS", "HS578T_BREAST","HS766T_PANCREAS", "HUPT3_PANCREAS",  "KMRC1_KIDNEY", "KMRC2_KIDNEY","KP2_PANCREAS", "KP4_PANCREAS", "KPL1_BREAST",  "L33_PANCREAS", "MCF7_BREAST", "MDAMB157_BREAST",   "MDAMB175VII_BREAST","MDAMB231_BREAST",   "MDAMB361_BREAST",   "MDAMB415_BREAST",  "MDAMB436_BREAST",   "MDAMB453_BREAST",   "MDAMB468_BREAST",   "MIAPACA2_PANCREAS", "NCIH292_LUNG","NCIH460_LUNG", "OSRC2_KIDNEY", "PANC0327_PANCREAS", "PANC0813_PANCREAS", "PANC1005_PANCREAS", "PSN1_PANCREAS", "QGP1_PANCREAS", "SKRC20_KIDNEY", "SLR20_KIDNEY", "SLR21_KIDNEY", "SLR23_KIDNEY", "SLR24_KIDNEY", "SLR25_KIDNEY", "SLR26_KIDNEY", "SNU349_KIDNEY", "SU8686_PANCREAS",   "SUIT2_PANCREAS",    "T47D_BREAST",  "TUHR10TKB_KIDNEY",  "TUHR14TKB_KIDNEY", "TUHR4TKB_KIDNEY", "UACC812_BREAST", "UOK101_KIDNEY","ZR751_BREAST", "ZR7530_BREAST")
# 
# test_cl <- c('BT474_BREAST', 'MCF7_BREAST', 'MFM223_BREAST', 'MDMMB231_BREAST', 'MM468_BREAST', 'T47D_BREAST',
#              'HCC1599_BREAST', 'HCC1143_BREAST', 'MDMMB486_BREAST', 'HCC38_BREAST', 'HCC70_BREAST',
#              'NCIH292', 'NCI-H345', 'NCI-H460',
#              '786-O', '769-P', 'A-498')

write.table(cl_list, sprintf('%s/supporting_files/processed/cl_list.csv', data_dir), row.names = F, col.names = T, sep=',')

rm(weights_raw, cl_list)

#### Manual Testing ####
# HCC1428_BREAST ends up as a training cell line, HCC38_BREAST does not
# I think this is because HCC1428_BREAST was in depmap dependancy data but HCC38_BREAST wasn't
# Check quickly here if this is the case (Note: my assumption is correct)
# Rerun process_dependencies.R with last line commented out to ensure dependencies and dependencies_melt dataframes have values

filtered_df <- dependencies %>% 
  #filter(cl == 'HCC1428_BREAST')
  filter(cl == 'HCC38_BREAST')

filtered_df 

filtered_df <- dependencies_melt %>% 
  #filter(cl == 'HCC1428_BREAST')
  filter(cl == 'HCC38_BREAST')

filtered_df

rm(list = c('dependencies_melt', 'raw_dependencies', 'filtered_df'))