
weights_raw <- read.csv(sprintf('%s/weights/base_weights.csv', data_dir))
weights_cls <- melt(weights_raw, id.vars = 'Name') %>%
  select(gene=Name, cl=variable, weight=value) %>%
  select(cell_line = cl) %>%
  distinct() %>%
  unlist() %>%
  as.vector()

weights_cls

dep_cls <- dependencies %>%
  select(cl) %>%
  distinct() %>%
  unlist() %>%
  as.vector()

training_cls <- data.frame(cl=weights_cls[weights_cls %in% dep_cls], train=1)
testing_cls <- data.frame(cl=weights_cls[!weights_cls %in% dep_cls], train=0)
cl_list <- rbind(training_cls, testing_cls) %>%
  mutate(disease = gsub('^[^_]*_', '\\1', cl)) %>%
  arrange(disease)

cl_list

### a few manual tweeks!

cl_list <- cl_list %>%
  filter(cl != 'NCIH460_LUNG')


# training_cl <-  c("X769P_KIDNEY", "X786O_KIDNEY", "A498_KIDNEY", "A704_KIDNEY", "ACHN_KIDNEY", "ASPC1_PANCREAS", "AU565_BREAST", "BFTC909_KIDNEY", "BT20_BREAST",  "BT474_BREAST","BXPC3_PANCREAS", "CAKI2_KIDNEY", "CAL120_BREAST","CAL51_BREAST", "CAL54_KIDNEY","CFPAC1_PANCREAS", "DANG_PANCREAS","EFM19_BREAST", "HCC1143_BREAST",  "HCC1187_BREAST",    "HCC1395_BREAST", "HCC1428_BREAST",    "HCC1500_BREAST", "HCC1569_BREAST","HCC1599_BREAST", "HCC1806_BREAST", "HCC1937_BREAST", "HCC1954_BREAST", "HCC202_BREAST","HCC2218_BREAST", "HCC38_BREAST", "HCC70_BREAST", "HEKTE_KIDNEY", "HPAC_PANCREAS","HPAFII_PANCREAS", "HS578T_BREAST","HS766T_PANCREAS", "HUPT3_PANCREAS",  "KMRC1_KIDNEY", "KMRC2_KIDNEY","KP2_PANCREAS", "KP4_PANCREAS", "KPL1_BREAST",  "L33_PANCREAS", "MCF7_BREAST", "MDAMB157_BREAST",   "MDAMB175VII_BREAST","MDAMB231_BREAST",   "MDAMB361_BREAST",   "MDAMB415_BREAST",  "MDAMB436_BREAST",   "MDAMB453_BREAST",   "MDAMB468_BREAST",   "MIAPACA2_PANCREAS", "NCIH292_LUNG","NCIH460_LUNG", "OSRC2_KIDNEY", "PANC0327_PANCREAS", "PANC0813_PANCREAS", "PANC1005_PANCREAS", "PSN1_PANCREAS", "QGP1_PANCREAS", "SKRC20_KIDNEY", "SLR20_KIDNEY", "SLR21_KIDNEY", "SLR23_KIDNEY", "SLR24_KIDNEY", "SLR25_KIDNEY", "SLR26_KIDNEY", "SNU349_KIDNEY", "SU8686_PANCREAS",   "SUIT2_PANCREAS",    "T47D_BREAST",  "TUHR10TKB_KIDNEY",  "TUHR14TKB_KIDNEY", "TUHR4TKB_KIDNEY", "UACC812_BREAST", "UOK101_KIDNEY","ZR751_BREAST", "ZR7530_BREAST")
# 
# test_cl <- c('BT474_BREAST', 'MCF7_BREAST', 'MFM223_BREAST', 'MDMMB231_BREAST', 'MM468_BREAST', 'T47D_BREAST',
#              'HCC1599_BREAST', 'HCC1143_BREAST', 'MDMMB486_BREAST', 'HCC38_BREAST', 'HCC70_BREAST',
#              'NCIH292', 'NCI-H345', 'NCI-H460',
#              '786-O', '769-P', 'A-498')

write.table(cl_list, sprintf('%s/supporting_files/processed/cl_list.csv', data_dir), row.names = F, col.names = T, sep=',')

rm(weights_raw, cl_list)


