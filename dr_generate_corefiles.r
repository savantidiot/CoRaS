#### This script takes PanCan33 TSV file (2GB), patient clinical TSV file (20MB), Nienke's 165 drugs data and generates the following:

### RData object with four pre-processed files
# 1 - BRCA expression data with mRNA expression of 1097 cancer patient samples
# 2 - patient clinical variables for survival analysis - based on "days_to_last_followup", "days_to_death"
# 3 - drug gene signatures for 165 drugs from Nienke's small molecule list
# 4 - subtype data for 1097 BRCA patients from a separate file because subtype data is not availabe in the original patient clinical variables - from GDC legacy archive: https://portal.gdc.cancer.gov/legacy-archive/files/735bc5ff-86d1-421a-8693-6e6f92055563

# 1) brca_cancer - BRCA expression data with mRNA expression of 1097 cancer patient samples + 1 column for gene_IDs
# 2) brca_days - clinical variables (vital status, time in days) of 1097 patients for survival analysis
# 3) gene_array - contains list of genes stored as list object for 165 drugs
# 4) gene_index - contains the number associated with drug in HMS LINCS to link back during analysis
# 5) subtype_cancer - contains molecular subtype (ER+, HER2+, PR+, Triplenegative) for 1097 brca patients

## Source files
# PanCan33: gene expression matrix normalized for 33 tissue types (TCGA data): https://www.synapse.org/#!Synapse:syn497636
# Clinical data of patients whose samples are part of PanCan33: https://www.synapse.org/#!Synapse:syn4983466.1
# TCGA barcodes and study abbreviations (GDC)
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
# Clinical data of BRCA patient samples used in Cox analysis (GDC): https://portal.gdc.cancer.gov/legacy-archive/files/735bc5ff-86d1-421a-8693-6e6f92055563
# HMS LINCS small molecule drugs resource: http://www.smallmoleculesuite.org/apps/hms_small_mol/
# Subtype clinical characteristics from new clinical data from GDC legacy archive: https://portal.gdc.cancer.gov/legacy-archive/files/735bc5ff-86d1-421a-8693-6e6f92055563



# install.packages("caret", dependencies = c("Depends", "Suggests"))
# source("https://bioconductor.org/biocLite.R")
# biocLite("survcomp") # for concordance.index
require(data.table) # for fread
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)

## set source directory
sourcedir <- './sourcedata/'


## 1- Store BRCA expression data with mRNA expression of 1097 cancer patient samples

# Get TCGA and clinical data from raw files
clin_file = paste0(sourcedir,'clinical_PANCAN_patient_with_followup.tsv')
tcga_file = paste0(sourcedir,'EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv')
clin_data <- fread(clin_file, header=T, sep='\t', check.names = FALSE, na.strings = c("NA", "[Not Available]","[Not Applicable]"))
tcga_data <- fread(tcga_file, header=T, sep='\t', check.names = FALSE)

# Get subset of brca data using values from clinical data barcodes
brca_clin_data <- clin_data %>% filter(acronym=="BRCA")  # select only BRCA patients from clinical data set
brca_cols_tcga <- which( substr(colnames(tcga_data), 1, 12) %in%  brca_clin_data$bcr_patient_barcode )  # get index of columns in TCGA data that are BRCA related
brca_data <- tcga_data %>% select(brca_cols_tcga)  # select columns based on index for BRCA samples
rownames(brca_data) <- tcga_data$gene_id  # assign gene names to brca data

# Identify tumor samples from normal samples for classification and store as factor
norm_codes <- c("11","12","13","14","15","16","17","18","19")  # In the TCGA barcode, characters 14-15 indicate whether tissue is normal or cancerous.  Characters "11"-"19" indicate normal samples
brca_normal_cols <- which( substr(colnames(brca_data), 14, 15) %in% norm_codes )  # capture index of normal samples
brca_data_cancer <- brca_data %>% select(-brca_normal_cols)

# Tissue samples contain a few duplicates from the same patient -- we will remove these duplicates to retain one sample per patient
brca_tiss_dup <- names(which(table( substr(colnames(brca_data_cancer), 1, 12) ) > 1))  # get counts of substring; search for counts >1; get the column names of those indexes
brca_tiss_dup <- grep(paste0("\\b",brca_tiss_dup,"\\b",collapse="|"), colnames(brca_data_cancer), ignore.case = T)  # gets all duplicates
brca_tiss_dup <- brca_tiss_dup[c(TRUE,FALSE)] # removes alternate indexes as the duplicates are adjacent to each other
brca_data_cancer <- brca_data_cancer %>% select(-brca_tiss_dup)  # remove duplicates from tissue sample data table

## Finally save the de-duplicated and subset BRCA cancer expression data table along with gene IDs to RData object
brca_cancer <- brca_data_cancer %>% mutate(gene_name = tcga_data$gene_id)


## 2 - Store patient clinical variables for survival analysis - based on "days_to_last_followup", "days_to_death"

# Filter patients vital status ("Alive" vs. "Dead") for the de-duplicated sample data
brca_vital <- brca_clin_data %>% filter( bcr_patient_barcode %in% substr(colnames(brca_data_cancer), 1, 12) ) %>% select(c(bcr_patient_barcode,vital_status))
brca_vital <- factor(brca_vital$vital_status)  # make vital status a factor type vector

# Filter patients total survival in days (for patients in class 'Alive', take value from "days_to_last_followup"; for patients in class 'Alive', take value from "days_to_death")
brca_days <- brca_clin_data %>% filter( bcr_patient_barcode %in% substr(colnames(brca_data_cancer), 1, 12) ) %>% select(c(bcr_patient_barcode,as.numeric(days_to_last_followup), as.numeric(days_to_death)))
brca_days <- brca_days %>% mutate(days = 0)
# Assign "days_to_death" value to column "days" if status is 'Dead', else assign value in "days_to_last_followup"
brca_days[which(brca_vital=="Dead"), "days"] <- brca_days[which(brca_vital=="Dead"), "days_to_death"]
brca_days[which(brca_vital=="Alive"), "days"] <- brca_days[which(brca_vital=="Alive"), "days_to_last_followup"]


# Final pre-processing and checks on data - remove negative days by making them equal to zero (for survival analysis) and change NA values to zero (similar to lost to follow up)
days_na <- which(is.na(brca_days$days))
ifelse(length(days_na) > 0, brca_days[days_na, "days"] <- brca_days[days_na,"days_to_last_followup"], brca_days <- brca_days)

days_lessthanzero <- which(brca_days$days < 0)
ifelse(length(days_lessthanzero) > 0, brca_days[days_lessthanzero, "days"] <- 0, brca_days <- brca_days)

## Add vital status to the data table
brca_days$days <- as.numeric(brca_days$days)
brca_days <- brca_days %>% mutate(status = brca_vital)

# cox model requires censoring status to be 0 or 1; here we censor if Alive (0) as we're modeling surivival outcomes
brca_days$status <- ifelse(brca_days$status =="Alive",0,1)

## 3 - Store drug-gene signatures to a data table
## Objects stored include:
# gene_array contains list of genes stored as list object for 165 drugs
# gene_index contains the number associated with drug in HMS LINCS to link back during analysis

# Parses a .gmt file and puts it into the list format
n_file <- paste0(sourcedir,'Nienke_10genes.gmt')

# Get gene list for drug from Nienke's file
n_data <- fread(n_file, header=F, sep='\t', check.names = FALSE, nrows=165, fill = TRUE)
gene_array <- list() # gene signatures associated with each of the 165 compounds
gene_index <- vector() # compound ID
for (r in 1:nrow(n_data)){
  gene_temp <- as.matrix(n_data[r,3:ncol(n_data),with=F])
  gene_temp <- gene_temp[ which(gene_temp != "") ]
  g_index <- strsplit(as.character(n_data[r,1]),"=")[[1]][2] # split the first term, i.e., link ("http://lincs.hms.harvard.edu/db/?search=10001") to obtain ID
  gene_array[r] <- list(gene_temp)
  gene_index <- c(gene_index, g_index)
}


## 3 - Store subtype data to a data table
## Get new clinical data from GDC archive as molecular subtype info not available in previous clinical file
subtype_file <- paste0(sourcedir,'nationwidechildrens.org_clinical_patient_brca.txt')
subtype_data <- fread(subtype_file, header=T, sep='\t', check.names = FALSE)
subtype_data <- subtype_data[c(-1,-2),] ## first two rows are header lablel repeats
subtype_data <- subtype_data %>% mutate(mol_type = case_when( (er_status_by_ihc=='Positive' & pr_status_by_ihc=='Positive' & her2_status_by_ihc=='Positive') ~ "ER+PR+HER2+",
                                                                        (er_status_by_ihc=='Positive' & pr_status_by_ihc=='Positive') ~ "ER+PR+",
                                                                        (er_status_by_ihc=='Positive' & her2_status_by_ihc=='Positive') ~ "ER+Her2+",
                                                                        (pr_status_by_ihc=='Positive' & her2_status_by_ihc=='Positive') ~ "PR+Her2+",
                                                                        (er_status_by_ihc=='Positive') ~ "ER+",
                                                                        (pr_status_by_ihc=='Positive') ~ "PR+",
                                                                        (her2_status_by_ihc=='Positive') ~ "Her2+",
                                                                        TRUE ~ "TripleNegative") )

# Finally, select only required columns for use in analysis
subtype_cancer <- subtype_data %>% select(bcr_patient_uuid, bcr_patient_barcode, er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc,mol_type)


## Finally, store the processed data objects to an RData object: 
# 1) brca_cancer - BRCA expression data with mRNA expression of 1097 cancer patient samples + 1 column for gene_IDs
# 2) brca_days - clinical variables (vital status, time in days) of 1097 patients for survival analysis
# 3) gene_array - contains list of genes stored as list object for 165 drugs
# 4) gene_index - contains the number associated with drug in HMS LINCS to link back during analysis
# 5) subtype_cancer - contains molecular subtype (ER+, HER2+, PR+, Triplenegative) for 1097 brca patients

## set destination directory (store all generated data here)
destdir <- './generateddata/'

save(brca_cancer, brca_days, gene_array, gene_index, subtype_cancer, file = paste0(destdir,"drugrepoData.RData"))

