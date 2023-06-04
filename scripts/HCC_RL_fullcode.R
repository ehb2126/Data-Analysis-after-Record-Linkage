# example R version of original record linkage 
# 1. Load packages
#install.packages("RecordLinkage")
library(RecordLinkage) 
library(dplyr)
library(ggplot2)

# 2. Read in full HCC data and prepare files to link
# openly available data on Washington State Health Care Workers
# source: https://data.wa.gov/Health/Health-Care-Provider-Credential-Data/qxh8-f4bddf 
#df <- read.csv("Health_Care_Provider_Credential_Data.csv") # (08/23/2022 version)
#df <- df %>% filter(FirstIssueDate <= 20220228) # remove records that were probably not in the data analyzed by Dr. Flaxman
#df <- df %>% mutate(LastName = tolower(LastName), FirstName = tolower(FirstName),
#                    MiddleName = tolower(MiddleName))
#df$MiddleName <- ifelse(is.na(df$MiddleName) | df$MiddleName == "", "missing", df$MiddleName)
#df$LastInitial <- substr(df$LastName, 1, 1)

# create file A and file B (individual files to link)
#file_a <- df[df$CredentialType == 'Registered Nurse License',]
#file_b <- df[df$CredentialType == 'Registered Nurse Temporary Practice Permit',]
#row.names(file_a) <- seq(1,nrow(file_a), by = 1)
#row.names(file_b) <- seq(1,nrow(file_b))

## full analysis files
#file_a <- file_a %>% select(LastName, FirstName, MiddleName, BirthYear, FirstIssueDate, LastInitial)
#file_b <- file_b %>% select(LastName, FirstName, MiddleName, BirthYear, FirstIssueDate, LastInitial)
# nrow(file_a) - 292862
# nrow(file_b) - 84045

file_a <- read.csv("hcc_file_a.csv")
file_b <- read.csv("hcc_file_b.csv")

head(file_a)
head(file_b)

# 3. compare records (obtain comparison score vectors)
# obtain match score for comparison pairs after blocking
compare_files_full <- compare.linkage(file_a, file_b,
                                      blockfld = c("BirthYear", "LastInitial"),
                                      strcmp = c("LastName", "FirstName", "MiddleName"),
                                      exclude = c("FirstIssueDate", "BirthYear", "LastInitial"),
                                      strcmpfun = jarowinkler)

comparison_vectors_full <- compare_files_full$pairs %>% select(id1, id2, LastName, FirstName, MiddleName)

head(comparison_vectors_full, 4) 

# 4. Classify links 
# plot distribution of match scores
par(mfrow=c(1,3))
hist(comparison_vectors_full[,1], xlab = "Last Name Score", main = "", cex.lab=1.4, cex.axis=1.4)
hist(comparison_vectors_full[,2], xlab = "First Name Score", main = "", cex.lab=1.4, cex.axis=1.4)
hist(comparison_vectors_full[,3], xlab = "Middle Name Score", main = "", cex.lab=1.4, cex.axis=1.4)
par(mfrow=c(1,1))

# select a threshold from the histogram (as in full data analysis)
threshold = 0.85

# classify link candidates
links_full <- compare_files_full$pairs %>% filter(LastName >= threshold & FirstName >= threshold
                                                  & MiddleName >= threshold)
nrow(links_full) 

# 5. Obtain linked data set
### obtain individual data frames with matched records
ida_full <- file_a[links_full$id1,] %>% mutate("id.a" = links_full$id1) %>% 
  select(id.a, LastName, FirstName, MiddleName, BirthYear, FirstIssueDate)
idb_full <- file_b[links_demo$id2,] %>% mutate("id.b" = links_demo$id2) %>% 
  select(id.b, LastName_TP = LastName, FirstName_TP = FirstName, 
         MiddleName_TP = MiddleName, BirthYear_TP = BirthYear, 
         FirstIssueDate_TP = FirstIssueDate)
### find linked data set
lds_full <- cbind(ida_full, idb_full)

head(lds_full, 4) 