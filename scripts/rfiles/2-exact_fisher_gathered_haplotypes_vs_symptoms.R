################################################### 2- Exact Fisher tests - Gathered haplotypes vs symptoms ###################################################

# Here we do the same as in stats "1- Exact Fisher tests - haplotypes vs symptoms", except that we gather haplotypes I, VIII, IX and OTHER.
# Thus, we have haplotypes II, III and OTHER in the contengency tables.
# Note that we keep only the rows where genome coverage > 80%. The filter is made by the join with df_subgroups.

# Set working directory:
current_path = ".../scripts/rfiles/" # TO ADAPT!!!
setwd(current_path) 

# Data loading:
data <- read.csv2("../../data/clinical_data.csv", header=TRUE, sep=",", na.strings=c(""))
df <- data.frame(data)

# Remove lines where variant == NA:
df <- df[!is.na(df$Haplotype.group.simplified), ] 

# Subgroups data loading:
# Contains only sample ID with genome coverage > 80%
subgroups <- read.csv2("../../data/phate_haplotype_clusters.csv", header=TRUE, sep=",", na.strings=c(""))
df_subgroups <- data.frame(subgroups)

# Join df and df_subgroups according to Sample ID column (merge makes a inner join by default):
joined_df <- merge(df, df_subgroups, by.x="Sample.ID", by.y="Sample_ID", sort = TRUE)

# Gather haplotypes I, VIII, IX and OTHER:
gathered_haplotypes <- as.character(rep(c(0), length(joined_df$Haplotype.group.simplified)))
for (i in 1:length(joined_df$Haplotype.group.simplified)){
  if (joined_df$Haplotype.group.simplified[i] == "II" || joined_df$Haplotype.group.simplified[i] == "III"){
     gathered_haplotypes[i] <- joined_df$Haplotype.group.simplified[i]
  }
  else{
    gathered_haplotypes[i] <- "OTHER"
  }
}

# Save the results of the statistical tests in a file.txt:
cat("Results exact Fisher tests  - gathered haplotypes vs symptoms", file = "2-results_fisher_gathered_haplotypes.txt")
# add 2 newlines
cat("\n\n", file = "2-results_fisher_gathered_haplotypes.txt", append = TRUE)

# Check presence of another symptom. Make the distinction between sore throat and the other symptoms:
check_other_symptom <- as.character(rep(c(0), length(joined_df$Other.symptom)))
sore_throat <- as.character(rep(c(0), length(joined_df$Other.symptom)))
for (i in 1:length(joined_df$Other.symptom)){
  if (is.na(joined_df$Other.symptom[i])){
    check_other_symptom[i] = NA
    sore_throat[i] = NA
  } else if (joined_df$Other.symptom[i] == "No"){
    check_other_symptom[i] = "No"
    sore_throat[i] = "No"
  } else{
    check_other_symptom[i] = "Yes"
    if (grepl(joined_df$Other.symptom[i], "mal_de_gorge")){
      check_other_symptom[i] = "No"
      sore_throat[i] = "Yes"
    } else{
      check_other_symptom[i] = "Yes"
      sore_throat[i] = "No"
    }
  }
}

# Add column "sore_throat", "check_other_symptom" and "gathered_haplotypes" to joined_df:
joined_df <- cbind(joined_df, sore_throat, check_other_symptom, gathered_haplotypes)

# Variables to test:
variables_to_test <- c("Employee", "Fever", "Cough", "Fatigue", "Headache", "Myalgia", "Nausea_vomiting_diarrhea", "Dyspnea", "Confusion", "sore_throat", "check_other_symptom", "Co-morbidity", "Deceased")
# Column indexes of variables to test:
col_idx <- c(19) # "Employee" column index
col_idx <- c(col_idx, c(22:29)) # Add columns "Fever" to "Confusion"
col_idx <- c(col_idx, 41, 42, 37, 39)  # Add columns "sore_throat", "check_other_symptom", "Co-morbidity" and "Deceased" 

j = 1 # Index used to find test name

# For each variable to test:
for (i in col_idx){
	# Create test name:
	test_name <- paste("Check independence Gathered haplotypes (simplified group)/", variables_to_test[j], sep="") 
	# Delete rows where current variable to test == NA:
	reduced_joined_df <- joined_df[!is.na(joined_df[i]), ]
	# Create the contingency table:
	CT <- table(as.factor(reduced_joined_df$gathered_haplotypes), as.factor(reduced_joined_df[,i]))
	# Exact Fisher test:
	test <- fisher.test(CT) 
	# Export test outputs:
	cat(test_name, "\n", file = "2-results_fisher_gathered_haplotypes.txt", append = TRUE)
	capture.output(CT, file = "2-results_fisher_gathered_haplotypes.txt", append = TRUE)
	capture.output(test, file = "2-results_fisher_gathered_haplotypes.txt", append = TRUE)
	# Add 2 newlines:
	cat("\n --------------------------------------------------------------------------------- \n", file = "2-results_fisher_gathered_haplotypes.txt", append = TRUE)
	# Increment index to find test name
	j = j + 1
}

# End of "2-results_fisher_gathered_haplotypes.txt" :
cat("END", file = "2-results_fisher_gathered_haplotypes.txt", append = TRUE)

