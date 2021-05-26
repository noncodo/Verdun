################################################### 4- Exact Fisher tests - Symptoms vs deceased ###################################################

# Set working directory:
current_path = ".../scripts/rfiles/" # TO ADAPT!!!
setwd(current_path) 

# Data loading:
data <- read.csv2("../../data/clinical_data.csv", header=TRUE, sep=",", na.strings=c(""))
df <- data.frame(data)

# Remove lines where Deceased == NA:
df <- df[!is.na(df$Deceased), ] 

# Save the results of the statistical tests in a file.txt:
cat("Results exact Fisher tests  - Symptoms vs Deceased", file = "4-results_fisher_symptoms_vs_deceased.txt")
# add 2 newlines
cat("\n\n", file = "4-results_fisher_symptoms_vs_deceased.txt", append = TRUE)

# Check presence of another symptom. Make the distinction between sore throat and the other symptoms:
check_other_symptom <- as.character(rep(c(0), length(df$Other.symptom)))
sore_throat <- as.character(rep(c(0), length(df$Other.symptom)))
for (i in 1:length(df$Other.symptom)){
  if (is.na(df$Other.symptom[i])){
    check_other_symptom[i] = NA
    sore_throat[i] = NA
  } else if (df$Other.symptom[i] == "No"){
    check_other_symptom[i] = "No"
    sore_throat[i] = "No"
  } else{
    check_other_symptom[i] = "Yes"
    if (grepl(df$Other.symptom[i], "mal_de_gorge")){
      check_other_symptom[i] = "No"
      sore_throat[i] = "Yes"
    } else{
      check_other_symptom[i] = "Yes"
      sore_throat[i] = "No"
    }
  }
}

# Add column "sore_throat" and "check_other_symptom" to df:
df <- cbind(df, sore_throat, check_other_symptom)

# Variables to test:
variables_to_test <- c("Fever", "Cough", "Fatigue", "Headache", "Myalgia", "Nausea_vomiting_diarrhea", "Dyspnea", "Confusion", "sore_throat", "check_other_symptom", "Co-morbidity")
# Column indexes of variables to test:
col_idx <- c(22:29) # Add columns "Fever" to "Confusion"
col_idx <- c(col_idx, 40, 41, 37)  # Add columns "sore_throat", "check_other_symptom" and "Co-morbidity" 

j = 1 # Index used to find test name

# For each variable to test:
for (i in col_idx){
	# Create test name:
	test_name <- paste("Check independence ", variables_to_test[j],  "/Deceased", sep="") 
	# Delete rows where current variable to test == NA:
	reduced_df <- df[!is.na(df[i]), ]
	# Create the contingency table:
	CT <- table(as.factor(reduced_df[,i]), as.factor(reduced_df$Deceased))
	# Exact Fisher test:
	test <- fisher.test(CT) 
	# Export test outputs:
	cat(test_name, "\n", file = "4-results_fisher_symptoms_vs_deceased.txt", append = TRUE)
	capture.output(CT, file = "4-results_fisher_symptoms_vs_deceased.txt", append = TRUE)
	capture.output(test, file = "4-results_fisher_symptoms_vs_deceased.txt", append = TRUE)
	# Add 2 newlines:
	cat("\n --------------------------------------------------------------------------------- \n", file = "4-results_fisher_symptoms_vs_deceased.txt", append = TRUE)
	# Increment index to find test name
	j = j + 1
}

# End of "4-results_fisher_symptoms_vs_deceased.txt" :
cat("END", file = "4-results_fisher_symptoms_vs_deceased.txt", append = TRUE)

