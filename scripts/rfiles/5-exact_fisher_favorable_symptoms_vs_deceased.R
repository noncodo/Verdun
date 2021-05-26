################################################### 5- Exact Fisher tests - Favorable symptoms vs Deceased ###################################################

# Set working directory:
current_path = ".../scripts/rfiles/" # TO ADAPT!!!
setwd(current_path) 

# Data loading:
data <- read.csv2("../../data/clinical_data.csv", header=TRUE, sep=",", na.strings=c(""))
df <- data.frame(data)

# Remove lines where Deceased == NA:
df <- df[!is.na(df$Deceased), ] 

# Save the results of the statistical tests in a file.txt:
cat("Results exact Fisher tests  - Favorable symptoms vs Deceased", file = "5-results_fisher_favorable_symptoms_vs_deceased.txt")
# add 2 newlines
cat("\n\n", file = "5-results_fisher_favorable_symptoms_vs_deceased.txt", append = TRUE)

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

# Check presence of a favorable symptom (headache, myalgia or sore throat):
favorable_symptoms <- as.character(rep(c(0), length(df$Patient.ID)))
for (i in 1:length(df$Patient.ID)){
  if (!is.na(df$Headache[i]) && !is.na(df$Myalgia[i]) && !is.na(df$sore_throat[i])){
    if (df$Headache[i] == "Yes" || df$Myalgia[i] == "Yes" || df$sore_throat[i] == "Yes"){
      favorable_symptoms[i] = "Yes"
    } else{
      favorable_symptoms[i] = "No"
    }
  }  else{
    favorable_symptoms[i] = NA
  }
}

# Add column "favorable_symptoms" to df:
df <- cbind(df, favorable_symptoms)

# Create test name:
test_name <- paste("Check independence Favorable symptoms (Headache, myalgia, sore throat)/Deceased", sep="") 
# Delete rows where favorable_symptoms == NA:
reduced_df <- df[!is.na(df$favorable_symptoms), ]
# Create the contingency table:
CT <- table(as.factor(reduced_df$favorable_symptoms), as.factor(reduced_df$Deceased))
# Exact Fisher test:
test <- fisher.test(CT) 
# Export test outputs:
cat(test_name, "\n", file = "5-results_fisher_favorable_symptoms_vs_deceased.txt", append = TRUE)
capture.output(CT, file = "5-results_fisher_favorable_symptoms_vs_deceased.txt", append = TRUE)
capture.output(test, file = "5-results_fisher_favorable_symptoms_vs_deceased.txt", append = TRUE)
# Add 2 newlines:
cat("\n --------------------------------------------------------------------------------- \n", file = "5-results_fisher_favorable_symptoms_vs_deceased.txt", append = TRUE)



# End of "5-results_fisher_favorable_symptoms_vs_deceased.txt" :
cat("END", file = "5-results_fisher_favorable_symptoms_vs_deceased.txt", append = TRUE)

