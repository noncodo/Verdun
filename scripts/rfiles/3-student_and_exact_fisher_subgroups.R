################################################### 3- Student and exact Fisher tests - Subgroups ###################################################

# Here we do not consider the variable "haplotype (simplified group)" but the variable "subgroup".
# Thus, we have subgroups (haplotypes) IIa, IIb, III and OTHER in the contengency tables.
# Note that we keep only the rows where genome coverage > 80%. The filter is made by the join with df_subgroups.

# Index:
#### 1) Compare mean age of patients from subgroup IIa with mean age of other patients
#### 2) Compare mean age of patients from subgroup IIb with mean age of other patients
#### 3) Check independence between employee and all subgroups
#### 4) Check independence between employee and subgroup IIa
#### 5) Check independence between employee and subgroup IIb
#### 6) Check independence between subgroups and symptoms
#### 7) Check independence between subgroup IIa and symptoms
#### 8) Check independence between subgroup IIb and symptoms

# Set working directory:
current_path = ".../scripts/rfiles/" # TO ADAPT!!!
setwd(current_path) 

# Data loading:
data <- read.csv2("../../data/clinical_data.csv", header=TRUE, sep=",", na.strings=c(""))
df <- data.frame(data)

# Subgroups data loading:
subgroups <- read.csv2("../../data/phate_haplotype_clusters.csv", header=TRUE, sep=",", na.strings=c(""))
df_subgroups <- data.frame(subgroups)

# Join df and df_subgroups according to Sample ID column (merge makes a inner join by default):
joined_df <- merge(df, df_subgroups, by.x="Sample.ID", by.y="Sample_ID", sort = TRUE)

# Save the results of the statistical tests in a file.txt:
cat("Results exact Fisher tests - subgroups", file = "3-results_student_and_exact_fisher_subgroups.txt")
# add 2 newlines
cat("\n\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)

# Check presence of another symptom in joined_df. Make the distinction between sore throat and the other symptoms:
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

# Create a new vector to check presence of subgroup IIb in joined_df:
check_2b <- as.character(rep(c(0), length(joined_df$Subgroups)))
for (i in 1:length(joined_df$Subgroups)){
  if (joined_df$Subgroups[i] == "IIb"){
    check_2b[i] = "Yes"
  } else{
    check_2b[i] = "No"
  }
}

# Create a new vector to check presence of subgroup IIa in joined_df:
check_2a <- as.character(rep(c(0), length(joined_df$Subgroups)))
for (i in 1:length(joined_df$Subgroups)){
  if (joined_df$Subgroups[i] == "IIa"){
    check_2a[i] = "Yes"
  } else{
    check_2a[i] = "No"
  }
}

# Add columns sore_throat, check_other_symptom, check_2a and check_2b to joined_df :
joined_df <- cbind(joined_df, check_2a, check_2b, sore_throat, check_other_symptom)



#### 1) Compare mean age of patients from subgroup IIa with mean age of other patients:
ages_2a <- c()
ages_non_2a <- c()
for (i in 1:dim(joined_df)[1]){
  if (joined_df$Subgroups[i] == "IIa"){
    ages_2a <- c(ages_2a, joined_df$Age[i])
  } else{
    ages_non_2a <- c(ages_non_2a, joined_df$Age[i])
  }
}
# Check normality and homoscedasticity (required for T-test):
hist(ages_2a)
hist(ages_non_2a)
var.test(ages_2a, ages_non_2a, alternative = "two.sided") 
# Export test outputs:
cat("T-test age patients with IIa vs age other patients\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
cat("T test age patients IIa vs age other patients IMPOSSIBLE because normality and homoscedasticity not verified!\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
# Add 2 newlines:
cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)



#### 2) Compare mean age of patients from subgroup IIb with mean age of other patients:
ages_2b <- c()
ages_non_2b <- c()
for (i in 1:dim(joined_df)[1]){
  if (joined_df$Subgroups[i] == "IIb"){
    ages_2b <- c(ages_2b, joined_df$Age[i])
  } else{
    ages_non_2b <- c(ages_non_2b, joined_df$Age[i])
  }
}
# Check normality and homoscedasticity (required for T-test):
hist(ages_2b)
hist(ages_non_2b)
var.test(ages_2b, ages_non_2b, alternative = "two.sided") #p-value = 0.6946 > 5% --> we can not reject H0 (true ratio of variances = 1)
# T-test:
ttest_2b_age <- t.test(ages_2b, ages_non_2b, paired = FALSE)
# Export test outputs:
cat("T-test age patients with IIb vs age other patients\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
capture.output(ttest_2b_age, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
# Add 2 newlines:
cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)



#### 3) Check independence between employee and all subgroups:
# Delete rows where $Employee == NA:
reduced_joined_df <- joined_df[!is.na(joined_df$Employee), ]
# Create the contingency table:
CT <- table(as.factor(reduced_joined_df$Employee), as.factor(reduced_joined_df$Subgroups))
# Exact Fisher test:
test <- fisher.test(CT) 
# Export test outputs:
cat("Check independence Employee/Subgroups\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
capture.output(CT, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
capture.output(test, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
# Add 2 newlines:
cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)



#### 4) Check independence between employee and subgroup IIa:
# Create the contingency table:
CT <- table(as.factor(reduced_joined_df$Employee), as.factor(reduced_joined_df$check_2a))
# Exact Fisher test:
test <- fisher.test(CT) 
# Export test outputs:
cat("Check independence Employee/Subgroup IIa\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
capture.output(CT, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
capture.output(test, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
# Add 2 newlines:
cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)



#### 5) Check independence between employee and subgroup IIb:
# Create the contingency table:
CT <- table(as.factor(reduced_joined_df$Employee), as.factor(reduced_joined_df$check_2b))
# Exact Fisher test:
test <- fisher.test(CT) 
# Export test outputs:
cat("Check independence Employee/Subgroup IIb\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
capture.output(CT, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
capture.output(test, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
# Add 2 newlines:
cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)



#### 6) Check independence between subgroups and symptoms:
# Variables to test:
variables_to_test <- c("Fever", "Cough", "Fatigue", "Headache", "Myalgia", "Nausea_vomiting_diarrhea", "Dyspnea", "Confusion", "sore_throat", "check_other_symptom", "Co-morbidity", "Deceased")
# Column indexes of variables to test:
col_idx <- c(22:29) # Columns "Fever" to "Confusion"
col_idx <- c(col_idx, 43, 44, 37, 39)  # Add columns "sore_throat", "check_other_symptom", "Co-morbidity" and "Deceased" 
j = 1 # Index used to find test name
# For each variable to test:
for (i in col_idx){
  # Create test name:
  test_name <- paste("Check independence Subgroups/", variables_to_test[j], sep="") 
  # Delete rows where current variable to test == NA:
  reduced_joined_df <- joined_df[!is.na(joined_df[i]), ]
  # Create the contingency table:
  CT <- table(as.factor(reduced_joined_df$Subgroups), as.factor(reduced_joined_df[,i]))
  # Exact Fisher test:
  test <- fisher.test(CT) 
  # Export test outputs:
  cat(test_name, "\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  capture.output(CT, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  capture.output(test, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  # Add 2 newlines:
  cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  # Increment index to find test name
  j = j + 1
} 



#### 7) Check independence between subgroup IIa and symptoms:
j = 1 # Index used to find test name
# For each variable to test:
for (i in col_idx){
  # Create test name:
  test_name <- paste("Check independence Subgroup IIa/", variables_to_test[j], sep="") 
  # Delete rows where current variable to test == NA:
  reduced_joined_df <- joined_df[!is.na(joined_df[i]), ]
  # Create the contingency table:
  CT <- table(as.factor(reduced_joined_df$check_2a), as.factor(reduced_joined_df[,i]))
  # Exact Fisher test:
  test <- fisher.test(CT) 
  # Export test outputs:
  cat(test_name, "\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  capture.output(CT, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  capture.output(test, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  # Add 2 newlines:
  cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  #Increment index to find test name
  j = j + 1
} 



#### 8) Check independence between subgroup IIb and symptoms:
j = 1 # Index used to find test name
# For each variable to test:
for (i in col_idx){
  # Create test name:
  test_name <- paste("Check independence Subgroup IIb/", variables_to_test[j], sep="") 
  # Delete rows where current variable to test == NA:
  reduced_joined_df <- joined_df[!is.na(joined_df[i]), ]
  # Create the contingency table:
  CT <- table(as.factor(reduced_joined_df$check_2b), as.factor(reduced_joined_df[,i]))
  # Exact Fisher test:
  test <- fisher.test(CT) 
  # Export test outputs:
  cat(test_name, "\n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  capture.output(CT, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  capture.output(test, file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  # Add 2 newlines:
  cat("\n --------------------------------------------------------------------------------- \n", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
  # Increment index to find test name
  j = j + 1
}



# End of "3-results_student_and_exact_fisher_subgroups.txt" :
cat("END", file = "3-results_student_and_exact_fisher_subgroups.txt", append = TRUE)
