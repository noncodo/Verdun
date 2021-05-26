################################################### 0- Check genome coverage ###################################################

# Set working directory:
current_path = ".../scripts/rfiles/" # TO ADAPT!!!
setwd(current_path) 

# Data loading:
data <- read.csv2("../../data/clinical_data.csv", header=TRUE, sep=",", na.strings=c(""))
df <- data.frame(data)

# Subgroups data loading:
subgroups <- read.csv2("../../data/phate_haplotype_clusters.csv", header=TRUE, sep=",", na.strings=c(""))
df_subgroups <- data.frame(subgroups)

# Find samples that exist only in df but not in df_subgroups:
removed_sample_id <- setdiff(df$Sample.ID, df_subgroups$Sample_ID) 

# Split df in two dataframes: 
# - df1: dataframe where samples are in df_subgroups too
# - df2: dataframe where samples are not in df_subgroups (hypothesis: genome coverage < 80%)
# Create the `%notin%` operator:
`%notin%` <- Negate(`%in%`)
df1 <- df[df$Sample.ID %notin% removed_sample_id,]
df2 <- df[df$Sample.ID %in% removed_sample_id,]
# Convert $Genome.coverage.p.c in numeric:
g1 <- as.numeric(as.character(df1$Genome.coverage.p.c.))
g2 <- as.numeric(as.character(df2$Genome.coverage.p.c.))
# Remove NA from g2:
g2 <- g2[!is.na(g2)] 
# Compute min and max:
min1 <- min(g1)
max1 <- max (g1)
min2 <- min(g2)
max2 <- max(g2)
results <- cat("min1 = ", min1, " ; max1 = ", max1, " ; min2 = ", min2, " ; max2 = ", max2)


