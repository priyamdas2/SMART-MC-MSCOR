# 1 = B-cell depletion (rituximab + ocrelizumab)  (original code = 1)
# 2 = glatiramer acetate                          (original code = 2)
# 3 = Interferon-beta                             (original code = 3)
# 4 = dimethyl fumarate                           (original code = 4)
# 5 = natalizumab                                 (original code = 5)
# 6 = S1P modulators (fingolimod + teriflunomide) (original code = 6,7)
# 7 = Aggressive / Legacy therapies (cyclophosphamide, mitoxantrone, alemtuzumab) (original code = 8, 9, 10)


setwd("U:/SMART-MC/Real Data Analysis/Real Data")
data <- read.csv("MS treatment sequences.csv", header = FALSE)

# Apply collapsing logic on second column
data$V2 <- ifelse(data$V2 %in% c(6,7), 6,
                  ifelse(data$V2 %in% c(8,9,10), 7, data$V2))

# View first few rows to verify
head(data)

# Write the modified data to a new file (optional)
write.table(data, "MS_treatment_sequences_collapsed.csv", 
            sep = ",", row.names = FALSE, col.names = FALSE)