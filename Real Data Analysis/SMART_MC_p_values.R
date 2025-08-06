rm(list=ls())

################################################################################
### p-values
################################################################################

setwd("U:/SMART-MC/Real Data Analysis")

coeff_raw <- read.csv("Output_CoeffSummary.csv", header = FALSE)
se_raw <- read.csv("Output_SeSummary.csv", header = FALSE)

# Convert to numeric matrix in case they're read as characters
coeff <- as.matrix(sapply(coeff_raw[,5:9], as.numeric))
se <- as.matrix(sapply(se_raw[,5:9], as.numeric))

# Compute z-scores
z_scores <- coeff / se

# Compute two-sided p-values
p_values <- round(2 * (1 - pnorm(abs(z_scores))), 3)

p_values_v2 <- cbind(coeff_raw[, 1:2], p_values) 

write.csv(p_values_v2, "Output_PValue_Summary.csv", row.names = FALSE)