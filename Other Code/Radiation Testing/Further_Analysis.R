# Read the text file into a data frame
time_means_data <- read.delim("C:/Users/wzt0020/Box/HERT_Box/Radiation Testing/RadTestData/time_means_data.txt", header=FALSE)

# Assign meaningful column names
colnames(time_means_data) <- c("Time", "Channel1", "Channel2", "Channel3", "Channel4", "Channel5", "Channel6", "Channel7", "Channel8")

# Perform linear regression for each channel
linear_models <- lapply(2:ncol(time_means_data), 
                        function(i) lm(time_means_data[, i] ~ time_means_data$Time))

# Create matrices for slope and intercept-related values
slope_matrix <- matrix(nrow = ncol(time_means_data) - 1, ncol = 4)
intercept_matrix <- matrix(nrow = ncol(time_means_data) - 1, ncol = 4)

# Extract R-squared values from each model
r_squared_values <- sapply(linear_models, function(model) summary(model)$r.squared)

# Populate matrices
for (i in 2:ncol(time_means_data)) {
  # Get the linear model
  model <- linear_models[[i - 1]]
  
  # Extract slope-related values
  slope_matrix[i - 1, 1] <- coef(model)[2]  # Slope coefficient
  slope_matrix[i - 1, 2] <- sd(resid(model))  # Standard deviation
  slope_matrix[i - 1, 3:4] <- confint(model)[2, ]  # Confidence interval bounds
  
  # Extract intercept-related values
  intercept_matrix[i - 1, 1] <- coef(model)[1]  # Intercept coefficient
  intercept_matrix[i - 1, 2] <- sd(resid(model))  # Standard deviation
  intercept_matrix[i - 1, 3:4] <- confint(model)[1, ]  # Confidence interval bounds
}

# Assign column names
colnames(slope_matrix) <- c("Slope", "Std_Dev", "Lower_Bound", "Upper_Bound")
colnames(intercept_matrix) <- c("Intercept", "Std_Dev", "Lower_Bound", "Upper_Bound")

# Save matrices to files
write.table(slope_matrix, "slope_matrix.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(intercept_matrix, "intercept_matrix.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


