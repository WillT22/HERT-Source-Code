# Your distance and radiation rate data
distance <- c(0.5, 1, 2, 3)
radiation <- c(0.57031, 0.14243, 0.037307, 0.01711)

# Create a data frame with your data
data <- data.frame(distance, radiation)

# Fit a power regression model
model <- nls(radiation ~ a * distance^b, data = data, start = list(a = 1, b = -2))
summary(model)
RSS.p <- sum(residuals(model)^2)  # Residual sum of squares
TSS <- sum((radiation - mean(radiation))^2)  # Total sum of squares
1 - (RSS.p/TSS)  # R-squared measure

# Extract the coefficients
a <- coef(model)["a"]
b <- coef(model)["b"]

# Print the coefficients
cat("Coefficient 'a':", a, "\n")
cat("Coefficient 'b':", b, "\n")

# Calculate the covariance matrix
cov_matrix <- vcov(model)

# Get the confidence intervals for coefficients
conf_intervals <- confint(model)

# Extract the confidence intervals for 'a' and 'b'
conf_interval_a <- conf_intervals["a", ]
conf_interval_b <- conf_intervals["b", ]

# Print the confidence intervals
cat("Confidence interval for 'a':", conf_interval_a, "\n")
cat("Confidence interval for 'b':", conf_interval_b, "\n")
