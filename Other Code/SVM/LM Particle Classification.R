# Libraries
library(MASS)      # For logistic regression model fitting

# Set the file path (replace "your_file.txt" with your actual file path)
electron_file <- "C:\\Users\\William Teague\\Box\\HERT_Box\\Data\\Aggregate Electron Data 1-100.txt"
electron_file <- "C:\\Users\\wzt0020\\Box\\HERT_Box\\Data\\Aggregate Electron Data 1-100.txt"
proton_file   <- "C:\\Users\\William Teague\\Box\\HERT_Box\\Data\\Aggregate Proton Data 1-50.txt"
proton_file   <- "C:\\Users\\wzt0020\\Box\\HERT_Box\\Data\\Aggregate Proton Data 1-50.txt"

# Read the data
imported_electron_data <- read.table(electron_file, sep="", skip = 1)
imported_proton_data   <- read.table(proton_file, sep="", skip = 1)

electron_data <- imported_electron_data
proton_data <- imported_proton_data

# Assign header labels
column_names <- c("Einc", "Detector1", "Detector2", "Detector3",
                  "Detector4","Detector5","Detector6",
                  "Detector7","Detector8","Detector9") 
colnames(electron_data) <- column_names
colnames(proton_data)   <- column_names

# Create a new column named "Particle Type" with "electron" for values of 1 in column 1
electron_data$Particle_Type <- 1
proton_data$Particle_Type <- 0

# Rearrange data to particle type, deposited energy
electron_data <- electron_data[ , c(ncol(electron_data),1:(ncol(electron_data)-1))]
proton_data <- proton_data[ , c(ncol(proton_data),1:(ncol(proton_data)-1))]

# Combine data into training, verification, and test groups (80/10/10 split)
combined_data    <- rbind(electron_data,proton_data)[,-2]
set.seed(123)
random_indices <- sample(nrow(combined_data), size = nrow(combined_data), replace = FALSE)
training_data <- combined_data[random_indices[1:800000], ]
verification_data <- combined_data[random_indices[800001:900000], ]
test_data         <- combined_data[random_indices[900001:1000000], ]

# Combine 10% of data
set.seed(1)
training_data_sample     <- training_data[sample(1:nrow(training_data), 
                                                 size = 100000,
                                                 replace = FALSE),] 
# Combine 1% of data
set.seed(2)
training_data_sample2     <- training_data[sample(1:nrow(training_data), 
                                                  size = 10000,
                                                  replace = FALSE),]

# Subset training data where Detector1 > 0
filtered_data <- training_data_sample2[training_data_sample2$Detector1 > 0, c(1,2,3,5,10)]


# Fit the logistic regression model
class_model <- glm(Particle_Type ~ .,data = training_data, family = binomial)
class_model_simplified <- glm(Particle_Type ~ Detector1+Detector2+Detector4+Detector9,
                   data = training_data, family = binomial)

# Extract coefficients (feature importance in logistic regression)
cm_coef <- coef(class_model)
cms_coef <- coef(class_model_simplified)
# View coefficients
summary(class_model)
summary(class_model_simplified)

# Predict class labels for new data
predictions <- as.numeric(predict(class_model, test_data)>0)
table(predict=factor(predictions, levels = c(0, 1), labels = c("Proton", "Electron")), 
      truth=factor(test_data[, 1], levels = c(0, 1), labels = c("Proton", "Electron")))

predictions <- as.numeric(predict(class_model_simplified, test_data)>0)
table(predict=factor(predictions, levels = c(0, 1), labels = c("Proton", "Electron")), 
      truth=factor(test_data[, 1], levels = c(0, 1), labels = c("Proton", "Electron")))

predictions_test <- as.numeric((cm_coef[1] 
                                + cm_coef[2] * test_data$Detector1
                                + cm_coef[3] * test_data$Detector2
                                + cm_coef[4] * test_data$Detector3
                                + cm_coef[5] * test_data$Detector4
                                + cm_coef[6] * test_data$Detector5
                                + cm_coef[7] * test_data$Detector6
                                + cm_coef[8] * test_data$Detector7
                                + cm_coef[9] * test_data$Detector8
                                + cm_coef[10] * test_data$Detector9)>0)
table(predict=factor(predictions_test, levels = c(0, 1), labels = c("Proton", "Electron")), 
      truth= factor(test_data[, 1], levels = c(0, 1), labels = c("Proton", "Electron")))

test_sample    <- test_data[sample(1:nrow(test_data), size = 10, replace = FALSE),]

predictions_test_sample <- as.numeric((cm_coef[1] 
                                + cm_coef[2] * test_sample$Detector1
                                + cm_coef[3] * test_sample$Detector2
                                + cm_coef[4] * test_sample$Detector3
                                + cm_coef[5] * test_sample$Detector4
                                + cm_coef[6] * test_sample$Detector5
                                + cm_coef[7] * test_sample$Detector6
                                + cm_coef[8] * test_sample$Detector7
                                + cm_coef[9] * test_sample$Detector8
                                + cm_coef[10] * test_sample$Detector9)>0)

table(predict=factor(predictions_test_sample, levels = c(0, 1), labels = c("Proton", "Electron")), 
      truth= factor(test_sample[, 1], levels = c(0, 1), labels = c("Proton", "Electron")))

