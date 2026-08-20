# Prepare packages for SVM
library(e1071)
library(caret)
library(ggplot2)
library(GGally)
library(ggthemes)
library(colorBlindness)
library(pracma)
library(cowplot)
library(magick)

setwd("C:\\Users\\wzt0020\\Box\\HERT_Box\\Particle Classification")

# Set the file path (replace "your_file.txt" with your actual file path)
electron_file <- "D:\\HERT_Drive\\Matlab Main\\Result\\Electron_FS\\Aggregate Data\\Aggregate_Electron_FS_Data_new.txt"
proton_file   <- "D:\\HERT_Drive\\Matlab Main\\Result\\Proton_FS\\Aggregate Data\\Aggregate_Proton_FS_Data_new.txt"

# Read the data
imported_electron_data <- read.table(electron_file, sep="", skip = 1)
imported_proton_data   <- read.table(proton_file, sep="", skip = 1)

# Assign header labels
column_names <- c("E_Inc",
                  "Detector1", "Detector2", "Detector3",
                  "Detector4","Detector5","Detector6",
                  "Detector7","Detector8","Detector9")
colnames(imported_electron_data) <- column_names
colnames(imported_proton_data)   <- column_names

# Create a new column named "Particle Type" with "electron" for values of 1 in column 1
imported_electron_data$Particle_Type <- as.factor("Electron")
imported_proton_data$Particle_Type <- as.factor("Proton")

# Apply the condition to the full data sets first.
electron_data_filtered <- imported_electron_data[imported_electron_data$Detector1 >= 0.1, ]
proton_data_filtered   <- imported_proton_data[imported_proton_data$Detector1 >= 0.1, ]
rm(imported_electron_data)
rm(imported_proton_data)

# Sum detectors 7 and 8 after the data is filtered.
electron_data_filtered$Detector7_8_sum <- electron_data_filtered$Detector7 + electron_data_filtered$Detector8
proton_data_filtered$Detector7_8_sum <- proton_data_filtered$Detector7 + proton_data_filtered$Detector8

cols_to_modify <- names(electron_data_filtered)[!names(electron_data_filtered) %in% c("E_Inc", "Detector7", "Detector8", "Particle_Type")]
electron_data_filtered[, cols_to_modify][electron_data_filtered[, cols_to_modify] < 0.1] <- 0
proton_data_filtered[, cols_to_modify][proton_data_filtered[, cols_to_modify] < 0.1] <- 0

# Set the number of rows you want for each data set.
num_training_rows <- 400000
num_validation_rows <- 50000
num_test_rows <- 50000
num_points<- (num_training_rows+num_validation_rows+num_test_rows)*2

# Randomly sample row indices from the filtered data strictly without overlap
num_electron_rows <- nrow(electron_data_filtered)
num_proton_rows   <- nrow(proton_data_filtered)

total_requested <- num_training_rows + num_validation_rows + num_test_rows
max_available <- min(num_electron_rows, num_proton_rows)

if (total_requested > max_available) {
  cat(sprintf("Warning: Only %d rows available per particle. Scaling down splits.\n", max_available))
  actual_train <- floor(max_available * (num_training_rows / total_requested))
  actual_val   <- floor(max_available * (num_validation_rows / total_requested))
  actual_test  <- max_available - actual_train - actual_val
} else {
  actual_train <- num_training_rows
  actual_val   <- num_validation_rows
  actual_test  <- num_test_rows
}

# Generate Mutually Exclusive Indices
set.seed(42)
all_e_indices <- sample(1:num_electron_rows, max_available, replace = FALSE)
e_training_indices   <- all_e_indices[1:actual_train]
e_validation_indices <- all_e_indices[(actual_train + 1):(actual_train + actual_val)]
e_test_indices       <- all_e_indices[(actual_train + actual_val + 1):(max_available)]

set.seed(43)
all_p_indices <- sample(1:num_proton_rows, max_available, replace = FALSE)
p_training_indices   <- all_p_indices[1:actual_train]
p_validation_indices <- all_p_indices[(actual_train + 1):(actual_train + actual_val)]
p_test_indices       <- all_p_indices[(actual_train + actual_val + 1):(max_available)]

training_data   <- rbind(electron_data_filtered[e_training_indices, ], proton_data_filtered[p_training_indices, ])
validation_data <- rbind(electron_data_filtered[e_validation_indices, ], proton_data_filtered[p_validation_indices, ])
test_data       <- rbind(electron_data_filtered[e_test_indices, ], proton_data_filtered[p_test_indices, ])

# Save each data frame to a separate text file (.csv)
write.csv(training_data, "training_data.csv", row.names = FALSE)
write.csv(validation_data, "validation_data.csv", row.names = FALSE)
write.csv(test_data, "test_data.csv", row.names = FALSE)

training_data$Particle_Type <- factor(training_data$Particle_Type, levels = c("Electron", "Proton"))
validation_data$Particle_Type <- factor(validation_data$Particle_Type, levels = c("Electron", "Proton"))
test_data$Particle_Type <- factor(test_data$Particle_Type, levels = c("Electron", "Proton"))

#Identify ideal HERT test data
HERT_data_electrons <- test_data[test_data$Particle_Type=='Electron'&test_data$E_Inc>=0.6&test_data$E_Inc<=7.5,]
HERT_data_protons <- test_data[test_data$Particle_Type=='Proton'&test_data$E_Inc>=14&test_data$E_Inc<=70,]
HERT_data <- rbind(HERT_data_electrons, HERT_data_protons)
HERT_data <- test_data

# Combine 10% of data
set.seed(1)
training_data_sample    <- training_data[sample(1:nrow(training_data), 
                                                size = actual_train*0.1,
                                                replace = FALSE),] 
# Combine 1% of data
set.seed(2)
training_data_sample2     <- training_data[sample(1:nrow(training_data), 
                                                  size = actual_train*0.01,
                                                  replace = FALSE),]

# Start Dual Logger
sink("svm_model_output_R.txt", split=TRUE)

### Linear SVM ### SV = 3012 
# Create SVM with optimized cost
svm_linear <- svm(Particle_Type ~ Detector1+Detector2+Detector3+Detector4+Detector5+Detector6+Detector7_8_sum+Detector9,
                  data = training_data_sample2, 
                  kernel = "linear", scale=FALSE, cost = 10,
                  class.weights = c("Electron" = 1,"Proton" = 1))
# Print a summary of the model
summary(svm_linear)
# Get the coefficients of the hyperplane
linear_hp_coefs <- coef(svm_linear)
print(linear_hp_coefs)
weights <- abs(linear_hp_coefs/linear_hp_coefs[[1]])
print(weights)
# Get the support vectors as a data frame
linear_sv <- svm_linear$SV
# Predict to validate model
linear_predictions <- predict(svm_linear, HERT_data)
# Print predictions
prop.table(table(predict = linear_predictions, truth = HERT_data[, 11]),2)*100

# Test the hyperplane
linear_hp_test <- factor((linear_hp_coefs[1]
                          + linear_hp_coefs[2] * HERT_data$Detector1[HERT_data$Detector1>0.1]
                          + linear_hp_coefs[3] * HERT_data$Detector2[HERT_data$Detector1>0.1]
                          + linear_hp_coefs[4] * HERT_data$Detector3[HERT_data$Detector1>0.1]
                          + linear_hp_coefs[5] * HERT_data$Detector4[HERT_data$Detector1>0.1]
                          + linear_hp_coefs[6] * HERT_data$Detector5[HERT_data$Detector1>0.1]
                          + linear_hp_coefs[7] * HERT_data$Detector6[HERT_data$Detector1>0.1]
                          + linear_hp_coefs[8] * HERT_data$Detector7_8_sum[HERT_data$Detector1>0.1]
                          + linear_hp_coefs[9] * HERT_data$Detector9[HERT_data$Detector1>0.1])>0,
                         levels = c(TRUE, FALSE), labels = c("Electron", "Proton"))
prop.table(table(predict = linear_hp_test, truth = HERT_data[HERT_data$Detector1 > 0.1, 11]), 2) * 100

### Polynomial SVM ###
# Create SVM with optimized cost
svm_poly <- svm(Particle_Type ~ Detector1+Detector2+Detector3+Detector4+Detector5+Detector6+Detector7_8_sum+Detector9,
                data = training_data, 
                kernel = "polynomial", scale=FALSE, cost = 10)
# Print a summary of the model
summary(svm_poly)
# Get the support vectors as a data frame
poly_sv <- svm_poly$SV
# Predict to validate model
poly_predictions <- predict(svm_poly, HERT_data)
# Print predictions (0: proton, 1: electron)
prop.table(table(predict = poly_predictions, truth = HERT_data[, 11]),2) * 100


### Radial SVM ###
# Create SVM with optimized cost
svm_radial <- svm(Particle_Type ~ Detector1+Detector2+Detector3+Detector4+Detector5+Detector6+Detector7_8_sum+Detector9,
                  data = training_data, 
                  kernel = "radial", scale=FALSE, cost = 10,
                  class.weights = c("Electron" = 1,"Proton" = 1))
# Print a summary of the model
summary(svm_radial)
# Get the support vectors as a data frame
radial_sv <- svm_radial$SV
# Predict to validate model
radial_predictions <- predict(svm_radial, HERT_data)
# Print predictions (0: proton, 1: electron)
prop.table(table(predict = radial_predictions, truth = HERT_data[, 11]),2)*100

### Simplified Linear SVM ###
svm_linearsi <- svm(Particle_Type ~ Detector1+Detector2,
                    data = training_data_sample2, 
                    kernel = "linear", scale=FALSE, cost = 10)
# Print a summary of the model
summary(svm_linearsi)
# Get the coefficients of the hyperplane
linearsi_hp_coefs <- coef(svm_linearsi)
print(linearsi_hp_coefs)
# Get the support vectors as a data frame
linearsi_sv <- svm_linearsi$SV
# Test model
linearsi_test <- predict(svm_linearsi, HERT_data)
# Print test results (0: proton, 1: electron)
prop.table(table(predict = linearsi_test, truth = HERT_data[, 11]),2) * 100

# Test the simplified hyperplane
linearsi_hp_test <- factor((linearsi_hp_coefs[1]
                            + linearsi_hp_coefs[2] * HERT_data$Detector1
                            + linearsi_hp_coefs[3] * HERT_data$Detector2)>0,
                           levels = c(TRUE, FALSE), labels = c("Electron", "Proton"))
prop.table(table(predict=linearsi_hp_test, truth= HERT_data[,11]),2) * 100

### Slant and logic equations from Khoo 2022 for REPTile-2 ###
REPTile2_data_electrons <- electron_data_filtered[electron_data_filtered$Particle_Type=='Electron'&electron_data_filtered$E_Inc>=0.3&electron_data_filtered$E_Inc<=4,]
REPTile2_data_protons <- proton_data_filtered[proton_data_filtered$Particle_Type=='Proton'&proton_data_filtered$E_Inc>=6.7&proton_data_filtered$E_Inc<=35,]
REPTile2_data <- rbind(REPTile2_data_electrons, REPTile2_data_protons)
REPTile2_data <- test_data

slant_eq_D12 <- (REPTile2_data$Detector1/2.8 + REPTile2_data$Detector2/4.2)>1
slant_eq_D34 <- (REPTile2_data$Detector3/13.5 + REPTile2_data$Detector2/30)>1
rng_p <- (REPTile2_data$Detector1>0.1 &  slant_eq_D12 & (!REPTile2_data$Detector4>0.1 | slant_eq_D12) & sum(REPTile2_data[,1:9]<=35))
pen_p <- (REPTile2_data$Detector1>0.1 &  slant_eq_D12 & (REPTile2_data$Detector4>0.1 | !slant_eq_D12) & sum(REPTile2_data[,1:9]<=35))
rng_e <- (REPTile2_data$Detector1>0.1 & !slant_eq_D12 & !REPTile2_data$Detector4>0.1 & sum(REPTile2_data[,1:9]<=4))
pen_e <- (REPTile2_data$Detector1>0.1 & !slant_eq_D12 & REPTile2_data$Detector4>0.1 & sum(REPTile2_data[,1:9]<=4))

khoo_rngetab <- table(predict=factor(rng_e[REPTile2_data$Detector1>0.1],
                                     levels = c(TRUE, FALSE), labels = c("Electron", "Rejected RNG_E")),
                      truth= REPTile2_data[REPTile2_data$Detector1>0.1,11])
prop.table(khoo_rngetab, 2) * 100
khoo_penetab <- table(predict=factor(pen_e[REPTile2_data$Detector1>0.1],
                                     levels = c(TRUE, FALSE), labels = c("Electron", "Rejected PEN_E")),
                      truth= REPTile2_data[REPTile2_data$Detector1>0.1,11])
prop.table(khoo_penetab, 2) * 100

khoo_etab <- table(predict=factor(rng_e | pen_e,
                                  levels = c(TRUE, FALSE), labels = c("Electron", "Rejected Electron")),
                   truth= REPTile2_data[,11])
prop.table(khoo_etab, 2) * 100

khoo_rngptab <- table(predict=factor(rng_p,
                                     levels = c(TRUE, FALSE), labels = c("Proton", "Rejected RNG_P")),
                      truth= REPTile2_data[,11])
prop.table(khoo_rngptab, 2) * 100
khoo_penptab <- table(predict=factor(pen_p,
                                     levels = c(TRUE, FALSE), labels = c("Proton", "Rejected PEN_P")),
                      truth= REPTile2_data[,11])
prop.table(khoo_penptab, 2) * 100
khoo_ptab <- table(predict=factor(rng_p | pen_p,
                                  levels = c(TRUE, FALSE), labels = c("Proton", "Rejected Proton")),
                   truth= REPTile2_data[,11])
prop.table(khoo_ptab, 2) * 100

REPTile2_logic <- rbind(khoo_etab, khoo_ptab)
REPTile2_logic <- (sweep(REPTile2_logic, 2, colSums(REPTile2_logic), FUN = "/")*100 
                   + sweep(REPTile2_logic, 2, colSums(REPTile2_logic), FUN = "/")*100)
print(REPTile2_logic)

### Logic Equations from Baker 2013 for REPT ###
REPT_data_electrons <- test_data[test_data$Particle_Type=='Electron'&test_data$E_Inc>=1.6&test_data$E_Inc<=18.9,]
REPT_data_protons <- proton_data_filtered[proton_data_filtered$Particle_Type=='Proton'&proton_data_filtered$E_Inc>=18&proton_data_filtered$E_Inc<=75,]
REPT_data <- rbind(REPT_data_electrons, REPT_data_protons)
REPT_data <- test_data
# Creating initial logic functions
Rxy <- function(x,y){
  Rxy_result <- rep(0, length(REPT_data$Detector1))
  for (l in x:y){
    Rxy_result <- Rxy_result + REPT_data[[l]]
  }
  return(Rxy_result)
}
Rbarexy <- function(x,y){
  Rbarexy_result <- rep(TRUE, length(REPT_data$Detector1))
  for (l in x:y){
    Rbarexy_result <- Rbarexy_result & REPT_data[[l]] < 0.4
  }
  return(Rbarexy_result)
}
Rbarpxy <- function(x,y){
  Rbarpxy_result <- rep(TRUE, length(REPT_data$Detector1))
  for (l in x:y){
    Rbarpxy_result <- Rbarpxy_result & REPT_data[[l]] < 0.5
  }
  return(Rbarpxy_result)
}
# Electron Logic Equations
EL1  <- (REPT_data$Detector1>=1.0 & REPT_data$Detector1<=1.2 & REPT_data$Detector2<=1.5
         & Rxy(1,2) >= 1.1 & Rxy(1,2)<=1.2 & Rbarexy(3,9))
EL2  <- (REPT_data$Detector1>=0.4 & REPT_data$Detector2>=0.4 & Rxy(1,2)>=1.3 
         & Rxy(1,2)<=1.7 & Rbarexy(3,9))
EL3  <- (REPT_data$Detector1>=0.4 & REPT_data$Detector2>=0.4 & Rxy(1,4)>=1.85 
         & Rxy(1,4)<=2.25 & Rbarexy(5,9))
EL4  <- (REPT_data$Detector1>=0.4 & REPT_data$Detector2>=0.4 & Rxy(1,4)>=2.65 
         & Rxy(1,4)<=2.95 & Rbarexy(5,9))
EL5  <- (REPT_data$Detector1>=0.4 & Rxy(2,4)>=0.4 & Rxy(1,6)>=3.35 & Rxy(1,6)<=3.95 & Rbarexy(7,9))
EL6  <- (REPT_data$Detector1>=0.4 & Rxy(2,6)>=0.4 & Rxy(1,8)>=4.4 & Rxy(1,8)<=5.0 & REPT_data$Detector9<0.4)
EL7  <- (REPT_data$Detector1>=0.4 & REPT_data$Detector1<=2.0 & REPT_data$Detector2>=0.4 
         & REPT_data$Detector2<=2.0 & Rxy(3,6)>=0.4 & Rxy(1,8)>=5.5 & Rxy(1,8)<=6.25 & REPT_data$Detector9<0.4)
EL8  <- (REPT_data$Detector1>=0.4 & REPT_data$Detector2>=0.4 & REPT_data$Detector2<=1.0 
         & Rxy(3,6)>=2.4 & Rxy(3,9)>=5.75 & Rxy(3,9)<=6.6)
EL9  <- (REPT_data$Detector1>=0.4 & REPT_data$Detector2>=0.4 & REPT_data$Detector2<=1.0
         & Rxy(3,4)>=0.4 & Rxy(3,4)<=2.0 & Rxy(5,6)>=0.4 & Rxy(7,9)>=0.4 
         & Rxy(3,9)>=8.0 & Rxy(3,9)<=9.0)
EL10 <- (REPT_data$Detector1>=0.4 & REPT_data$Detector2>=0.4 & Rxy(3,4)>=0.4 
         & Rxy(5,6)>=0.4 & Rxy(7,8)>=0.4 & REPT_data$Detector9>=0.1 
         & Rxy(3,9)>=10.3 & Rxy(3,9)<=12.5)
EL11 <- (REPT_data$Detector1>=0.4 & REPT_data$Detector1<=1.0 & REPT_data$Detector2>=0.4 
         & Rxy(3,4)>=0.4 & Rxy(5,9)>=0.4 & Rxy(7,9)>=11)
EL12 <- (REPT_data$Detector1>=0.4 & REPT_data$Detector1<=1.0 & REPT_data$Detector2>=0.4
         & REPT_data$Detector2<=1.0 & Rxy(3,4)>=0.4 & Rxy(3,4)<=1.5 
         & Rxy(5,9) >= 0.4 & Rxy(7,9)>=15)

ELOGIC <- data.frame(EL1, EL2, EL3, EL4, EL5, EL6, EL7, EL8, EL9, EL10, EL11, EL12)
Erow_predict <- factor(rowSums(ELOGIC[REPT_data$Detector1>0.4,])>0,
                       levels = c(TRUE, FALSE), 
                       labels = c("Electron","Rejected Electron"))

Enproblem_rows <- sum(rowSums(ELOGIC[REPT_data$Detector1>0.4,]) > 1)
ratio_Eproblem_rows <- Enproblem_rows/length(REPT_data$Detector1>0.4)
print(ratio_Eproblem_rows*100)

REPT_etab <- table(predict=factor(Erow_predict,
                                  levels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton"),
                                  labels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton")),
                   truth= factor(REPT_data[REPT_data$Detector1>0.4,11],
                                 levels = c("Electron", "Proton"),
                                 labels = c("Electron", "Proton")))
prop.table(REPT_etab, 2) * 100

# Proton Logic Equations
PL1  <- (REPT_data$Detector1>8.2 & REPT_data$Detector1<16 & REPT_data$Detector2<6.5
         & Rxy(1,2)>8.2 & Rxy(1,2)<18 & Rxy(3,9)<0.5 & Rbarpxy(3,9))
PL2  <- (REPT_data$Detector1>5.4 & REPT_data$Detector1<12.2 & REPT_data$Detector2>5.0
         & REPT_data$Detector2<16.9 & Rxy(3,4)>0.1 & Rxy(3,4)<11 & Rxy(1,4)>15.9 
         & Rxy(1,4)<25.7 & Rxy(5,9)<0.5 & Rbarpxy(5,9))
PL3  <- (REPT_data$Detector1>4 & REPT_data$Detector1<7 & REPT_data$Detector2>4
         & REPT_data$Detector2<9.5 & Rxy(3,4)>10 & Rxy(5,6)<12.5
         & Rxy(1,6)>24 & Rxy(1,6)<35.5 & Rxy(7,9)<0.5 & Rbarpxy(7,9))
PL4  <- (REPT_data$Detector1>3.1 & REPT_data$Detector1<4.9 & REPT_data$Detector2>3.2
         & REPT_data$Detector2<5.7 & Rxy(3,4)>7.6 & Rxy(3,4)<16.8
         & Rxy(5,6)>9.2 & Rxy(5,6)<24 & Rxy(7,8)<23 
         & REPT_data$Detector9<4.1 & Rxy(5,9)>11.5 & Rxy(5,9)<33.0)
PL5  <- (REPT_data$Detector1>2.2 & REPT_data$Detector1<4 & REPT_data$Detector2>1.9
         & REPT_data$Detector2<4.2 & Rxy(3,4)>5.5 & Rxy(3,4)<12.5
         & Rxy(5,6)>5.8 & Rxy(5,6)<12.5 & Rxy(7,8)>7 & Rxy(7,8)<22.7 
         & REPT_data$Detector9>1 & REPT_data$Detector9<13 & Rxy(7,9)>5 & Rxy(7,9)<45)
PL6  <- (REPT_data$Detector1>1.5 & REPT_data$Detector1<3.3 & REPT_data$Detector2>1.0
         & REPT_data$Detector2<3.3 & Rxy(3,4)>4.1 & Rxy(3,4)<6.5 & Rxy(5,6)>4.5
         & Rxy(5,6)<7.2 & Rxy(7,8)>4.8 & Rxy(7,8)<8.0 & REPT_data$Detector9>2.0
         & REPT_data$Detector9<8.5 & Rxy(1,6)>11 & Rxy(1,6)<22 & Rxy(1,9)<65)
PL7  <- (REPT_data$Detector1>1.4 & REPT_data$Detector1<2.5 & REPT_data$Detector2>1.4
         & REPT_data$Detector2<2.8 & Rxy(3,4)>3.4 & Rxy(3,4)<5.4 & Rxy(5,6)>3.4 
         & Rxy(5,6)<5.9 & Rxy(7,8)>3.5 & Rxy(7,8)<6 & Rxy(1,9)>10 & Rxy(1,9)<45)
PL8  <- (REPT_data$Detector1>0.8 & REPT_data$Detector1<3.0 & REPT_data$Detector2>0.8
         & REPT_data$Detector2<3.0 & Rxy(3,4)>2.5 & Rxy(3,4)<5 & Rxy(5,6)>2.5 
         & Rxy(5,6)<5.5 & Rxy(7,8)>2.5 & Rxy(7,8)<5.5 & REPT_data$Detector9>1 
         & REPT_data$Detector9<6 & Rxy(1,9)<8 & Rxy(1,9)<32.0)

PLOGIC <- data.frame(PL1, PL2, PL3, PL4, PL5, PL6, PL7, PL8)
Prow_predict <- factor(rowSums(PLOGIC[REPT_data$Detector1>0.4,])>0,
                       levels = c(TRUE, FALSE), 
                       labels = c("Proton","Rejected Proton"))

Pnproblem_rows <- sum(rowSums(PLOGIC[REPT_data$Detector1>0.4,]) > 1)
ratio_Pproblem_rows <- Pnproblem_rows/length(REPT_data$Detector1)
print(ratio_Pproblem_rows*100)

REPT_ptab <- table(predict=factor(Prow_predict,
                                  levels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton"),
                                  labels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton")),
                   truth= factor(REPT_data[REPT_data$Detector1>0.4,11],
                                 levels = c("Electron", "Proton"),
                                 labels = c("Electron", "Proton")))
prop.table(REPT_ptab, 2) * 100

REPT_logic <- (sweep(REPT_etab, 2, colSums(REPT_etab), FUN = "/")*100 
               + sweep(REPT_ptab, 2, colSums(REPT_ptab), FUN = "/")*100)
prop.table(REPT_logic, 2) * 200

# Stop Dual Logger
sink()


### Plot Hyperplane ###
test_data_plot <- test_data
columns_to_keep <- c("Detector1", "Detector2", "Detector3", "Detector4", "Detector5", "Detector6", "Detector7_8_sum", "Detector9", "Particle_Type")
test_data_plot <- test_data_plot[, columns_to_keep]

# Define binning edges globally so find_max_2d_density can reference them
nbins <- 100
edges <<- logspace(log10(0.1), log10(max(test_data_plot[,1:8], na.rm = TRUE)), nbins+1)

find_max_2d_density <- function(data, particle_type, edges, cols_to_use) {
  subset_data <- data[data$Particle_Type == particle_type, cols_to_use]
  max_density <- 0
  for (i in 1:(length(cols_to_use) - 1)) {
    for (j in (i + 1):length(cols_to_use)) {
      bin_x <- cut(subset_data[, i], breaks = edges)
      bin_y <- cut(subset_data[, j], breaks = edges)
      pair_table <- table(bin_x, bin_y)
      pair_max_frac <- max(pair_table) / sum(pair_table)
      if (pair_max_frac > max_density) {
        max_density <- pair_max_frac
      }
    }
  }
  return(ceiling(max_density * 2000) / 2000)
}

columns = c(1:8)
matrix_size <- 5 * length(columns)

# Calculate e_max and p_max BEFORE plotting to prevent ggpairs from crashing
e_max <<- find_max_2d_density(test_data_plot, "Electron", edges, columns)
p_max <<- find_max_2d_density(test_data_plot, "Proton", edges, columns)

binning_fn <- function(data){
  ecounts <- lapply(c(1:8), function(x) {
    ebins <- cut(data[data$Particle_Type == "Electron",x], breaks = edges)
    data.frame(table(ebins))
  })
  ecounts_new <- matrix(0, nrow = nbins, ncol = length(ecounts))
  max_ecounts <- matrix(0, nrow = 1, ncol = length(ecounts))
  for (i in 1:length(ecounts)){
    ecounts_new[,i] <- ecounts[[i]][,2]
    max_ecounts[i] <- max(ecounts_new[,i])
  }
  
  pcounts <- lapply(c(1:8), function(x) {
    pbins <- cut(data[data$Particle_Type == "Proton",x], breaks = edges)
    data.frame(table(pbins))
  })
  pcounts_new <- matrix(0, nrow = nbins, ncol = length(pcounts))
  max_pcounts <- matrix(0, nrow = 1, ncol = length(pcounts))
  for (i in 1:length(pcounts)){
    pcounts_new[,i] <- pcounts[[i]][,2]
    max_pcounts[i] <- max(pcounts_new[,i])
  }
  ep_ratio <- max(max_ecounts)/max(max_pcounts)
  return(list(edges = edges, ecounts_new = ecounts_new, pcounts_new = pcounts_new, ep_ratio = ep_ratio))
}

density_fn_electron <- function(data, mapping, pt, ...) {
  ggplot(data = data[data$Particle_Type == "Electron",], mapping = mapping) +
    stat_bin_2d(aes(fill = after_stat(count/sum(count))), breaks = edges) +
    scale_fill_gradient(low = "white", high = "red2", limits = c(0, e_max), na.value = "white") +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 16)) +
    theme_few()
}

density_fn_proton <- function(data, mapping, pt, ...) {
  ggplot(data = data[data$Particle_Type == "Proton",], mapping = mapping) +
    stat_bin_2d(aes(fill = after_stat(count/sum(count))), breaks = edges) +
    scale_fill_gradient(low = "white", high = "blue2", limits = c(0, p_max), na.value = "white") +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 16)) +
    theme_few()
}

diag_fn <- function(data, mapping, pt, ...) {
  mapping_index <- which(names(data) == quo_name(mapping[[1]]))
  returned_data <- binning_fn(data)
  edges       = returned_data[[1]]
  ep_ratio = returned_data[[4]]
  midpoints = matrix(0, nrow = length(edges)-1, ncol = 1)
  for (i in 1:length(edges)-1){
    midpoints[i] = (edges[i+1]+edges[i])/2
  }
  ecounts_new <- returned_data[[2]]
  pcounts_new <- returned_data[[3]]
  plot_data <- data.frame(midpoints, ecounts_new[,mapping_index], pcounts_new[,mapping_index])
  colnames(plot_data) <- c("midpoints","ecounts", "pcounts")
  current_weight <- sprintf("%.3f", weights[mapping_index+1])
  
  ggplot(data = plot_data, aes(midpoints)) +
    geom_line(aes(y = ecounts/sum(ecounts)*100), colour = 'red2', lwd = 2) +
    geom_line(aes(y = pcounts/sum(pcounts)*100), colour = 'blue2', lwd = 2) + 
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 16)) +
    theme_few() +
    annotate("text", x = Inf, y = Inf, label = "Y-Axis: Particle Density (%)  ", 
             hjust = 1, vjust = 2, size = 8, fontface = "italic", color = "darkgray") +
    annotate("text", x = Inf, y = Inf, label = paste("Weight:", current_weight," "), 
             hjust = 1, vjust = 4, size = 8, fontface = "bold")
}

# Generate your main ggpairs plot
my_plot <- ggpairs(test_data_plot,
                   columns = columns,
                   columnLabels = c("Detector 1", "Detector 2", "Detector 3", "Detector 4", "Detector 5", "Detector 6", "Detector 7&8", "Detector 9"),
                   xlab = "Deposited Energy (MeV)",
                   ylab = "Deposited Energy (MeV)",
                   upper = list(continuous = density_fn_proton),
                   diag = list(continuous = diag_fn),
                   lower = list(continuous = density_fn_electron)) +
  theme(axis.text = element_text(size = 40), 
        axis.title = element_text(size = 45, face = "bold"),
        strip.text = element_text(size = 40, color = "black"))

electron_palette <- colorRampPalette(c("white", "red2"))
electron_colors <- electron_palette(100)
proton_palette <- colorRampPalette(c("white", "blue2"))
proton_colors <- proton_palette(100)

plot_elec <- ggplot(data.frame(x=1, y=1, z=c(0, e_max*100)), aes(x, y, fill=z)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = electron_colors, 
    name = "Electron\nDensity\n(%)\n",
    guide = guide_colorbar(barheight = unit(matrix_size*4, "lines"), barwidth = unit(2, "lines"))
  ) +
  theme(legend.position = "right", legend.text = element_text(size = 30), legend.title = element_text(size = 30), legend.margin = margin(1, 2, 2, 1, "cm"))
leg_elec <- get_legend(plot_elec)

plot_prot <- ggplot(data.frame(x=1, y=1, z=c(0, p_max*100)), aes(x, y, fill=z)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = proton_colors, 
    name = "Proton\nDensity\n(%)\n",
    guide = guide_colorbar(barheight = unit(matrix_size*4, "lines"), barwidth = unit(2, "lines"))
  ) +
  theme(legend.position = "right", legend.text = element_text(size = 30), legend.title = element_text(size = 30), legend.margin = margin(1, 2, 2, 1, "cm"))
leg_prot <- get_legend(plot_prot)

combined_legends <- plot_grid(leg_elec, leg_prot, ncol = 2)

final_plot <- plot_grid(
  ggmatrix_gtable(my_plot), 
  NULL, 
  combined_legends, 
  ncol = 3,
  rel_widths = c(1, 0.02, 0.1) 
)

ggsave("square_pairs_plot_with_vertical_legends.png", 
       plot = final_plot, 
       width = matrix_size + (matrix_size * 0.30),
       height = matrix_size, 
       units = "in", 
       dpi = 300,
       limitsize = FALSE)

massive_plot <- image_read("square_pairs_plot_with_vertical_legends.png")
resized_plot <- image_scale(massive_plot, "1920x")
image_write(resized_plot, "square_pairs_1920px.png")