# Prepare packages for SVM
library(e1071)
library(caret)
library(ggplot2)
library(GGally)
library(ggthemes)
library(colorBlindness)
library(pracma)

setwd("C:\\Users\\wzt0020\\Box\\HERT_Box\\Particle Classification")

# Set the file path (replace "your_file.txt" with your actual file path)
electron_file <- "E:\\HERT_Drive\\Matlab Main\\Result\\Electron_FS\\Aggregate Data\\Aggregate Electron_FS Data.txt"
proton_file   <- "E:\\HERT_Drive\\Matlab Main\\Result\\Proton_FS\\Aggregate Data\\Aggregate Proton_FS Data.txt"

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

# Randomly sample row indices from the filtered data.
# Set a seed for reproducibility.
num_electron_rows <- nrow(electron_data_filtered)
set.seed(42)
training_indices   <- sample(1:num_electron_rows, num_training_rows, replace = FALSE)
set.seed(43)
validation_indices <- sample(1:num_electron_rows, num_validation_rows, replace = FALSE)
set.seed(44)
test_indices       <- sample(1:num_electron_rows, num_test_rows, replace = FALSE)

# Select the random rows from the filtered data sets.
num_proton_rows   <- nrow(proton_data_filtered)
training_data   <- rbind(electron_data_filtered[training_indices, ], proton_data_filtered[training_indices, ])
validation_data <- rbind(electron_data_filtered[validation_indices, ], proton_data_filtered[validation_indices, ])
test_data       <- rbind(electron_data_filtered[test_indices, ], proton_data_filtered[test_indices, ])

# Save each data frame to a separate text file (.csv)
write.csv(training_data, "training_data.csv", row.names = FALSE)
write.csv(validation_data, "validation_data.csv", row.names = FALSE)
write.csv(test_data, "test_data.csv", row.names = FALSE)

rm(electron_data_filtered)
rm(proton_data_filtered)

# Read in the saved data frames.
# training_data  <- read.csv("training_data.csv")
# validation_data <- read.csv("validation_data.csv")
# test_data       <- read.csv("test_data.csv")

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
                                                 size = num_training_rows*0.1,
                                                 replace = FALSE),] 
# Combine 1% of data
set.seed(2)
training_data_sample2     <- training_data[sample(1:nrow(training_data), 
                                                  size = num_training_rows*0.01,
                                                  replace = FALSE),]

### Linear SVM ### SV = 3012 
# Tune SVM to find optimized cost
#tuner.lin=tune(svm,Particle_Type ~ .,
#           data=validation_data, kernel="linear",
#           ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100)))
#summary(tuner.lin)

# Create SVM with optimized cost
svm_linear <- svm(Particle_Type ~ Detector1+Detector2+Detector3+Detector4+Detector5+Detector6+Detector7_8_sum+Detector9,
                  data = training_data, 
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
prop.table(table(predict=linear_hp_test, truth=HERT_data[,11]), 2) * 100

# electron_misclassified_mask <- (linear_predictions == 'Proton') & (HERT_data$Particle_Type == "Electron")
# electron_misclassified_indices <- which(electron_misclassified_mask)
# hist(HERT_data[electron_misclassified_mask,1],
#      seq(from = 0, to = 8, by = 0.1),
#      main = '',
#      xlab = "Electron Incident Energy (MeV)",
#      ylab = "Misclassification as Protons",
#      col = "lightblue",
#      border = "black")
# 
# proton_misclassified_mask <- (linear_predictions == 'Electron') & (HERT_data$Particle_Type == "Proton")
# proton_misclassified_indices <- which(proton_misclassified_mask)
# hist(HERT_data[proton_misclassified_mask,1],
#      seq(from = 0, to = 200, by = 1),
#      main = '',
#      xlab = "Proton Incident Energy (MeV)",
#      ylab = "Misclassification as Electrons",
#      col = "lightblue",
#      border = "black")
# hist(HERT_data[proton_misclassified_mask,2],
#      seq(from = 0, to = 5, by = 0.1),
#      main = '',
#      xlab = "Proton D1 Energy (MeV)",
#      ylab = "Misclassification as Electrons",
#      col = "lightblue",
#      border = "black")
# hist(rowSums(HERT_data[proton_misclassified_mask,3:10]),
#      breaks = seq(from = 0, to = 11, by = 1),
#      main = '',
#      xlab = "Proton D2-D9 Energy (MeV)",
#      ylab = "Misclassification as Electrons",
#      col = "lightblue",
#      border = "black")
# hist(HERT_data[proton_misclassified_mask,][rowSums(HERT_data[proton_misclassified_mask,3:10])>0,1],
#      breaks = seq(from = 0, to = 200, by = 1),
#      main = '',
#      xlab = "Proton D2-D9 Energy (MeV)",
#      ylab = "Misclassification as Electrons",
#      col = "lightblue",
#      border = "black")

### Polynomial SVM ###
# Tune SVM to find optimized cost
#tuner.poly=tune(svm,Particle_Type ~ Detector1+Detector2+Detector3+Detector4+Detector5+Detector6+Detector7_8_sum+Detector9,,
#               data=validation_data, kernel="polynomial",
#               ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100,1000)))
#summary(tuner.poly)

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
# Tune SVM to find optimized cost
#tuner.rad=tune(svm,Particle_Type ~ Detector1+Detector2+Detector3+Detector4+Detector5+Detector6+Detector7_8_sum+Detector9,
#               data=validation_data, kernel="radial",
#               ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100,1000)))
#summary(tuner.rad)

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

# proton_misclassified_mask <- (pen_e[REPTile2_data$Detector1 > 0.1] == TRUE) & (REPTile2_data[REPTile2_data$Detector1 > 0.1,]$Particle_Type == "Proton")
# proton_misclassified_indices <- which(proton_misclassified_mask)
# hist(REPTile2_data[proton_misclassified_mask,1],
#      seq(from = 0, to = 200, by = 1),
#      main = '',
#      xlab = "Proton Incident Energy (MeV)",
#      ylab = "Misclassification as Penetrating Electrons",
#      col = "lightblue",
#      border = "black")

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
# Electron Logic Equations (with veto constraints)
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
Etotal_true <- sum(ELOGIC[REPT_data$Detector1>0.4,])
Erow_predict <- factor(rowSums(ELOGIC[REPT_data$Detector1>0.4,])>0,
                       levels = c(TRUE, FALSE), 
                       labels = c("Electron","Rejected Electron"))
# How many particles meet more than one logic condition?
Enproblem_rows <- sum(rowSums(ELOGIC[REPT_data$Detector1>0.4,]) > 1)
Eproblem_rows <- which(rowSums(ELOGIC[REPT_data$Detector1>0.4,]) > 1)
ratio_Eproblem_rows <- Enproblem_rows/length(REPT_data$Detector1>0.4)
print(ratio_Eproblem_rows*100)
REPT_etab <- table(predict=factor(Erow_predict,
                                  levels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton"),
                                  labels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton")),
                   truth= factor(REPT_data[REPT_data$Detector1>0.4,11],
                                 levels = c("Electron", "Proton"),
                                 labels = c("Electron", "Proton")))
prop.table(REPT_etab, 2) * 100

# proton_misclassified_mask <- (rowSums(ELOGIC[REPT_data$Detector1>0.4,])>0) & (REPT_data[REPT_data$Detector1 > 0.4,]$Particle_Type == "Electron")
# proton_misclassified_indices <- which(proton_misclassified_mask)
# hist(REPT_data[proton_misclassified_mask,1],
#      seq(from = 0, to = 200, by = 1),
#      main = '',
#      xlab = "Proton Incident Energy (MeV)",
#      ylab = "Misclassification as Electrons",
#      col = "lightblue",
#      border = "black")

# Proton Logic Equations (with veto constraints)
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
Ptotal_true <- sum(PLOGIC[REPT_data$Detector1>0.4,])
Prow_predict <- factor(rowSums(PLOGIC[REPT_data$Detector1>0.4,])>0,
                       levels = c(TRUE, FALSE), 
                       labels = c("Proton","Rejected Proton"))
# How many particles meet more than one logic condition?
Pnproblem_rows <- sum(rowSums(PLOGIC[REPT_data$Detector1>0.4,]) > 1)
Pproblem_rows <- which(rowSums(PLOGIC[REPT_data$Detector1>0.4,]) > 1)
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
print(ratio_Eproblem_rows*100)
print(ratio_Pproblem_rows*100)


### Plot Hyperplane ###
# Plot Edep v Einc for each detector, color coding particles
par(mar = c(5, 6, 4, 2) + 0.1)
plot(test_data$E_Inc, test_data$Detector9, 
     pch = 20,cex = 0.25, # Solid circle for all points
     #xlim = c(0,20), ylim = c(0,20),
     xlab = 'Incident Energy (MeV)', ylab = 'E_dep on Detector 9 (MeV)',
     col = ifelse(test_data$Particle_Type == "Electron", "red","blue"),  # Color based on particle type
     cex.axis = 2.0, # Adjust axis tick label size
     cex.lab = 2.0)  # Adjust axis label size

legend("topright", # Position the legend in the top left corner.
       legend = c("Electron", "Proton"), # Labels for the legend.
       col = c("red", "blue"), # Colors for the legend.
       pch = 20, # Use the same symbol as in the plot.
       bty = "n", # Do not draw a box around the legend.
       cex = 2.0) # Adjust the size of the legend text.

# abline(a = -linear_hp_coefs[1]/linear_hp_coefs[4], 
#        b = -linear_hp_coefs[2]/linear_hp_coefs[4], lwd = 2, lty = 2)
# legend("topright", legend = c("Linear Model", "SVM"), lty = c(2,3), pch = 20)

test_data_plot <- test_data
columns_to_keep <- c("Detector1", "Detector2", "Detector3", "Detector4", "Detector5", "Detector6", "Detector7_8_sum", "Detector9", "Particle_Type")
test_data_plot <- test_data_plot[, columns_to_keep]

pairs(~ Detector1+Detector2+Detector3+Detector4+Detector5
      +Detector6+Detector7_8_sum+Detector9, data=test_data_plot,
      panel=function(x,y){
        points(x,y,col = ifelse(test_data_plot$Particle_Type == "Electron", "red","blue"),
               xlim = c(0,20), ylim = c(0,20),
               pch = 20, cex = 0.25 # Solid circle for all points
        )
        for (i in 1:8){
          if (sum(x-test_data_plot[i])==0){
            for (j in 1:8){
              if (sum(y-test_data_plot[j])==0){
                abline(a = -linear_hp_coefs[1]/linear_hp_coefs[j+1],
                       b = -linear_hp_coefs[i+1]/linear_hp_coefs[j+1], lwd = 3, lty = 2)
              }
            }
          }
        }
      })

output_col_names <- names(test_data_plot)[1:(ncol(test_data_plot) - 1)]
test_data_means_electrons <- matrix(0, nrow = length(output_col_names), ncol = length(output_col_names))
test_data_means_protons <- matrix(0, nrow = length(output_col_names), ncol = length(output_col_names))
colnames(test_data_means_electrons) <- output_col_names
colnames(test_data_means_protons) <- output_col_names
test_data_medians_electrons <- test_data_means_electrons
test_data_medians_protons <- test_data_means_protons
# Loop through each detector column.
for (i in 1:length(output_col_names)) {
  for (j in 1:i) {
    # Filter the data based on particle type.
    electrons <- test_data_plot[test_data_plot$Particle_Type == "Electron", ]
    protons <- test_data_plot[test_data_plot$Particle_Type == "Proton", ]
    
    if ((i+1) < length(output_col_names)) {
      # Filter to find rows where current column is non-zero
      # AND all higher columns are zero.
      mask_electrons <- (electrons[, j] != 0) & (rowSums(electrons[, (i + 1):length(output_col_names)]) == 0)
      mask_protons <- (protons[, j] != 0) & (rowSums(protons[, (i + 1):length(output_col_names)]) == 0)
    } else if ((i+1) == length(output_col_names)){
      mask_electrons <- (electrons[, j] != 0) & (electrons[, length(output_col_names)] == 0)
      mask_protons <- (protons[, j] != 0) & (protons[, length(output_col_names)] == 0)
    } else if (i == length(output_col_names)){
      mask_electrons <- (electrons[, j] != 0)
      mask_protons <- (protons[, j] != 0)
    }
    # Calculate the mean and assign it to the correct row and column.
    test_data_means_electrons[i, j] <- mean(electrons[mask_electrons, j])
    test_data_means_protons[i, j] <- mean(protons[mask_protons, j])
    # Calculate the median and assign it to the correct row and column.
    test_data_medians_electrons[i, j] <- median(electrons[mask_electrons, j])
    test_data_medians_protons[i, j] <- median(protons[mask_protons, j])
  }
}
print(test_data_means_electrons)
print(test_data_means_protons)
print(test_data_medians_electrons)
print(test_data_medians_protons)

binning_fn <- function(data){
  nbins <- 100
  #finding bin edges
  edges = logspace(log10(0.1), log10(max(data[,1:8])), nbins+1)
  
  # Binning electron and proton counts
  ecounts <- lapply(c(1:8), function(x) {
    ebins <- cut(data[data$Particle_Type == "Electron",x], 
                 breaks = edges)
    data.frame(table(ebins))
  })
  ebins_midpoints <- matrix(0, nrow = nbins, ncol = length(ecounts))
  ecounts_new <- matrix(0, nrow = nbins, ncol = length(ecounts))
  max_ecounts <- matrix(0, nrow = 1, ncol = length(ecounts))
  colnames(ecounts_new) <- c("Detector1", "Detector2", "Detector3", "Detector4", "Detector5", "Detector6", "Detector7_8_sum", "Detector9")
  colnames(max_ecounts) <- c("Detector1", "Detector2", "Detector3", "Detector4", "Detector5", "Detector6", "Detector7_8_sum", "Detector9")
  for (i in 1:length(ecounts)){
    ecounts_new[,i] <- ecounts[[i]][,2]
    max_ecounts[i] <- max(ecounts_new[,i])
  }
  
  pcounts <- lapply(c(1:8), function(x) {
    pbins <- cut(data[data$Particle_Type == "Proton",x], 
                 breaks = edges)
    data.frame(table(pbins))
  })
  pbins_midpoints <- matrix(0, nrow = nbins, ncol = length(pcounts))
  pcounts_new <- matrix(0, nrow = nbins, ncol = length(pcounts))
  max_pcounts <- matrix(0, nrow = 1, ncol = length(pcounts))
  colnames(pcounts_new) <- c("Detector1", "Detector2", "Detector3", "Detector4", "Detector5", "Detector6", "Detector7_8_sum", "Detector9")
  colnames(max_pcounts) <- c("Detector1", "Detector2", "Detector3", "Detector4", "Detector5", "Detector6", "Detector7_8_sum", "Detector9")
  for (i in 1:length(pcounts)){
    pcounts_new[,i] <- pcounts[[i]][,2]
    max_pcounts[i] <- max(pcounts_new[,i])
  }
  
  ep_ratio <- max(max_ecounts)/max(max_pcounts)
  
  return(list(edges = edges, ecounts_new = ecounts_new, pcounts_new = pcounts_new, ep_ratio = ep_ratio))
}

density_fn_electron <- function(data, mapping, pt, ...) {
  # Access variable indices based on mapping
  map_i <- which(names(data) == quo_name(mapping[[1]]))  # Index of x variable
  map_j <- which(names(data) == quo_name(mapping[[2]]))  # Index of y variable
  
  max_detect <- max(map_i,map_j)
  temp_boundary_guess <- tcounts_min_E
  
  ggplot(data = data[data$Particle_Type == "Electron",], mapping = mapping) +
    stat_bin_2d(aes(fill = after_stat(count/num_test_rows)), breaks = edges, limits=c(0,e_max)) +
    scale_fill_gradient(low = "yellow", high = "red1") +
    geom_abline(intercept = -sum(c(linear_hp_coefs[1],temp_boundary_guess[-c(map_j,map_i)]*linear_hp_coefs[-c(1,map_j+1,map_i+1)])) / linear_hp_coefs[map_j + 1],
                slope = -linear_hp_coefs[map_i + 1] / linear_hp_coefs[map_j + 1], lwd = 2, lty = 1) +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 16)) +
    theme_few()
}

density_fn_proton <- function(data, mapping, pt, ...) {
  # Access variable indices based on mapping
  map_i <- which(names(data) == quo_name(mapping[[1]]))  # Index of x variable
  map_j <- which(names(data) == quo_name(mapping[[2]]))  # Index of y variable
  
  max_detect <- max(map_i,map_j)
  temp_boundary_guess <- tcounts_min_E
  
  ggplot(data = data[data$Particle_Type == "Proton",], mapping = mapping) +
    stat_bin_2d(aes(fill = after_stat(count/num_test_rows)), breaks = edges, limits=c(0,p_max)) +
    scale_fill_gradient(low = "#00E2FF", high = "blue1") +
    geom_abline(intercept = -sum(c(linear_hp_coefs[1],temp_boundary_guess[-c(map_j,map_i)]*linear_hp_coefs[-c(1,map_j+1,map_i+1)])) / linear_hp_coefs[map_j + 1],
                slope = -linear_hp_coefs[map_i + 1] / linear_hp_coefs[map_j + 1], lwd = 2, lty = 1) +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 16)) +
    theme_few()
}

diag_fn <- function(data, mapping, pt, ...) {
  # Access variable indices based on mapping
  mapping_index <- which(names(data) == quo_name(mapping[[1]]))  # Index of x variable
  
  #use binning function to bin data
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
  
  ggplot(data = plot_data, aes(midpoints)) +
    geom_line(aes(y = ecounts/sum(ecounts)), colour = 'red1', lwd = 2) +
    geom_line(aes(y = pcounts/sum(pcounts)), colour = 'blue1', lwd = 2) + 
    geom_line(aes(y = pcounts/sum(pcounts)+ecounts/sum(ecounts)), colour = 'black', lwd = 2) +
    geom_vline(xintercept=ecounts_max_E[mapping_index], col='red1', linetype='dashed', lwd = 1) +
    geom_vline(xintercept=pcounts_max_E[mapping_index], col='blue1', linetype='dashed', lwd = 1) +
    geom_vline(xintercept=tcounts_min_E[mapping_index], col='black', linetype='dashed', lwd = 1) +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.05)) +
    theme_few()
}

returned_data <- binning_fn(test_data_plot)
edges <- returned_data[[1]]
ecounts_new <- returned_data[[2]]
ecounts_max_index <- apply(ecounts_new, 2, which.max)
ecounts_max_E <- (edges[ecounts_max_index+1]+edges[ecounts_max_index])/2
e_max <- ceiling(max(ecounts_new/num_test_rows) * 100) / 100
#e_max <- max(ecounts_new/num_test_rows)
pcounts_new <- returned_data[[3]]
pcounts_max_index <- apply(pcounts_new, 2, which.max)
pcounts_max_E <- (edges[pcounts_max_index+1]+edges[pcounts_max_index])/2
p_max <- ceiling(max(pcounts_new/num_test_rows) * 100) / 100
tcounts <- ecounts_new+pcounts_new
tcounts_min_index <- apply(tcounts[ecounts_max_index:pcounts_max_index,], 2, which.min)+ecounts_max_index
tcounts_min_E <- (edges[tcounts_min_index+1]+edges[tcounts_min_index])/2
#p_max <- max(pcounts_new/num_test_rows)

ggpairs(test_data_plot,
        columns = c(1:8),
        #columns = c(1, 2, 6, 8),
        #columns = c(1, 2),
        upper = list(continuous = density_fn_proton),
        diag = list(continuous = diag_fn),
        lower = list(continuous = density_fn_electron)) +
  theme(axis.text = element_text(size = 28), 
        strip.text = element_text(size = 28, color = "black"),
        panel.spacing=grid::unit(1.2,"lines"))

# Creating color scales
electron_palette <- colorRampPalette(c("yellow", "red1"))
electron_colors <- electron_palette(100)

proton_palette <- colorRampPalette(c("#00E2FF", "blue1"))
proton_colors <- proton_palette(100)

v <- ggplot(faithful[1:100,], aes(waiting, eruptions, fill = (seq(e_max, 0, -e_max/99)))) +
  geom_tile()
v + scale_fill_gradientn(colours = electron_colors, name = "Electron Density     ", breaks = c(0, 0.02,0.04,0.06, 0.08)) + 
  theme(legend.key.width  = unit(5, "lines"), legend.position = "bottom",
        legend.text = element_text(size = 26), legend.title = element_text(size = 26))
v <- ggplot(faithful[1:100,], aes(waiting, eruptions, fill = (seq(p_max, 0, -p_max/99)))) +
  geom_tile()
v + scale_fill_gradientn(colours = proton_colors, name = "Proton Density     ", breaks = c(0, 0.02,0.04,0.06)) + 
  theme(legend.key.width  = unit(5/(e_max/p_max), "lines"), legend.position = "bottom",
        legend.text = element_text(size = 26), legend.title = element_text(size = 26))
