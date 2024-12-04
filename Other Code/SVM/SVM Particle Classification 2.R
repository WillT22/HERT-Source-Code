# Prepare packages for SVM
library(e1071)
library(caret)
library(ggplot2)
library(GGally)
library(ggthemes)
library(colorBlindness)

# Set the file path (replace "your_file.txt" with your actual file path)
electron_file <- "C:\\Users\\Will\\Box\\HERT_Box\\Data\\Aggregate Electron Data 1-100.txt"
electron_file <- "C:\\Users\\wzt0020\\Box\\HERT_Box\\Data\\Aggregate Electron Data 1-100.txt"
proton_file   <- "C:\\Users\\Will\\Box\\HERT_Box\\Data\\Aggregate Proton Data 1-50.txt"
proton_file   <- "C:\\Users\\wzt0020\\Box\\HERT_Box\\Data\\Aggregate Proton Data 1-50.txt"

# Read the data
imported_electron_data <- read.table(electron_file, sep="", skip = 1)
imported_proton_data   <- read.table(proton_file, sep="", skip = 1)

electron_data <- imported_electron_data[-1]
proton_data <- imported_proton_data[-1]

# Assign header labels
column_names <- c("Detector1", "Detector2", "Detector3",
                  "Detector4","Detector5","Detector6",
                  "Detector7","Detector8","Detector9") 
colnames(electron_data) <- column_names
colnames(proton_data)   <- column_names

# Create a new column named "Particle Type" with "electron" for values of 1 in column 1
electron_data$Particle_Type <- as.factor("Electron")
proton_data$Particle_Type <- as.factor("Proton")

# Combine data into training, validation, and test groups (80/10/10 split)
training_data     <- rbind(electron_data[1:440000, ],proton_data[1:440000, ])
validation_data <- rbind(electron_data[440001:450000, ],proton_data[440001:450000, ])
test_data         <- rbind(electron_data[450001:500000, ],proton_data[450001:500000, ])

training_data <- training_data[training_data$Detector1>0.1,];
validation_data <- validation_data[validation_data$Detector1>0.1,];
test_data <- test_data[test_data$Detector1>0.1,];

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

### Linear SVM ### SV = 3012 
# Tune SVM to find optimized cost
tuner.lin=tune(svm,Particle_Type ~ .,
           data=validation_data, kernel="linear",
           ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100)))
summary(tuner.lin)

# Create SVM with optimized cost
svm_linear <- svm(Particle_Type ~ .,
                  data = training_data_sample, 
                  kernel = "linear", scale=FALSE, cost = 10,
                  class.weights = c("Electron" = 1,"Proton" = 4))

# Print a summary of the model
summary(svm_linear)

# Get the coefficients of the hyperplane
linear_hp_coefs <- -coef(svm_linear)
print(linear_hp_coefs)

# Get the support vectors as a data frame
linear_sv <- svm_linear$SV

# Predict to validate model
linear_predictions <- predict(svm_linear, test_data[test_data$Detector1>0.1,1:9])

# Print predictions
prop.table(table(predict = linear_predictions, truth = test_data[test_data$Detector1>0.1, 10]),2)*100

linear_hp_test <- factor((linear_hp_coefs[1]
                            + linear_hp_coefs[2] * test_data$Detector1[test_data$Detector1>0.1]
                            + linear_hp_coefs[3] * test_data$Detector2[test_data$Detector1>0.1]
                            + linear_hp_coefs[4] * test_data$Detector3[test_data$Detector1>0.1]
                            + linear_hp_coefs[5] * test_data$Detector4[test_data$Detector1>0.1]
                            + linear_hp_coefs[6] * test_data$Detector5[test_data$Detector1>0.1]
                            + linear_hp_coefs[7] * test_data$Detector6[test_data$Detector1>0.1]
                            + linear_hp_coefs[8] * test_data$Detector7[test_data$Detector1>0.1]
                            + linear_hp_coefs[9] * test_data$Detector8[test_data$Detector1>0.1]
                            + linear_hp_coefs[10] * test_data$Detector9[test_data$Detector1>0.1])>0,
                           levels = c(TRUE, FALSE), labels = c("Electron", "Proton"))

prop.table(table(predict=linear_hp_test, truth=test_data[test_data$Detector1>0.1,10]), 2) * 100

### Polynomial SVM ###
# Tune SVM to find optimized cost
tuner.poly=tune(svm,Particle_Type ~ .,
               data=validation_data, kernel="polynomial",
               ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100,1000)))
summary(tuner.poly)

# Create SVM with optimized cost
svm_poly <- svm(Particle_Type ~ .,
                data = training_data_sample, 
                kernel = "polynomial", scale=FALSE, cost = 10)

# Print a summary of the model
summary(svm_poly)

# Get the support vectors as a data frame
poly_sv <- svm_poly$SV

# Predict to validate model
poly_predictions <- predict(svm_poly, test_data)

# Print predictions (0: proton, 1: electron)
table(predict = poly_predictions, truth = test_data[, 10])


### Radial SVM ###
# Tune SVM to find optimized cost
tuner.rad=tune(svm,Particle_Type ~ .,
               data=validation_data, kernel="radial",
               ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100,1000)))
summary(tuner.rad)

# Create SVM with optimized cost
svm_radial <- svm(Particle_Type ~ .,
                  data = training_data_sample, 
                  kernel = "radial", scale=FALSE, cost = 10,
                  class.weights = c("Electron" = 1,"Proton" = 4))

# Print a summary of the model
summary(svm_radial)

# Get the support vectors as a data frame
radial_sv <- svm_radial$SV

# Predict to validate model
radial_predictions <- predict(svm_radial, test_data[test_data$Detector1>0.1,])

# Print predictions (0: proton, 1: electron)
prop.table(table(predict = radial_predictions, truth = test_data[test_data$Detector1>0.1, 10]),2)*100


### Sigmoid SVM ###
# Tune SVM to find optimized cost
tuner.sig=tune(svm,Particle_Type ~ .,
               data=validation_data, kernel="sigmoid",
               ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100,1000)))
summary(tuner.sig)

# Create SVM with optimized cost
svm_sigmoid <- svm(Particle_Type ~ .,
                  data = training_data_sample, 
                  kernel = "sigmoid", scale=FALSE, cost = 0.1)

# Print a summary of the model
summary(svm_sigmoid)

# Get the support vectors as a data frame
sigmoid_sv <- svm_sigmoid$SV

# Predict to validate model
sigmoid_predictions <- predict(svm_sigmoid, test_data)

# Print predictions (0: proton, 1: electron)
table(predict = sigmoid_predictions, truth = test_data[, 10])


### Simplified Linear SVM ###
svm_linearsi <- svm(Particle_Type ~ Detector1+Detector2+Detector4+Detector9,
                    data = training_data_sample, 
                    kernel = "linear", scale=FALSE, cost = 10)

# Print a summary of the model
summary(svm_linearsi)

# Get the coefficients of the hyperplane
linearsi_hp_coefs <- -coef(svm_linearsi)
print(linearsi_hp_coefs)

# Get the support vectors as a data frame
linearsi_sv <- svm_linearsi$SV

# Test model
linearsi_test <- predict(svm_linearsi, test_data)

# Print test results (0: proton, 1: electron)
table(predict = linearsi_test, truth = test_data[, 10])

linearsi_hp_test <- factor((linearsi_hp_coefs[1]
                                + linearsi_hp_coefs[2] * test_data$Detector1
                                + linearsi_hp_coefs[3] * test_data$Detector2
                                + linearsi_hp_coefs[4] * test_data$Detector4
                                + linearsi_hp_coefs[5] * test_data$Detector9)>0,
                             levels = c(TRUE, FALSE), labels = c("Electron", "Proton"))

table(predict=linearsi_hp_test, truth= test_data[,10])

### Slant and logic equations from Khoo 2022 for REPTile-2 ###
slant_eq_D12 <- (test_data$Detector1/2.8 + test_data$Detector2/4.2)>1
slant_eq_D34 <- (test_data$Detector3/13.5 + test_data$Detector2/30)>1
rng_p <- (test_data$Detector1>0.1 &  slant_eq_D12 & (!test_data$Detector4>0.1 | slant_eq_D12) & sum(test_data[,1:9]<=35))
pen_p <- (test_data$Detector1>0.1 &  slant_eq_D12 & (test_data$Detector4>0.1 | !slant_eq_D12) & sum(test_data[,1:9]<=35))
rng_e <- (test_data$Detector1>0.1 & !slant_eq_D12 & !test_data$Detector4>0.1 & sum(test_data[,1:9]<=4))
pen_e <- (test_data$Detector1>0.1 & !slant_eq_D12 & test_data$Detector4>0.1 & sum(test_data[,1:9]<=4))

khoo_rngetab <- table(predict=factor(rng_e[test_data$Detector1>0.1],
                   levels = c(TRUE, FALSE), labels = c("Electron", "Rejected RNG_E")),
                   truth= test_data[test_data$Detector1>0.1,10])
khoo_penetab <- table(predict=factor(pen_e[test_data$Detector1>0.1],
                                     levels = c(TRUE, FALSE), labels = c("Electron", "Rejected PEN_E")),
                      truth= test_data[test_data$Detector1>0.1,10])
khoo_etab <- table(predict=factor(rng_e[test_data$Detector1>0.1] | pen_e[test_data$Detector1>0.1],
                                     levels = c(TRUE, FALSE), labels = c("Electron", "Rejected Electron")),
                      truth= test_data[test_data$Detector1>0.1,10])
prop.table(khoo_etab, 2) * 100

khoo_rngptab <- table(predict=factor(rng_p[test_data$Detector1>0.1],
                                     levels = c(TRUE, FALSE), labels = c("Proton", "Rejected RNG_P")),
                      truth= test_data[test_data$Detector1>0.1,10])
khoo_penptab <- table(predict=factor(pen_p[test_data$Detector1>0.1],
                                     levels = c(TRUE, FALSE), labels = c("Proton", "Rejected PEN_P")),
                      truth= test_data[test_data$Detector1>0.1,10])
khoo_ptab <- table(predict=factor(rng_p[test_data$Detector1>0.1] | pen_p[test_data$Detector1>0.1],
                                  levels = c(TRUE, FALSE), labels = c("Proton", "Rejected Proton")),
                   truth= test_data[test_data$Detector1>0.1,10])
prop.table(khoo_ptab, 2) * 100

REPTile2_logic <- rbind(khoo_etab, khoo_ptab)
REPTile2_logic <- (sweep(REPTile2_logic, 2, colSums(REPTile2_logic), FUN = "/")*100 
               + sweep(REPTile2_logic, 2, colSums(REPTile2_logic), FUN = "/")*100)
print(REPTile2_logic)

### Logic Equations from Baker 2013 for REPT ###
# Creating initial logic functions
Rxy <- function(x,y){
  Rxy_result <- rep(0, length(test_data$Detector1))
  for (l in x:y){
    Rxy_result <- Rxy_result + test_data[[l]]
  }
  return(Rxy_result)
}
Rbarexy <- function(x,y){
  Rbarexy_result <- rep(TRUE, length(test_data$Detector1))
  for (l in x:y){
    Rbarexy_result <- Rbarexy_result & test_data[[l]] < 0.4
  }
  return(Rbarexy_result)
}
Rbarpxy <- function(x,y){
  Rbarpxy_result <- rep(TRUE, length(test_data$Detector1))
  for (l in x:y){
    Rbarpxy_result <- Rbarpxy_result & test_data[[l]] < 0.5
  }
  return(Rbarpxy_result)
}
# Electron Logic Equations (with veto constraints)
EL1  <- (test_data$Detector1>=1.0 & test_data$Detector1<=1.2 & test_data$Detector2<=1.5
         & Rxy(1,2) >= 1.1 & Rxy(1,2)<=1.2 & Rbarexy(3,9))
EL2  <- (test_data$Detector1>=0.4 & test_data$Detector2>=0.4 & Rxy(1,2)>=1.3 
         & Rxy(1,2)<=1.7 & Rbarexy(3,9))
EL3  <- (test_data$Detector1>=0.4 & test_data$Detector2>=0.4 & Rxy(1,4)>=1.85 
         & Rxy(1,4)<=2.25 & Rbarexy(5,9))
EL4  <- (test_data$Detector1>=0.4 & test_data$Detector2>=0.4 & Rxy(1,4)>=2.65 
         & Rxy(1,4)<=2.95 & Rbarexy(5,9))
EL5  <- (test_data$Detector1>=0.4 & Rxy(2,4)>=0.4 & Rxy(1,6)>=3.35 & Rxy(1,6)<=3.95 & Rbarexy(7,9))
EL6  <- (test_data$Detector1>=0.4 & Rxy(2,6)>=0.4 & Rxy(1,8)>=4.4 & Rxy(1,8)<=5.0 & test_data$Detector9<0.4)
EL7  <- (test_data$Detector1>=0.4 & test_data$Detector1<=2.0 & test_data$Detector2>=0.4 
         & test_data$Detector2<=2.0 & Rxy(3,6)>=0.4 & Rxy(1,8)>=5.5 & Rxy(1,8)<=6.25 & test_data$Detector9<0.4)
EL8  <- (test_data$Detector1>=0.4 & test_data$Detector2>=0.4 & test_data$Detector2<=1.0 
         & Rxy(3,6)>=2.4 & Rxy(3,9)>=5.75 & Rxy(3,9)<=6.6)
EL9  <- (test_data$Detector1>=0.4 & test_data$Detector2>=0.4 & test_data$Detector2<=1.0
         & Rxy(3,4)>=0.4 & Rxy(3,4)<=2.0 & Rxy(5,6)>=0.4 & Rxy(7,9)>=0.4 
         & Rxy(3,9)>=8.0 & Rxy(3,9)<=9.0)
EL10 <- (test_data$Detector1>=0.4 & test_data$Detector2>=0.4 & Rxy(3,4)>=0.4 
         & Rxy(5,6)>=0.4 & Rxy(7,8)>=0.4 & test_data$Detector9>=0.1 
         & Rxy(3,9)>=10.3 & Rxy(3,9)<=12.5)
EL11 <- (test_data$Detector1>=0.4 & test_data$Detector1<=1.0 & test_data$Detector2>=0.4 
         & Rxy(3,4)>=0.4 & Rxy(5,9)>=0.4 & Rxy(7,9)>=11)
EL12 <- (test_data$Detector1>=0.4 & test_data$Detector1<=1.0 & test_data$Detector2>=0.4
         & test_data$Detector2<=1.0 & Rxy(3,4)>=0.4 & Rxy(3,4)<=1.5 
         & Rxy(5,9) >= 0.4 & Rxy(7,9)>=15)

ELOGIC <- data.frame(EL1, EL2, EL3, EL4, EL5, EL6, EL7, EL8, EL9, EL10, EL11, EL12)
Etotal_true <- sum(ELOGIC[test_data$Detector1>0.4,])
Erow_predict <- factor(rowSums(ELOGIC[test_data$Detector1>0.4,])>0,
                       levels = c(TRUE, FALSE), 
                       labels = c("Electron","Rejected Electron"))
Enproblem_rows <- sum(rowSums(ELOGIC[test_data$Detector1>0.4,]) > 1)
Eproblem_rows <- which(rowSums(ELOGIC[test_data$Detector1>0.4,]) > 1)
ratio_Eproblem_rows <- Enproblem_rows/length(test_data$Detector1>0.4)
print(ratio_Eproblem_rows*100)
REPT_etab <- table(predict=factor(Erow_predict,
                                  levels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton"),
                                  labels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton")),
                   truth= factor(test_data[test_data$Detector1>0.4,10],
                                 levels = c("Electron", "Proton"),
                                 labels = c("Electron", "Proton")))
print(REPT_etab)

# Proton Logic Equations (with veto constraints)
PL1  <- (test_data$Detector1>8.2 & test_data$Detector1<16 & test_data$Detector2<6.5
         & Rxy(1,2)>8.2 & Rxy(1,2)<18 & Rxy(3,9)<0.5 & Rbarpxy(3,9))
PL2  <- (test_data$Detector1>5.4 & test_data$Detector1<12.2 & test_data$Detector2>5.0
         & test_data$Detector2<16.9 & Rxy(3,4)>0.1 & Rxy(3,4)<11 & Rxy(1,4)>15.9 
         & Rxy(1,4)<25.7 & Rxy(5,9)<0.5 & Rbarpxy(5,9))
PL3  <- (test_data$Detector1>4 & test_data$Detector1<7 & test_data$Detector2>4
         & test_data$Detector2<9.5 & Rxy(3,4)>10 & Rxy(5,6)<12.5
         & Rxy(1,6)>24 & Rxy(1,6)<35.5 & Rxy(7,9)<0.5 & Rbarpxy(7,9))
PL4  <- (test_data$Detector1>3.1 & test_data$Detector1<4.9 & test_data$Detector2>3.2
         & test_data$Detector2<5.7 & Rxy(3,4)>7.6 & Rxy(3,4)<16.8
         & Rxy(5,6)>9.2 & Rxy(5,6)<24 & Rxy(7,8)<23 
         & test_data$Detector9<4.1 & Rxy(5,9)>11.5 & Rxy(5,9)<33.0)
PL5  <- (test_data$Detector1>2.2 & test_data$Detector1<4 & test_data$Detector2>1.9
         & test_data$Detector2<4.2 & Rxy(3,4)>5.5 & Rxy(3,4)<12.5
         & Rxy(5,6)>5.8 & Rxy(5,6)<12.5 & Rxy(7,8)>7 & Rxy(7,8)<22.7 
         & test_data$Detector9>1 & test_data$Detector9<13 & Rxy(7,9)>5 & Rxy(7,9)<45)
PL6  <- (test_data$Detector1>1.5 & test_data$Detector1<3.3 & test_data$Detector2>1.0
         & test_data$Detector2<3.3 & Rxy(3,4)>4.1 & Rxy(3,4)<6.5 & Rxy(5,6)>4.5
         & Rxy(5,6)<7.2 & Rxy(7,8)>4.8 & Rxy(7,8)<8.0 & test_data$Detector9>2.0
         & test_data$Detector9<8.5 & Rxy(1,6)>11 & Rxy(1,6)<22 & Rxy(1,9)<65)
PL7  <- (test_data$Detector1>1.4 & test_data$Detector1<2.5 & test_data$Detector2>1.4
         & test_data$Detector2<2.8 & Rxy(3,4)>3.4 & Rxy(3,4)<5.4 & Rxy(5,6)>3.4 
         & Rxy(5,6)<5.9 & Rxy(7,8)>3.5 & Rxy(7,8)<6 & Rxy(1,9)>10 & Rxy(1,9)<45)
PL8  <- (test_data$Detector1>0.8 & test_data$Detector1<3.0 & test_data$Detector2>0.8
         & test_data$Detector2<3.0 & Rxy(3,4)>2.5 & Rxy(3,4)<5 & Rxy(5,6)>2.5 
         & Rxy(5,6)<5.5 & Rxy(7,8)>2.5 & Rxy(7,8)<5.5 & test_data$Detector9>1 
         & test_data$Detector9<6 & Rxy(1,9)<8 & Rxy(1,9)<32.0)

PLOGIC <- data.frame(PL1, PL2, PL3, PL4, PL5, PL6, PL7, PL8)
Ptotal_true <- sum(PLOGIC[test_data$Detector1>0.4,])
Prow_predict <- factor(rowSums(PLOGIC[test_data$Detector1>0.4,])>0,
                       levels = c(TRUE, FALSE), 
                       labels = c("Proton","Rejected Proton"))
Pnproblem_rows <- sum(rowSums(PLOGIC[test_data$Detector1>0.4,]) > 1)
Pproblem_rows <- which(rowSums(PLOGIC[test_data$Detector1>0.4,]) > 1)
ratio_Pproblem_rows <- Pnproblem_rows/length(test_data$Detector1)
print(ratio_Pproblem_rows*100)
REPT_ptab <- table(predict=factor(Prow_predict,
                                  levels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton"),
                                  labels = c("Electron","Rejected Electron",  "Proton", "Rejected Proton")),
                   truth= factor(test_data[test_data$Detector1>0.4,10],
                                 levels = c("Electron", "Proton"),
                                 labels = c("Electron", "Proton")))
print(REPT_ptab)

REPT_logic <- (sweep(REPT_etab, 2, colSums(REPT_etab), FUN = "/")*100 
              + sweep(REPT_ptab, 2, colSums(REPT_ptab), FUN = "/")*100)
print(REPT_logic)
print(ratio_Eproblem_rows*100)
print(ratio_Pproblem_rows*100)

### K-fold Cross-Validation ###
# create folds
set.seed(123)
k_folds = 10
folds = createFolds(training_data_sample$Particle_Type, k = k_folds)

# k-fold cross validation for linear kernel
lin_accuracy <- vector(length = k_folds)
cm_lin <- matrix(0, nrow = 2, ncol = 2)

for (i in 1:k_folds) {
  # Subset training and validation data based on indices
  train_fold <- training_data_sample[folds[[i]], ]
  test_fold  <- training_data_sample[-folds[[i]], ]
  
  # Train the SVM model on the training data (replace with your desired kernel)
  svm_model <- svm(Particle_Type ~ .,
                   data = train_fold, 
                   kernel = "linear", scale=FALSE, cost = 10)
  
  # Make predictions on the test data
  k_pred <- predict(svm_model, test_fold)
  
  # Calculate confusion matrix for the current fold
  fold_cm <- table(predict = k_pred, truth = test_fold$Particle_Type)
  lin_accuracy[i] = ((fold_cm[1,1] + fold_cm[2,2]) 
                     / (fold_cm[1,1] + fold_cm[2,2] + fold_cm[1,2] + fold_cm[2,1]))
  
  # Update the overall confusion matrix by adding elements
  cm_lin <- cm_lin + fold_cm
}
lin_ave_acc <- mean(lin_accuracy)
print(paste("Average Accuracy across Folds:", lin_ave_acc))
cm_lin


# k-fold cross validation for polynomial kernel
poly_accuracy <- vector(length = k_folds)
cm_poly <- matrix(0, nrow = 2, ncol = 2)

for (i in 1:k_folds) {
  # Subset training and validation data based on indices
  train_fold <- training_data_sample[folds[[i]], ]
  test_fold  <- training_data_sample[-folds[[i]], ]
  
  # Train the SVM model on the training data (replace with your desired kernel)
  svm_model <- svm(Particle_Type ~ .,
                   data = train_fold, 
                   kernel = "polynomial", scale=FALSE, cost = 10)
  
  # Make predictions on the test data
  k_pred <- predict(svm_model, test_fold)
  
  # Calculate confusion matrix for the current fold
  fold_cm <- table(predict = k_pred, truth = test_fold$Particle_Type)
  poly_accuracy[i] = ((fold_cm[1,1] + fold_cm[2,2]) 
                     / (fold_cm[1,1] + fold_cm[2,2] + fold_cm[1,2] + fold_cm[2,1]))
  
  # Update the overall confusion matrix by adding elements
  cm_poly <- cm_poly + fold_cm
}
poly_ave_acc <- mean(poly_accuracy)
print(paste("Average Accuracy across Folds:", poly_ave_acc))
cm_poly


# k-fold cross validation for radial kernel
rad_accuracy <- vector(length = k_folds)
cm_rad <- matrix(0, nrow = 2, ncol = 2)

for (i in 1:k_folds) {
  # Subset training and validation data based on indices
  train_fold <- training_data_sample[folds[[i]], ]
  test_fold  <- training_data_sample[-folds[[i]], ]
  
  # Train the SVM model on the training data (replace with your desired kernel)
  svm_model <- svm(Particle_Type ~ .,
                   data = train_fold, 
                   kernel = "radial", scale=FALSE, cost = 10)
  
  # Make predictions on the test data
  k_pred <- predict(svm_model, test_fold)
  
  # Calculate confusion matrix for the current fold
  fold_cm <- table(predict = k_pred, truth = test_fold$Particle_Type)
  rad_accuracy[i] = ((fold_cm[1,1] + fold_cm[2,2]) 
                 / (fold_cm[1,1] + fold_cm[2,2] + fold_cm[1,2] + fold_cm[2,1]))
  
  # Update the overall confusion matrix by adding elements
  cm_rad <- cm_rad + fold_cm
}
rad_ave_acc <- mean(rad_accuracy)
print(paste("Average Accuracy across Folds:", rad_ave_acc))
cm_rad


# k-fold cross validation for simplified linear kernel
linsi_accuracy <- vector(length = k_folds)
cm_linsi <- matrix(0, nrow = 2, ncol = 2)

for (i in 1:k_folds) {
  # Subset training and validation data based on indices
  train_fold <- training_data_sample[folds[[i]], ]
  test_fold  <- training_data_sample[-folds[[i]], ]
  
  # Train the SVM model on the training data (replace with your desired kernel)
  svm_model <- svm(Particle_Type ~ Detector1+Detector2+Detector4+Detector9,
                   data = train_fold, 
                   kernel = "linear", scale=FALSE, cost = 10)
  
  # Make predictions on the test data
  k_pred <- predict(svm_model, test_fold)
  
  # Calculate confusion matrix for the current fold
  fold_cm <- table(predict = k_pred, truth = test_fold$Particle_Type)
  linsi_accuracy[i] = ((fold_cm[1,1] + fold_cm[2,2]) 
                       / (fold_cm[1,1] + fold_cm[2,2] + fold_cm[1,2] + fold_cm[2,1]))
  
  # Update the overall confusion matrix by adding elements
  cm_linsi <- cm_linsi + fold_cm
}
linsi_ave_acc <- mean(linsi_accuracy)
print(paste("Average Accuracy across Folds:", linsi_ave_acc))
cm_linsi


### Plot Hyperplane ###
# Plot Edep v Einc for each detector, color coding particles
plot(training_data_sample$Detector2, training_data_sample$Detector4, 
     pch = 20,cex = 0.25, # Solid circle for all points
     xlim = c(0,20), ylim = c(0,20),
     xlab = 'Detector 2', ylab = 'Detector 4',
     col = ifelse(training_data_sample$Particle_Type == "Electron", "red","blue"))  # Color based on particle type

abline(a = -linear_hp_coefs[1]/linear_hp_coefs[4], 
       b = -linear_hp_coefs[2]/linear_hp_coefs[4], lwd = 2, lty = 2)

# Add a legend
legend("topright", legend = c("Linear Model", "SVM"), lty = c(2,3), pch = 20)


linearsi_hp_coefs <- c(linearsi_hp_coefs[1], linearsi_hp_coefs[2],linearsi_hp_coefs[3],0,linearsi_hp_coefs[4],0,0,0,0,linearsi_hp_coefs[5])

pairs(~ Detector1+Detector2+Detector3+Detector4+Detector5
      +Detector6+Detector7+Detector8+Detector9, data=training_data_sample2,
      panel=function(x,y){
      points(x,y,col = ifelse(training_data_sample2$Particle_Type == "Electron", "red","blue"),
                              xlim = c(0,20), ylim = c(0,20),
                              pch = 20, cex = 0.25 # Solid circle for all points
      )
      for (i in 1:9){
        if (sum(x-training_data_sample2[i])==0){
          for (j in 1:9){
            if (sum(y-training_data_sample2[j])==0){
              abline(a = -linear_hp_coefs[1]/linear_hp_coefs[j+1], 
                     b = -linear_hp_coefs[i+1]/linear_hp_coefs[j+1], lwd = 3, lty = 2)
            }
          }
        }
      }
})

density_fn <- function(data, mapping, ...) {
  # Access variable indices based on mapping
  i <- which(names(data) == quo_name(mapping[[1]]))  # Index of x variable
  j <- which(names(data) == quo_name(mapping[[2]]))  # Index of y variable
  
  ggplot(data = data, mapping = mapping) +
    stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, h = c(1,1) ) +
    scale_fill_gradientn(colours=rev(c(rainbow(100)[1:70],"#00abff"))) +
    geom_abline(intercept = -linear_hp_coefs[1]/linear_hp_coefs[j+1], 
                slope = -linear_hp_coefs[i+1]/linear_hp_coefs[j+1], lwd = 1, lty = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 20))
}

binning_fn <- function(data){
  #finding bin edges
  edges = seq(from=0,to=max(apply(data[,1:9], 2, max)),length.out = 101)
  
  # Binning electron and proton counts
  ecounts <- lapply(c(1:9), function(x) {
    ebins <- cut(data[data$Particle_Type == "Electron",x], 
                 breaks = edges)
    data.frame(table(ebins))
  })
    
  ebins_midpoints <- matrix(0, nrow = 100, ncol = length(ecounts))
  ecounts_new <- matrix(0, nrow = 100, ncol = length(ecounts))
  max_ecounts <- matrix(0, nrow = 1, ncol = length(ecounts))
  colnames(ecounts_new) <- column_names
  colnames(max_ecounts) <- column_names
  for (i in 1:length(ecounts)){
    ecounts_new[,i] <- ecounts[[i]][,2]
    max_ecounts[i] <- max(ecounts_new[,i])
  }
  
  pcounts <- lapply(c(1:9), function(x) {
    pbins <- cut(data[data$Particle_Type == "Proton",x], 
                 breaks = edges)
    data.frame(table(pbins))
  })
  
  pcounts_new <- matrix(0, nrow = 100, ncol = length(pcounts))
  max_pcounts <- matrix(0, nrow = 1, ncol = length(pcounts))
  colnames(pcounts_new) <- column_names
  colnames(max_pcounts) <- column_names
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
  
  # Create color mapping
  returned_data <- binning_fn(data)
  ep_ratio = returned_data[[4]]
  d1_emax <- max(returned_data[[2]][,1])
  di_emax <- max(returned_data[[2]][,map_i])
  d_eratio <- di_emax/d1_emax
  
  e_colors <- c(rainbow(10000)[(1+floor((1-d_eratio)*10000/ep_ratio)):floor(10000/ep_ratio)],NA)

  ggplot(data = data[data$Particle_Type == "Electron",], mapping = mapping) +
    stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, h = c(1, 1)) +
    scale_fill_gradientn(colours = rev(e_colors), name = "Electron Density") +
    geom_abline(intercept = -linear_hp_coefs[1] / linear_hp_coefs[map_j + 1],
                slope = -linear_hp_coefs[map_i + 1] / linear_hp_coefs[map_j + 1], lwd = 2, lty = 1) +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 20)) +
    theme_few()
}

density_fn_proton <- function(data, mapping, pt, ...) {
  # Access variable indices based on mapping
  map_i <- which(names(data) == quo_name(mapping[[1]]))  # Index of x variable
  map_j <- which(names(data) == quo_name(mapping[[2]]))  # Index of y variable
  
  # Create color mapping
  returned_data <- binning_fn(data)
  ep_ratio = returned_data[[4]]
  d1_pmax <- max(returned_data[[3]][,1])
  di_pmax <- max(returned_data[[3]][,map_j])
  d_pratio <- di_pmax/d1_pmax
  
  p_colors <- c(rainbow(10000)[(4001+ceiling(10000/ep_ratio+(1-d_pratio)*ep_ratio)):
                                (4001+(ceiling(10000/ep_ratio+ep_ratio)))],NA)
  
  ggplot(data = data[data$Particle_Type == "Proton",], mapping = mapping) +
    stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE, h = c(1, 1)) +
    scale_fill_gradientn(colours = rev(p_colors), name = "Proton Density") +
    geom_abline(intercept = -linear_hp_coefs[1] / linear_hp_coefs[map_j + 1],
                slope = -linear_hp_coefs[map_i + 1] / linear_hp_coefs[map_j + 1], lwd = 2, lty = 1) +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 20)) +
    theme_few()
}

diag_fn <- function(data, mapping, pt, ...) {
  # Access variable indices based on mapping
  mapping_index <- which(names(data) == quo_name(mapping[[1]]))  # Index of x variable
  
  #use binning function to bin data
  returned_data <- binning_fn(data)
  edges       = returned_data[[1]]
  midpoints = matrix(0, nrow = length(edges)-1, ncol = 1)
  for (i in 1:length(edges)-1){
    midpoints[i] = (edges[i+1]+edges[i])/2
  }
  ecounts_new = returned_data[[2]]
  pcounts_new = returned_data[[3]]
  
  plot_data <- data.frame(midpoints, ecounts_new[,mapping_index], pcounts_new[,mapping_index])
  colnames(plot_data) <- c("midpoints","ecounts", "pcounts")
  
  ggplot(data = plot_data, aes(midpoints)) +
    geom_line(aes(y = ecounts/sum(ecounts), colour = "blue"), lwd = 2) + 
    geom_line(aes(y = pcounts/sum(pcounts), colour = "red"), lwd = 2) +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, 0.4)) +
    theme_few()
}

ggpairs(training_data_sample2[training_data_sample2$Detector1 > 0.1,],
        #columns = c(1:9),
        columns = c(1, 2, 4, 9),
        # Mapping function with additional argument for particle type
        #mapping = aes(color = Particle_Type),
        upper = list(continuous = density_fn_proton),
        diag = list(continuous = diag_fn),
        lower = list(continuous = density_fn_electron)) +
  theme(axis.text = element_text(size = 28), 
        strip.text = element_text(size = 28, color = "black"),
        panel.spacing=grid::unit(1.2,"lines"))

# Creating color scales
returned_data <- binning_fn(training_data_sample2[training_data_sample2$Detector1 > 0.1,])
ep_ratio = returned_data[[4]]
e_colors_length <- floor(10000/ep_ratio)+1
p_colors_start <- 4001+ceiling(10000/ep_ratio)
p_colors_end <- 4001+ceiling(10000/ep_ratio)+e_colors_length

e_colors <- c(rainbow(10000)[1:floor(10000/ep_ratio)],NA)

p_colors <- c(NA,rainbow(10000)[(4001+ceiling(10000/ep_ratio)):
                              (4001+ceiling(10000/ep_ratio)+e_colors_length)])

v <- ggplot(faithful[1:100,], aes(waiting, eruptions, fill = (returned_data[[2]][,1]/sum(ecounts_new[,1])))) +
  geom_tile()
v + scale_fill_gradientn(colours = rev(e_colors), name = "Electron Density   ") + 
  theme(legend.key.width  = unit(5, "lines"), legend.position = "bottom",
        legend.text = element_text(size = 26), legend.title = element_text(size = 26))

v <- ggplot(faithful[1:100,], aes(waiting, eruptions, fill = (returned_data[[2]][,1]/sum(ecounts_new[,1])))) +
  geom_tile()
v + scale_fill_gradientn(colours = p_colors, name = "Proton Density   ") + 
  theme(legend.key.width  = unit(5, "lines"), legend.position = "bottom",
        legend.text = element_text(size = 26), legend.title = element_text(size = 26))
