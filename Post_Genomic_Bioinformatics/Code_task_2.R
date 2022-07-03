
# Clear Environment
rm(list = ls())

# Load Classyfire and R.Matlab
library('classyfire')
library('R.matlab')

# Read in the Faecal control vs CD samples
bgw <- readMat("BWG_FA_CDvCTRL.mat")

# Extract the data we need from bgw

# Data matrix with the input from GCMS
X <- bgw$XTIC
# The output vector containing the disease state
y <- bgw$CLASS

# copy sample names from SAM to row names of X
rownames(X) <- as.character(unlist(bgw$SAM))

########################################

## Task 1 ## Faecal ##

# Get a feel for the data set
summary(X)
dim(X)

# Get the retention times (RT)
RT <- bgw$RT
i = 1 # set index number to sample you want to plot
plot(RT, X[i,], type="l", xlab="RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(X)[i]))


## Task 2 ##

# Data Pre-Proccessing

# Perform PCA on the data
Xpca <- prcomp(X, scale= TRUE)

# Look at the PCA results
summary(Xpca)

ssss <- summary(Xpca)
TablePCA <- data.frame(t(ssss$importance))
TablePCA <- TablePCA$Proportion.of.Variance
TablePCA$Standard.deviation <- NULL
TablePCA$Cumulative.Proportion <- NULL


barplot(TablePCA$Proportion.of.Variance,
        ylab = 'Proportion of variance',
        xlab = 'Principle components',
        names.arg = rownames(TablePCA))


# Plot scatter plot of PCA1 vs PCA2
xscores <- Xpca$x
my_cols <- c('blue', 'green')
plot(xscores[,1], xscores[,2], pch=4, col=my_cols[y])

# Retain first 4 PCAs
PCA1_4 <- xscores[,1:4]

## Task 3 ##



# Building and testing machine learning models

# Create a model - bootstrap 100 - his paper
ml <- cfBuild(PCA1_4, inputClass = y, cpus = 4)

# Check to see how many ensaembles should be used ~ 49
ggEnsTrend(ml)

#Redo with better number of ensembles
mlens <- cfBuild(PCA1_4, inputClass = y, cpus = 4, bootNum = 100, ensNum = 49)

# Check performace of model
getConfMatr(mlens)

# Test accuracy 
mean(mlens$testAcc) # 79.85% accurate
mean(mlens$trainAcc)



# Determine is any samples were particularly difficult
zzz <- mlens$missNames
new.pp <- unlist(zzz,recursive=FALSE)
xxx <- as.data.frame(table(new.pp))
rownames(xxx) <- xxx$new.pp
xxx$new.pp <- NULL
xxx <- t(xxx)

# Increase margin size
par(mar=c(11,4,4,4))
barplot(xxx, las=2, ylab = 'Frequency') # Plot miss-identified samples

par(mfrow = c(1,2))

ggEnsHist(mlens, density = FALSE, percentiles=TRUE, median=TRUE, mean = TRUE)

ggPermHist(permute, density=FALSE, percentiles = TRUE, mean = TRUE, median = TRUE)

# Which condition was it better at identifying
ggClassPred(mlens, position = "stack", displayAll = TRUE, showText = TRUE)

# Perumtation testing
permute <- cfPermute(PCA1_4,
                     y,
                     bootNum = 100,
                     ensNum = 49,
                     permNum = 100,
                     parallel = TRUE, 
                     cpus = 4)



# Histogram of perumation accuracies
ggPermHist(permute, position = 'stack', displayALL = TRUE, showtext = TRUE)

# Histogram of test accuracies
hist(testacc)
# Better one
Accuracies <- testacc
ggplot() + aes(Accuracies)+ geom_histogram(binwidth=4, colour="black", fill="white")

# Box blot comparing permutations to tests
par(mfrow = c(1,1))
testacc <- mlens$testAcc
permacc <- permute$avgAcc

max_length <- max(c(length(testacc), length(permacc)))    # Find out maximum length
max_length                


modelstat <- data.frame(Ensemble = c(testacc,             
                            rep(NA, max_length - length(testacc))),
                   Permutation = c(permacc,
                            rep(NA, max_length - length(permacc))))
                                           
mmdf <- data.frame(Accuracy = c(modelstat[,2], modelstat[1:49, 1]),
                 Type = rep(c("Permutation","Ensemble"), times = c(100,49)))

# Basic boxplot 
boxplot(Accuracy~Type,
        data=mmdf,
        col = c('cyan3','firebrick4'),
        notch = TRUE,
        staplewex = TRUE,
        varwidth = TRUE,
        horizontal = TRUE)

# Comparison of permutation and test means

mean(testacc)
mean(permacc)

#######################################################
########################################################

### Blood samples ###


# Read in the Faecal control vs CD samples
bloodgw <- readMat("BWG_BL_CDvCTRL.mat")

# Extract the data we need from bgw

# Data matrix with the input from GCMS
bX <- bloodgw$XTIC
# The output vector containing the disease state
by <- bloodgw$CLASS

# copy sample names from SAM to row names of bX
rownames(bX) <- as.character(unlist(bloodgw$SAM))

########################################

## Task 1 ##

# Get a feel for the data set
summary(bX)
dim(bX)

# Get the retention times (RT)
RT <- bloodgw$RT
i = 1 # set index number to sample you want to plot
plot(RT, bX[i,], type="l", xlab="RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(bX)[i]))


## Task 2 ##

# Data Pre-Proccessing

# Perform PCA on the data
bXpca <- prcomp(bX, scale= TRUE)

# Look at the PCA results
summary(bXpca)
plot(bXpca)

# Plot scatter plot of PCA1 vs PCA2
bxscores <- bXpca$x
my_cols <- c('blue', 'green')
plot(bxscores[,1], bxscores[,2], pch=4, col=my_cols[y])

# Retain first 4 PCAs
bPCA1_3 <- bxscores[,1:3]

## Task 3 ##

# Building and testing machine learning models

# Create a model - bootstrap 100 - his paper
bml <- cfBuild(bPCA1_3, inputClass = by, cpus = 4)

# Check to see how many ensembles should be used ~ 78
ggEnsTrend(bml)

#Redo with better number of ensembles
bmlens <- cfBuild(bPCA1_3, inputClass = by, cpus = 4, bootNum = 100, ensNum = 78)

# Check performace of model
getConfMatr(bmlens)

# Test accuracy 
mean(bmlens$testAcc) # 53.61% accurate
mean(bmlens$trainAcc)


#######################################################
########################################################

### Breath samples ###


# Read in the Faecal control vs CD samples
breathgw <- readMat("BWG_BR_CDvCTRL.mat")

# Extract the data we need from bgw

# Data matrix with the input from GCMS
brX <- breathgw$XTIC
# The output vector containing the disease state
bry <- breathgw$CLASS

# copy sample names from SAM to row names of bX
rownames(brX) <- as.character(unlist(breathgw$SAM))

########################################

## Task 1 ##

# Get a feel for the data set
summary(brX)
dim(brX)

# Get the retention times (RT)
RT <- breathgw$RT
i = 1 # set index number to sample you want to plot
plot(RT, brX[i,], type="l", xlab="RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(brX)[i]))


## Task 2 ##

# Data Pre-Proccessing

# Perform PCA on the data
brXpca <- prcomp(brX, scale= TRUE)

# Look at the PCA results
summary(brXpca)
plot(brXpca)

# Plot scatter plot of PCA1 vs PCA2
brxscores <- brXpca$x
my_cols <- c('blue', 'green')
plot(brxscores[,1], brxscores[,2], pch=4, col=my_cols[y])

# Retain first 4 PCAs
brPCA1_2 <- brxscores[,1:2]

## Task 3 ##

# Building and testing machine learning models

# Create a model - bootstrap 100 - his paper
brml <- cfBuild(brPCA1_2, inputClass = bry, cpus = 4)

# Check to see how many ensembles should be used ~ 70
ggEnsTrend(brml)

#Redo with better number of ensembles
brmlens <- cfBuild(brPCA1_2, inputClass = bry, cpus = 4, bootNum = 100, ensNum = 70)

# Check performace of model
getConfMatr(brmlens)

# Test accuracy 
mean(brmlens$testAcc) # 58.27% accurate
mean(brmlens$trainAcc)

#######################################################
########################################################

### Urine samples ###


# Read in the Faecal control vs CD samples
Urinegw <- readMat("BWG_UR_CDvCTRL.mat")

# Extract the data we need from bgw

# Data matrix with the input from GCMS
uX <- Urinegw$XTIC
# The output vector containing the disease state
uy <- Urinegw$CLASS

# copy sample names from SAM to row names of bX
rownames(uX) <- as.character(unlist(Urinegw$SAM))

########################################

## Task 1 ##

# Get a feel for the data set
summary(uX)
dim(uX)

# Get the retention times (RT)
RT <- Urinegw$RT
i = 1 # set index number to sample you want to plot
plot(RT, uX[i,], type="l", xlab="RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(uX)[i]))


## Task 2 ##

# Data Pre-Proccessing

# Perform PCA on the data
uXpca <- prcomp(uX, scale= TRUE)

# Look at the PCA results
summary(uXpca)
plot(uXpca)

# Plot scatter plot of PCA1 vs PCA2
uxscores <- uXpca$x
my_cols <- c('blue', 'green')
plot(uxscores[,1], uxscores[,2], pch=4, col=my_cols[uy])

# Retain first 3 PCAs
uPCA1_3 <- uxscores[,1:3]

## Task 3 ##

# Building and testing machine learning models

# Create a model - bootstrap 100 - his paper
uml <- cfBuild(uPCA1_3, inputClass = uy, cpus = 4)

# Check to see how many ensembles should be used ~ 75
ggEnsTrend(uml)

#Redo with better number of ensembles
umlens <- cfBuild(uPCA1_3, inputClass = uy, cpus = 4, bootNum = 100, ensNum = 75)

# Check performace of model
getConfMatr(umlens)

# Test accuracy 
mean(umlens$testAcc) # 61% accurate
mean(umlens$trainAcc)


#######################################################
########################################################

## Random Forest ##

# Install caret
install.packages("caret", dependencies = TRUE)
install.packages("randomForest")

# Load packages
library(caret)
library(randomForest)


# Variables that we are using - first 4 PCs
forest_PCA1_4 <- as.data.frame(xscores[,1:4])

# Class of each sample
forest_y <- y

# Add class column to variables data frame 
forest_PCA1_4$Diagnosis <- forest_y


# Convert diagnosis to a factor
forest_PCA1_4$Diagnosis <- as.factor(forest_PCA1_4$Diagnosis)

#Set a random seed
set.seed(998)

# Model 
forest_model <- train(Diagnosis~PC1+PC2+PC3+PC4, 
                      data = forest_PCA1_4,
                      method = 'rf',
                      trControl=trainControl(method = 'cv'),
                      number = 10,
                      repeats = 10)

plot(forest_model)

forest_model

citation('caret')

par(mfrow = c(2,2))

ggEnsHist(mlens, density = TRUE, percentiles=TRUE, mean=TRUE)
ggEnsHist(bmlens, density = TRUE, percentiles=TRUE, mean=TRUE)
ggEnsHist(brmlens, density = TRUE, percentiles=TRUE, mean=TRUE)
ggEnsHist(umlens, density = TRUE, percentiles=TRUE, mean=TRUE)


