# Load the data into a data frame
data1 <- read.table('Experiment1_data_set 11 .txt')

# Check the data frame and variables
str(data1)
summary(data1)

# Change Treatment category to a factor
data1$Probiotic <- as.factor(data1$Probiotic)

# Exploratory data analysis
par(mfrow=c(1,3))

plot(data1$Days, data1$AMP)
# Appears to show a negative correlation between severity and AMP

plot(data1$Probiotic, data1$AMP)
# Patients taking the probiotic appear to have higher AMP levels

plot(data1$Days, data1$AMP, col = data1$Probiotic)
# Appears to be different slope between the variables suggesting an interaction

# Fitting a model - ANCOVA - with interaction
model <- lm(data1$AMP ~ data1$Days*data1$Probiotic, data = data1)

# Model diagnostic plots
par(mfrow=c(2,2))
plot(model)

#Normal distribution of errors#
#Normal Q-Q plot - closely follows line of unity
# definitely no systematic deviation from line of unity

#Homogeneous variance#
#No real heteroscedacity in residuals vs fitted

#Linearity#
#no real curve in residuals vs fitted 

# Looking at the summary table of the model
summary(model)
#Slope placebo = -1.2673, intercept = 70.4220 
#Probiotic slope = 0.11, intercept = 68.51 

# Test interaction for significance with a deletion test
drop1(model, test='F')
# Significant interaction p = 3.215e-06, F-value = 25.28

#Plot data with a fitted model lines added for each group
par(mfrow=c(1,1))
plot(data1$Days, data1$AMP, col = data1$Probiotic, pch=18, ylim=c(20,120), xlab='Number of days with symptoms pre-treatment', ylab = 'Antimicrobial peptides (mg/ml)')
legend(15,120,legend=c('Placebo', 'Probiotic'), col=1:length(data1$Probiotic),pch=18, bty='n')
abline(model)
abline(a=model$coefficients[1]+model$coefficients[3],b=model$coefficients[2]+model$coefficients[4], col='red' )


## Data set 2 ##

# Load data
data2 <- read.table('Experiment2_data_set 11 .txt')

# Check the data frame
str(data2)
summary(data2)

# Change variables to factors
data2$Diet <- as.factor(data2$Diet)
data2$Probiotic <- as.factor(data2$Probiotic)

#Explore the data
par(mfrow=c(1,1))
boxplot(data2$Shannon.diversity ~ data2$Diet + data2$Probiotic)

# Fit a model - Two-way ANOVA
# Set reference variables
data2$Probiotic <- relevel(data2$Probiotic, ref = 'No')
data2$Diet <- relevel(data2$Diet, ref = 'Normal')
model2 <- lm(data2$Shannon.diversity ~ data2$Diet * data2$Probiotic, data = data2)

# Diagnostic plots
par(mfrow=c(2,2))
plot(model2)

#Normal distribution of errors#
#Normal Q-Q plot - closely follows line of unity
# definitely no systematic deviation from line of unity

#Homogeneous variance#
#no heteroscedacity in residuals vs fitted

#Linearity#
#no curve at all in residuals vs fitted 

#Interpret results
# Summary Table
summary(model2)
#Anova table
anova(model2)
#Signififcant interaction

#Plot data 
par(mfrow=c(1,1))
boxplot(data2$Shannon.diversity ~ data2$Diet + data2$Probiotic,  
        data = data2, xlab='Treatment Type', ylab='Shannon Diversity', 
        names= c('Placebo', 'Placebo', 'Probiotic', 'Probiotic'), 
        col= c('green','orange'))
legend('bottomright',legend = c('High Fibre','Normal'), pch = 15, col = c('green','orange'), title = 'Diet Type', bty = 'y')

# Make a better plot using ggplot
# Install ggpubr so that the ggboxplot and ggplot functions can be called
install.packages('ggpubr')
library('ggpubr')
#Plot a boxplot
ggboxplot(data2, x = 'Probiotic', y = 'Shannon.diversity', color = 'Diet',
          palette = c("#E69F00", "#009E73"), ylab = 'Shannon Diversity') + 
          theme_gray()


#Interaction diagram with ggpplot

#Install tidyverse so that the summarise function can be used 
#Used to produce a table which will be used for the interaction plot
install.packages("tidyverse")
library(tidyverse)
#Create a mean table
data2mean <- data2 %>% group_by(Probiotic,Diet) %>% summarise(Shannon_Diversity = mean(Shannon.diversity))

#Plot interaction diagram
ggplot(data2mean, 
       aes(x = Probiotic, y = Shannon_Diversity, colour = Diet, group = Diet) ) +
  geom_point(size = 3) + geom_line() + scale_color_manual(values = c("#E69F00", "#009E73"))

#Citations

citation("tidyverse")

