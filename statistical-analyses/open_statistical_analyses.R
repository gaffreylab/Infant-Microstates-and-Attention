################################################################################
##################################### Hello ####################################
################################################################################

# Author: Armen Bagdasarov
# Institution: Duke University
# Purpose: Primary stat analyses for attention & EEG microstates manuscript

################################################################################

# Load packages
library(lm.beta)
library(ggplot2)
library(jtools)
library(MASS)
library(olsrr)
library(stats)
library(car)
library(splithalfr)

################################################################################

# Set working directory
setwd("XXX")

################################################################################

# Correlations to determine covariate for regression models

# Read data
data <- read.csv("open_complete_data_frame.csv")

# Run correlations
cor.test(data$m4_duration, data$age_months)
cor.test(data$m4_duration, data$eeg_time_recoded)
cor.test(data$m4_duration, data$sex_recoded)
cor.test(data$m4_duration, data$et_time_seconds)

################################################################################

# Regression model 1

# Dependent Variable: Microstate 4 duration
# Independent Variable: Child-led joint attention rate
# Covariates: Age & time of EEG data collection

# Read data
data <- read.csv("open_complete_data_frame.csv")

# Create subset of data
data_subset <- subset(data, select = c(m4_duration, 
                                       age_months,
                                       eeg_time_recoded,
                                       et_ja_rate_child_led))

# Identify multivariate outliers
# Outlier detection via minimum covariance determinant (MCD)
output75 <- cov.mcd(data_subset, quantile.used = nrow(data_subset)*.75)
mhmcd75 <- mahalanobis(data_subset, output75$center, output75$cov)
alpha <- .001
cutoff <- (qchisq(p = 1 - alpha, df = ncol(data_subset)))
names_outlier_MCD75 <- which(mhmcd75 > cutoff)
  # No outliers

# Step 1: Model 1a (n = 43)
model_1a <- lm(m4_duration ~ age_months + eeg_time_recoded, 
               data = data_subset)
summary(model_1a)
lm.beta(model_1a)
confint(model_1a, level = 0.95)

# Compare to intercept-only model
model_0 <- lm(m4_duration ~ 1, data = data_subset)
anova(model_0, model_1a)

# Step 2: Model 1b (n = 43)
model_1b <- lm(m4_duration ~ age_months + eeg_time_recoded + 
                 et_ja_rate_child_led, data = data_subset)
summary(model_1b)
lm.beta(model_1b)
confint(model_1b, level = 0.95)

# Compare models 1a and 1b
anova(model_1a, model_1b)

################################################################################

# Regression model 2

# Dependent Variable: Microstate 4 duration
# Independent Variable: Caregiver-led joint attention rate
# Covariates: Age & time of EEG data collection

# Read data
data <- read.csv("open_complete_data_frame.csv")

# Create subset of data
data_subset <- subset(data, select = c(m4_duration, 
                                       age_months,
                                       eeg_time_recoded,
                                       et_ja_rate_caregiver_led))

# Identify multivariate outliers
# Outlier detection via minimum covariance determinant (MCD)
output75 <- cov.mcd(data_subset, quantile.used = nrow(data_subset)*.75)
mhmcd75 <- mahalanobis(data_subset, output75$center, output75$cov)
alpha <- .001
cutoff <- (qchisq(p = 1 - alpha, df = ncol(data_subset)))
names_outlier_MCD75 <- which(mhmcd75 > cutoff)

# Remove MCD outliers as identified (n = 4)
data_reduced <- data_subset[-c(names_outlier_MCD75), ]
  # Use this data frame moving forward

# Step 1: Model 2a (n = 39)
model_2a <- lm(m4_duration ~ age_months + eeg_time_recoded, 
               data = data_reduced)
summary(model_2a)
lm.beta(model_2a)
confint(model_2a, level = 0.95)

# Compare to intercept-only model
model_0 <- lm(m4_duration ~ 1, data = data_reduced)
anova(model_0, model_2a)

# Step 2: Model 2b (n = 39)
model_2b <- lm(m4_duration ~ age_months + eeg_time_recoded + 
                 et_ja_rate_caregiver_led, data = data_reduced)
summary(model_2b)
lm.beta(model_2b)
confint(model_2b, level = 0.95)

# Compare models 2a and 2b
anova(model_2a, model_2b)

################################################################################

# Regression model 3

# Dependent Variable: Microstate 4 duration
# Independent Variable: Infant attention shifts rate
# Covariates: Age & time of EEG data collection

# Read data
data <- read.csv("open_complete_data_frame.csv")

# Create subset of data
data_subset <- subset(data, select = c(m4_duration, 
                                       age_months,
                                       eeg_time_recoded,
                                       et_shifts_rate))

# Identify multivariate outliers
# Outlier detection via minimum covariance determinant (MCD)
output75 <- cov.mcd(data_subset, quantile.used = nrow(data_subset)*.75)
mhmcd75 <- mahalanobis(data_subset, output75$center, output75$cov)
alpha <- .001
cutoff <- (qchisq(p = 1 - alpha, df = ncol(data_subset)))
names_outlier_MCD75 <- which(mhmcd75 > cutoff)

# Remove MCD outliers as identified (n = 2)
data_reduced <- data_subset[-c(names_outlier_MCD75), ]
  # Use this data frame moving forward

# Step 1: Model 3a (n = 41)
model_3a <- lm(m4_duration ~ age_months + eeg_time_recoded, 
               data = data_reduced)
summary(model_3a)
lm.beta(model_3a)
confint(model_3a, level = 0.95)

# Compare to intercept-only model
model_0 <- lm(m4_duration ~ 1, data = data_reduced)
anova(model_0, model_3a)

# Step 2: Model 3b (n = 41)
model_3b <- lm(m4_duration ~ age_months + eeg_time_recoded + et_shifts_rate, 
               data = data_reduced)
summary(model_3b)
lm.beta(model_3b)
confint(model_3b, level = 0.95)

# Compare models 3a and 3b
anova(model_3a, model_3b)

################################################################################

# Regression Model 4

# Dependent Variable: Microstate 4 duration
# Independent Variable: Infant sustained attention duration
# Covariates: Age & time of EEG data collection

# Read data
data <- read.csv("open_complete_data_frame.csv")

# Create subset of data
data_subset <- subset(data, select = c(m4_duration, 
                                       age_months,
                                       eeg_time_recoded,
                                       et_sa_duration))

# Identify multivariate outliers
# Outlier detection via minimum covariance determinant (MCD)
output75 <- cov.mcd(data_subset, quantile.used = nrow(data_subset)*.75)
mhmcd75 <- mahalanobis(data_subset, output75$center, output75$cov)
alpha <- .001
cutoff <- (qchisq(p = 1 - alpha, df = ncol(data_subset)))
names_outlier_MCD75 <- which(mhmcd75 > cutoff)

# Remove MCD outliers as identified (n = 2)
data_reduced <- data_subset[-c(names_outlier_MCD75), ]
  # Use this data frame moving forward

# Step 1: Model 4a (n = 39)
model_4a <- lm(m4_duration ~ age_months + eeg_time_recoded, 
               data = data_reduced)
summary(model_4a)
lm.beta(model_4a)
confint(model_4a, level = 0.95)

# Compare to intercept-only model
model_0 <- lm(m4_duration ~ 1, data = data_reduced)
anova(model_0, model_4a)

# Step 2: Model 4b (n = 39)
model_4b <- lm(m4_duration ~ age_months + eeg_time_recoded + et_sa_duration, 
               data = data_reduced)
summary(model_4b)
lm.beta(model_4b)
confint(model_4b, level = 0.95)

# Compare models 4a and 4b
anova(model_4a, model_4b)

################################################################################

# Comparison for multiple corrections 

# Original p values
p = c(0.001581, # Model 1
      0.0295, # Model 2
      0.0003457, # Model 3
      0.0007458) # Model 4

# Perform correction
p.adjust(p, method = "bonferroni")

################################################################################

# Models without outliers removed can be run by setting data = data instead of 
# data = data_reduced

################################################################################

# Infant vs. caregiver-led joint attention t-test

# Read data
data <- read.csv("open_complete_data_frame.csv")

# Create subset of data
data_subset <- subset(data, select = c(et_ja_rate_child_led, 
                                       et_ja_rate_caregiver_led))

# Run paired-samples t-test
t.test(data_subset$et_ja_rate_child_led, data_subset$et_ja_rate_caregiver_led, 
       paired = TRUE)

################################################################################

# Split-half reliability of microstate 4 duration

# Read data
data <- read.csv("open_complete_data_frame.csv")

# Get Spearman-Brown split-half reliability coefficient 
spearman_brown(data$m4_duration_odd_split_half, data$m4_duration_even_split_half)

################################################################################
##################################### Done #####################################
################################################################################