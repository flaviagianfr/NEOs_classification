# Retrieval of the libraries necessary for the correct import of the data, manipulation and data visualization

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(tidymodels)
library(leaps)
library(glmnet)
library(pROC)
library(rsample)
library(correlation)
library(DataExplorer)
library(knitr)
library(corrplot)
library(regclass)
library(rsample)
library(corrplot)
library(outliers)
library(dplyr)
library(caret)
library(class)

# Data Uploading and Pre-processing

nasa_orig <- read_delim("nasa.csv", delim = ",")

head(nasa_orig,5)

names(nasa_orig)<- gsub("\\s","_",names(nasa_orig))

anyNA(nasa_orig) 

toremove<- c("Neo_Reference_ID","Name", "Est_Dia_in_M(min)", 
             "Est_Dia_in_M(max)", "Est_Dia_in_Miles(min)",
             "Est_Dia_in_Miles(max)", "Est_Dia_in_Feet(min)",
             "Est_Dia_in_Feet(max)", "Close_Approach_Date",
             "Epoch_Date_Close_Approach", 
             "Relative_Velocity_km_per_sec",
             "Miles_per_hour", "Miss_Dist.(lunar)",
             "Miss_Dist.(kilometers)",
             "Miss_Dist.(miles)", "Orbiting_Body", "Orbit_ID", 
             "Orbit_Determination_Date", "Epoch_Osculation",
             "Equinox", "Est_Dia_in_KM(min)")  

nasa<- nasa_orig %>%
  select(-all_of(toremove))

colnames(nasa)[2]<- "Est_Dia_in_KM_max"
colnames(nasa)[4]<- "Miss_Dist_Astronomical"

summary(nasa)

nasa$Orbit_Uncertainity<- nasa$Orbit_Uncertainity %>%
  as.factor() %>%
  as.numeric()

nasa$Perihelion_Distance<- nasa$Perihelion_Distance %>%
  as.factor() %>%
  as.numeric()

nasa <- nasa %>% mutate(Hazardous = ifelse(Hazardous == TRUE, 1, 0))

nasa$Hazardous<- nasa$Hazardous %>%
  as.factor() 

prop.table(table(nasa$Hazardous))

# Train and test split

set.seed(0607)

split <- initial_split(nasa, prop = 0.75) 
train <- training(split)
test <- testing(split)

prop.table(table(train$Hazardous))
prop.table(table(test$Hazardous))

# Data analysis with some correlations plot

Haz<- train$Hazardous

# Absolute Magnitude 
a<- ggplot(train, aes(x = Absolute_Magnitude, fill = Haz)) + 
  geom_density(alpha = 0.4)+
  ggtitle("Absolute Magnitude - Density Plot") +
  xlab("Absolute Magnitude")

# Orbit Uncertainity 
b<- ggplot(train, aes(x = Orbit_Uncertainity, fill = Haz)) +
  geom_bar(aes(y = ..prop..),position = "dodge") +
  ggtitle("Orbit Uncertainity - Bar Plot") +
  xlab("Orbit Uncertainity")

# Minimum Orbit Intersection
c<- ggplot(train, aes(x = Minimum_Orbit_Intersection, fill = Haz)) +
  geom_density(alpha = 0.4) +
  xlim(0, 0.5) +
  ggtitle("Minimum Orbit Intersection - Density Plot") +
  xlab("Minimum Orbit Intersection")

# Eccentricity
d<- ggplot(train, aes(x = Eccentricity, fill = Haz)) +
  geom_density(alpha = 0.4)+
  ggtitle("Eccentricity - Density Plot") +
  xlab("Eccentricity")

# Perihelion Distance
e<- ggplot(train, aes(x = Perihelion_Distance, fill = Haz)) +
  geom_density(alpha = 0.4) +
  ggtitle("Perihelion Distance - Density Plot") +
  xlab("Perihelion Distance")

# Relative Velocity
f<- ggplot(train, aes(x = Relative_Velocity_km_per_hr, fill = Haz)) +
  geom_density(alpha = 0.4) +
  ggtitle("Relative Velocity - Density Plot") +
  xlab("Relative Velocity")

grid.arrange(a, b, c, d, e, f, nrow = 2)

corrplot(cor(train[,-19]),
         method = "number",
         diag = FALSE,
         tl.cex = 0.4,
         number.cex = 0.5,
         tl.col = "black")

# Partial correlations

correlation(train, partial = TRUE)

a1<- ggplot(data = train, aes(Jupiter_Tisserand_Invariant, Semi_Major_Axis)) + 
  geom_jitter(color = "dark green") +
  xlab("Jupiter Tisserand Invariant") +
  ylab("Semi-Major Axis")

a2<- ggplot(data = train, aes(Mean_Motion, Orbital_Period)) +
  geom_jitter(color = "dark green") +
  xlab("Mean Motion") +
  ylab("Orbital Period")

a3<- ggplot(data = train, aes(Orbital_Period, Semi_Major_Axis)) +
  geom_jitter(color = "dark green")  +
  xlab("Orbital Period") + 
  ylab("Semi-Major Axis")

a4<- ggplot(data = train, aes(Absolute_Magnitude, Est_Dia_in_KM_max))+
  geom_jitter(color = "dark green") +
  xlab("Absolute Magnitude") +
  ylab("Estimated Diameter")

a5<- ggplot(data = train, aes(as.factor(Orbit_Uncertainity), Absolute_Magnitude, 
                              fill = as.factor(Orbit_Uncertainity))) +
  geom_boxplot() +
  labs(x = "Orbit Uncertainity", y = " Absolute Magnitude") + 
  theme(legend.position = "none") 

grid.arrange(a1, a2, a3, a4, nrow = 2)
a5

# Some useful definition

library(ROSE)
library(caret)
library(cowplot)
library(leaps)

# We use the function "ovun.sample" from the "ROSE" package that
# creates possibly balanced samples by random over-sampling
# minority examples and under-sampling majority examples.

train_balanced<- ovun.sample(Hazardous~., data = train, 
                             method = "both", p = 0.5, 
                             N = 3515, seed = 1)$data


# Calling the "data" part of the output we obtain the resulting
# new dataset.

# We define three different thresholds for the classification task:
# 0.4, 0.5, 0.6

threshold4<- 0.4
threshold5<- 0.5
threshold6<- 0.6

# Simple Logistic Regression

# Model definition:

glm_compl<- glm(data = train,
                Hazardous ~ .,
                family = "binomial")

# We compute the reference level R-Squared

s<- summary(glm_compl)
r2<- 1 - (s$deviance/s$null.deviance)

1/(1-r2)

# Using the VIF function and comparing the obtained values with the 
# computed quantity:
# (The process is done iteratively where we delete one variable at time)

VIF(glm_compl)

glm_compl<- glm(data = train,
                Hazardous ~.-Aphelion_Dist-Semi_Major_Axis-
                  Jupiter_Tisserand_Invariant-
                  Eccentricity-Mean_Motion-Est_Dia_in_KM_max,
                family = "binomial")

# Observation of the model summary:

summary(glm_compl)

# Computing the predictions with the model on the test set: 

pred_glm_compl<- predict(glm_compl, test, type = "response")

# Converting the prediction in {0,1} according to the chosen threshold:

pred_glm_compl_04<- ifelse(pred_glm_compl > threshold4, 1, 0)
pred_glm_compl_05<- ifelse(pred_glm_compl > threshold5, 1, 0)
pred_glm_compl_06<- ifelse(pred_glm_compl > threshold6, 1, 0)

# Confusion matrix with different thresholds

table(test$Hazardous, pred_glm_compl_04)
mean(pred_glm_compl_04!=test$Hazardous)

table(test$Hazardous, pred_glm_compl_05)
mean(pred_glm_compl_05!=test$Hazardous)

table(test$Hazardous, pred_glm_compl_06)
mean(pred_glm_compl_06!=test$Hazardous)

# Logistic Model with Stepwise Selection

library(leaps)
library(MASS)

# Here we don't re-apply the VIF method because we start from the
# previous result.

glm_compl<- glm(data = train,
                Hazardous ~.-Aphelion_Dist-Semi_Major_Axis-
                  Jupiter_Tisserand_Invariant-
                  Eccentricity-Mean_Motion-Est_Dia_in_KM_max,
                family = "binomial")

# Application of the Stepwise method, specifying that we consider
# both the forward and the backward directions. We consider as 
# reference metric the Akaike Information Criterion: 

glm_compl_step <- stepAIC(glm_compl, direction = "both", 
                          trace = FALSE)

# Observation of the model summary:

summary(glm_compl_step)


# Computing the predictions with the model on the test set:

pred_glm_compl_step = predict(glm_compl_step, test, type = "response")


# Converting the predictions in {0,1} according to the chosen threshold:

pred_glm_compl_step_04 = ifelse(pred_glm_compl_step > threshold4, 1, 0)
pred_glm_compl_step_05 = ifelse(pred_glm_compl_step > threshold5, 1, 0)
pred_glm_compl_step_06 = ifelse(pred_glm_compl_step > threshold6, 1, 0)

# Confusion matrix with different thresholds

table(test$Hazardous, pred_glm_compl_step_04)
mean(pred_glm_compl_step_04!=test$Hazardous)

table(test$Hazardous, pred_glm_compl_step_05)
mean(pred_glm_compl_step_05!=test$Hazardous)

table(test$Hazardous, pred_glm_compl_step_06)
mean(pred_glm_compl_step_06!=test$Hazardous)

# Logistic Model with Balancing (over-sampling and under-sampling)

# Model definition:

glm_bal<- glm(data = train_balanced,
              Hazardous ~ .,
              family = "binomial")

# We compute the reference level R-Squared

s<- summary(glm_bal)
r2<- 1 - (s$deviance/s$null.deviance)

1/(1-r2)

# Using the VIF function and comparing the obtained values with the 
# computed quantity:
# (The process is done iteratively where we delete one variable at time)

VIF(glm_bal)


glm_bal<- glm(data = train_balanced,
              Hazardous ~.-Aphelion_Dist-Semi_Major_Axis-
                Jupiter_Tisserand_Invariant-Eccentricity-
                Est_Dia_in_KM_max-Mean_Motion,
              family = "binomial")

# Observation of the model summary:

summary(glm_bal)


# Computing the predictions with the model on the test set:

pred_glm_bal<- predict(glm_bal, test, type = "response")


# Converting the predictions in {0,1} according to the chosen threshold:

pred_glm_bal_04<- ifelse(pred_glm_bal > threshold4, 1, 0)
pred_glm_bal_05<- ifelse(pred_glm_bal > threshold5, 1, 0)
pred_glm_bal_06<- ifelse(pred_glm_bal > threshold6, 1, 0)

# Confusion matrix with different thresholds

table(test$Hazardous, pred_glm_bal_04)
mean(pred_glm_bal_04!=test$Hazardous)

table(test$Hazardous, pred_glm_bal_05)
mean(pred_glm_bal_05!=test$Hazardous)

table(test$Hazardous, pred_glm_bal_06)
mean(pred_glm_bal_06!=test$Hazardous)

# Logistic Model with Balancing (weights)

# Definition of the weights

w<- rep(1, nrow(train))

sum(train$Hazardous==0)/sum(train$Hazardous==1)

w[train$Hazardous == 1]<- 5

# Model definition:

glm_weighted<- glm(data = train,
                   Hazardous ~ .,
                   family = "binomial", weights = w)

# We compute the reference level R-Squared

s<- summary(glm_weighted)
r2<- 1 - (s$deviance/s$null.deviance)

1/(1-r2)

# Using the VIF function and comparing the obtained values with the 
# computed quantity:
# (The process is done iteratively where we delete one variable at time)

VIF(glm_weighted)


glm_weighted<- glm(data = train_balanced,
                   Hazardous ~.-Aphelion_Dist-Semi_Major_Axis-
                     Jupiter_Tisserand_Invariant-Eccentricity-
                     Mean_Motion-Est_Dia_in_KM_max,
                   family = "binomial", weights = w)

# Observation of the model summary:

summary(glm_weighted)

# Computing the predictions with the model on the test set:

pred_glm_weighted<- predict(glm_weighted, test, type = "response")

# Converting the predictions in {0,1} according to the chosen threshold:

pred_glm_weighted_04<- ifelse(pred_glm_weighted > threshold4, 1, 0)
pred_glm_weighted_05<- ifelse(pred_glm_weighted > threshold5, 1, 0)
pred_glm_weighted_06<- ifelse(pred_glm_weighted > threshold6, 1, 0)

# Confusion matrix with different thresholds

table(test$Hazardous, pred_glm_weighted_04)
mean(pred_glm_weighted_04!=test$Hazardous)

table(test$Hazardous, pred_glm_weighted_05)
mean(pred_glm_weighted_05!=test$Hazardous)

table(test$Hazardous, pred_glm_weighted_06)
mean(pred_glm_weighted_06!=test$Hazardous)

# Logistic Model with Balancing and Stepsize

# Application of the Stepwise method :

# (Here we don't re-apply the VIF method because we start from the
# previous result.)

#We start from the glm_bal model where we have already applied the
# VIF process, then we proceed with the Stepwise method.

glm_bal_step<- stepAIC(glm_bal, direction = "both", 
                       trace = FALSE)

# Observing the summary of the model:

summary(glm_bal_step)


# Computing the predictions with the model on the test set:

pred_glm_bal_step<- predict(glm_bal_step, test, type = "response")


# Converting the predictions in {0,1} according to the chosen threshold:

pred_glm_bal_step_04 = ifelse(pred_glm_bal_step > threshold4, 1, 0)
pred_glm_bal_step_05 = ifelse(pred_glm_bal_step > threshold5, 1, 0)
pred_glm_bal_step_06 = ifelse(pred_glm_bal_step > threshold6, 1, 0)

# Confusion matrix with different thresholds

table(test$Hazardous, pred_glm_bal_step_04)
mean(pred_glm_bal_step_04!=test$Hazardous)

table(test$Hazardous, pred_glm_bal_step_05)
mean(pred_glm_bal_step_05!=test$Hazardous)

table(test$Hazardous, pred_glm_bal_step_06)
mean(pred_glm_bal_step_06!=test$Hazardous)

Haz_test<- test$Hazardous

# Comparison between correct classification and classification of models

ggplot(test, aes(x = Absolute_Magnitude,
                 y = Minimum_Orbit_Intersection,
                 color = Haz_test)) +
  geom_point()+
  labs(x = "Absolute Magnitude",
       y = "Minimum Orbit Intersection",
       color = "Hazardous")  +
  theme(legend.position = c(0.8, 0.8))

a<- ggplot(test, aes(x = Absolute_Magnitude,
                     y = Minimum_Orbit_Intersection,
                     color = as.factor(pred_glm_compl_04))) +
  geom_point()+
  labs(x = "Absolute Magnitude",
       y = "Minimum Orbit Intersection",
       color = "Hazardous",
       title = "Simple GLM : 0.4")  +
  theme(legend.position = c(0.8, 0.8))



b<- ggplot(test, aes(x = Absolute_Magnitude,
                     y = Minimum_Orbit_Intersection,
                     color = as.factor(pred_glm_compl_step_04))) +
  geom_point()+
  labs(x = "Absolute Magnitude",
       y = "Minimum Orbit Intersection",
       color = "Hazardous",
       title = "GLM with Stepwise : 0.4")  +
  theme(legend.position = c(0.8, 0.8))



c<- ggplot(test, aes(x = Absolute_Magnitude,
                     y = Minimum_Orbit_Intersection,
                     color = as.factor(pred_glm_bal_06))) +
  geom_point()+
  labs(x = "Absolute Magnitude",
       y = "Minimum Orbit Intersection",
       color = "Hazardous",
       title = "GLM with Balanced dataset : 0.6")  +
  theme(legend.position = c(0.8, 0.8))


d<- ggplot(test, aes(x = Absolute_Magnitude,
                     y = Minimum_Orbit_Intersection,
                     color = as.factor(pred_glm_bal_step_06))) +
  geom_point()+
  labs(x = "Absolute Magnitude",
       y = "Minimum Orbit Intersection",
       color = "Hazardous",
       title = "GLM with Balanced dataset and Stepwise : 0.6")  +
  theme(legend.position = c(0.8, 0.8))


grid.arrange(a, b, c, d, nrow = 2)

# Logistic Regression curve

# We compare the results obtained with the four different models, plotting now an estimation of the logistic curve using the predictions given by the models:

predicted_data<- data.frame(prob.of.Haz = pred_glm_compl, Haz = test$Hazardous)
predicted_data<- predicted_data[order(predicted_data$prob.of.Haz, decreasing = FALSE),]
predicted_data$rank<-  1:nrow(predicted_data)

a<- ggplot(data = predicted_data, aes(x = rank, y = prob.of.Haz)) +
  geom_point(aes(color = as.factor(Haz)), alpha = 1, shape = 1, stroke = 1) +
  xlab("Index")+
  ylab("Predicted probability")+
  ggtitle("Estimated Logistic Curve - Simple GLM")


predicted_data<- data.frame(prob.of.Haz = pred_glm_compl_step, Haz = test$Hazardous)
predicted_data<- predicted_data[order(predicted_data$prob.of.Haz, decreasing = FALSE),]
predicted_data$rank<- 1:nrow(predicted_data)

b<- ggplot(data = predicted_data, aes(x = rank, y = prob.of.Haz)) +
  geom_point(aes(color = as.factor(Haz)), alpha = 1, shape = 1, stroke = 1) +
  xlab("Index")+
  ylab("Predicted probability")+
  ggtitle("Estimated Logistic Curve - GLM with Stepwise")



predicted_data<- data.frame(prob.of.Haz = pred_glm_bal, Haz = test$Hazardous)
predicted_data<- predicted_data[order(predicted_data$prob.of.Haz, decreasing = FALSE),]
predicted_data$rank<- 1:nrow(predicted_data)

c<- ggplot(data = predicted_data, aes(x = rank, y = prob.of.Haz)) +
  geom_point(aes(color = as.factor(Haz)), alpha = 1, shape = 1, stroke = 1) +
  xlab("Index")+
  ylab("Predicted probability")+
  ggtitle("Estimated Logistic Curve - GLM with Balanced data")



predicted_data<- data.frame(prob.of.Haz = pred_glm_bal_step, Haz = test$Hazardous)
predicted_data<- predicted_data[order(predicted_data$prob.of.Haz, decreasing = FALSE),]
predicted_data$rank<- 1:nrow(predicted_data)

d<- ggplot(data = predicted_data, aes(x = rank, y = prob.of.Haz)) +
  geom_point(aes(color = as.factor(Haz)), alpha = 1, shape = 1, stroke = 1) +
  xlab("Index")+
  ylab("Predicted probability")+
  ggtitle("Estimated Logistic Curve - GLM with Balanced data and Stepwise")


library(gridExtra)
grid.arrange(a, b, c, d, nrow = 2)

# Discriminant Analysis (requirements to check)

# Normality
# We apply the Shapiro - Wilks test on each covariate, considering the two different classes:

shapiro.test(train_balanced$Absolute_Magnitude[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Absolute_Magnitude[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Est_Dia_in_KM_max[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Est_Dia_in_KM_max[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Relative_Velocity_km_per_hr[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Relative_Velocity_km_per_hr[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Miss_Dist_Astronomical[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Miss_Dist_Astronomical[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Orbit_Uncertainity[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Orbit_Uncertainity[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Minimum_Orbit_Intersection[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Minimum_Orbit_Intersection[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Jupiter_Tisserand_Invariant[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Jupiter_Tisserand_Invariant[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Eccentricity[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Eccentricity[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Semi_Major_Axis[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Semi_Major_Axis[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Inclination[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Inclination[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Asc_Node_Longitude[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Asc_Node_Longitude[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Orbital_Period[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Orbital_Period[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Perihelion_Distance[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Perihelion_Distance[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Perihelion_Arg[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Perihelion_Arg[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Mean_Anomaly[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Mean_Anomaly[train_balanced$Hazardous==1])
shapiro.test(train_balanced$Mean_Motion[train_balanced$Hazardous==0])
shapiro.test(train_balanced$Mean_Motion[train_balanced$Hazardous==1])

#Outliers

# We consider the previous output given by the summary of the dataset.

# We report here some examples of outliers removal. Then we will 
# define the new training dataset without these examples.

# Absolute Magnitude

g1<- ggplot(data = train, aes(y = Absolute_Magnitude,fill = 2)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16,
               outlier.size = 2)+
  theme(legend.position="none") +
  ylab("Absolute Magnitude")

# We look for the presence of outliers:

chisq.out.test(train$Absolute_Magnitude)

# Removal of the found outlier:

which(train$Absolute_Magnitude == 11.16) # 256


# We do the same:

g2<- ggplot(data = train, aes(y = Est_Dia_in_KM_max,fill = 2)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16,
               outlier.size = 2)+
  theme(legend.position="none") +
  ylab("Estimated Diameter")

chisq.out.test(train$Est_Dia_in_KM_max)

which(train$Est_Dia_in_KM_max == 34.836938254) # 256


g3<- ggplot(data = train, aes(y = Orbital_Period,fill = 2)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16,
               outlier.size = 2)+
  theme(legend.position="none") +
  ylab("Orbital Period")

chisq.out.test(train$Orbital_Period)

which(train$Orbital_Period == 2912.0220196159) # 1556

grid.arrange(g1, g2, g3, nrow = 1)

train<- train[-c(256,1556),]

# Simple Linear Discriminant Analysis LDA

library(MASS)

# Model definition:

lda_compl<- lda(Hazardous ~ . - Aphelion_Dist - Semi_Major_Axis - 
                  Jupiter_Tisserand_Invariant - Eccentricity - Mean_Motion - 
                  Est_Dia_in_KM_max, family = "binomial", data = train)

# Observing the model summary:

lda_compl

# Computing predictions:

pred_lda_compl<- predict(lda_compl, test, type = "response") # threshold: 0.5
post_lda_compl<- pred_lda_compl$posterior

# Converting the predictions in {0,1} according to the chosen threshold:

pred_lda_compl_04<- as.factor(ifelse(post_lda_compl[,2] > threshold4, 1, 0))
pred_lda_compl_05<- pred_lda_compl$class
pred_lda_compl_06<- as.factor(ifelse(post_lda_compl[,2] > threshold6, 1, 0))

# Confusion matrix with different thresholds

table(test$Hazardous, pred_lda_compl_04)
mean(pred_lda_compl_04!=test$Hazardous)

table(test$Hazardous, pred_lda_compl_05)
mean(pred_lda_compl_05!=test$Hazardous)

table(test$Hazardous, pred_lda_compl_06)
mean(pred_lda_compl_06!=test$Hazardous)

pred_lda_compl_02<- as.factor(ifelse(post_lda_compl[,2] > 0.2, 1, 0))
pred_lda_compl_03<- as.factor(ifelse(post_lda_compl[,2] > 0.3, 1, 0))

table(test$Hazardous, pred_lda_compl_02)
mean(pred_lda_compl_02!=test$Hazardous)

table(test$Hazardous, pred_lda_compl_03)
mean(pred_lda_compl_03!=test$Hazardous)

# We use now the information given by:
# - x: linear combination of the variables that better describe the examples
# - class: assigned class

ldahist(pred_lda_compl$x[,1], g = pred_lda_compl$class, col = 2)

# Linear Discriminant Analysis LDA with balancing

# We consider the previous output given by the summary of the dataset.

# We report here some examples of outliers removal. Then we will 
# define the new training dataset without these examples.

# Absolute Magnitude

g4<-ggplot(data = train_balanced, aes(y = Relative_Velocity_km_per_hr,fill = 2)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16,
               outlier.size = 2)+
  theme(legend.position="none") +
  ylab("Relative Velocity")

# We look for the presence of outliers:

chisq.out.test(train_balanced$Relative_Velocity_km_per_hr)

# Removal of the found outlier:

which(train_balanced$Relative_Velocity_km_per_hr >= 160681.487851189) # 2313 2580


# We do the same:

g5<- ggplot(data = train_balanced, aes(y = Est_Dia_in_KM_max,fill = 2)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16,
               outlier.size = 2)+
  theme(legend.position="none") +
  ylab("Estimated Diameter")

chisq.out.test(train_balanced$Est_Dia_in_KM_max)

which(train_balanced$Est_Dia_in_KM_max == 34.836938254) # 1301 1703


g6<- ggplot(data = train_balanced, aes(y = Inclination,fill = 2)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16,
               outlier.size = 2)+
  theme(legend.position="none") +
  ylab("Inclination")

chisq.out.test(train_balanced$Inclination)

which(train_balanced$Inclination >= 75.406666841) # 3044

grid.arrange(g4, g5, g6, nrow = 1)

train_balanced<- train_balanced[-c(1301, 1703, 2313, 2580, 3044),]

# Model definition starting from the previous glm_bal model:

lda_bal<- lda(data = train_balanced,
              Hazardous ~.-Aphelion_Dist-Semi_Major_Axis-
                Jupiter_Tisserand_Invariant-Eccentricity-
                Est_Dia_in_KM_max-Mean_Motion,
              family = "binomial")

# Observing the model summary:

lda_bal

# Computing the predictions with the model on the test set:

pred_lda_bal<- predict(lda_bal, test, type = "response")
post_lda_bal<- pred_lda_bal$posterior

# Converting the predictions in {0,1} according to the chosen threshold:

pred_lda_bal_03<- as.factor(ifelse(post_lda_bal[,2] > 0.3, 1, 0))
pred_lda_bal_04<- as.factor(ifelse(post_lda_bal[,2] > threshold4, 1, 0))
pred_lda_bal_05<- pred_lda_bal$class
pred_lda_bal_06<- as.factor(ifelse(post_lda_bal[,2] > threshold6, 1, 0))

# Confusion matrix with different thresholds

table(test$Hazardous, pred_lda_bal_03)
mean(pred_lda_bal_03!=test$Hazardous)

table(test$Hazardous, pred_lda_bal_04)
mean(pred_lda_bal_04!=test$Hazardous)

table(test$Hazardous, pred_lda_bal_05)
mean(pred_lda_bal_05!=test$Hazardous)

table(test$Hazardous, pred_lda_bal_06)
mean(pred_lda_bal_06!=test$Hazardous)

ldahist(pred_lda_bal$x[,1], g = pred_lda_bal$class, col = 2)



# Simple Quadratic Discriminant Analysis QDA

# Model definition starting from the previous glm_compl:

qda_compl<- qda(Hazardous ~ . - Aphelion_Dist - Semi_Major_Axis - 
                  Jupiter_Tisserand_Invariant - Eccentricity - Mean_Motion - 
                  Est_Dia_in_KM_max- Orbit_Uncertainity, 
                family = "binomial", data = train)


# Computing predictions:

pred_qda_compl<- predict(qda_compl, test, type = "response") # threshold: 0.5
post_qda_compl<- pred_qda_compl$posterior


# Converting the predictions in {0,1} according to the chosen threshold:

pred_qda_compl_03<- as.factor(ifelse(post_qda_compl[,2] > 0.3, 1, 0))
pred_qda_compl_04<- as.factor(ifelse(post_qda_compl[,2] > threshold4, 1, 0))
pred_qda_compl_05<- pred_qda_compl$class
pred_qda_compl_06<- as.factor(ifelse(post_qda_compl[,2] > threshold6, 1, 0))

# Confusion matrix with different thresholds

table(test$Hazardous, pred_qda_compl_03)
mean(pred_qda_compl_03!=test$Hazardous)

table(test$Hazardous, pred_qda_compl_04)
mean(pred_qda_compl_04!=test$Hazardous)

table(test$Hazardous, pred_qda_compl_05)
mean(pred_qda_compl_05!=test$Hazardous)

table(test$Hazardous, pred_qda_compl_06)
mean(pred_qda_compl_06!=test$Hazardous)

# Quadratin Discriminant Analysis with balancing

# Model definition starting from the previous glm_bal model:

qda_bal<- qda(data = train_balanced,
              Hazardous ~.-Aphelion_Dist-Semi_Major_Axis-
                Jupiter_Tisserand_Invariant-Eccentricity-
                Est_Dia_in_KM_max-Mean_Motion- Orbit_Uncertainity,
              family = "binomial")

# Observing the model summary:

qda_bal

# Computing the predictions with the model on the test set:

pred_qda_bal<- predict(qda_bal, test, type = "response")
post_qda_bal<- pred_qda_bal$posterior

# Converting the predictions in {0,1} according to the chosen threshold:

pred_qda_bal_03<- as.factor(ifelse(post_qda_bal[,2] > 0.3, 1, 0))
pred_qda_bal_04<- as.factor(ifelse(post_qda_bal[,2] > threshold4, 1, 0))
pred_qda_bal_05<- pred_qda_bal$class
pred_qda_bal_06<- as.factor(ifelse(post_qda_bal[,2] > threshold6, 1, 0))

# Confusion matrix with different thresholds

table(test$Hazardous, pred_qda_bal_03)
mean(pred_qda_bal_03!=test$Hazardous)

table(test$Hazardous, pred_qda_bal_04)
mean(pred_qda_bal_04!=test$Hazardous)

table(test$Hazardous, pred_qda_bal_05)
mean(pred_qda_bal_05!=test$Hazardous)

table(test$Hazardous, pred_qda_bal_06)
mean(pred_qda_bal_06!=test$Hazardous)

# Ridge Regression

# We use here the balanced training dataset because we have seen that, generally, we
# obtain better results.

train_bal_mat<- as.matrix(train_balanced[,-19])
test_mat<- as.matrix(test[,-19])

# We look for the best value for lambda in order to define the best model:

# We use cross validation glmnet

ridge_cv <- cv.glmnet(train_bal_mat, train_balanced$Hazardous,
                      alpha = 0, family = "binomial", type.measure = "class")

plot(ridge_cv)

# We identify th best lambda value

lambda_opt_ridge <- ridge_cv$lambda.min
lambda_opt_ridge

# We compute predictions with the ridge model using the best value for 
# lambda that we have obtained

pred_ridge<- predict(ridge_cv, test_mat, type = "class", s = lambda_opt_ridge)

table(test$Hazardous, pred_ridge)

# Lasso Regression

# We look for the best value for lambda in order to define the best model:

# We use cross validation glmnet

lasso_cv <- cv.glmnet(train_bal_mat, train_balanced$Hazardous, 
                      alpha = 1, family = "binomial", type.measure = "class")

plot(lasso_cv)

# We identify th best lambda value

lambda_opt_lasso <- lasso_cv$lambda.min
lambda_opt_lasso

# We compute predictions with the lasso model using the best value for 
# lambda that we have obtained

pred_lasso<- predict(lasso_cv, test_mat, type = "class", s = lambda_opt_lasso)

table(test$Hazardous, pred_lasso)

# K-Nearest Neighbors 

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# We normalize the columns

nasa_n <- as.data.frame(lapply(nasa[,-19], min_max_norm))

# We re-add the Hazardous features to re-complete the dataset

nasa_n$Hazardous <- nasa$Hazardous

# We re-do the split obtaining the same differentiation as before

set.seed(0607)

split_n <- initial_split(nasa_n, prop = 0.75) 
train_n <- training(split_n)
test_n <- testing(split_n)

# We balance again our dataset

train_n_balanced<- ovun.sample(Hazardous~., data = train_n, 
                               method = "under", p = 1/5.35,
                               seed = 1)$data

# We implement a first KNN model removing the variables that show collinearity

library(class)

# We look now for the best value of the parameter k

kmax <- 100

test_error <- numeric(kmax)

# For each possible value of k we consider the obtained accuracy of the model

for (k in 1:kmax) {
  knn_pred <- as.factor(knn(train_n_balanced[,-c(2,7,8,9,15,18,19)],
                            test_n[,-c(2,7,8,9,15,18,19)], 
                            cl = train_n_balanced$Hazardous, k = k))
  
  cm <- confusionMatrix(data = knn_pred, reference = as.factor(test_n$Hazardous))
  
  test_error[k] <- 1 - cm$overall[1]
}

# We took the minimum value of the error

k_min <- which.min(test_error)

# We compute now the prediction with the value of k that gives us the minimum error

knn<- knn(train_n_balanced[,-c(2,7,8,9,15,18,19)], test_n[,-c(2,7,8,9,15,18,19)], 
          cl = train_n_balanced$Hazardous, k = k_min)

knn_pred_min <- as.factor(knn)

# Confusion matrix for KNN on the test set

tab<- table(test_n$Hazardous, knn)
tab
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
accuracy(tab)

cm <- confusionMatrix(data = knn_pred_min, reference = as.factor(test_n$Hazardous))

ggplot(data.frame(test_error), aes(x = 1:kmax, y = test_error)) +
  geom_line(colour="blue") +
  geom_point(colour="blue") +
  xlab("K (#neighbors)") + ylab("Test error") +
  ggtitle(paste0("Best value of K = ", k_min,
                 " (minimal error = ",
                 format((test_error[k_min])*100, digits = 4), "%)"))

ggplot(test, aes(x = Absolute_Magnitude,
                 y = Minimum_Orbit_Intersection,
                 color = as.factor(knn))) +
  geom_point()+
  labs(x = "Absolute Magnitude",
       y = "Minimum Orbit Intersection",
       color = "Hazardous",
       title = "KNN with K = 9")  +
  theme(legend.position = c(0.8, 0.8))

# Final comparisons

# Best GLM model

glm_best<- glm(data = train_balanced,
               Hazardous ~ Absolute_Magnitude+Minimum_Orbit_Intersection+
                 Orbit_Uncertainity+Orbital_Period,
               family = "binomial")

pred_glm_best<- predict(glm_best, test, type = "response")

# Best LDA model

lda_best<- lda(data = train_balanced,
               Hazardous ~Absolute_Magnitude+Miss_Dist_Astronomical+
                 Orbit_Uncertainity+Minimum_Orbit_Intersection,
               family = "binomial")

pred_lda_best<- predict(lda_best, test, type = "response")
post_lda_best<- pred_lda_best$posterior

# Best QDA model

qda_best<- qda(Hazardous ~ . - Aphelion_Dist - Semi_Major_Axis - 
                 Jupiter_Tisserand_Invariant - Eccentricity - Mean_Motion - 
                 Est_Dia_in_KM_max- Orbit_Uncertainity, 
               family = "binomial", data = train)

pred_qda_best<- predict(qda_best, test, type = "response")
post_qda_best<- pred_qda_best$posterior

# Best Ridge model

ridge_best<- glmnet(train_bal_mat, train_balanced$Hazardous, 
                    alpha = 0, family = "binomial", lambda = lambda_opt_ridge)

pred_ridge_best<- predict(ridge_best, test_mat, type = "response", s = lambda_opt_ridge)

# Best Lasso model

lasso_best<- glmnet(train_bal_mat, train_balanced$Hazardous, 
                    alpha = 0, family = "binomial", lambda = lambda_opt_lasso)

pred_lasso_best<- predict(lasso_best, test_mat, type = "response", s = lambda_opt_lasso)

prediction <- tibble(truth = as.factor(test$Hazardous))
prediction <- prediction %>% mutate(pred = as.numeric(pred_glm_best))%>%
  mutate(model= "GLM")%>% 
  add_row(truth = as.factor(test$Hazardous), pred = as.numeric(post_qda_best[,2]), 
          model= "LDA")%>%
  add_row(truth = as.factor(test$Hazardous), pred = as.numeric(post_qda_best[,2]), 
          model= "QDA")%>%
  add_row(truth = as.factor(test$Hazardous), pred = as.numeric(pred_ridge_best), 
          model= "Ridge")%>%
  add_row(truth = as.factor(test$Hazardous), pred = as.numeric(pred_lasso_best), 
          model= "Lasso")


roc <- prediction %>% group_by(model) %>% 
  roc_curve(truth, pred, event_level = "second") %>% 
  autoplot()
roc

auc(test$Hazardous, pred_glm_best)
auc(test$Hazardous, post_lda_best[,2])
auc(test$Hazardous, post_qda_best[,2])
auc(test$Hazardous, pred_ridge_best)
auc(test$Hazardous, pred_lasso_best)

e<- ggplot(test, aes(x = Absolute_Magnitude, y =
                       Minimum_Orbit_Intersection, color = as.factor(ifelse(pred_glm_best > 0.5, 1, 0)))) +
  geom_point()+ labs(x = "Absolute Magnitude", y = "Minimum Orbit
Intersection", color = "Hazardous", title = "GLM") + theme(legend.position = c(0.8, 0.8))

f<- ggplot(test, aes(x = Absolute_Magnitude, y =
                       Minimum_Orbit_Intersection, color =
                       as.factor(pred_qda_best$class))) + geom_point()+ labs(x = "Absolute
Magnitude", y = "Minimum Orbit Intersection", color = "Hazardous", title
                                                                             = "LDA") + theme(legend.position =
                                                                                                c(0.8, 0.8))

g<- ggplot(test, aes(x = Absolute_Magnitude, y =
                       Minimum_Orbit_Intersection, color = as.factor(pred_qda_best$class))) +
  geom_point()+ labs(x = "Absolute Magnitude", y = "Minimum Orbit
Intersection", color = "Hazardous", title = "QDA") + theme(legend.position = c(0.8, 0.8))

h<- ggplot(test, aes(x = Absolute_Magnitude, y =
                       Minimum_Orbit_Intersection, color = as.factor(ifelse(pred_ridge_best > 0.5, 1, 0)))) +
  geom_point()+ labs(x = "Absolute Magnitude", y = "Minimum Orbit
Intersection", color = "Hazardous", title = "Ridge") + theme(legend.position = c(0.8,
                                                                                 0.8))

i<- ggplot(test, aes(x = Absolute_Magnitude, y =
                       Minimum_Orbit_Intersection, color = as.factor(ifelse(pred_lasso_best > 0.5, 1, 0)))) +
  geom_point()+ labs(x = "Absolute Magnitude", y = "Minimum Orbit
Intersection", color = "Hazardous", title = "Lasso") + theme(legend.position = c(0.8,
                                                                                 0.8))

grid.arrange(e,f,g, h, i, nrow = 2)