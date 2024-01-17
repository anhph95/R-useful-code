# Import data
data <- read.csv(file.choose())
# data <- read.csv(choose.files())
data

# Normality
# Sharpio-Wilk test
shapiro.test(data$Sepal.Length)

# Homoscedasticity
library(car)
leveneTest(data$Sepal.Length~data$Species)

### Statistical tests
## T.test
# Method 1
setosa <- data[data$Species=='setosa',] # Filter for setosa
virginica <- data[data$Species=='virginica',]
t.test(setosa$Sepal.Length,virginica$Sepal.Length,alternative='two.sided',paired=FALSE,var.equal=FALSE)

# Method 2
library(dplyr)
subdata <- data[data$Species %in% c('setosa','virginica'),] # Filter for both setosa and virginica
t.test(subdata$Sepal.Length~subdata$Species,alternative='two.sided',paired=FALSE,var.equal=FALSE)

## ANOVA
# One-way ANOVA
anova1 <- aov(data$Sepal.Length~data$Species)
summary(anova1)
TukeyHSD(anova1)

# Two-way ANOVA
data$Site <- rep(c('A','B','C'),50)
anova2 <- aov(data$Sepal.Length~data$Species*data$Site)
summary(anova2)

# Linear regression
linear <- lm(data$Sepal.Length~data$Petal.Width)
summary(linear)

# Multiple linear regression
mlinear <- lm(data$Sepal.Length~data$Petal.Width+data$Sepal.Width+data$Petal.Length)
summary(mlinear)


## Data visulization
library(ggplot2)
# Histogram
ggplot(data,aes(x=Sepal.Length,fill=Species))+
  geom_histogram()+
  geom_density(alpha=0)+
  theme_classic()

# Box plot
ggplot(data, aes(x=Species,y=Sepal.Length,fill=Species))+
  geom_boxplot()+
  theme_classic()+
  theme(text = element_text(size = 20))

# Bar plot
ggplot(data, aes(x=Species,y=Sepal.Length,fill=Species))+
  stat_summary(fun=mean,geom='bar')+
  stat_summary(fun=mean,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               geom='errorbar')+
  theme_classic()+
  theme(text = element_text(size = 20))

# Scatter plot
ggplot(data, aes(x=Sepal.Length,y=Petal.Width))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  theme(text = element_text(size = 20))
