# Libraries
library(ggplot2)
library(dplyr)

## One factor
# Example data
data(iris)

# Calculate means and standard errors
summary_data <- iris %>%
  group_by(Species) %>%
  summarise(your_mean_column_name = mean(Sepal.Length), your_sd_column_name = sd(Sepal.Length))

# Plot
ggplot(summary_data, aes(x = Species, y = your_mean_column_name)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", fill = "green") +
  geom_errorbar(aes(ymin = your_mean_column_name - your_sd_column_name, ymax = your_mean_column_name + your_sd_column_name), width = .2, position = position_dodge(.9)) +
  theme_minimal() +
  labs(title = "Example bar graph", x = "Species", y = "Sepal Length")

## Two factors
# Add another factor to the dataset
iris$SecondFactor <- sample(c('Low', 'Medium', 'High'), nrow(iris), replace = TRUE)

# Calculate means and standard errors
summary_data <- iris %>%
  group_by(Species, SecondFactor) %>%
  summarise(your_mean_column_name = mean(Sepal.Length), your_sd_column_name = sd(Sepal.Length),.groups = 'drop')

# Create the bar plot with error bars
ggplot(summary_data, aes(x = Species, y = your_mean_column_name, fill = SecondFactor)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = your_mean_column_name - your_sd_column_name, ymax = your_mean_column_name + your_sd_column_name), width = .2, position = position_dodge(0.9)) +
  theme_minimal() +
  labs(title = "Example bar graph", x = "Species", y = "Sepal Length")

