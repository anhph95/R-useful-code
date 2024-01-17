library(multcompView)
library(ggplot2)
library(dplyr)

# Function to generate compact letter display
compare_mean <- function(dataset,x,y){
  # list of cld
  cld_list <- levels(dataset[[x]])
  cld <- setNames(rep(NA, length(cld_list)), cld_list)
  # Test for normality
  shapiro <- by(dataset[[y]],dataset[[x]],shapiro.test)
  p_val <- sapply(shapiro, function(x) x$p.value)
    # Parametric test
  test <- aov(dataset[[y]]~dataset[[x]])
  if (summary(test)[[1]][1, 5] < 0.05){
    posthoc <- as.data.frame(TukeyHSD(test)$dataset)
    Diff <- posthoc$`p adj` < 0.05
    names(Diff) <- row.names(posthoc) 
    cld <- multcompLetters(Diff)$Letters
  }
  # use factor groups as x coordinate
  x_cor <- as.numeric(factor(names(cld),ordered = TRUE))
  # Calculate y coordinate
  y_cor <- (max(dataset[[y]],na.rm=T)-min(dataset[[y]],na.rm=T))*0.2
  cld <- data.frame(cld,x_cor,y_cor)
  return(list(shapiro=p_val,test=test,posthoc=posthoc,cld=cld))
}

# Box plot
mybox <- function(dataset,x,y,rev=F){
  # Color
  p <- ggplot(dataset,aes(.data[[x]],.data[[y]]))
  # Setting base theme
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 axis.text = element_text(size=10),axis.title = element_text(size=12,face="bold"),
                 plot.margin = margin(0.5,2,0.5,0.5,"cm"),
                 legend.position = "none")
  # Drawing whisker box
  p <- p + geom_boxplot(aes(fill=.data[[x]]),lwd=0.5)
  # Add CLD
  if (rev){
    p <- p + geom_text(data=compare_mean(dataset,x,y)$cld,aes(x=x_cor,y=min(dataset[[y]],na.rm=T)-y_cor,label=cld))
  } else {
    p <- p + geom_text(data=compare_mean(dataset,x,y)$cld,aes(x=x_cor,y=max(dataset[[y]],na.rm=T)+y_cor,label=cld))
  }
  return(p)
}