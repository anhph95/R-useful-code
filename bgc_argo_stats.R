# Habitat paper stats
# Clean up
rm(list=ls())

library(dplyr)
library(tidyr)
library(readxl)
library(rstatix)
library(multcompView)
library(purrr)
# Functions
cld <- function(posthoc,hbt=hbt){
  temp<-data.frame(Habitat=hbt,group=length(hbt):1)
  pval <- posthoc$p.adj
  names(pval) <- paste(posthoc$group1,posthoc$group2,sep='-')
  mycld <- multcompLetters3('Habitat','group',pval,temp)
  group <- setNames(list(unique(posthoc$Def_variable)),'Def_Variable')
  return(as.data.frame(c(group,mycld$Letters)))
}

append_list_as_row <- function(df, lst) {
  # Create a complete row with NA for missing columns
  complete_row <- setNames(vector("list", length(df)), names(df))
  for (col in names(lst)) {
    complete_row[[col]] <- lst[[col]]
  }
  return(rbind(df, as.data.frame(t(complete_row), stringsAsFactors = FALSE)))
}

# BGC data
bgc <- read.csv(file.choose())
bgc$Habitat <- factor(bgc$Cluster,levels=c(1,3,2,4,5,6),ordered = TRUE)
#output <- "Habitat_defining_variable_stats.txt"
hbt <- factor(c(1,3,2,4,5,6),levels=c(1,3,2,4,5,6),ordered=T)
# Long format
def_var <- c('SST','SSS','MLD_2','DCM','NAI')
bgc_long <- bgc %>% pivot_longer(cols = all_of(def_var), names_to = "Def_variable", values_to = "Def_value")
bgc_long$Def_variable <- factor(bgc_long$Def_variable,levels=def_var,ordered=T)
#capture.output(cat('Habitat defining variable comparisons',sep="\n"),file=output)

# Normality
sw_test <- bgc_long %>% group_by(Def_variable,Habitat) %>% shapiro_test(Def_value) %>% as.data.frame()
sw_test <- sw_test %>% mutate(Sig. = case_when(p<=0.0001~"***",p<=0.001~"**",p<=0.01~"*",p<=0.05~".",TRUE~""))
#capture.output(cat("\n", 'Table 1. Shapiro-Wilk test for normality of habitat types',sep="\n"),file=output,append=T)
#capture.output(sw_test,file=output,append=T)
#capture.output(cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1",sep="\n"),file=output,append=T)

# Test of homogeneity
lv_test <- bgc_long %>%  group_by(Def_variable) %>% levene_test(Def_value ~ Habitat) %>% as.data.frame()
lv_test <- lv_test %>% mutate(Sig. = case_when(p<=0.0001~"***",p<=0.001~"**",p<=0.01~"*",p<=0.05~".",TRUE~""))
#capture.output(cat("\n", 'Table 2. Levene test for homogeneity of habitat types',sep="\n"),file=output,append=T)
#capture.output(lv_test,file=output,append=T)
#capture.output(cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1",sep="\n"),file=output,append=T)

# Welch's ANOVA
welch <- bgc_long %>% group_by(Def_variable) %>% welch_anova_test(Def_value~Habitat) %>% as.data.frame() %>% select(-'.y.')
welch <- welch %>% mutate(Sig. = case_when(p<=0.0001~"***",p<=0.001~"**",p<=0.01~"*",p<=0.05~".",TRUE~""))
#capture.output(cat("\n", "Table 3. Welch's ANOVA test for mean comparison of habitat types",sep="\n"),file=output,append=T)
#capture.output(welch,file=output,append=T)
#capture.output(cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",sep="\n"),file=output,append=T)

# Games Howell post hoc
gh_posthoc <- bgc_long %>% group_by(Def_variable) %>% games_howell_test(Def_value~Habitat) %>% as.data.frame() %>% select(-'.y.')
# Create cld table
posthoc_cld <- list()
for (i in def_var){
  result<-gh_posthoc %>% filter(Def_variable==i) %>% cld(hbt=hbt)
  posthoc_cld[[i]] <- result
}
posthoc_cld <- do.call(rbind, posthoc_cld)
rownames(posthoc_cld) <- NULL
#capture.output(cat("\n", "Table 4. Games-Howell post hoc test for mean comparison of habitat types",sep="\n"),file=output,append=T)
#capture.output(posthoc_cld,file=output,append=T)

# Kruskal-Wallis
kw <- bgc_long %>% group_by(Def_variable) %>% kruskal_test(Def_value~Habitat) %>% as.data.frame() %>% select(-'.y.')
kw <- kw %>% mutate(Sig. = case_when(p<=0.0001~"***",p<=0.001~"**",p<=0.01~"*",p<=0.05~".",TRUE~""))
dunn <- bgc_long %>% group_by(Def_variable) %>% dunn_test(Def_value~Habitat,p.adjust.method='hochberg') %>% as.data.frame() %>% select(-'.y.')
posthoc_cld <- list()
for (i in def_var){
  result<- dunn %>% filter(Def_variable==i) %>% cld(hbt=hbt)
  posthoc_cld[[i]] <- result
}
posthoc_cld <- do.call(rbind, posthoc_cld)
rownames(posthoc_cld) <- NULL

# Comparison among cruise
welch_cruise <- bgc_long %>%
  group_by(Habitat, Def_variable, Cruise) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_split(Def_variable,Habitat) %>%
  map_dfr(.f = function(df) {
    if(n_distinct(df$Cruise) > 1) {
      welch_test_result <- df %>%
        welch_anova_test(Def_value ~ Cruise) %>%
        as.data.frame() %>%
        select(-'.y.')
      # Add Habitat and Def_variable and rearrange columns
      welch_test_result <- cbind(Def_variable = unique(df$Def_variable),Habitat = unique(df$Habitat), welch_test_result)
    } else {
      # Return only Habitat and Def_variable for groups that don't have enough data
      data.frame(Def_variable = unique(df$Def_variable),Habitat = unique(df$Habitat))
    }
  })
welch_cruise <- welch_cruise %>% mutate(Sig. = case_when(p<0.001~"***",p<0.01~"**",p<0.05~"*",TRUE~""))

# Games-Howell for cruises
cruise_list <- c('KN197','MV1110','AT2104','EN614','EN640','M174')
posthoc_table <- setNames(data.frame(matrix(ncol = length(cruise_list), nrow = 0)), cruise_list)
complete_row <- setNames(vector("list", length(df)), names(df))
df <- welch_cruise
for(i in seq_along(df$Def_variable)){
  if (!is.na(df$p[i]) && df$p[i]<0.05){
    test <- bgc_long %>% filter(Def_variable==df$Def_variable[i],Habitat==df$Habitat[i])
    count <- table(test$Cruise)
    test <- test %>% filter(!(test$Cruise %in% names(count)[count<=1]))
    posthoc <- test %>% games_howell_test(Def_value~Cruise)
    pval <- posthoc$p.adj
    names(pval) <- paste(posthoc$group1,posthoc$group2,sep='-')
    cld <- as.list(multcompLetters(pval)$Letters)
    posthoc_table <- append_list_as_row(posthoc_table, cld)
  } else {
    posthoc_table <- append_list_as_row(posthoc_table, complete_row)
  }
}
posthoc_table[posthoc_table=='NULL'] <- ""
welch_cruise <- cbind(welch_cruise,posthoc_table)
options(width=1000)
capture.output(cat("\n", "Table 5. Welch ANOVA and Games-Howell post hoc for mean comparison of cruise",sep="\n"),file=output,append=T)
capture.output(welch_cruise,file=output,append=T)
capture.output(cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1",sep="\n"),file=output,append=T)
options(width=80)


test <- bgc_long %>%
  group_by(Habitat, Def_variable, Cruise) %>%
  filter(n() > 1) %>% group_by(Habitat) %>% mean(na.rm=T)

test <- bgc_long %>% filter(Def_variable=='SSS') %>% group_by(Habitat,Cruise) %>%
  summarise(across(all_of(c('Def_value')), 
                   ~ifelse(is.na(mean(., na.rm = TRUE)) && is.na(sd(., na.rm = TRUE)), 
                           NA, 
                           paste(round(mean(., na.rm = TRUE), 4), 
                                 "±", 
                                 round(sd(., na.rm = TRUE), 4)))))
test <- bgc_long %>% filter(Def_variable=='SSS',Habitat=='WPM')
ggplot(test,aes(x=Cruise,y=Def_value)) + geom_boxplot(notch = T)


?Chisquare
test <- bgc %>% filter(Cruise=='EN614') %>% group_by(Habitat) %>%
  summarise(across(all_of(c('SST','SSS')), 
                   ~ifelse(is.na(mean(., na.rm = TRUE)) && is.na(sd(., na.rm = TRUE)), 
                           NA, 
                           paste(round(mean(., na.rm = TRUE), 2), 
                                 "±", 
                                 round(sd(., na.rm = TRUE),2)))))
