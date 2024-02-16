# Clean workspace
rm(list=ls())
dev.off()

# Load library
library(readxl)
library(gridExtra)
library(tidyverse)
library(multcompView)
library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)
library(rstatix)
compare_mean <- function(dataset,x,y){
  # Formula
  myformula <- as.formula(paste(y,x,sep='~'))
  # Compare means test
  test <- dataset %>% kruskal_test(myformula) %>% as.data.frame() %>% select(-'.y.')
  test <- test %>% mutate(Sig. = case_when(p<=0.0001~"***",p<=0.001~"**",p<=0.01~"*",p<=0.05~".",TRUE~""))
  # Posthoc test
  posthoc <- dataset %>% dunn_test(myformula,p.adjust.method='bonferroni') %>% as.data.frame() %>% select(-'.y.')
  # Get CLD
  mygroup <- levels(dataset[[x]]) # Get levels of group
  temp<-data.frame(group=factor(mygroup,levels=mygroup,ordered=T),rank=length(mygroup):1) # Generate rank of group
  pval <- posthoc$p.adj # Get post hoc p-val
  names(pval) <- paste(posthoc$group1,posthoc$group2,sep='-') # Add name
  cld <- multcompLetters3('group','rank',pval,temp)$Letters # Extract cld
  # use factor groups as x coordinate
  x_cor <-1:length(cld)
  # Calculate y coordinate
  y_cor <- (max(dataset[[y]],na.rm=T)-min(dataset[[y]],na.rm=T))*0.2
  cld <- data.frame(cld,x_cor,y_cor)
  return(list(test=test,posthoc=posthoc,cld=cld))
}

# Box plot
mybox <- function(dataset,x,y,rev=F){
  # Color
  p <- ggplot(dataset,aes(.data[[x]],.data[[y]]))
  # Setting base theme
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 axis.text = element_text(size=12),axis.title = element_text(size=12),
                 plot.margin = margin(0.5,2,0.5,0.5,"cm"),
                 legend.position = "none")
  # Drawing whisker box
  p <- p + geom_boxplot(aes(fill=.data[[x]]),color='black',lwd=0.5)
  # Change border and fill colors
  p <- p + scale_fill_manual(values=mycolor)
  # Add CLD
  if (rev){
    p <- p + geom_text(data=compare_mean(dataset,x,y)$cld,aes(x=x_cor,y=min(dataset[[y]],na.rm=T)-y_cor,label=cld),size=6)
  } else {
    p <- p + geom_text(data=compare_mean(dataset,x,y)$cld,aes(x=x_cor,y=max(dataset[[y]],na.rm=T)+y_cor,label=cld),size=6)
  }
  return(p)
}

# generate_distinct_colors <- function(n) {
#   set.seed(123)
#   colors <- vector("character", length = n)
#   for (i in 1:n) {
#     hue <- (i -1) * (360 / n)  # Equally spaced hues
#     saturation <- 0.8  # Adjust saturation and luminance as desired
#     luminance <- 0.6
#     rgb <- grDevices::hcl(h = hue, c = saturation * 100, l = luminance * 100)
#     #color <- grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
#     colors[i] <- rgb
#   }
#   return(sample(colors))  # Randomly shuffle the colors
# }


# Load data
mydata <- read_csv(choose.files())
## Single factor analysis
mydata$Habitat <- factor(mydata$Cluster, levels=c(1,2,3,4,5,6),labels=c('TPI','TP','TGE','TGWS','TGSS','TGSF'),ordered = TRUE)
mycolor <- c('#df1f3e', '#dfcf1f', '#1fdf1f', '#1fd0df', '#3f1fdf', '#df1faf')
p1 <- mybox(mydata,'Habitat','MLD_2',rev = T)+
  ylab(expression(paste("Z"[MLD], " (m)")))+
  scale_y_reverse(breaks=seq(0,200,50))
p2 <- mybox(mydata,'Habitat','DCM',rev=T)+
  ylab(expression(paste("Z"[DCM], " (m)")))+
  scale_y_reverse(breaks=seq(0,200,50))
p3 <- mybox(mydata,'Habitat','SST')+
  ylab("SST (°C)")
p4 <- mybox(mydata,'Habitat','SSS')+
  ylab("SSS (psu)")
p5 <- mybox(mydata,'Habitat','NAI')+
  ylab("NAI")+
  scale_y_continuous(breaks=seq(-200,0,50))
p6 <- mybox(mydata,'Habitat','DNC',rev=T)+
  ylab("DNC")+
  scale_y_reverse(breaks=seq(0,200,25))
# Combining plots
grid.arrange(p3,p4,p1,p2,p5,ncol=2)

test <- compare_mean(mydata,'Habitat','SST')
########################################################
df <- filter(mydata,Habitat %in% c('TPI','TP','TGE','TGWS','TGSS','TGSF'))
hbt <- df$Habitat
plot(df$SSS,df$SST,pch=21,bg=mycolor[hbt])


table(mydata$Habitat,mydata$Province)


########################################################
df<-mydata
df <- na.omit(df)
pft <- select(df,'DIATO','DINO','HAPTO','GREEN','PROKAR','PROCHLO')
hbt <- df$Habitat
pft_chi <- decostand(pft,method='log',na.rm=TRUE)
nmds <- metaMDS(pft_chi,distance ="bray")
nmds$stress

nmds$points[,1]<- -nmds$points[,1]
# Plot & edit ordination
par(mar = c(5,5,2,2)) # Set plot margin 
par(mfrow=c(1,1))
plot(nmds$points,xlab="NMDS1",ylab="NMDS2",pch=21,cex=2,bg=mycolor[hbt],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,scaling=3,correlation=T,xlim=c(-0.25,0.45),ylim=c(-0.35,0.35),type='n')
abline(h=0,lty="dotted",lwd=2)
abline(v=0,lty="dotted",lwd=2)
ordiellipse(nmds$points, group = hbt,col=mycolor,draw='polygon',alpha=175,kind='sd',lwd=0.05,label=F)
ctd <- select(df,'SST','SSS','MLD_2','DCM','NAI')
fit <- envfit(nmds,ctd,permutations=999,na.rm=TRUE)
arrows(0,0,fit$vectors$arrows[,1]*0.25,fit$vectors$arrows[,2]*0.25,length=0.1,col='red',lwd=2)
lab <- ordiArrowTextXY(fit$vectors$arrows,rescale=FALSE)
text(lab[,1]*0.26,lab[,2]*0.26,labels=rownames(fit$vectors$arrows),cex=2,font=2,col='red')
#ordispider(nmds$points, group = hbt,col=mycolor)
score <- scores(nmds,display = 'species',choices=c(1,2))
score[,1] <- -score[,1]
points(score,pch="x",cex=1.5)
lab = c('Diatoms','Dinoflagellates','Haptophytes','Green algae','Prokaryotes','Prochlo')
text(score*1.1,cex=1.5,label=lab,font=3,adj=c(0,0.5))
# corr <- cor(pft_chi,scores(nmds,choices=c(1,2),display="sites"),method="pearson")
# arrows(0,0,corr[,1],corr[,2],length=0.1,col='black',lwd=4)
# lab <- ordiArrowTextXY(corr,rescale=FALSE)
# text(lab[,1],lab[,2],labels=rownames(corr),cex=0.75,col='black')
pft_dist = vegdist(pft_chi,method="bray")
adonis2(pft_dist~hbt)
pairwise.adonis(pft_dist,hbt)
# corr2 <- cor(ctd,scores(nmds,choices=c(1,2),display="sites"),method="pearson")
# arrows(0,0,corr2[,1],corr2[,2],length=0.1,col='red',lwd=4)
pft_cap <- cca(pft_chi~hbt,method='bray')


pft_cap <- capscale(pft_chi~ctd$SST+ctd$SSS+ctd$MLD_2+ctd$DCM+ctd$NAI,distance="bray",sqrt.dist=TRUE,comm=TRUE)
ev <- as.data.frame(summary(eigenvals(pft_cap,model="all"))) # Extract eigenvals of unconstrained model
ev_con <- as.data.frame(summary(eigenvals(pft_cap,model="constrained"))) # Extract eigenvals of constrained model

# Plot & edit CAP
par(mar = c(5,5,2,2)) # Set plot margin 
par(mfrow=c(1,2))
plot(pft_cap,type='n',
     xlab=bquote("dbRDA1 ("*.(round(ev[2,1]*100,1))*"% of total variation,"~.(round(ev_con[2,1]*100,1))*"% of fitted variation)"),
     ylab=bquote("dbRDA2 ("*.(round(ev[2,2]*100,1))*"% of total variation,"~.(round(ev_con[2,2]*100,1))*"% of fitted variation)"),
     cex.lab=1,cex.axis=1,
)
ordiellipse(pft_cap, group = hbt,col=mycolor,draw='polygon',alpha=127)
scrs <- scores(pft_cap,display="wa",choices=c(1,2))
scrs[,1] <- scrs[,1]
scrs <- jitter(scrs,factor=2,amount=0.25)
points(scrs,bg=mycolor[hbt],pch=21,cex=1.5,lwd=2)
legend("topright",legend=levels(hbt),col=co,pch=20)

orditorp(pft_cap,display='species')
# Calculate Pearson correlation of species scores with ordination axes
corr <- cor(pft_chi,scores(pft_cap,choices=c(1,2),display="sites",scaling=3),method="pearson")

plot.default(corr,type="n",xlim=range(-1.5,1.5),ylim=range(-1.5,1.5),xaxt='n',yaxt='n',
             xlab=paste("Correlation with dbRDA1"),
             ylab=paste("Correlation with dbRDA2"),
             cex.lab=1,cex.axis=1)
axis(1, at = seq(-1,1,by=0.5),cex.axis=1)
axis(2, at = seq(-1,1,by=0.5),cex.axis=1)
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
arrows(0,0,corr[,1],corr[,2],length=0.1,col='black',lwd=2)
lab <- ordiArrowTextXY(corr,rescale=FALSE)
text(lab[,1],lab[,2],labels=rownames(corr),cex=1,col='black')

df <- pft_chi
df$Habitat <- hbt
p1 <- mybox(df,'Habitat','DIATO')+ylab(expression(paste("Diatoms (log"[10],"(", mg~m^-3, "))")))
p2 <- mybox(df,'Habitat','DINO')+ylab(expression(paste("Dinophytes (log"[10],"(", mg~m^-3, "))")))
p3 <- mybox(df,'Habitat','HAPTO')+ylab(expression(paste("Haptophytes (log"[10],"(", mg~m^-3, "))")))
p4 <- mybox(df,'Habitat','GREEN')+ylab(expression(paste("Green algae (log"[10],"(", mg~m^-3, "))")))
p5 <- mybox(df,'Habitat','PROKAR')+ylab(expression(paste("Prokaryotes (log"[10],"(", mg~m^-3, "))")))
p6 <- mybox(df,'Habitat','PROCHLO')+ylab(expression(paste("Prochlorococcus (log"[10],"(", mg~m^-3, "))")))
# Combining plots
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)


subdata <- mydata#subset(mydata,CHL_profile=='DBM')
p1 <- mybox(subdata,'Habitat','CHL_Zeu') + ylab(expression(paste(group(langle, CHLA, rangle)[Zeu], " (",mg~m^-2, ")")))
p3 <- mybox(subdata,'Habitat','CHL_Zpd') + ylab(expression(paste(group(langle, CHLA, rangle)[Zpd], " (",mg~m^-2, ")")))
p5 <- mybox(subdata,'Habitat','CHL_Zlow') + ylab(expression(paste(group(langle, CHLA, rangle)[Zlow], " (",mg~m^-2, ")")))
p2 <- mybox(subdata,'Habitat','CHL_mean_Zeu') + ylab(expression(paste("[" * CHLA * "]"[Zeu], " (",mg~m^-3, ")")))
p4 <- mybox(subdata,'Habitat','CHL_mean_Zpd') + ylab(expression(paste("[" * CHLA * "]"[Zpd], " (",mg~m^-3, ")")))
p6 <- mybox(subdata,'Habitat','CHL_mean_Zlow') + ylab(expression(paste("[" * CHLA * "]"[Zlow], " (",mg~m^-3, ")")))
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
grid.arrange(p1,p3,p5,ncol=1)

p1 <- mybox(mydata,'Habitat','CHL') + ylab(expression(paste("CHL_merged ( ", mg~m^-3, ")")))
p2 <- mybox(subdata,'Habitat','CHL_surface') + ylab(expression(paste("CHL_Argo ( ", mg~m^-3, ")")))
p3 <- mybox(subdata,'Habitat','CHL_sat') + ylab(expression(paste("CHL_weighted_mean ( ", mg~m^-3, ")")))
p4 <- mybox(subdata,'Habitat','CHL_mean_Zeu') + ylab(expression(paste("CHL_Zeu_mean ( ", mg~m^-3, ")")))
p5 <- mybox(subdata,'Habitat','CHL_mean_Zpd') + ylab(expression(paste("CHL_Zpd_mean ( ", mg~m^-3, ")")))
p6 <- mybox(subdata,'Habitat','CHL_mean_Zlow') + ylab(expression(paste("CHL_Zlow_mean ( ", mg~m^-3, ")")))
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)

p1 <- mybox(subdata,'Habitat','DCM',rev=T) + scale_y_reverse(breaks=seq(0,200,25))  
p2 <- mybox(subdata,'Habitat','Zeu',rev=T) + scale_y_reverse(breaks=seq(0,200,25))
p3 <- mybox(subdata,'Habitat','Zpd',rev=T) + scale_y_reverse(breaks=seq(0,200,10))
p4 <- mybox(subdata,'Habitat','Zlow',rev=T) + scale_y_reverse(breaks=seq(0,500,100))
grid.arrange(p1,p2,p3,p4,ncol=2)

test <- filter(subdata,!(is.na(CHL))&(CHL<=0.35)&(CHL_mean_Zeu<=0.35)) 
p1 <- mybox(test,'Habitat','CHL_mean_Zeu') + ylab(expression(paste("[" * CHLA * "]"[Zeu], " (",mg~m^-3, ")")))
p2 <- mybox(test,'Habitat','CHL') + ylab(expression(paste("[" * CHLA * "]"[Sat], " (",mg~m^-3, ")")))
grid.arrange(p1,p2,ncol=1)


test <- subdata
test$sat_log <- log(subdata$CHL)
test$zeu_log <- log(subdata$CHL_mean_Zeu)
p1 <- mybox(test,'Habitat','zeu_log') + ylab(expression(paste("[" * CHLA * "]"[Zeu], " (",mg~m^-3, ")")))
p2 <- mybox(test,'Habitat','sat_log') + ylab(expression(paste("[" * CHLA * "]"[Sat], " (",mg~m^-3, ")")))
grid.arrange(p1,p2,ncol=1)

x='CHL_Zpdmean'
y='CHL'
df=subdata
linear<-lm(df[[y]]~df[[x]])
plot(df[[x]],df[[y]],pch=21,cex=2,bg=mycolor[df$Habitat],xlab=x,ylab=y)
abline(linear)
coefficients <- round(coef(linear),4)
# Construct and print the linear equation
cat("y =", coefficients[2], "* x +", coefficients[1], "\n")
cat("R-squared=", round(summary(linear)$r.squared,4), "\n")


library(tidyr)
df_long <- mydata %>% pivot_longer(cols = all_of(c('DCM','Zpd','Zeu','Zlow','MLD_2','DNC')), names_to = "Depth_name", values_to = "Depths")
df_long$Depth_name <- factor(df_long$Depth_name,levels=c('Zpd','Zeu','Zlow','DCM','MLD_2','DNC'),ordered=T)
ggplot(df_long, aes(x = Habitat, y = Depths, fill = Depth_name)) +
  geom_boxplot(notch=T) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size=12),axis.title = element_text(size=12),
          plot.margin = margin(0.5,2,0.5,0.5,"cm"))+
  scale_y_reverse(breaks=seq(0,400,50))+
  scale_fill_manual(values = c('yellow','gold','goldenrod3','green','gray','blue'))

library(ggpattern)
df_long <- mydata %>% pivot_longer(cols = all_of(c('DCM','Zpd','Zeu','Zlow','MLD_2','DNC')), names_to = "Depth_name", values_to = "Depths")
df_long$Depth_name <- factor(df_long$Depth_name,levels=c('Zpd','Zeu','Zlow','DCM','MLD_2','DNC'),ordered=T)
df_long <- df_long %>% mutate(temp = sub("\\s.*", "", Habitat))
p <- ggplot(df_long, aes(x = Habitat, y = Depths)) +
  geom_boxplot_pattern(aes(pattern=Depth_name,pattern_density=Depth_name,fill=Depth_name,pattern_size=Depth_name),outlier.size = 0.75,lwd=0.5,
                       notch=T, pattern_spacing=0.02,pattern_fill='black',width=0.5,position=position_dodge(width=0.7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=12),axis.title = element_text(size=12),
        plot.margin = margin(0.5,2,0.5,0.5,"cm"))+
  theme(legend.key.size = unit(2, 'cm'),legend.text=element_text(size=12),legend.key=element_rect(fill='white'))+
  scale_y_reverse(breaks=seq(0,300,50),limits=c(300,0))+
  scale_pattern_manual(values=c('stripe','stripe', 'stripe','stripe','wave','wave'))+
  scale_pattern_density_manual(values=c(0,0,0.2,0.2,0.2,0.2))+
  scale_pattern_size_manual(values=c(0.1,0.1,0.1,0.1,0.1,0.1))+
  scale_fill_manual(values=c('white','grey','white','grey','white','grey'))+
  guides(pattern=guide_legend(ncol=6))+
  theme(legend.position = "None")+
  ylab('Depth (m)')+xlab('Habitats')

library(cowplot)
library(grid)
leg <- get_legend(p)  
grid.newpage()
grid.draw(leg)
p

def_var <- factor(c('SST','SSS','MLD_2','DCM','NAI'),ordered=T)
def_var <- factor(c('CHL_Zeu','CHL_Zpd','CHL_Zlow','CHL_mean_Zeu','CHL_mean_Zpd','CHL_mean_Zlow'),ordered=T)
def_var <- factor(c('Zeu','Zpd','Zlow','CHL_max'),ordered=T)
def_var <- factor(c('CHL','DIATO','DINO','HAPTO','GREEN','PROKAR','PROCHLO'),ordered=T)
df_long <- mydata %>% pivot_longer(cols = all_of(def_var), names_to = "Def_variable", values_to = "Def_value")
result <- df_long %>% filter(!is.na(Def_value)) %>%
  group_by(Def_variable,Habitat) %>%
  summarize(
    n = n(),
    Median = round(median(Def_value),2),
    IQR = round(IQR(Def_value),2),
    `25th` = round(quantile(Def_value, probs = 0.25),2),
    `75th` = round(quantile(Def_value, probs = 0.75),2)
  ) %>% t() %>% as.data.frame()

write.csv(result,"summary_pft.csv")

install.packages("patternplot")
library(patternplot)
group<-df_long$Depth_name
y<-df_long$Depths
x<-df_long$Habitat

par(mfrow=c(1,1))
par(mar = c(5,5,5,5)) 
#groups = c(rep('low',500),rep('high',600))
ordiplot(nmds$points,type='n',cex.lab=2,cex.axis=2,choices=c(1,2),xlim=c(-3,5),ylim=c(-2,2))
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
#ordiellipse(nmds,groups=groups,label=F,lwd=2,col='red')
# score <- scores(nmds,display = 'sites',choices=c(1,2))
# points(score, bg=co[hbt],
#     pch=21, cex=2, main = "Pigment NMDS",
#     xlab="NMDS1",ylab="NMDS2")
#ordiellipse(nmds, group = hbt,label = T,col=co,draw='polygon',font=1,cex=1.5,alpha=170,kind='sd')
#ordiellipse(nmds, group = hbt,label = T,show.groups=c(2,8,9,10),col=co,draw='polygon',font=1,cex=1.5,alpha=170,kind='sd')
ordiellipse(nmds, group = hbt,label = F,col=mycolor,draw='polygon',font=1,cex=1.5,alpha=170,kind='sd')

fit <- envfit(nmds,ctd,permutations=999,na.rm=TRUE)
arrows(0,0,fit$vectors$arrows[,1],fit$vectors$arrows[,2],length=0.1,col='red',lwd=5)
lab <- ordiArrowTextXY(fit$vectors$arrows,rescale=FALSE)
text(lab[,1]*1.1,lab[,2]*1.1,labels=rownames(fit$vectors$arrows),cex=2,font=2,col='red')


score <- scores(nmds,display = 'species',choices=c(1,2))
points(score,pch="x",cex=2)
lab = c('Diatoms','Dinoflagellates','Haptophytes','Greenalgae','Prokaryotes','Prochlorococcus')
text(score*1.2,cex=1.5,label=lab,font=3,adj=c(0,0.5))
#legend("topright",legend=levels(hbt),col=co,pch=20)

def_var <- factor(c('DIATO','DINO','HAPTO','GREEN','PROKAR'),ordered=T)
df_long <- mydata %>% pivot_longer(cols = all_of(def_var), names_to = "Def_variable", values_to = "Def_value") 
df_long <- df_long %>%  na.omit(subset=Def_value) %>% group_by(Habitat) %>% mutate(Percentage = Def_value / sum(Def_value) * 100)
ggplot(df_long, aes(fill=Def_variable, y=Percentage, x=Habitat)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  labs(fill="Species", y="Count", x="Site", title="Species Composition")

test <- mydata[!(mydata$Cluster %in% c(1,2)),]
plot(test$PROKAR,test$PROCHLO)
abline(a = 0, b = 0.5, col = "red")


library(dplyr)
library(reshape)
install.packages("reshape2")
species_data <- pft
species_data$Site <- rownames(species_data)

# Melting the data for ggplot
long_data <- species_data %>% pivot_longer(
    cols = -Site, 
    names_to = "Species", 
    values_to = "Count"
)

# Joining NMDS coordinates
nmds_coords <- as.data.frame(scores(nmds$points))
colnames(nmds_coords) <- c("NMDS1", "NMDS2")
nmds_coords$Site <- rownames(nmds_coords)
nmds_coords$Habitat <- hbt
long_data <- left_join(long_data, nmds_coords, by = "Site")

# Plotting
ggplot(long_data, aes(x = NMDS1, y = NMDS2, size = Count, color = Species,fill=Habitat)) +
  geom_point(alpha = 0.75,pch=21) +
  theme_minimal() +
  labs(size = "Abundance", color = "Species", 
       x = "NMDS Axis 1", y = "NMDS Axis 2", 
       title = "Segmented Bubble Plot on NMDS")+
  scale_fill_manual(values=mycolor)
