#R script for Trait resolution simulations - summary analysis and plots  - Kohli and Jarzyna 2021 Global Ecology and Biogeography
#Author: Brooks Kohli
#Date Created: 19 Nov 2019
#Date last modified: 2 February 2021
#############################################


require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)


#####
#MEDIANS
####

mAll <- read.csv("AppendixS3_MedianSESs.csv") #this file is all the Median SES values across all simulation runs from the main text.  It is available as Supplemental Information.

#split convergence and divergence processes into separate files
mAll.conv <- mAll[mAll$proc=="Conv",]
mAll.conv

mAll.div <- mAll[mAll$proc=="Div",]
mAll.div


######
##Wilcoxon test estimates
####


#Read in  the appended file, disregarding the first header column (NNN.1 columns are p-values and are not needed for following analyses)
wAll <- read.csv("AppendixS4_WilcoxEstsPvals.csv", skip = 1)   #this file is all the Median SES values across all simulation runs from the main text.  It is available as Supplemental Information.

###Note that these were calculated without changing INf to NAs.  It is needed to plot Medians, but for tests, the Inf need to be reflected to get most reliable results.


#split convergence and divergence processes into separate files
wAll.conv <- wAll[wAll$proc=="Conv",]
wAll.conv

wAll.div <- wAll[wAll$proc=="Div",]
wAll.div

#############################################################################################

#summary stats reported in main text results - proportion of median SESs that shifted significantly
summary(wAll.conv)
summary(wAll.conv[wAll.conv$coarseness!="Rounded",])
summary(wAll.conv[wAll.conv$coarseness=="Cat02" | wAll.conv$coarseness=="Cat04",])

summary(wAll.div)
summary(wAll.div[wAll.div$coarseness!="Rounded",])
summary(wAll.div[wAll.div$coarseness=="Cat02" | wAll.div$coarseness=="Cat04",])


###conv
#count for proportions
wAll.conv.cat <- wAll.conv[wAll.conv$coarseness!="Rounded",]  ##16
summary(wAll.conv.cat)  

wAll.conv.24 <- wAll.conv[wAll.conv$coarseness=="Cat02" | wAll.conv$coarseness=="Cat04",]  ##8
summary(wAll.conv.24)

##number of tests with signif shifts from continuous traits (out of 16)
#***relies on the combined WilcoxEstsPvals file, with estimates calculated with Inf values***
nrow(wAll.conv.cat[wAll.conv.cat$MNND.1 <= 0.05, ])  
nrow(wAll.conv.cat[wAll.conv.cat$FDis.1 <= 0.05, ])  
nrow(wAll.conv.cat[wAll.conv.cat$MPD.1 <= 0.05, ])   
nrow(wAll.conv.cat[wAll.conv.cat$Rao.1 <= 0.05, ])   

#proportions (/16)
nrow(wAll.conv.cat[wAll.conv.cat$MNND.1 <= 0.05, ])/16  
nrow(wAll.conv.cat[wAll.conv.cat$FDis.1 <= 0.05, ])/16  
nrow(wAll.conv.cat[wAll.conv.cat$MPD.1 <= 0.05, ])/16   
nrow(wAll.conv.cat[wAll.conv.cat$Rao.1 <= 0.05, ])/16   

#proportions (/20  ==  including rounded contin)
nrow(wAll.conv[wAll.conv$MNND.1 <= 0.05, ])/20  
nrow(wAll.conv[wAll.conv$FDis.1 <= 0.05, ])/20  
nrow(wAll.conv[wAll.conv$MPD.1 <= 0.05, ])/20   
nrow(wAll.conv[wAll.conv$Rao.1 <= 0.05, ])/20  

#proportions (/8  == just 2 and 4 state categ)
nrow(wAll.conv.24[wAll.conv.24$MNND.1 <= 0.05, ])/8  
nrow(wAll.conv.24[wAll.conv.24$FDis.1 <= 0.05, ])/8  
nrow(wAll.conv.24[wAll.conv.24$MPD.1 <= 0.05, ])/8   
nrow(wAll.conv.24[wAll.conv.24$Rao.1 <= 0.05, ])/8  


###div
#count for proportions
wAll.div.cat <- wAll.div[wAll.div$coarseness!="Rounded",]
summary(wAll.div.cat)

wAll.div.24 <- wAll.div[wAll.div$coarseness=="Cat02" | wAll.div$coarseness=="Cat04",]
summary(wAll.div.24)

##number of tests with signif shifts from continuous traits (out of 16)

nrow(wAll.div.cat[wAll.div.cat$MNND.1 <= 0.05, ]) 
nrow(wAll.div.cat[wAll.div.cat$FDis.1 <= 0.05, ]) 
nrow(wAll.div.cat[wAll.div.cat$MPD.1 <= 0.05, ])  
nrow(wAll.div.cat[wAll.div.cat$Rao.1 <= 0.05, ])  

#proportions (/16)
nrow(wAll.div.cat[wAll.div.cat$MNND.1 <= 0.05, ])/16  
nrow(wAll.div.cat[wAll.div.cat$FDis.1 <= 0.05, ])/16  
nrow(wAll.div.cat[wAll.div.cat$MPD.1 <= 0.05, ])/16   
nrow(wAll.div.cat[wAll.div.cat$Rao.1 <= 0.05, ])/16   

#proportions (/20  ==  including rounded contin)
nrow(wAll.div[wAll.div$MNND.1 <= 0.05, ])/20  
nrow(wAll.div[wAll.div$FDis.1 <= 0.05, ])/20  
nrow(wAll.div[wAll.div$MPD.1 <= 0.05, ])/20   
nrow(wAll.div[wAll.div$Rao.1 <= 0.05, ])/20   

#proportions (/8  == just 2 and 4 state categ)
nrow(wAll.div.24[wAll.div.24$MNND.1 <= 0.05, ])/8  
nrow(wAll.div.24[wAll.div.24$FDis.1 <= 0.05, ])/8 
nrow(wAll.div.24[wAll.div.24$MPD.1 <= 0.05, ])/8   
nrow(wAll.div.24[wAll.div.24$Rao.1 <= 0.05, ])/8   


##div and conv combined
#proportions (/40 == all resolution levels)
nrow(wAll[wAll$MNND.1 <= 0.05, ])/40  
nrow(wAll[wAll$FDis.1 <= 0.05, ])/40  
nrow(wAll[wAll$MPD.1 <= 0.05, ])/40   
nrow(wAll[wAll$Rao.1 <= 0.05, ])/40   

#proportions (/16  == just 2 and 4 state categ)
wAll.24 <- wAll[wAll$coarseness=="Cat02" | wAll$coarseness=="Cat04",]
summary(wAll.24)

nrow(wAll.24[wAll.24$MNND.1 <= 0.05, ])/16  
nrow(wAll.24[wAll.24$FDis.1 <= 0.05, ])/16
nrow(wAll.24[wAll.24$MPD.1 <= 0.05, ])/16   
nrow(wAll.24[wAll.24$Rao.1 <= 0.05, ])/16   

###################################################################################################

####ANOVA####    
#**do this with the Wilcoxon estimates (median shifts). 
#from Wilcoxon test documentation: Note that in the two-sample case the estimator for the difference in location parameters does not estimate the difference in medians (a common misconception) but rather the median of the difference between a sample from x and a sample from y.

#changing Inf to NA, to allow ANOVA to proceed (only occurs for Convergence of Cat02 MNND w estimates)
wAll.conv[wAll.conv==Inf]<-NA
wAll.conv

###Testing Median shift estimates (wilcoxon test outputs (w))

#convergence results
aov.mnnd.conv.w <- aov(MNND ~ (coarseness + p)^2, data=wAll.conv) 
anova(aov.mnnd.conv.w)

aov.fdis.conv.w <- aov(FDis ~ (coarseness + p)^2, data=wAll.conv)
anova(aov.fdis.conv.w)

aov.mpd.conv.w <- aov(MPD ~ (coarseness + p)^2, data=wAll.conv)
anova(aov.mpd.conv.w)

aov.rao.conv.w <- aov(Rao ~ (coarseness + p)^2, data=wAll.conv)
anova(aov.rao.conv.w)


#divergence results
aov.mnnd.div.w <- aov(MNND ~ (coarseness + p)^2, data=wAll.div)
anova(aov.mnnd.div.w)

aov.fdis.div.w <- aov(FDis ~ (coarseness + p)^2, data=wAll.div)
anova(aov.fdis.div.w)

aov.mpd.div.w <- aov(MPD ~ (coarseness + p)^2, data=wAll.div)
anova(aov.mpd.div.w)

aov.rao.div.w <- aov(Rao ~ (coarseness + p)^2, data=wAll.div)
anova(aov.rao.div.w)

###############################################################


###############
#####summarizing regression results - w - for each index and process
####


#MNND
mnnd.conv.w <- as.data.frame(anova(aov.mnnd.conv.w))
mnnd.conv.w$index <- "MNND"
mnnd.conv.w$proc <- "Conv"
mnnd.conv.w$ConvRelSS <- mnnd.conv.w$`Sum Sq`/(sum(mnnd.conv.w$`Sum Sq`))
mnnd.conv.w

mnnd.div.w <- as.data.frame(anova(aov.mnnd.div.w))
mnnd.div.w$index <- "MNND"
mnnd.div.w$proc <- "Div"
mnnd.div.w$DivRelSS <- mnnd.div.w$`Sum Sq`/(sum(mnnd.div.w$`Sum Sq`))
mnnd.div.w

RegResults.mnnd <- cbind.data.frame(rownames(mnnd.conv.w), mnnd.conv.w$Df, mnnd.conv.w$ConvRelSS, mnnd.div.w$DivRelSS)
RegResults.mnnd
colnames(RegResults.mnnd) <- c("Variable", "Df", "SES_MNND_Conv", "SES_MNND_Div")
RegResults.mnnd


#FDis
fdis.conv.w <- as.data.frame(anova(aov.fdis.conv.w))
fdis.conv.w$index <- "fdis"
fdis.conv.w$proc <- "Conv"
fdis.conv.w$ConvRelSS <- fdis.conv.w$`Sum Sq`/(sum(fdis.conv.w$`Sum Sq`))
fdis.conv.w

fdis.div.w <- as.data.frame(anova(aov.fdis.div.w))
fdis.div.w$index <- "fdis"
fdis.div.w$proc <- "Div"
fdis.div.w$DivRelSS <- fdis.div.w$`Sum Sq`/(sum(fdis.div.w$`Sum Sq`))
fdis.div.w

RegResults.fdis <- cbind.data.frame(rownames(fdis.conv.w), fdis.conv.w$Df, fdis.conv.w$ConvRelSS, fdis.div.w$DivRelSS)
RegResults.fdis
colnames(RegResults.fdis) <- c("Variable", "Df", "SES_FDis_Conv", "SES_FDis_Div")
RegResults.fdis


#MPD
mpd.conv.w <- as.data.frame(anova(aov.mpd.conv.w))
mpd.conv.w$index <- "mpd"
mpd.conv.w$proc <- "Conv"
mpd.conv.w$ConvRelSS <- mpd.conv.w$`Sum Sq`/(sum(mpd.conv.w$`Sum Sq`))
mpd.conv.w

mpd.div.w <- as.data.frame(anova(aov.mpd.div.w))
mpd.div.w$index <- "mpd"
mpd.div.w$proc <- "Div"
mpd.div.w$DivRelSS <- mpd.div.w$`Sum Sq`/(sum(mpd.div.w$`Sum Sq`))
mpd.div.w

RegResults.mpd <- cbind.data.frame(rownames(mpd.conv.w), mpd.conv.w$Df, mpd.conv.w$ConvRelSS, mpd.div.w$DivRelSS)
RegResults.mpd
colnames(RegResults.mpd) <- c("Variable", "Df", "SES_mpd_Conv", "SES_mpd_Div")
RegResults.mpd


#Rao         
rao.conv.w <- as.data.frame(anova(aov.rao.conv.w))
rao.conv.w$index <- "rao"
rao.conv.w$proc <- "Conv"
rao.conv.w$ConvRelSS <- rao.conv.w$`Sum Sq`/(sum(rao.conv.w$`Sum Sq`))
rao.conv.w

rao.div.w <- as.data.frame(anova(aov.rao.div.w))
rao.div.w$index <- "rao"
rao.div.w$proc <- "Div"
rao.div.w$DivRelSS <- rao.div.w$`Sum Sq`/(sum(rao.div.w$`Sum Sq`))
rao.div.w

RegResults.rao <- cbind.data.frame(rownames(rao.conv.w), rao.conv.w$Df, rao.conv.w$ConvRelSS, rao.div.w$DivRelSS)
RegResults.rao
colnames(RegResults.rao) <- c("Variable", "Df", "SES_rao_Conv", "SES_rao_Div")
RegResults.rao


###COMPILE to a single table of relative sum of squares
RegResults.All.w <- cbind(RegResults.mnnd, RegResults.fdis[,3:4], RegResults.mpd[,3:4], RegResults.rao[,3:4])

RegResults.All.w ###This is Table S1 reproduced in the Supplementary Materials of our paper (prior to reformatting columns, row names, etc)

################################################################################











###############################################################################

######
###PLOTTING
#####

###For plotting, convert all to factor except index SES values
mAll_plots <- mAll
str(mAll_plots)
mAll_plots$p <- as.factor(mAll_plots$p)
str(mAll_plots)
summary(mAll_plots)

#create more informative facet labels
facetlabProc <- c(Conv = "Convergence", Div = "Divergence")
facetlabp <-  c("50" = "50 sp.", "100" = "100 sp.", "500" = "500 sp.", "1000" = "1000 sp.")
facetlabInd <- c("MNND" = "MNND", "MPD" = "MPD", "FDis" = "FDis", "Rao" = "Rao")


################################
# faceted line graphs   

###Re-parse the data with dplyr to have the metric name as a column.

#######Produces the overall summary figure in main text (Figure 4).


##
mAll_long <- gather(mAll_plots, "index", "value", 3:6) ###check that all cols are here
str(mAll_long)
head(mAll_long)

#make index a factor and order the way I want it displayed
mAll_long$index = factor(mAll_long$index, levels=c('MNND','MPD','FDis', 'Rao'))
str(mAll_long)

#create new column for splitting out the lines by concatenating
mAll_long$proc_ind <- paste(mAll_long$proc, mAll_long$index, sep="_")
mAll_long$proc_ind <- as.factor(mAll_long$proc_ind)
str(mAll_long)
summary(mAll_long)

mAll_long$proc_p <- paste(mAll_long$proc, mAll_long$p, sep="_")
mAll_long$proc_p <- as.factor(mAll_long$proc_p)
str(mAll_long)
summary(mAll_long)

mAll_long$proc_coarse <- paste(mAll_long$proc, mAll_long$coarseness, sep="_")
mAll_long$proc_coarse <- as.factor(mAll_long$proc_coarse)
str(mAll_long)
summary(mAll_long)

mAll_long$ind_coarse <- paste(mAll_long$index, mAll_long$coarseness, sep="_")
mAll_long$ind_coarse <- as.factor(mAll_long$ind_coarse)
levels(mAll_long$ind_coarse) <- c('MNND_TrContin', 'MNND_Rounded', 'MNND_Cat16', 'MNND_Cat08', 'MNND_Cat04', 'MNND_Cat02', 'FDis_TrContin', 'FDis_Rounded', 'FDis_Cat16', 'FDis_Cat08', 'FDis_Cat04', 'FDis_Cat02', 'MPD_TrContin', 'MPD_Rounded', 'MPD_Cat16', 'MPD_Cat08', 'MPD_Cat04', 'MPD_Cat02', 'Rao_TrContin', 'Rao_Rounded', 'Rao_Cat16', 'Rao_Cat08', 'Rao_Cat04', 'Rao_Cat02')
str(mAll_long)
summary(mAll_long)

##########################

#Main summary of simulation results

SummPlotAll <- ggplot(mAll_long, aes(x=rev(coarseness), y=value, group=proc_p, color=proc, fill=coarseness)) + #####change groupings etc
  geom_line(aes(linetype=p)) +
  scale_linetype_discrete(name='Pool size') + #change to pool size?
  
  scale_x_discrete(labels=c("Con", "ConR", "16", "8", "4", "2"), limits = rev(levels(mAll_long$coarseness))) +
  #scale_size_manual(values=c(.25,1.5)) +
  scale_color_manual(values = c("#D55E00", "#0072B2", "gray"), labels=c('Convergent', 'Divergent'), guide=guide_legend(reverse=T), name='Pattern') +
  geom_point(shape=21, size=2) +
  expand_limits(y=0) +
  
facet_wrap( ~index, ncol = 4, scales = "fixed") +
  ##to unlink row scales:   , scales = "free_y"
  scale_color_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill", labels=c("2-state categorical", "4-state categorical", "8-state categorical", "16-state categorical", "Rounded (ConR)", "Continuous (Con)"), guide=guide_legend(reverse = T), name= 'Trait resolution') +
  
  # scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  
  
  ###AXES###
  #axis labels
  ylab("Functional structure index (SES)") +
  xlab("Trait resolution") +
  
  
  #axis scale adjustments to make axes always uniform
  #scale_y_continuous(limits = c(-25, NA)) +
  #coord_cartesian(ylim= c(-20,20)) +
 
  
  ###BACKGROUND and ADD-INS###
  #removes background but keeps the border lines
  theme_bw() + #simple, neat, black and white theme
  
  #removes gridlines, but keeps border
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #adjust axis text
  
  #formatting axis titles: 
  theme(axis.title.x= element_text(face="bold", colour="black", size=16),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=60,hjust=1, vjust=1, size=14,  colour="black")) + #adjusts size and color of y axis title
  theme(axis.title.y = element_text(face="bold", colour="black", size=16), #adjusts angles, size, etc of labels for y axis
        axis.text.y  = element_text(angle=0, size=14, colour="black")) +
  
  #adjust facet chars
  theme(strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14, angle = 0, hjust = 0),
        strip.background = element_rect(colour=NA, fill=NA)) +
  
  #superimpose line for Zero
  geom_hline(aes(yintercept=0), color="gray", linetype="dashed", size=.5) +
  
  #legend
  theme (legend.position = "right")


plot(SummPlotAll)

########################################################################
