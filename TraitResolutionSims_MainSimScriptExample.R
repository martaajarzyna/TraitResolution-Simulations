##R script for Trait resolution simulations - Kohli and Jarzyna 2021, Global Ecology and Biogeography
#Simulations via Botta-Dukat and Czucz 2016 process-based assembly algorithms
#Created by: BA Kohli 
#last modified: 2 February 2021

##############################################################################
#load required packages

require(dplyr)
require(FD)
require(picante)
require(tidyr)
require(ggplot2)
require(cowplot)

####################################################


###Botta-Dukat and Czucz 2016 MEE Appendix S6
##modified 22 Oct 2020 By BAK to make the default sigma.b = 0.05 (was 0.03) to match the values used from the original publication.

##################################################
#                                                #
#           Community simulation                 #
#                                                #
##################################################
#
# Input parameters:
#
# S = number of species in the regional species pool
# m = probability of colonization from meta-community
# n = number of local communities
# J = number of individuals in a local community
# sigma = tolerance width (equal for all species)
#         has to be positive
#         lower values means more specialist species
#         sigma=Inf means that species are maximally generalist,
#               thus abiotic conditions
#                   do not influence their abundance
#         sigma=0 would mean that species are maximally specialist,
#                 they can occur at only one point of the 
#                 environmental gradient(s)
# sigma.b =  width of competition kernel
#            sigma.b=0 means no interspecific competition (no effect
#                       of trait B on competition)
#            sigma.b=Inf leads to equally strong inter- and 
#                        intraspecific competition 
# If both sigma and sigma.b equal to Inf, species are neutral,
# community composition influenced by random drift only
# b0 = probability of birth without competition 
# K = carrying capacity
# xrange =  the range of the environmental gradient, along which the 
#           simulated sites lie
# distrib = parameter influencing the shape of distribution. It
#           should be positive!
#             distrib<1 U-shaped distribution
#             distrib=1 uniform  distribution
#             distrib>1 bell-shaped distribution
# correl = correlation between traits, no ceoorelation if correl=0 
# (default)
# rand.seed = seed for random number generation
#               the default NULL initialize the random number 
#               genarator using current time
# sim.length = length of the simulation
#
#
#Output:
#
#  List of parameters +
#          Y = plot-by-species matrix of abundances
#          trait.env, trait.compet, trait.neutr = three vector of 
#          trait values

# S=200; m=0.1; n=50; J=300; sigma=0.05; sigma.b=0.03; b0=1;
# xrange=.8;distrib=1;
# correl=0;rand.seed=NULL; sim.length=100; K=200
# correl=-.6
# sigma = 0.05; sigma.b = 0.25; S = 200; J = 450; n = 50; distrib = 1; correl = 0; m = 0.1; b0 = 1; sim.length = 100; K = 200; xrange = 0.8; n.random = 1000; sig.level = 0.05
# rand.seed=NULL
traitsimul<-function(S=200, n=50, J=300, sigma=0.05, sigma.b=0.05, 
                     m=0.1, b0=1, xrange=.8, distrib=1, correl=0, 
                     rand.seed=NULL, sim.length=100, K=200, ...) 
{
  # the position of the sites along the environmental gradient:
  #  the xrange long central part of the gradient, sampled at equidistant points 
  x <- seq((1-xrange)/2,1-(1-xrange)/2,xrange/(n-1))
  set.seed(rand.seed)
  cat("Generating species pool... \n")
  trait.a <- rbeta(S, distrib, distrib) # values of trait1 (related to tolerance)
  tmp <- rbeta(S, distrib, distrib) # values of trait2 (related to resource use)
  trait.b <- switch(sign(correl)+2, 
                    abs(correl)*(1-trait.a)+(1-abs(correl))*tmp,  # correl: negative 
                    tmp,                                          # correl=0 (no correlation) 
                    correl*trait.a+(1-correl)*tmp)                # correl: posistive
  # correl values of +/-.6 result in cor(trait.a,trait.b) of ~.8 for distrib=1 using this algorithm
  trait.c <- rbeta(S, distrib, distrib) # values of trait3 (neutral)
  
  dist.b <- as.matrix(dist(trait.b))
  dist.a <- as.matrix(dist(trait.a))
  compet <- matrix(0,S,S)
  if (sigma.b==0) diag(compet) <- 1
  if (sigma.b==Inf) compet <- matrix(1,S,S)
  if ((sigma.b>0) & (sigma.b<Inf)) compet <- exp(-dist.b^2/sigma.b)
  
  
  Y<-matrix(NA,n,S) # species abundances
  off.spring<-vector()
  X<-matrix(rep(x,S),ncol=S)  #position along the gradient
  A<-t(matrix(rep(trait.a,n),ncol=n))
  
  survive <- if (sigma<Inf) pmax(exp(-((X-A)^2)/sigma)-0.01,0) else matrix(0.99,nrow(X),ncol(X))
  
  cat("Generating starting community composition...\n")
  for (i in 1:n) Y[i,]<-table(c(sample(1:S,J,replace=T,prob=survive[i,]),seq(1,S)))-1
  
  cat("Community assembly...\n")
  pb <- txtProgressBar (min = 0, max = sim.length, char = ".", width = 45, style = 3)
  
  # epoch=1; j=1
  for (epoch in 1:sim.length) {
    for (j in 1:J) {
      seed<-matrix(0,nrow=n,ncol=S)
      for (i in 1:n) {
        death<-sample(1:S,1,prob=Y[i,])
        Y[i,death] <- Y[i,death]-1
        NE <- compet %*% Y[i,]
        birth.limit <- b0*(K-NE)/K
        birth.limit[birth.limit<0] <- 0
        occurrence <- (Y[i,]>0)
        seed[i,occurrence] <- rbinom(sum(as.numeric(occurrence)),
                                     Y[i,occurrence],birth.limit[occurrence])
      }
      off.spring <- matrix(rbinom(n*S,size=seed,prob=(1-m)),nrow=n,ncol=S)
      seed <- seed-off.spring
      p <- matrix(1/(n-1),nrow=n,ncol=n)
      diag(p) <- 0
      
      for (i in 1:n)
        for (k in 1:S)
          if (seed[i,k]>0) off.spring[,k] <- off.spring[,k] + 
        rmultinom(1, size=seed[i,k], prob=p[,i])
      
      for (i in 1:n) {
        if (sum(off.spring[i,]*survive[i,])>0) {
          birth <- sample(1:S,1,prob=off.spring[i,]*survive[i,])
        } else {
          birth <- sample(1:S,1,prob=as.numeric(Y[i,]>0))
        }
        Y[i,birth] <- Y[i,birth]+1
      }
    }
    setTxtProgressBar(pb, epoch)  
  }
  res <- list(S=S, m=m, n=n, J=J, sigma=sigma, sigma.b=sigma.b, b0=b0, x=x, 
              distrib=distrib, correl=correl, rand.seed=rand.seed,
              sim.length=sim.length,K=K, Y=Y,trait.env=trait.a,
              trait.compet=trait.b,trait.neutr=trait.c)
  close(pb)
  return(res)
}


#################################################################
##############################################################



#############################################################################################
###De Bello etal 2016 Oecologia - 'melodic' function for correcting errors in how abundance-weighted MPD and Rao's Q are calculated in existing packages (picante and dbFD)####
###########################################################################################

##MEan DIssimilarity Components##
#samp:  community matrix; sites in lines, species in columns
#dis:   dissimilarity matrix
#type:  "both" for results with abundances weighted and non-weighted
#       "abundance" for results with abundances weighted
#       "presence" for results with abundances non-weighted

melodic <- function(samp,dis,type="both"){
  if(class(samp)!="matrix"){samp <- as.matrix(samp)}
  if(class(dis)!="matrix"){dis <- as.matrix(dis)}
  if(is.null(colnames(samp)) | is.null(colnames(dis)) ){
    stop("Both samp and dis must have colnames.\n")
  }
  N<-dim(samp)[1]
  melodic<-list()
  if (type=="both"){
    melodic$abundance<-list()
    melodic$abundance$mpd<-melodic$abundance$rao<-melodic$abundance$simpson<-numeric(N)
    melodic$presence<-list()
    melodic$presence$mpd<-melodic$presence$rao<-melodic$presence$simpson<-numeric(N)
  }
  if (type=="abundance"){ 
    melodic$abundance<-list()
    melodic$abundance$mpd<-melodic$abundance$rao<-melodic$abundance$simpson<-numeric(N)
  }
  if (type=="presence"){ 
    melodic$presence<-list()
    melodic$presence$mpd<-melodic$presence$rao<-melodic$presence$simpson<-numeric(N)
  }
  for (i in 1:N){
    sppInSample<-names(samp[i,samp[i, ]>0])
    melodic$richness[i]<-rowSums(samp>0)[i]
    if (length(sppInSample)>1){
      sample.dis<-dis[sppInSample,sppInSample]
      abund.w<-numeric(length(sppInSample))
      if (type=="both" | type=="abundance"){
        abund.w <- samp[i , sppInSample] / sum(samp[i , sppInSample])
        sample.weights <- outer(abund.w , abund.w)
        melodic$abundance$mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
        melodic$abundance$rao[i] <- sum(sample.weights * sample.dis)
        melodic$abundance$simpson[i] <- sum(2*sample.weights[lower.tri(sample.weights)])
      } 	
      if (type=="both" | type=="presence"){
        abund.nw <- rep(1 , length(sppInSample)) / length(sppInSample)
        sample.weights <- outer(abund.nw , abund.nw)
        melodic$presence$mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
        melodic$presence$rao[i] <- sum(sample.weights * sample.dis)
        melodic$presence$simpson[i] <- sum(2*sample.weights[lower.tri(sample.weights)])
      }	
    }	else {
      if (type=="both" | type=="abundance"){
        melodic$abundance$mpd[i] <- NA
        melodic$abundance$rao[i] <- melodic$abundance$simpson[i] <-0
      }
      if (type=="both" | type=="presence"){
        melodic$presence$mpd[i] <- NA
        melodic$presence$rao[i] <- melodic$presence$simpson[i] <-0
      }
    }
  }  	
  out<-melodic
  return(out)	
}
#########################################################################
#########################################################################

###########
##fdisp2 function - Collin VanBuren modified fdisp to allow it to proceed even when communities have no species or when a species does not occur in any of the communities
###this solution is needed because randomized communities were being constructed as such.
###########

fdisp2 <-   function (d, a, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (!is.matrix(a)) 
    stop("'a' must be a matrix.")
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  #if (any(abun.sum == 0)) 
  #  stop("At least one community has zero-sum abundances (no species).", 
  #       "\n")
  abun.sum2 <- apply(a, 2, sum)
  #if (any(abun.sum2 == 0)) 
  #  stop("At least one species does not occur in any community (zero total abundance across all communities).", 
  #       "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}
#########################################################################################
#########################################################################



##############################################################################
##BAK code


###OVERALL STEPS
#1. assemble 100 local communities from a regional pool with 1 trait each relating to environmental filtering, symmetric competition, and neutral processes

#2. coarsen traits at multiple benchmarks for entire pool

#3.for each local community calculate functional structure indices for each trait resolution level

#4.Compare each nonrandom community index value to the null distribution of the 1000 randomized communities values (SES). 

#5.Calculate the Skewness of distributions

####Then Repeat and change starting regional pool size 

##############################################################

##############################################################

###
##Part 1 : Simultaneously create regional pools and assemble local communities
###

#set pool size - choose from 50, 100, 500, 1000 
p <- 100 

################################################################################
#create process-based communities via 'traitsimul'
#vary regional pool size (S)
#

#regional pool (S) = 100, number of end communities (n) = 100, all else default values of Botta-Dukat and Czucz 2016.

#baseline scenario: environmental filtering along a gradient with competition
outcomm <- traitsimul(S=100, n=100) #sigma and sigma.b = 0.05 defaults, traits are drawn from uniform distributions and not correlated

Y <- as.data.frame(outcomm$Y) # plot-by-species matrix of abundances

#check matrices
rowSums(Y) #total abundance per site (J; all equal)
colSums(Y) #species with 0 don't occur in any local comm

#compile traits into single dataframe
Tr <- as.data.frame(cbind(outcomm$trait.env, outcomm$trait.compet, outcomm$trait.neutr)) #species-by-trait matrix
colnames(Tr) <- c("TrEF", "TrComp", "TrNeut")
row.names(Tr) <- colnames(Y)
Tr$sp <- row.names(Tr)


###
#remove species not occuring in any comms
ToRemove <- select_if(Y, colSums(Y) == 0)
tRemove <- as.data.frame(t(ToRemove)) #transpose to enable row matching
#row names are the species not found in any of the local communities after simulation
tRemove$sp <- row.names(tRemove)

Y2 <- select_if(Y, colSums(Y) != 0 ) #retain only present species in the abundance matrix

Tr2 <- anti_join(Tr, tRemove, by = "sp")#match the trait file accordingly, retain only present species


###############################################################################################


#################
##Part 2. Coarsen all traits at multiple resolutions
###########

###starting with the continuous traits for a given Pool size
continTr_100sp <- Tr2

##first level of coarsening: rounding (maintains fine structure; small loss of information)
rounded_100sp <- continTr_100sp %>% mutate_if(is.numeric, round, 2)

##successive levels of coarsening via categorical binning
##preliminary analyses revealed it unnessary to do 64 and 32-group categorizations
cuts <- c(16, 8, 4, 2) #using multiples of 2 maintains bin boundaries and a perfectly nested coarsening.
matrices <- vector("list", 4)
for (i in cuts) {matrices[[i]] <- continTr_100sp %>% mutate_if(is.numeric, cut, breaks=i, labels = F)} 
matrices

#save them separately
categ2_100sp <- matrices[[2]]
categ4_100sp <- matrices[[4]]
categ8_100sp <- matrices[[8]]
categ16_100sp <- matrices[[16]]

#row names
row.names(rounded_100sp)<-row.names(continTr_100sp)
row.names(categ16_100sp)<-row.names(continTr_100sp)
row.names(categ8_100sp)<-row.names(continTr_100sp)
row.names(categ4_100sp)<-row.names(continTr_100sp)
row.names(categ2_100sp)<-row.names(continTr_100sp)

##################################

#####LISTS OF FINAL MATRICES - from highest to lowest resolution#####

#combine all resolution matrices into a list per Reg Pool Size

ResMatrices100sp <- list(continTr_100sp, rounded_100sp, categ16_100sp, categ8_100sp, categ4_100sp, categ2_100sp)

############################################################################

########
## Parts 3 and 4. Calculate functional structure indices and SES values at each trait resolution level 
############

##Create species distance matrices via gower distance (Pavoine 2009) so that the same dissimilarity measure is used for all iterations.

#isolate each trait - functional distance matrix with just the single associated trait to the process

#functions to create distance matrix from species pool trait matrices
distmat.ef <- function(trts){trts2 <- as.data.frame(trts[1])
Dist.G.ef <-gowdis(trts2) ###create gower dissim matrix
#DistMat.continG.ef <- as.matrix(Dist.G.ef)  #convert to matrix form
}

distmat.comp <- function(trts){trts2 <- as.data.frame(trts[2])
Dist.G.comp <-gowdis(trts2) ###create gower dissim matrix
#DistMat.continG.comp <- as.matrix(Dist.G.comp)  #convert to matrix form
}

####################

#Do for each sp X trait pool per resolution - loop through trait resolution sets to create each distance matrix
AllDists.ef <- lapply(ResMatrices100sp, distmat.ef)
AllDists.comp <- lapply(ResMatrices100sp, distmat.comp)

###################################################################

#convert Species abundances to matrix for use in some functions below
Y2.matrix <- as.matrix(Y2)

#####################################################################


#######
#calculate indices
####

###FDis###

####EF/convergence sub-function
FDis.func.ef <- function(x){
  fdisp(x, randomizeMatrix(Y2.matrix, null.model = "frequency"))$FDis  #null model choice based on Gotzenberger etal 2016
}

####Competition/divergence sub-function
FDis.func.comp <- function(x){
  fdisp2(x, randomizeMatrix(Y2.matrix, null.model = "richness"))$FDis  #null model choice based on Gotzenberger etal 2016
}

#the main functions for calculating observed and null index values
###EF
SES.FDis.func.ef <- function(x){
  obs.null.output <- cbind(fdisp(x, Y2.matrix)$FDis, replicate(999, FDis.func.ef(x)))
  obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
  ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / apply(obs.null.output, 1, sd)
  p.value <- apply(cbind(obs.null.output[,1], obs.null.output), MARGIN = 1, rank)[1,] / 1000
  ses.fdis.out <- cbind.data.frame(obs.null.output[,1], ses.value, obs.rank, p.value)
  
}

###Comp
SES.FDis.func.comp <- function(x){
  obs.null.output <- cbind(fdisp2(x, Y2.matrix)$FDis, replicate(999, FDis.func.comp(x)))
  obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
  ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / apply(obs.null.output, 1, sd)
  p.value <- apply(cbind(obs.null.output[,1], obs.null.output), MARGIN = 1, rank)[1,] / 1000
  ses.fdis.out <- cbind.data.frame(obs.null.output[,1], ses.value, obs.rank, p.value)
  
}

#summarize observed and SES values
###EF/conv
FDis.ef.list <- lapply(AllDists.ef, SES.FDis.func.ef) #with post processing and SES calcs
FDis.ef <- data.frame(t(matrix(unlist(FDis.ef.list), nrow=length(FDis.ef.list), byrow=T)))
colnames(FDis.ef) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (ses.value column)
SES.FDis.ef <- FDis.ef[101:200,]
SES.FDis.ef$proc <- c(rep("Conv", 100)) #add a column to ID the assembly process used to create communities (covergence, divergence)

###Comp/Div
FDis.comp.list <- lapply(AllDists.comp, SES.FDis.func.comp) #with post processing and SES calcs
FDis.comp <- data.frame(t(matrix(unlist(FDis.comp.list), nrow=length(FDis.comp.list), byrow=T)))
colnames(FDis.comp) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (ses.value column)
SES.FDis.comp <- FDis.comp[101:200,]
SES.FDis.comp$proc <- c(rep("Div", 100)) #add a column to ID the assembly process used to create communities (covergence, divergence)

#join the ses outputs together
SES.FDis.ALL <- full_join(SES.FDis.ef, SES.FDis.comp,  all = TRUE)


##########################


###Rao's Q, quadratic entropy###          

#calculated with the melodic function of De Bello etal 2016

####EF/convergence sub-function
Rao.func.ef <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "frequency")),  type = "abundance")$abundance$rao  #null model choice based on Gotzenberger etal 2016
}

####Competition/divergence sub-function
Rao.func.comp <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "richness")),  type = "abundance")$abundance$rao  #null model choice based on Gotzenberger etal 2016
}

#the main functions for calculating observed and null index values
###EF
SES.Rao.func.ef <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$rao, replicate(999, Rao.func.ef(x)))
  obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
  ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / apply(obs.null.output, 1, sd)
  p.value <- apply(cbind(obs.null.output[,1], obs.null.output), MARGIN = 1, rank)[1,] / 1000
  mel.Rao.out <- cbind.data.frame(obs.null.output[,1], ses.value, obs.rank, p.value)
  
}

###Comp
SES.Rao.func.comp <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$rao, replicate(999, Rao.func.comp(x)))
  obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
  ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / apply(obs.null.output, 1, sd)
  p.value <- apply(cbind(obs.null.output[,1], obs.null.output), MARGIN = 1, rank)[1,] / 1000
  ses.Rao.out <- cbind.data.frame(obs.null.output[,1], ses.value, obs.rank, p.value)
  
}

#summarize observed and SES values
###EF/conv
Rao.ef.list <- lapply(AllDists.ef, SES.Rao.func.ef) #with post processing and SES calcs
Rao.ef <- data.frame(t(matrix(unlist(Rao.ef.list), nrow=length(Rao.ef.list), byrow=T)))
colnames(Rao.ef) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (ses.value column)
SES.Rao.ef <- Rao.ef[101:200,]
SES.Rao.ef$proc <- c(rep("Conv", 100)) #add a column to ID the assembly process used to create communities (covergence, divergence)

###Comp/Div
Rao.comp.list <- lapply(AllDists.comp, SES.Rao.func.comp) #with post processing and SES calcs
Rao.comp <- data.frame(t(matrix(unlist(Rao.comp.list), nrow=length(Rao.comp.list), byrow=T)))
colnames(Rao.comp) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (ses.value column)
SES.Rao.comp <- Rao.comp[101:200,]
SES.Rao.comp$proc <- c(rep("Div", 100)) #add a column to ID the assembly process used to create communities (covergence, divergence)

#join the ses outputs together
SES.Rao.ALL <- full_join(SES.Rao.ef, SES.Rao.comp,  all = TRUE)


##########################


###SES.MNND###

#observed and ses values
#EF/conv
MNTD.list.ef <- lapply(AllDists.ef, ses.mntd, samp=Y2, null.model = "frequency", abundance.weighted = T, runs = 999, iterations = 1000)
MNTD.ef <- data.frame(t(matrix(unlist(MNTD.list.ef), nrow=length(MNTD.list.ef), byrow=T)))
colnames(MNTD.ef) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (mntd.obs.z column)
SES.MNTD.ef <- MNTD.ef[501:600, ]
SES.MNTD.ef$proc <- c(rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)

#Comp/divergence
MNTD.list.comp <- lapply(AllDists.comp, ses.mntd, samp=Y2, null.model = "richness", abundance.weighted = T, runs = 999, iterations = 1000)
MNTD.comp <- data.frame(t(matrix(unlist(MNTD.list.comp), nrow=length(MNTD.list.comp), byrow=T)))
colnames(MNTD.comp) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (mntd.obs.z column)
SES.MNTD.comp <- MNTD.comp[501:600, ]
SES.MNTD.comp$proc <- c(rep("Div", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)

#join the ses outputs together
SES.MNTD.ALL <- full_join(SES.MNTD.ef, SES.MNTD.comp,  all = TRUE)


#######################

###SES.MPD###

#calculated with melodic function of de Bello etal 2016 (picante mpd actually calculates Raos Q, but with some incovenient peculiarities)

####EF/convergence sub-function
MPD.func.ef <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "frequency")),  type = "abundance")$abundance$mpd  #null model choice based on Gotzenberger etal 2016
}

####Competition/divergence sub-function
MPD.func.comp <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "richness")),  type = "abundance")$abundance$mpd  #null model choice based on Gotzenberger etal 2016
}

#the main functions for calculating observed and null index values
###EF
SES.MPD.func.ef <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$mpd, replicate(999, MPD.func.ef(x)))
  obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
  ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / apply(obs.null.output, 1, sd)
  p.value <- apply(cbind(obs.null.output[,1], obs.null.output), MARGIN = 1, rank)[1,] / 1000
  mel.MPD.out <- cbind.data.frame(obs.null.output[,1], ses.value, obs.rank, p.value)
  
}

###Comp
SES.MPD.func.comp <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$mpd, replicate(999, MPD.func.comp(x)))
  obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
  ses.value <- (obs.null.output[,1] - apply(obs.null.output, 1, mean)) / apply(obs.null.output, 1, sd)
  p.value <- apply(cbind(obs.null.output[,1], obs.null.output), MARGIN = 1, rank)[1,] / 1000
  ses.MPD.out <- cbind.data.frame(obs.null.output[,1], ses.value, obs.rank, p.value)
  
}

#summarize observed and SES values
###EF/conv
MPD.ef.list <- lapply(AllDists.ef, SES.MPD.func.ef) #with post processing and SES calcs
MPD.ef <- data.frame(t(matrix(unlist(MPD.ef.list), nrow=length(MPD.ef.list), byrow=T)))
colnames(MPD.ef) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (ses.value column)
SES.MPD.ef <- MPD.ef[101:200,]
SES.MPD.ef$proc <- c(rep("Conv", 100)) #add a column to ID the assembly process used to create communities (covergence, divergence)

###Comp/Div
MPD.comp.list <- lapply(AllDists.comp, SES.MPD.func.comp) #with post processing and SES calcs
MPD.comp <- data.frame(t(matrix(unlist(MPD.comp.list), nrow=length(MPD.comp.list), byrow=T)))
colnames(MPD.comp) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
#save out only the SES values (ses.value column)
SES.MPD.comp <- MPD.comp[101:200,]
SES.MPD.comp$proc <- c(rep("Div", 100)) #add a column to ID the assembly process used to create communities (covergence, divergence)

#join the ses outputs together
SES.MPD.ALL <- full_join(SES.MPD.ef, SES.MPD.comp,  all = TRUE)


########################################################################################
#to allow certain calculations and plotting to proceed normally, may need to ignore values of Inf for some indices
SES.MNTD.ALL[SES.MNTD.ALL==Inf]<-NA


#########
#PLOT
####

#switch from wide to long format 

#FDis
SES.FDis.ALL.long <- gather(SES.FDis.ALL, "coarseness", "value", 1:6)

#MNTD
SES.MNTD.ALL.long <- gather(SES.MNTD.ALL, "coarseness", "value", 1:6)

#Rao's Q
SES.Rao.ALL.long <- gather(SES.Rao.ALL, "coarseness", "value", 1:6)

#MPD
SES.MPD.ALL.long <- gather(SES.MPD.ALL, "coarseness", "value", 1:6)

#######################

#find the median of each group
med.MNTD <- group_by(SES.MNTD.ALL.long, proc, coarseness) %>% summarise(med=median(value, na.rm = T))

med.MPD <- group_by(SES.MPD.ALL.long, proc, coarseness) %>% summarise(med=median(value, na.rm = T))

med.Rao <- group_by(SES.Rao.ALL.long, proc, coarseness) %>% summarise(med=median(value, na.rm = T))

med.FDis <- group_by(SES.FDis.ALL.long, proc, coarseness) %>% summarise(med=median(value, na.rm = T))

#Combine ALL
medians.all <- cbind.data.frame(med.MNTD, med.MPD$med, med.Rao$med, med.FDis$med, rep(p, length(med.MNTD)))  
colnames(medians.all)<- c("proc", "coarseness", "MNND", "MPD", "Rao", "FDis", "p") 

#######################


###Clean plots

# Overlay with transparent density plot
F1 <- ggplot(SES.MNTD.ALL.long, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.MNTD, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
  scale_shape_manual(values = c(24, 24))+
  
  ###AXES###
  #axis labels
  xlab("SES-MNND") +
  ylab("Density") +
  
  #formatting axis titles: 
  theme(axis.title.x= element_text(face="bold", colour="black", size=22),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=0,size=14, colour="black")) + #adjusts size and color of y axis title
  theme(axis.title.y = element_text(face="bold", colour="black", size=22), #adjusts angles, size, etc of labels for y axis
        axis.text.y  = element_text(angle=0, vjust=0.5, size=14, colour="black")) +
  
  #axis scale adjustments to make axes always uniform
  scale_x_continuous(limits = c(-5, 5)) +
  

  ###BACKGROUND and ADD-INS###
  #removes background but keeps the border lines
  theme_bw() + #simple, neat, black and white theme
  
  #removes gridlines, but keeps border
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  
  #superimpose line for Zero
  geom_vline(aes(xintercept=0), color="gray", linetype="dashed", size=.5) +
  
  #legend
  theme (legend.position = "none")

#######

# Overlay with transparent density plot
F2<-ggplot(SES.Rao.ALL.long, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.Rao, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
  scale_shape_manual(values = c(24, 24))+
  
  ###AXES###
  #axis labels
  xlab("SES-Rao") +
  ylab("Density") +
  
  #formatting axis titles: 
  theme(axis.title.x= element_text(face="bold", colour="black", size=22),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=0,size=14, colour="black")) + #adjusts size and color of y axis title
  theme(axis.title.y = element_text(face="bold", colour="black", size=22), #adjusts angles, size, etc of labels for y axis
        axis.text.y  = element_text(angle=0, vjust=0.5, size=14, colour="black")) +
  
  #axis scale adjustments to make axes always uniform
  scale_x_continuous(limits = c(-5, 5)) +
  
  
  ###BACKGROUND and ADD-INS###
  #removes background but keeps the border lines
  theme_bw() + #simple, neat, black and white theme
  
  #removes gridlines, but keeps border
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  
  #superimpose line for Zero
  geom_vline(aes(xintercept=0), color="gray", linetype="dashed", size=.5) +
  
  #legend
  theme (legend.position = "none")

#######

# Overlay with transparent density plot
F3<-ggplot(SES.FDis.ALL.long, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.FDis, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
  scale_shape_manual(values = c(24, 24))+
  
  ###AXES###
  #axis labels
  xlab("SES-FDis") +
  ylab("Density") +
  
  #formatting axis titles: 
  theme(axis.title.x= element_text(face="bold", colour="black", size=22),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=0,size=14, colour="black")) + #adjusts size and color of y axis title
  theme(axis.title.y = element_text(face="bold", colour="black", size=22), #adjusts angles, size, etc of labels for y axis
        axis.text.y  = element_text(angle=0, vjust=0.5, size=14, colour="black")) +
  
  #axis scale adjustments to make axes always uniform
  scale_x_continuous(limits = c(-5, 5)) +
  
  
  ###BACKGROUND and ADD-INS###
  #removes background but keeps the border lines
  theme_bw() + #simple, neat, black and white theme
  
  #removes gridlines, but keeps border
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  
  #superimpose line for Zero
  geom_vline(aes(xintercept=0), color="gray", linetype="dashed", size=.5) +
  
  #legend
  theme (legend.position = "none")

#########################################

# Overlay with transparent density plot
F4 <- ggplot(SES.MPD.ALL.long, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.MPD, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
  scale_shape_manual(values = c(24, 24))+
  
  ###AXES###
  #axis labels
  xlab("SES-MPD") +
  ylab("Density") +
  
  #formatting axis titles: 
  theme(axis.title.x= element_text(face="bold", colour="black", size=22),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=0,size=14, colour="black")) + #adjusts size and color of y axis title
  theme(axis.title.y = element_text(face="bold", colour="black", size=22), #adjusts angles, size, etc of labels for y axis
        axis.text.y  = element_text(angle=0, vjust=0.5, size=14, colour="black")) +
  
  #axis scale adjustments to make axes always uniform
  scale_x_continuous(limits = c(-5, 5)) +

  
  ###BACKGROUND and ADD-INS###
  #removes background but keeps the border lines
  theme_bw() + #simple, neat, black and white theme
  
  #removes gridlines, but keeps border
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  
  #superimpose line for Zero
  geom_vline(aes(xintercept=0), color="gray", linetype="dashed", size=.5) +
  
  #legend
  theme (legend.position = "none")

#######
#print(MFD.obs.plot.R)
library(cowplot)
#desired order: [MNND, MPD, FDis, Rao], 
indexPlots <- plot_grid(F1, F4, F3, F2, ncol = 2, nrow = 2, labels = "auto") + draw_figure_label(label = SimRunVar, position = "bottom.left")
indexPlots
ponents


#######################################################################################



#######
#Wilcoxon Signed Rank Tests of differences between mean SES of distributions
####

#1.MNTD
##divergence
wt.d.rd <- wilcox.test(SES.MNTD.ALL$Rounded[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) 

wt.d.16 <-wilcox.test(SES.MNTD.ALL$Cat16[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.MNTD.ALL$Cat08[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.MNTD.ALL$Cat04[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.MNTD.ALL$Cat02[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.MNTD.ALL$Rounded[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) 

wt.c.16 <- wilcox.test(SES.MNTD.ALL$Cat16[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.MNTD.ALL$Cat08[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.MNTD.ALL$Cat04[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.MNTD.ALL$Cat02[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_MNND <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

####################

#2.Rao
##divergence
wt.d.rd <- wilcox.test(SES.Rao.ALL$Rounded[SES.Rao.ALL$proc=="Div"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) 

wt.d.16 <-wilcox.test(SES.Rao.ALL$Cat16[SES.Rao.ALL$proc=="Div"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.Rao.ALL$Cat08[SES.Rao.ALL$proc=="Div"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.Rao.ALL$Cat04[SES.Rao.ALL$proc=="Div"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.Rao.ALL$Cat02[SES.Rao.ALL$proc=="Div"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.Rao.ALL$Rounded[SES.Rao.ALL$proc=="Conv"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) 

wt.c.16 <- wilcox.test(SES.Rao.ALL$Cat16[SES.Rao.ALL$proc=="Conv"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.Rao.ALL$Cat08[SES.Rao.ALL$proc=="Conv"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.Rao.ALL$Cat04[SES.Rao.ALL$proc=="Conv"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.Rao.ALL$Cat02[SES.Rao.ALL$proc=="Conv"], SES.Rao.ALL$TrContin[SES.Rao.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_Rao <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

########################
#3.FDis
##divergence
wt.d.rd <- wilcox.test(SES.FDis.ALL$Rounded[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.16 <-wilcox.test(SES.FDis.ALL$Cat16[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.FDis.ALL$Cat08[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.FDis.ALL$Cat04[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.FDis.ALL$Cat02[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.FDis.ALL$Rounded[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) 

wt.c.16 <- wilcox.test(SES.FDis.ALL$Cat16[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.FDis.ALL$Cat08[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.FDis.ALL$Cat04[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.FDis.ALL$Cat02[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_FDis <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

######################################################
####
#4.MPD
##divergence
wt.d.rd <- wilcox.test(SES.MPD.ALL$Rounded[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) 

wt.d.16 <-wilcox.test(SES.MPD.ALL$Cat16[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.MPD.ALL$Cat08[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.MPD.ALL$Cat04[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.MPD.ALL$Cat02[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.MPD.ALL$Rounded[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.16 <- wilcox.test(SES.MPD.ALL$Cat16[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.MPD.ALL$Cat08[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.MPD.ALL$Cat04[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.MPD.ALL$Cat02[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_MPD <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

####################
####################



##########################

#Combine ALL
WilcoxTests.all <- cbind.data.frame(c(rep("Div", 10), rep("Conv", 10)), c(rep("Rounded", 2), rep("Cat16", 2), rep("Cat08",2), rep("Cat04",2), rep("Cat02",2),rep("Rounded", 2), rep("Cat16", 2), rep("Cat08",2), rep("Cat04",2), rep("Cat02",2)) ,t(WilcoxTests_MNND), t(WilcoxTests_FDis), t(WilcoxTests_MPD), t(WilcoxTests_Rao), rep(p, length(WilcoxTests_MNND)))

colnames(WilcoxTests.all)<- c("proc", "coarseness", "MNND", "FDis", "MPD", "Rao", "p")

###############################################################################


#######
##Part 5. : Calculating Skewness of Null SES distributions
#####

#######
####functions and code for summarizing the skew of null randomizations to inform the reliability of SES values (see Botta-Dukat 2018 for reasoning and explanation)
########

#Created by: BA Kohli 
#last modified: 5 November 2020

library(moments)

#calcs skew value excluding first column, which is the observed values.
skew.func <- function(z){
  mean.skew <- mean(apply(z[,2:1000], MARGIN = 1, skewness), na.rm = T)
  sd.skew <- sd(apply(z[,2:1000], MARGIN = 1, skewness), na.rm = T)
  skew.out <- cbind.data.frame(mean.skew, sd.skew)
}


###EF
####EF/convergence sub-function
Rao.func.ef <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "frequency")),  type = "abundance")$abundance$rao  #null model choice based on Gotzenberger etal 2016
}

SES.Rao.func.ef.skew <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$rao, replicate(999, Rao.func.ef(x)))
}

Rao.ef.list.skew <- lapply(AllDists.ef, SES.Rao.func.ef.skew) #run it
skew.list <- lapply(Rao.ef.list.skew, skew.func)
Rao.ef.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
colnames(Rao.ef.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
Rao.ef.skew$value<- c("mean", "sd")
Rao.ef.skew$p <- p
Rao.ef.skew$index <- "Rao"
Rao.ef.skew$process <- "Conv"

####Competition/divergence sub-function
Rao.func.comp <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "richness")),  type = "abundance")$abundance$rao  #null model choice based on Gotzenberger etal 2016
}

SES.Rao.func.comp.skew <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$rao, replicate(999, Rao.func.comp(x)))
}

Rao.comp.list.skew <- lapply(AllDists.comp, SES.Rao.func.comp.skew) #run it
skew.list <- lapply(Rao.comp.list.skew, skew.func)
skew.list
Rao.comp.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
head(Rao.comp.skew)
colnames(Rao.comp.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
Rao.comp.skew$value<- c("mean", "sd")
Rao.comp.skew$p <- p
Rao.comp.skew$index <- "Rao"
Rao.comp.skew$process <- "Div"
head(Rao.comp.skew)

#join the ses outputs together
Rao.ALL.skew <- rbind.data.frame(Rao.ef.skew, Rao.comp.skew)
Rao.ALL.skew
summary(Rao.ALL.skew)
########
####
#######


###FDis###

####EF/convergence sub-function
FDis.func.ef <- function(x){
  fdisp(x, randomizeMatrix(Y2.matrix, null.model = "frequency"))$FDis  #null model choice based on Gotzenberger etal 2016
}

####Competition/divergence sub-function
FDis.func.comp <- function(x){
  fdisp2(x, randomizeMatrix(Y2.matrix, null.model = "richness"))$FDis  #null model choice based on Gotzenberger etal 2016
}

#the main functions for calculating observed and null index values
###EF
SES.FDis.func.ef.skew <- function(x){
  obs.null.output <- cbind(fdisp(x, Y2.matrix)$FDis, replicate(999, FDis.func.ef(x)))
}

###Comp
SES.FDis.func.comp.skew <- function(x){
  obs.null.output <- cbind(fdisp2(x, Y2.matrix)$FDis, replicate(999, FDis.func.comp(x)))
}

FDis.ef.list.skew <- lapply(AllDists.ef, SES.FDis.func.ef.skew) #run it
skew.list <- lapply(FDis.ef.list.skew, skew.func)
skew.list
FDis.ef.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
head(FDis.ef.skew)
colnames(FDis.ef.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
FDis.ef.skew$value<- c("mean", "sd")
FDis.ef.skew$p <- p
FDis.ef.skew$index <- "FDis"
FDis.ef.skew$process <- "Conv"
head(FDis.ef.skew)


FDis.comp.list.skew <- lapply(AllDists.comp, SES.FDis.func.comp.skew) #run it
skew.list <- lapply(FDis.comp.list.skew, skew.func)
skew.list
FDis.comp.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
head(FDis.comp.skew)
colnames(FDis.comp.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
FDis.comp.skew$value<- c("mean", "sd")
FDis.comp.skew$p <- p
FDis.comp.skew$index <- "FDis"
FDis.comp.skew$process <- "Div"
head(FDis.comp.skew)

#join the ses outputs together
FDis.ALL.skew <- rbind.data.frame(FDis.ef.skew, FDis.comp.skew)
FDis.ALL.skew
summary(FDis.ALL.skew)
#######
#####
#######


###MPD

####EF/convergence sub-function
MPD.func.ef <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "frequency")),  type = "abundance")$abundance$mpd  #null model choice based on Gotzenberger etal 2016
}

####Competition/divergence sub-function
MPD.func.comp <- function(x){
  melodic(x, samp = (randomizeMatrix(Y2.matrix, null.model = "richness")),  type = "abundance")$abundance$mpd  #null model choice based on Gotzenberger etal 2016
}

#the main functions for calculating observed and null index values
###EF
SES.MPD.func.ef.skew <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$mpd, replicate(999, MPD.func.ef(x)))
}

###Comp
SES.MPD.func.comp.skew <- function(x){
  obs.null.output <- cbind(melodic(x, samp = Y2.matrix, type = "abundance")$abundance$mpd, replicate(999, MPD.func.comp(x))) 
}

MPD.ef.list.skew <- lapply(AllDists.ef, SES.MPD.func.ef.skew) #run it
skew.list <- lapply(MPD.ef.list.skew, skew.func)
skew.list
MPD.ef.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
head(MPD.ef.skew)
colnames(MPD.ef.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
MPD.ef.skew$value<- c("mean", "sd")
MPD.ef.skew$p <- p
MPD.ef.skew$index <- "MPD"
MPD.ef.skew$process <- "Conv"
head(MPD.ef.skew)


MPD.comp.list.skew <- lapply(AllDists.comp, SES.MPD.func.comp.skew) #run it
skew.list <- lapply(MPD.comp.list.skew, skew.func)
skew.list
MPD.comp.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
head(MPD.comp.skew)
colnames(MPD.comp.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
MPD.comp.skew$value<- c("mean", "sd")
MPD.comp.skew$p <- p
MPD.comp.skew$index <- "MPD"
MPD.comp.skew$process <- "Div"
head(MPD.comp.skew)

#join the ses outputs together
MPD.ALL.skew <- rbind.data.frame(MPD.ef.skew, MPD.comp.skew)
MPD.ALL.skew
summary(MPD.ALL.skew)


#########
#####
#######

####MNTD 

###requires matrix form of distances for this index

#functions to create distance matrix from species pool trait matrices
distmat.ef <- function(trts){trts2 <- as.data.frame(trts[1])
Dist.G.ef <-gowdis(trts2) ###create gower dissim matrix
DistMat.continG.ef <- as.matrix(Dist.G.ef)  #convert to matrix form
}

distmat.comp <- function(trts){trts2 <- as.data.frame(trts[2])
Dist.G.comp <-gowdis(trts2) ###create gower dissim matrix
DistMat.continG.comp <- as.matrix(Dist.G.comp)  #convert to matrix form
}

####################

#Do for each sp X trait pool per resolution

#loop through trait resolution sets to create each distance matrix
AllDists.ef <- lapply(ResMatrices100sp, distmat.ef)
AllDists.comp <- lapply(ResMatrices100sp, distmat.comp)

#######

###EF
MNTD.func.ef <- function(x){
  mntd(x, samp = (randomizeMatrix(Y2.matrix, null.model = "frequency")), abundance.weighted = T)  #null model choice based on Gotzenberger etal 2016
}
SES.MNTD.func.ef.skew <- function(x){
  obs.null.output <- cbind(mntd(x, samp = Y2.matrix, abundance.weighted = T), replicate(999, MNTD.func.ef(x)))
}

####COMP
MNTD.func.comp <- function(x){
  mntd(x, samp = (randomizeMatrix(Y2.matrix, null.model = "richness")), abundance.weighted = T)  #null model choice based on Gotzenberger etal 2016
}
SES.MNTD.func.comp.skew <- function(x){
  obs.null.output <- cbind(mntd(x, samp = Y2.matrix, abundance.weighted = T), replicate(999, MNTD.func.comp(x)))
}

MNTD.ef.list.skew <- lapply(AllDists.ef, SES.MNTD.func.ef.skew) #run it
skew.list <- lapply(MNTD.ef.list.skew, skew.func)
skew.list
MNTD.ef.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
head(MNTD.ef.skew)
colnames(MNTD.ef.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
MNTD.ef.skew$value<- c("mean", "sd")
MNTD.ef.skew$p <- p
MNTD.ef.skew$index <- "MNTD"
MNTD.ef.skew$process <- "Conv"
head(MNTD.ef.skew)


MNTD.comp.list.skew <- lapply(AllDists.comp, SES.MNTD.func.comp.skew) #run it
skew.list <- lapply(MNTD.comp.list.skew, skew.func)
skew.list
MNTD.comp.skew <- data.frame(t(matrix(unlist(skew.list), nrow=length(skew.list), byrow=T)))
head(MNTD.comp.skew)
colnames(MNTD.comp.skew) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
MNTD.comp.skew$value<- c("mean", "sd")
MNTD.comp.skew$p <- p
MNTD.comp.skew$index <- "MNTD"
MNTD.comp.skew$process <- "Div"

head(MNTD.comp.skew)

#join the ses outputs together
MNTD.ALL.skew <- rbind.data.frame(MNTD.ef.skew, MNTD.comp.skew)
MNTD.ALL.skew
summary(MNTD.ALL.skew)  ####NOTE that for MNTD many of the SES values are NA, so the mean and sd of skew are strongly affected.  Most dramatic for Cat02, which does not invalidate the use of SES more broadly.

R1_Skew_p100 <- rbind.data.frame(MNTD.ALL.skew, MPD.ALL.skew, FDis.ALL.skew, Rao.ALL.skew)
R1_Skew_p100 <- R1_Skew_p100[with(R1_Skew_p100, order(value, index, process)), ]
R1_Skew_p100
####################################################################################