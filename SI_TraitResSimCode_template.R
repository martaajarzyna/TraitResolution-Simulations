##R script for Trait resolution simulations - Code for generating pools and calculating functional structure indices following the Supplementary Methods and Results described in Supplemental Information S1 and S2 from Kohli and Jarzyna 2021 Global Ecology and Biogeography.
#Created by: BA Kohli
#on: 18 Nov 2019
#last modified: 2 July 2020


####################################################################################################
#load required packages

require(dplyr)
require(FD)
require(picante)
require(vegan)
require(tidyr)
require(ggplot2)
require(psych)
require(cowplot)
require(ape)

####################################################

###############


#######################
###1. create universal regional species pool continuous traits (this ensures using the same traits and sp across runs)
############

#load simul.comms function of McPherson

#R code for simul.comms
#Created by: Mcpherson, Yeager and Baum 2018 Methods EE - Simulation tool for scrutinizing FD metrics
#downloaded by BAK: 11 Sept 2019
#last modified: Never


####################################################################################################
#load required packages
require(scales)
require(moments)

# ===== START OF R CODE =====
simul.comms<-function (s=c(5, 10, 15, 20, 25, 30, 35, 40), r=10, p=100,
                       t=3, tr.method="norm:1", tr.type="con", presence="random", abun.method="lnorm:1",
                       abundance="random", commin=NULL, commax=NULL, comnot=NULL,
                       comnot.type=NULL, dups=0) {
  # increase memory limit (4095Mb for 32-bit R; values up to 8Tb possible
  # for 64-bit R on 64-bit Windows)
  memory.limit(10000)
  # Load required packages and functions
  require(scales)
  require(moments)
  # function to fill list of abundance vectors based on user specifications
  fill.abun<-function(empty, comsize, presprob, abunprob, abundist, abunsd, ids){
    nn<-length(empty)/2
    set<-sample(1:length(empty), size=(comsize-ids), replace=FALSE, prob=presprob)
    if(ids>0){set<-c(set, nn+sample(set, size=dups, replace=FALSE, prob=NULL))}
    sorted<-set[order(abunprob[set])]
    if(abundist == "lnorm"){counts<-0.5+rlnorm(n=length(set), meanlog=0, sdlog=abunsd)}
    if(abundist == "norm"){counts<-0.5+rnorm(n=length(set), mean=0, sd=abunsd)}
    if(abundist == "unif"){counts<-0.5+runif(n=length(set), min=0, max=abunsd)}
    if(abundist == "fixed"){counts<-rep(abunsd, times=length(set))}
    empty[sorted]<-sort(counts)
    return(empty)
  } # end function fill.abun
  # function to ensure positive abundance values (strictly > 0)
  pos.abun<-function(vec){
    vmin<-min(vec)
    pres<-which(vec!=0)
    pmin<-min(abs(vec[pres]))
    nvec<-numeric(length=length(vec))
    nvec[pres]<-vec[pres]+abs(vmin)+pmin
    return(nvec)
  } # end function pos.abun
  # Error checks
  if(p < max(s)){
    stop("'p' must be greater than the maximum species richness level in 's'.")
  }
  if(length(tr.method)!=t && length(tr.method)>1){
    stop(paste("In tr.method supply either a single distribution to be used for all traits, ",
               "or as many distributions as there are traits in t.", sep=" "))
  }
  if(length(tr.type)!=t && length(tr.type)>1){
    stop(paste("In tr.type supply either a single type to be used for all traits, ",
               "or as many types as there are traits in t.", sep=" "))
  }
  if(!is.null(commin) && length(commin)<t){
    stop(paste("Vector commin must contain", t, "values.", sep=" "))
  }
  if(!is.null(commax) && length(commax)<t){
    stop(paste("Vector commax must contain", t, "values.", sep=" "))
  }
  if(!is.null(comnot) && length(comnot)<t){
    stop(paste("Option 'comnot' must be a list that contains", t, "elements.", sep=" "))
  }
  if(!is.null(comnot.type) && length(comnot.type)<t){
    stop(paste("Option 'comnot.type' must be a list that contains", t, "elements.", sep=" "))
  }
  if(dups>(0.5*min(s))){
    stop("The value in dups cannot exceed 50% of the minimum value specified in s!")
  }
  # assignments
  l.s <- length(s)
  r.rep <- rep(r, l.s)
  nb.sp <- rep(s, r.rep)
  if(is.null(presence)){presence<-"random"}
  presence<-match.arg(presence, c("random", "common", "rare"))
  if(is.null(abundance)){abundance<-"random"}
  abundance<-match.arg(abundance, c("random", "common", "rare"))
  if(is.null(commin)){commin<-rep(NA,t)}
  if(is.null(commax)){commax<-rep(NA,t)}
  if(is.null(comnot)){comnot<-as.list(rep(NA,t))}
  if(is.null(comnot.type)){comnot.type<-as.list(rep("single",t))}
  if(is.null(dups)){dups<-0}
  # check if trait distributions specified include standard deviations
  if(is.null(tr.method)){tr.method<-"norm"}
  stopper<-regexpr(":", tr.method, fixed=TRUE)
  tr.dist<-ifelse(stopper==-1, tr.method, substr(tr.method, start=0, stop=stopper-1))
  tr.dist <- match.arg(tr.dist, c("unif", "norm", "lnorm"), several.ok=TRUE)
  tr.sd<-as.numeric(ifelse(stopper==-1, 1, substr(tr.method, start=stopper+1,
                                                  stop=nchar(tr.method))))
  # check if trait type specified indicates number of levels desired for categorical
  # or ordinal variables
  if(is.null(tr.type)){tr.type<-"con"}
  var.type<-tr.type
  stopper<-regexpr(":", tr.type, fixed=TRUE)
  tr.type<-ifelse(stopper==-1, tr.type, substr(tr.type, start=0, stop=stopper-1))
  tr.type<-match.arg(tr.type, c("con", "ord", "cat"), several.ok=TRUE)
  ncat<-as.numeric(ifelse(stopper==-1, 10, substr(var.type, start=stopper+1,
                                                  stop=nchar(var.type))))
  if(length(ncat)==1){ncat<-rep(ncat, t)}
  # check if the abundance distribution specified includes a standard deviation
  if(is.null(abun.method)){abun.method<-"lnorm"}
  stopper<-regexpr(":", abun.method, fixed=TRUE)
  abun.dist<-ifelse(stopper==-1, abun.method, substr(abun.method, start=0, stop=stopper-1))
  abun.dist<-match.arg(abun.dist, c("unif", "norm", "lnorm", "fixed"), several.ok=FALSE)
  abun.sd<-as.numeric(ifelse(stopper==-1, 1, substr(abun.method, start=stopper+1,
                                                    stop=nchar(abun.method))))
  # Create a traits matrix for the p species in the species pool with t traits each
  # Amended from original simulation function in FD to allow different distributions
  # for each trait, including different standard deviations.
  # Also keep track of the likelihood (density) of trait values in traits.d
  traits <- matrix(NA, p, t)
  traits.d<-traits
  if(length(tr.dist)>1){
    for(i in 1:t){
      if(tr.dist[i]=="unif"){
        traits[,i]<-runif(p, min=0, max=tr.sd[i])
        traits.d[,i]<-dunif(traits[,i], min=0, max=tr.sd[i])
      }
      if(tr.dist[i]=="norm"){
        traits[,i]<-rnorm(p, mean=0, sd=tr.sd[i])
        traits.d[,i]<-dnorm(traits[,i], mean=0, sd=tr.sd[i])
      }
      if(tr.dist[i]=="lnorm"){
        traits[,i]<-rlnorm(p, meanlog=0, sdlog=tr.sd[i])
        traits.d[,i]<-dlnorm(traits[,i], meanlog=0, sdlog=tr.sd[i])
      }
    } # end for i
  } # end if
  if(length(tr.dist)==1 && tr.dist == "unif"){
    traits <- apply(traits, 2, function(n, min, max) runif(n=p, min=0, max=tr.sd))
    traits.d<-apply(traits, 2, dunif, min=0, max=tr.sd)
    tr.dist<-rep(tr.dist,t)
    tr.sd<-rep(tr.sd,t)
  }
  if(length(tr.dist)==1 && tr.dist == "norm") {
    traits <- apply(traits, 2, function(n,mean,sd) rnorm(n=p, mean=0, sd=tr.sd))
    traits.d<-apply(traits, 2, dnorm, mean=0, sd=tr.sd)
    tr.dist<-rep(tr.dist,t)
    tr.sd<-rep(tr.sd,t)
  }
  if(length(tr.dist)==1 && tr.dist == "lnorm"){
    traits <- apply(traits, 2, function(n,meanlog,sdlog) rlnorm(n=p, meanlog=0, sdlog=tr.sd))
    traits.d<-apply(traits, 2, dlnorm, meanlog=0, sdlog=tr.sd)
    tr.dist<-rep(tr.dist,t)
    tr.sd<-rep(tr.sd,t)
  }
  # Now convert traits to user-specified types. Needs to be done column by
  # column even when all traits are of the same type in case different columns
  # use different distributions or standard deviations.
  # For categorical and ordinal variables, compute mean of trait likelihood
  # values per categorical or ordinal value.
  traits.final<-as.data.frame(traits)
  traits.d.final<-as.data.frame(traits.d)
  if(length(tr.type)==1){tr.type<-rep(tr.type, t)}
  for(i in 1:t){
    #i<-2 # for testing
    if(tr.type[i]!="con"){
      if(tr.dist[i]!="unif"){
        # categorize
        #cuts<-c(-4, -3, -2, -1, 0, 1, 2, 3, 4)*tr.sd[i]
        # this fixes Number of categories at 10
        cuts<-seq(from=(-4*tr.sd[i]), to=(4*tr.sd[i]), length.out=ncat[i]-1)
        if(tr.dist[i]=="norm"){
          pass.one<-cut(traits.final[,i], breaks=cuts, labels=FALSE,
                        include.lowest=TRUE)
          pass.two<-ifelse(traits.final[,i]<cuts[1], 0, pass.one)
          pass.three<-ifelse(traits.final[,i]>=cuts[(ncat[i]-1)], ncat[i]-1,
                             pass.two)
        }else{
          pass.one<-cut(log(traits.final[,i]), breaks=cuts, labels=FALSE,
                        include.lowest=TRUE)
          pass.two<-ifelse(log(traits.final[,i])<cuts[1], 0, pass.one)
          pass.three<-ifelse(log(traits.final[,i])>=cuts[(ncat[i]-1)], ncat[i]-1,
                             pass.two)
        } # end if-else
        # compute mean likelihood for each categorized trait value
        group.means<-aggregate(traits.d[,i], by=list(pass.three), mean)
        dimnames(group.means)[[2]]<-c("level", "avg")
        for(g in 1:nrow(group.means)){
          traits.d.final[pass.three==group.means$level[g],i]<-
            group.means$avg[g]
        } # end for g
      }else{
        # categorize
        #cuts<-c(0, tr.sd[i]*c(1:10)/10) # this fixes number of categories at 10
        cuts<-seq(from=0, to=tr.sd[i], length.out=(ncat[i]+1))
        pass.three<-cut(traits.final[,i], breaks=cuts, labels=FALSE,
                        include.lowest=TRUE)
        # no need for mean trait likelihoods under a uniform distribution
      } # end if-else
      traits.final[,i]<-as.character(pass.three)
      if(tr.type[i]=="ord"){traits.final[,i]<-as.ordered(traits.final[,i])}
      if(tr.type[i]=="cat"){traits.final[,i]<-as.factor(traits.final[,i])}
    } # end if
  } # end for i
  # Now sample species from species pool for each community, based on
  # species richness specified in s, and repeated as specified in r
  # First, create a list of abundance vectors, with enough elements to have a
  # separate vector of length p for each combination of species richness and
  # replicate
  abun <- list(rep(0, p))
  if(dups>0){abun <- list(rep(0, 2*p))}
  abun <- rep(abun, r * l.s)
  # Next determine which species should be selected for non-zero abundances.
  # Check how communities should reflect distributions specified for the overall
  # species pool, and how abundances should be distributed among the species
  # selected.
  traits.prob<-apply(traits.d.final, 1, prod)
  if(presence=="random"){pres.prob<-rep(1,p)}
  if(presence=="common"){pres.prob<-traits.prob}
  if(presence=="rare"){pres.prob<-rescale(1/traits.prob)}
  if(abundance=="random"){abun.prob<-rep(1,p)}
  if(abundance=="common"){abun.prob<-traits.prob}
  if(abundance=="rare"){abun.prob<-rescale(1/traits.prob)}
  # Also enforce specified minima and maxima and excluded values
  if(sum(is.na(commin))<t || sum(is.na(commax))<t || sum(is.na(comnot))<t){
    for(i in 1:t){
      if(!is.na(commin[i])){pres.prob<-
        ifelse(as.numeric(as.character(traits.final[,i]))<commin[i], 0, pres.prob)}
      if(!is.na(commax[i])){pres.prob<-
        ifelse(as.numeric(as.character(traits.final[,i]))>commax[i], 0, pres.prob)}
      comnot.type[[i]]<-match.arg(comnot.type[[i]], c("single", "range"))
      if(comnot.type[[i]]=="single"){
        if(sum(!is.na(comnot[[i]]))>0){pres.prob<-
          ifelse(is.element(as.numeric(as.character(traits.final[,i])),
                            as.numeric(as.character(comnot[[i]]))), 0, pres.prob)}
      }else{
        if(sum(!is.na(comnot[[i]]))>0){
          istart<-min(comnot[[i]], na.rm=TRUE)
          istop<-max(comnot[[i]], na.rm=TRUE)
          pres.prob<-
            ifelse(findInterval(as.numeric(as.character(traits.final[,i])),
                                c(istart, istop), rightmost.closed=TRUE)==1, 0, pres.prob)
        }
      }
    }
  }
  if(sum(pres.prob!=0)<=max(s)){stop("Values in commin and/or commax and/or comnot are too
                                     restrictive!")}
  # adjust length of likelihood vectors if duplicate species are required
  if(dups>0){pres.prob<-c(pres.prob,rep(0,p))}
  if(dups>0){abun.prob<-rep(abun.prob,2)}
  # assign abundance values in each element of abun to the relevant number of
  # species picked from the species pool, taking into account duplicates and
  # how species and their abundances should reflect the pool's trait distributions.
  abun<-mapply(fill.abun, empty=abun, comsize=nb.sp, MoreArgs=list( presprob=pres.prob,
                                                                    abunprob=abun.prob, abundist=abun.dist, abunsd=abun.sd, ids=dups))
  # ensure that abundance values are positive (normal distribution can yield
  # negative values, which may be problematic in further analysis)
  if(min(abun)<0){abun<-apply(X=abun, MARGIN=2, FUN=pos.abun)}
  # transpose the abundance matrix for a site-by-species matrix with rows
  # equating to communities, and columns to species
  abun <- t(abun)
  # Now assign row and column names to the trait and community matrices
  names.tr<-paste("TR", c(1:t), sep="")
  names.sp<-paste("SP", c(1:p), sep="")
  if(dups>0){
    names.sp<-c(names.sp, paste("DP", c(1:p), sep=""))
    traits.final<-rbind(traits.final, traits.final)
    traits.prob<-rep(traits.prob,2)
  }
  names.com<-paste("COM", c(1:(r * l.s)), sep = "")
  traits.final<-cbind(traits.final, traits.prob)
  dimnames(traits.final)<-list(names.sp, c(names.tr, "LIKELIHOOD"))
  dimnames(abun)<-list(names.com, names.sp)
  # Finally, track min, max, range, skewness and kurtosis in trait values and
  # trait likelihood for the species pool and per community, plus user input values
  comstats<-as.data.frame(matrix(NA, nrow=1+r*l.s, ncol=26+t*13))
  dimnames(comstats)[[2]]<-c("Site", "SR_pool", "SR_site", "N_dups", "N_traits",
                             paste("T", rep(1:t, each=13), c("_dist", "_sd", "_type", "_levels", "_setmin", "_obsmin",
                                                             "_setmax", "_obsmax", "_range", "_unique.values", "_var", "_skewness", "_kurtosis"), sep=""),
                             "likeli_min", "likeli_max", "likeli_range", "likeli_var", "likeli_skewness", "likeli_kurtosis",
                             "presence", "abundance", "abun_dist", "abun_sd", "abun_min", "abun_max", "abun_var",
                             "abun_skewness", "abun_kurtosis", "tot_abun", "relabun_min", "relabun_max", "relabun_var",
                             "relabun_skewness", "relabun_kurtosis")
  comstats$Site<-c("pool", names.com)
  comstats$SR_pool<-p
  comstats$SR_site<-c(p, nb.sp)
  comstats$N_dups<-dups
  comstats$N_traits<-t
  comstats$presence<-presence
  comstats$abundance<-abundance
  comstats$abun_dist<-abun.dist
  comstats$abun_sd<-abun.sd
  for(i in 0:(r*l.s)){
    #i<-2 # for testing
    # deal with the special case of stats for the species pool
    if(i==0){
      indx<-c(1:nrow(traits.final))
      comstats$likeli_min[i+1]<-min(traits.prob[indx])
      comstats$likeli_max[i+1]<-max(traits.prob[indx])
      comstats$likeli_range[i+1]<-comstats$likeli_max[i+1]-comstats$likeli_min[i+1]
      comstats$likeli_var[i+1]<-var(traits.prob[indx])
      comstats$likeli_skewness[i+1]<-skewness(traits.prob[indx])
      comstats$likeli_kurtosis[i+1]<-kurtosis(traits.prob[indx])
      comstats$abun_min[i+1]<-NA
      comstats$abun_max[i+1]<-NA
      comstats$abun_var[i+1]<-NA
      comstats$abun_skewness[i+1]<-NA
      comstats$abun_kurtosis[i+1]<-NA
      comstats$tot_abun[i+1]<-NA
      comstats$relabun_min[i+1]<-NA
      comstats$relabun_max[i+1]<-NA
      comstats$relabun_var[i+1]<-NA
      comstats$relabun_skewness[i+1]<-NA
      comstats$relabun_kurtosis[i+1]<-NA
    }else{
      indx<-which(abun[i,]>0)
      comstats$likeli_min[i+1]<-min(traits.prob[indx])
      comstats$likeli_max[i+1]<-max(traits.prob[indx])
      comstats$likeli_range[i+1]<-comstats$likeli_max[i+1]-comstats$likeli_min[i+1]
      comstats$likeli_var[i+1]<-var(traits.prob[indx])
      comstats$likeli_skewness[i+1]<-skewness(traits.prob[indx])
      comstats$likeli_kurtosis[i+1]<-kurtosis(traits.prob[indx])
      comstats$abun_min[i+1]<-min(abun[i,indx])
      comstats$abun_max[i+1]<-max(abun[i,indx])
      comstats$abun_var[i+1]<-var(abun[i,indx])
      comstats$abun_skewness[i+1]<-skewness(abun[i,indx])
      comstats$abun_kurtosis[i+1]<-kurtosis(abun[i,indx])
      comstats$tot_abun[i+1]<-sum(abun[i,indx])
      comstats$relabun_min[i+1]<-min(abun[i,indx]/sum(abun[i,indx]))
      comstats$relabun_max[i+1]<-max(abun[i,indx]/sum(abun[i,indx]))
      comstats$relabun_var[i+1]<-var(abun[i,indx]/sum(abun[i,indx]))
      comstats$relabun_skewness[i+1]<-skewness(abun[i,indx]/sum(abun[i,indx]))
      comstats$relabun_kurtosis[i+1]<-kurtosis(abun[i,indx]/sum(abun[i,indx]))
    } # end if-else
    for(j in 1:t){
      #j<-1 # for testing
      names.col<-paste("T", rep(j,each=12), c("_dist", "_sd", "_type", "_levels",
                                              "_setmin", "_obsmin", "_setmax", "_obsmax", "_range", "_unique.values",
                                              "_var", "_skewness", "_kurtosis"), sep="")
      comstats[i+1,names.col[1]]<-tr.dist[j]
      comstats[i+1,names.col[2]]<-tr.sd[j]
      comstats[i+1,names.col[3]]<-tr.type[j]
      comstats[i+1,names.col[4]]<-ncat[j]
      comstats[i+1,names.col[5]]<-commin[j]
      comstats[i+1,names.col[7]]<-commax[j]
      # compute onbserved min, max, range, skewness and kurtosis per trait
      comstats[i+1, names.col[6]]<-min(as.numeric(as.character(traits.final[indx,j])))
      comstats[i+1, names.col[8]]<-max(as.numeric(as.character(traits.final[indx,j])))
      comstats[i+1, names.col[9]]<-comstats[i+1, names.col[8]]-comstats[i+1,
                                                                        names.col[6]]
      comstats[i+1, names.col[10]]<-
        length(unique(as.numeric(as.character(traits.final[indx,j]))))
      comstats[i+1, names.col[11]]<-var(as.numeric(as.character(traits.final[indx,j])))
      comstats[i+1, names.col[12]]<-
        skewness(as.numeric(as.character(traits.final[indx,j])))
      comstats[i+1, names.col[13]]<-
        kurtosis(as.numeric(as.character(traits.final[indx,j])))
    } # end for j
  } # end for i
  return(list(T=traits.final, A=abun, S=comstats))
  } # end function simul.comms


# ===== END OF simul.comms R CODE ===== 

###########################################################################

###
##create regional species and trait pools from a single universal starting pool simultaneously using simul.comms function
###

#s = num sp in comms
#r = num replicate comms 
#p = num sp in regional pool 
#t = num traits


RegPools <-simul.comms(p=1500, r=1, s=c(50, 100, 500, 1000), t=24, 
                       tr.method=c("lnorm", "norm", "unif", "lnorm", "norm", "unif", "lnorm", "norm", "unif", "lnorm", "norm", "unif", "lnorm", "norm", "unif", "lnorm", "norm", "unif", "lnorm", "norm", "unif", "lnorm", "norm", "unif"), #sampling distribution of traits
                       tr.type="con",
                       presence="random", #sp randomly selected from pool 
                       abun.method="fixed:1", #set to pres-abs
                       abundance="random",
                       dups=0)

#save these out to have the record of the original universal pool and parameters. Must use same universal pool for ALL runs for comparability and reproducability sake.
write.csv(RegPools$T, "RegPoolsTraits_demo.csv", row.names = T)
write.csv(RegPools$A, "RegPoolsPools_demo.csv", row.names = T)
write.csv(RegPools$S, "RegPoolsParams_demo.csv", row.names = T)

##################################################################################################


###########

###
##Load Kraft and Ackerly 2010 Trait-based community assembly simulation functions 
###


## Nathan Kraft
## University of California, Berkeley 
## and currently
## Biodiversity Research Centre, University of British Columbia
## nkraft@biodiversity.ubc.ca
##Last modified 1 Oct 2019 by Brooks Kohli to use gowdis to create dissimilarity matrices



## Supplemental code from
## Nathan J. B. Kraft and Davd D. Ackerly, 2010, Functional trait and phylogenetic tests of community assembly across spatial scales in an Amazonian forest, Ecological Monographs 80:401-422. 

## Functions here are based on the EVELYN community assembly model(which was originally written in JAVA) from:
## Kraft, N. J. B., W. K. Cornwell, C. O. Webb, and D. D. Ackerly. 2007. Trait evolution, community assembly, and the phylogenetic structure of ecological communities. American Naturalist 170:271-283.

## for the original verbal description of the competition algorithm used here and general inspiration see:  
## Colwell, R. K., and D. W. Winkler. 1984. A null model for null models in biogeography., Pages 344-359 in D. R. Strong, D. S. Simberloff, L. G. Abele, and A. B. Thistle, eds. Ecological communities: conceptual issues and the evidence. Princeton, NJ, Princeton University Press.

########################
### GENERAL FUNCTIONS ##
########################

#check for duplicate species names in a community vector
validate_community=function(community=community){
  community[,1]->spnames
  if(length(spnames)>length(unique(spnames))){
    print("warning- duplicate species names")
    return(FALSE)
  }
  
  return(TRUE)
}

#remove a species from the community (general function used in all models)
zap=function(target, community=community){
  #target can be either a species name (letters) or a column number(numeric).  Can also be a vector of names
  #defaults to looking for a vector called community with two columns (first for name and second for trait value)
  
  ##currently no error is given if some targets match species in the community and others do not- only matching species removed
  
  if(!is.numeric(target)){
    target<-which(target==community[,1])
    
    if(length(target)<1){
      print("target does not match anyone in community- can not zap")
      return(community)
    }
    
  }
  
  community<-community[-target,]
  return(community)
  
}

####################
## RANDOM ASSEMBLY #
####################

## Randomly removes species, weighted by abundance, from the species pool until a final_richness value is reached	
random_assembly=function(pool=pool, final_richness, abund=NULL){
  size<-nrow(pool)
  n_victim<-size-final_richness
  if(n_victim<1){
    print("pool has equal or lesser richness to final size- can't do random assembly")
    return(pool)
  }
  
  
  if(is.null(abund)){
    alive<-sample(1:size, final_richness)
    return(pool[alive,])
  }
  
  merge(pool, abund)->poolPlus
  
  alive<-sample((1:nrow(poolPlus)), final_richness, prob=poolPlus$abund)
  
  return(poolPlus[alive,c(1,2)])
  
  
  
}


######################
## COMPETITION MODEL #
######################

#find the most similar taxa within the community:
find_most_similar_taxa=function(community=community, speciesNamesCol1=TRUE){
  if(speciesNamesCol1!=TRUE){
    print("first col needs to be species names for now")
    return(community)
  }
  ##somewhat of an awkward fix here:
  labels<-community[,1]
  traits<-community[,-1]
  gowdis(traits)->m     ##on 1 Oct, changed from dist to gowdis function-BAK
  as.matrix(m)->m2
  which(m2==min(m), arr.ind=TRUE)->m3
  as.vector(m3[,1])->a
  as.vector(m3[,2])->b
  ##will contain multiples but this is good- more weight to a "sandwiched" taxa being removed
  taxa_indices<-c(a[which(a!=b)], b[which(a!=b)])
  return(labels[taxa_indices])
}


## find most similar species and remove it from a community	
one_round_competition=function(community=community){
  threatened_list<-find_most_similar_taxa(community)
  n<-length(threatened_list)
  if(n<2){
    print("couldn't find most similar taxa for competition")
    return(community)
  }
  return(zap(sample(threatened_list, 1), community=community))
}

## remove species via competition from the community until a specified number of taxa are left (nfinal)
compete_until=function(nfinal, community=community){
  nrow(community)->start
  if(nfinal>start){
    print("community is already smaller than target- can't compete")
    return(community)
  }
  
  if(nfinal<4){
    print("target for competition is too small- can't compete")
    return(community)
  }
  
  start-nfinal->to_kill
  
  for(i in 1:to_kill){
    one_round_competition(community)->community
    
  }
  
  if(nrow(community)!=nfinal){
    print('error- competition ended with incorrect number of species')
  }
  return(community)
  
}

######################
## Habitat Filtering #
######################

#Used to identify which species in a community is farthest from the trait optima used for habitat filtering- farthest species are removed first.  

find_farthest_from_optima=function(optima=optima, community=community, speciesNamesCol1=TRUE){
  if(speciesNamesCol1!=TRUE){
    print("first col needs to be species names for now")
    return(community)
  }
  
  if( (length(optima)!=(ncol(community)-1))){
    print("error- filtering optima is not right dimension for community")
    return(community)
  }
  
  ##cludgy- will generate warnings- not sure how to avoid
  paste(community[,1])->start
  labels<-c(start,"optima")
  traits<-suppressWarnings(rbind(community, optima)[,-1])
  
  gowdis(traits)->m   ##1 Oct, changed dist to gowdis function - BAK
  as.matrix(m)->m2
  ##get last row, which is the distance of each taxa from optimum
  m2[,nrow(m2)]->comp
  which(comp==max(comp))->index
  return(labels[index])
  
}


##identify farthest species from optima and remove it
one_round_filtering_optima=function(community=community, optima=optima){
  
  threatened_list<-find_farthest_from_optima(optima, community)
  n<-length(threatened_list)
  if(n<1){
    print("couldn't find farthest taxa for filtering")
    return(community)
  }
  return(zap(sample(threatened_list, 1), community=community))
  
}


##run filtering model until only nfinal number of taxa are left in the community
filter_until=function(nfinal, community=community, optima=optima){
  nrow(community)->start
  if(nfinal>start){
    print("community is already smaller than target- can't filter")
    return(community)
  }
  
  if(nfinal<4){
    print("target for filtering is too small- can't filter")
    return(community)
  }
  
  start-nfinal->to_kill
  
  for(i in 1:to_kill){
    one_round_filtering_optima(community, optima)->community
    
  }
  
  if(nrow(community)!=nfinal){
    print('error- filtering ended with incorrect number of species')
  }
  return(community)
  
}


#######END of Kraft and Ackerly functions####

#################################################################
##############################################################


##BAK code


###OVERALL STEPS
#1. create regional pools from universal pool created by using McPherson etal code above

#2. coarsen traits at multiple benchmarks for entire pool

#3. assemble (100) local comms with each process-based mechanism (Kraft and Ackerly code) using continuous data

#4. assemble (1000) random communities for use as null models

#5.for each local community from #3 and #4, calculate functional structure indices for each trait resolution level

#6.Compare each nonrandom community index value to the null distribution of the 1000 randomly assembled communities values from #5. =Calculate SES values. Summarize and plot.

#7.Repeat steps 3-6; change pool size, community size, number of traits, mixture of trait resolutions, etc 

##############################################################

##############################################################

###
##Part I : Modify and create analysis pools from the Universal pool 
###

##1. modify and create pools

#read in universal Trait and different sized regional pool matrices  generated by UnivPoolGen script
RegPoolsTr <- read.csv("RegPoolsTraits_demo.csv", row.names = 1)
RegPoolsPools <- read.csv("RegPoolsPools_demo.csv", row.names = 1)


# Universal RegPool Trait Values
dim(RegPoolsTr) #species and their trait info (sp X traits)

#logT the lognorm traits
RegPools1 <- as.data.frame(RegPoolsTr)
head(RegPools1)
#need a character column for Kraft code
RegPools1$Sp <- as.factor(c(1:1500))

TraitsLogT24 <- RegPools1 %>% transmute(Sp=Sp, TR1log = (log10(TR1)), TR2=TR2, TR3=TR3, TR4log = (log10(TR4)), TR5=TR5, TR6=TR6, TR7log = (log10(TR7)), TR8=TR8, TR9=TR9, TR10log = (log10(TR10)), TR11=TR11, TR12=TR12, TR13log = (log10(TR13)), TR14=TR14, TR15=TR15, TR16log = (log10(TR16)), TR17=TR17, TR18=TR18, TR19log = (log10(TR19)), TR20=TR20, TR21=TR21, TR22log = (log10(TR22)), TR23=TR23, TR24=TR24, LIKELIHOOD=LIKELIHOOD)

hist(TraitsLogT24$TR16log)
row.names(TraitsLogT24) <- row.names(RegPools1)


#Reg Pools of varying Size 
rowSums(RegPoolsPools) #size of regional pools in community matrix 

RegPool50 <- as.data.frame(t(RegPoolsPools[1,])) #added t when changed from list to read in df
colnames(RegPool50) <- "pool"
RegPool50$Sp <- as.factor(c(1:1500))  #NEED this to join with trait matrix, but not after
RegPool50 <- RegPool50[RegPool50$pool == 1,]
RegPool50$pool <- RegPool50$pool*2  #need to make these a number besides 1 for proper joining 
RegPool50

RegPool100 <- as.data.frame(t(RegPoolsPools[2,]))
colnames(RegPool100) <- "pool"
RegPool100$Sp <- as.factor(c(1:1500))
RegPool100 <- RegPool100[RegPool100$pool == 1,]
RegPool100$pool <- RegPool100$pool*2  #need to make these a number besides 1 for proper joining 
RegPool100

RegPool500 <- as.data.frame(t(RegPoolsPools[3,]))
colnames(RegPool500) <- "pool"
RegPool500$Sp <- as.factor(c(1:1500))
RegPool500 <- RegPool500[RegPool500$pool == 1,]
RegPool500$pool <- RegPool500$pool*2  #need to make these a number besides 1 for proper joining 
RegPool500

RegPool1000 <- as.data.frame(t(RegPoolsPools[4,]))
colnames(RegPool1000) <- "pool"
RegPool1000$Sp <- as.factor(c(1:1500))
RegPool1000 <- RegPool1000[RegPool1000$pool == 1,]
RegPool1000$pool <- RegPool1000$pool *2  #need to make these a number besides 1 for proper joining 
RegPool1000

#now join these with trait matrices 
Traits24Sp50 <- inner_join(TraitsLogT24, RegPool50, by = "Sp")
rownames(Traits24Sp50) <- rownames(RegPool50)
Traits24Sp50 <- Traits24Sp50[-27] #still contains Sp as col1 and Likelihoods as col26

Traits24Sp100 <- inner_join(TraitsLogT24, RegPool100, by = "Sp")
rownames(Traits24Sp100) <- rownames(RegPool100)
Traits24Sp100 <- Traits24Sp100[-27] #still contains Sp as col1 and Likelihoods as col26

Traits24Sp500 <- inner_join(TraitsLogT24, RegPool500, by = "Sp")
rownames(Traits24Sp500) <- rownames(RegPool500)
Traits24Sp500 <- Traits24Sp500[-27] #still contains Sp as col1 and Likelihoods as col26

Traits24Sp1000 <- inner_join(TraitsLogT24, RegPool1000, by = "Sp")
rownames(Traits24Sp1000) <- rownames(RegPool1000)
Traits24Sp1000 <- Traits24Sp1000[-27] #still contains Sp as col1 and Likelihoods as col26

#removing likelihood column
continTr_24tr50sp <- Traits24Sp50[-26]
continTr_24tr100sp <- Traits24Sp100[-26]  
continTr_24tr500sp <- Traits24Sp500[-26]
continTr_24tr1000sp <- Traits24Sp1000[-26]
#############################################################################################






###############################################################################################


###
##Part II : Example simulation for an iteration of a given pool size, community size, trait number
#####

#Conceptual trait resolution project - simulation steps
#author: Brooks Kohli
#Date Created: 18 Sept 2019
#Date last modified: 2 July 2020


#################
##2. Coarsen all traits at multiple resolutions
###########

###########################################
###starting with the continuous traits for a given Pool size...

##first level of coarsening: rounding (maintains fine structure; small loss of information)
rounded_24tr100sp <- continTr_24tr100sp %>% mutate_if(is.numeric, round, 2)

##successive levels of coarsening via categorical binning
##preliminary analyses revealed it unnessary to do 64 and 32-group categorizations
cuts <- c(16, 8, 4, 2) #using multiples of 2 maintains bin boundaries and a perfectly nested coarsening.
matrices <- vector("list", 4)
for (i in cuts) {matrices[[i]] <- continTr_24tr100sp %>% mutate_if(is.numeric, cut, breaks=i, labels = F)} 
matrices

#save them separately
categ2_24tr100sp <- matrices[[2]]
categ4_24tr100sp <- matrices[[4]]
categ8_24tr100sp <- matrices[[8]]
categ16_24tr100sp <- matrices[[16]]


#row names
row.names(rounded_24tr100sp)<-row.names(continTr_24tr100sp)
row.names(categ16_24tr100sp)<-row.names(continTr_24tr100sp)
row.names(categ8_24tr100sp)<-row.names(continTr_24tr100sp)
row.names(categ4_24tr100sp)<-row.names(continTr_24tr100sp)
row.names(categ2_24tr100sp)<-row.names(continTr_24tr100sp)

##################################

#subset various number of traits per resolution level
continTr_3tr100sp <- continTr_24tr100sp[1:4]
continTr_6tr100sp <- continTr_24tr100sp[1:7]
continTr_12tr100sp <- continTr_24tr100sp[1:13]

rounded_3tr100sp <- rounded_24tr100sp[1:4]
rounded_6tr100sp <- rounded_24tr100sp[1:7]
rounded_12tr100sp <- rounded_24tr100sp[1:13]

categ16_3tr100sp <- categ16_24tr100sp[1:4]
categ16_6tr100sp <- categ16_24tr100sp[1:7]
categ16_12tr100sp <- categ16_24tr100sp[1:13]

categ8_3tr100sp <- categ8_24tr100sp[1:4]
categ8_6tr100sp <- categ8_24tr100sp[1:7]
categ8_12tr100sp <- categ8_24tr100sp[1:13]

categ4_3tr100sp <- categ4_24tr100sp[1:4]
categ4_6tr100sp <- categ4_24tr100sp[1:7]
categ4_12tr100sp <- categ4_24tr100sp[1:13]

categ2_3tr100sp <- categ2_24tr100sp[1:4]
categ2_6tr100sp <- categ2_24tr100sp[1:7]
categ2_12tr100sp <- categ2_24tr100sp[1:13]

#####LISTS OF FINAL MATRICES - from highest to lowest resolution#####

#combine all resolution matrices into a list per Reg Pool Size
ResMatrices24tr100sp <- list(continTr_24tr100sp, rounded_24tr100sp, categ16_24tr100sp, categ8_24tr100sp, categ4_24tr100sp, categ2_24tr100sp)

ResMatrices12tr100sp <- list(continTr_12tr100sp, rounded_12tr100sp, categ16_12tr100sp, categ8_12tr100sp, categ4_12tr100sp, categ2_12tr100sp)

ResMatrices6tr100sp <- list(continTr_6tr100sp, rounded_6tr100sp, categ16_6tr100sp, categ8_6tr100sp, categ4_6tr100sp, categ2_6tr100sp)

ResMatrices3tr100sp <- list(continTr_3tr100sp, rounded_3tr100sp, categ16_3tr100sp, categ8_3tr100sp, categ4_3tr100sp, categ2_3tr100sp)

###########################################################################################



#############################################################
##3. assemble communities via random, filtering, and limiting similarity processes.
#######################

##call all functions in the Kraft R code
#need column of species names and trait values

#set pool size - choose from 50, 100, 500, 1000
p <- 100 

#set Species richness for communities - choose from p/3, p/4, p/5
nc <- round(((1/4)*p), 0)  #force it to round, otherwise can get mismatched numbers and warnings from the asssembly algorithms.

#number of traits to use - choose from 3, 6, 12, 24
tr <- 6

#which simulation run variant is this? = NumTraits_PoolSize_CommSize
SimRunVar <- "tr6p100nc1.4"  #will be used as plot label for clarity


#establish a baseline of all sp to re-join all to
AllSp <- RegPool100[-2]
AllSp <- as.data.frame(t(AllSp))
AllSp
##
###


######################################################################################
#create 100 competition-assembled communities
########

Comp100 <- vector("list", 100)
for (c in 1:100) {Comp100[[c]] <- compete_until(nc, continTr_6tr100sp)} 

#extract resulting community lists to combine with random communities

get.sp.filt <- function(comm, nc){comm2 <- comm[-1]
FiltComm <- as.data.frame(t(comm2))
FiltComm[7,] <- rep(1, nc)
FiltComm <- FiltComm[7,]
combined <- full_join(AllSp, FiltComm, all = TRUE)
cleaned <- combined[2,]
}


#run for all 100 competitively structured communities
complist <- lapply(Comp100, get.sp.filt)
CompComms <- data.frame(matrix(unlist(complist), nrow=length(complist), byrow=T))
colnames(CompComms) <- colnames(AllSp)
CompComms #this is the final set of 100 communities

################################################################################################

## for filtering use filter_until(). You need to specify what the optimum trait value is- in the 2010  Ecological Monographs paper Kraft and Ackerly randomly placed the optima in the environment.
## for multiple traits, the optima needs to be of the same number of dimensions
## ties are broken randomly so repeat calls will give slightly different answers in some cases

#####
#create environmentally filtered communities (100) with randomized optima
EF100 <- vector("list", 100)
for (e in 1:100) {EF100[[e]] <- filter_until(nc, continTr_6tr100sp, optima=c(
  #3 traits
  runif(1,min(continTr_6tr100sp[,2]),max(continTr_6tr100sp[,2])),
  runif(1,min(continTr_6tr100sp[,3]),max(continTr_6tr100sp[,3])),
  runif(1,min(continTr_6tr100sp[,4]),max(continTr_6tr100sp[,4])),
  #6 traits
  runif(1,min(continTr_6tr100sp[,5]),max(continTr_6tr100sp[,5])),
  runif(1,min(continTr_6tr100sp[,6]),max(continTr_6tr100sp[,6])), 
  runif(1,min(continTr_6tr100sp[,7]),max(continTr_6tr100sp[,7]))  ))}
#12 traits
#  runif(1,min(continTr_3tr100sp[,8]),max(continTr_3tr100sp[,8])),
#  runif(1,min(continTr_3tr100sp[,9]),max(continTr_3tr100sp[,9])),
#  runif(1,min(continTr_3tr100sp[,10]),max(continTr_3tr100sp[,10])),
#  runif(1,min(continTr_3tr100sp[,11]),max(continTr_3tr100sp[,11])),
#  runif(1,min(continTr_3tr100sp[,12]),max(continTr_3tr100sp[,12])),
#  runif(1,min(continTr_3tr100sp[,13]),max(continTr_3tr100sp[,13])),
#24 traits
#  runif(1,min(continTr_24tr100sp[,14]),max(continTr_24tr100sp[,14])),
#  runif(1,min(continTr_24tr100sp[,15]),max(continTr_24tr100sp[,15])),
#  runif(1,min(continTr_24tr100sp[,16]),max(continTr_24tr100sp[,16])),
#  runif(1,min(continTr_24tr100sp[,17]),max(continTr_24tr100sp[,17])),
#  runif(1,min(continTr_24tr100sp[,18]),max(continTr_24tr100sp[,18])),
#  runif(1,min(continTr_24tr100sp[,19]),max(continTr_24tr100sp[,19])),
#  runif(1,min(continTr_24tr100sp[,20]),max(continTr_24tr100sp[,20])),
#  runif(1,min(continTr_24tr100sp[,21]),max(continTr_24tr100sp[,21])),
# runif(1,min(continTr_24tr100sp[,22]),max(continTr_24tr100sp[,22])), 
# runif(1,min(continTr_24tr100sp[,23]),max(continTr_24tr100sp[,23])),
#  runif(1,min(continTr_24tr100sp[,24]),max(continTr_24tr100sp[,24])), 
# runif(1,min(continTr_24tr100sp[,25]),max(continTr_24tr100sp[,25]))  
#  ))} 


#run for all 100  environmental filter structured communities
eflist <- lapply(EF100, get.sp.filt)
EFComms <- data.frame(matrix(unlist(eflist), nrow=length(eflist), byrow=T))
colnames(EFComms) <- colnames(AllSp)
EFComms #this is the final set of 100 communities

#####################################################################



####################
###4. assemble random communities
###########

Rand1000 <- vector("list", 1000)
for (r in 1:1000) {Rand1000[[r]] <- random_assembly(continTr_6tr100sp, nc, abund = NULL)} 

#
randlist <- lapply(Rand1000, get.sp.filt)
RandComms <- data.frame(matrix(unlist(randlist), nrow=length(randlist), byrow=T))
colnames(RandComms) <- colnames(AllSp)



############################################################################

#combine random comms with the filtered comms

#competition
comp_rand <- full_join(RandComms, CompComms,  all = TRUE)
#add EF and Combine all together
AllComms <- full_join(comp_rand, EFComms, all = TRUE)
#change NAs to 0s
AllComms[is.na(AllComms)]<-0
#check that all sp are in at least 1 comm
range(colSums(AllComms))
#check occurrence of sp in each filtered dataset
colSums(AllComms[1001:1100,]) #competition
colSums(AllComms[1101:1200,]) #EF


################################################

##Create species distance matrices via gower distance (Pavoine 2009) so that the same dissimilarity measure is used for all iterations.


###HAVE TO REMOVE THE SPECIES NUMBER COLUMN
#function to create distance matrix from species pool trait matrices
distmat <- function(trts){trts2 <- trts[-1]
Dist.G <-gowdis(trts2) ###create gower dissim matrix
DistMat.continG <- as.matrix(Dist.G)  #convert to matrix form
}

####################

#Do for each sp X trait pool per resolution

#loop through trait resolution sets to create each distance matrix
AllDists <- lapply(ResMatrices6tr100sp, distmat)

###################################################################




########
###5. & 6. calculate functional structure indices at each trait resolution level
############

#######
#calculate indices
####


#function to calculate SES values 

ses.calc <- function(values){
  null.mean <- mean(values[1:1000])
  null.sd <- sd(values[1:1000])
  ses.value <- (values - null.mean)/ null.sd
}


###FDis###
#observed values
FDis.list <- lapply(AllDists, dbFD, a = AllComms, w.abun = FALSE, calc.FRic = FALSE, stand.FRic = FALSE, calc.FDiv = FALSE, calc.CWM = FALSE, print.pco = FALSE) #may take some time to complete

FDis.ALL <- data.frame(t(matrix(unlist(FDis.list), nrow=length(FDis.list), byrow=T)))

head(FDis.ALL)

colnames(FDis.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")

head(FDis.ALL) #this is the final set of index values calculated for all levels of resolution.

FDis.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(FDis.ALL)

#SES values
FDis.newlist <- list(FDis.list[[1]]$FDis, FDis.list[[2]]$FDis, FDis.list[[3]]$FDis, FDis.list[[4]]$FDis, FDis.list[[5]]$FDis, FDis.list[[6]]$FDis)

SES.FDis.list <- lapply(FDis.newlist, ses.calc)

SES.FDis.ALL <- data.frame(t(matrix(unlist(SES.FDis.list), nrow=length(SES.FDis.list), byrow=T)))

head(SES.FDis.ALL)

colnames(SES.FDis.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")

head(SES.FDis.ALL) 

SES.FDis.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)

head(SES.FDis.ALL)

##########################

###SES.MNND###

#observed values
MNTD.list <- lapply(AllDists, mntd, samp=AllComms)
MNTD.ALL <- data.frame(t(matrix(unlist(MNTD.list), nrow=length(MNTD.list), byrow=T)))
head(MNTD.ALL)
colnames(MNTD.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(MNTD.ALL) #this is the final set of index values calculated for all levels of resolution.
MNTD.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(MNTD.ALL)

#and do it for SES values
SES.MNTD.list <- lapply(MNTD.list, ses.calc)
SES.MNTD.ALL <- data.frame(t(matrix(unlist(SES.MNTD.list), nrow=length(SES.MNTD.list), byrow=T)))
head(SES.MNTD.ALL)
colnames(SES.MNTD.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(SES.MNTD.ALL) #this is the final set of index values calculated for all levels of resolution.
SES.MNTD.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(SES.MNTD.ALL)


#######################

##FD dendrogram sums (adapted from Aiba etal 2013 function 'petchey')

Fdendro <- function(traits, comms){
  tree <- hclust(gowdis(traits[-1]),method="average") #this uses UPGMA
  FD.out <- treedive(comms,tree)
}

#All resolution sets
FD.list <- lapply(ResMatrices6tr100sp, Fdendro, comms = AllComms)
FD.ALL <- data.frame(t(matrix(unlist(FD.list), nrow=length(FD.list), byrow=T)))
head(FD.ALL)
colnames(FD.ALL) <-c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(FD.ALL) #this is the final set of index values calculated for all levels of resolution.
FD.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(FD.ALL)

#SES
SES.FD.list <- lapply(FD.list, ses.calc)
SES.FD.ALL <- data.frame(t(matrix(unlist(SES.FD.list), nrow=length(SES.FD.list), byrow=T)))
head(SES.FD.ALL)
colnames(SES.FD.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(SES.FD.ALL) #this is the final set of index values calculated for all levels of resolution.
SES.FD.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(SES.FD.ALL)

##########################

###SES.MPD###

#observed values
MPD.list <- lapply(AllDists, mpd, samp=AllComms)
MPD.ALL <- data.frame(t(matrix(unlist(MPD.list), nrow=length(MPD.list), byrow=T)))
head(MPD.ALL)
colnames(MPD.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(MPD.ALL) #this is the final set of index values calculated for all levels of resolution.
MPD.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(MPD.ALL)

#and do it for SES values
SES.MPD.list <- lapply(MPD.list, ses.calc)
SES.MPD.ALL <- data.frame(t(matrix(unlist(SES.MPD.list), nrow=length(SES.MPD.list), byrow=T)))
head(SES.MPD.ALL)
colnames(SES.MPD.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(SES.MPD.ALL) #this is the final set of index values calculated for all levels of resolution.
SES.MPD.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(SES.MPD.ALL)



#######################################

#Functional distinctiveness (via dendrogram)  - evol.distinct in picante
AllCommsListOrig <- c(Rand1000, Comp100, EF100)  #combine all minimal species lists (SR = nc, not p) to calculate each tree separately  
AllCommsList <- lapply(AllCommsListOrig, t) #have to transpose to get sp names as columns
#create full functional dendrogram
#prune each tree
#calculate distinctiveness
#take geometric mean of sp values
FDphylo <- function(traits, comms){
  tree <- hclust(gowdis(traits[-1]),method="average") #this uses UPGMA
  tree.phylo <- as.phylo(tree) #as.phylo is called from ape
  pruned <- prune.sample(comms, tree.phylo)
  FDI.out <- evol.distinct(pruned, type = "fair.proportion")
  FDI.avg.out <- geometric.mean(FDI.out$w) #psych require
}


#All resolution levels
FDI.list <- vector("list", 6)
for (f in 1:6) {FDI.list[[f]] <- lapply(AllCommsList, FDphylo, traits = ResMatrices6tr100sp[[f]])} 
FDI.unlist <- lapply(FDI.list, unlist)
FDI.ALL <- data.frame(t(matrix(unlist(FDI.unlist), nrow=length(FDI.unlist), byrow=T))) 
head(FDI.ALL)
colnames(FDI.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(FDI.ALL) #this is the final set of index values calculated for all levels of resolution.
FDI.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(FDI.ALL)

#SES
SES.FDI.list <- lapply(FDI.unlist, ses.calc)
SES.FDI.ALL <- data.frame(t(matrix(unlist(SES.FDI.list), nrow=length(SES.FDI.list), byrow=T)))
head(SES.FDI.ALL)
colnames(SES.FDI.ALL) <- c("TrContin", "Rounded", "Cat16", "Cat08", "Cat04", "Cat02")
head(SES.FDI.ALL) #this is the final set of index values calculated for all levels of resolution.
SES.FDI.ALL$proc <- c(rep("Rdm", 1000), rep("Div", 100), rep("Conv", 100)) #add a column to ID the assembly process used to create communities (random, covergence, divergence)
head(SES.FDI.ALL)

########################################################################################

#write our all 1200 SES values
write.csv(SES.FDis.ALL, "SESvalues_tr6p100nc1.4_FDis.csv")
write.csv(SES.MNTD.ALL, "SESvalues_tr6p100nc1.4_MNTD.csv")
write.csv(SES.FD.ALL, "SESvalues_tr6p100nc1.4_FD.csv")
write.csv(SES.MPD.ALL, "SESvalues_tr6p100nc1.4_MPD.csv")
write.csv(SES.FDI.ALL, "SESvalues_tr6p100nc1.4_FDI.csv")



#########
#PLOT
####


#switch from wide to long format so that each col is now per row

#FDis
SES.FDis.ALL.long <- gather(SES.FDis.ALL, "coarseness", "value", 1:6)
#only ef and comp
SES.FDis.ALL.long.subset <- SES.FDis.ALL.long[SES.FDis.ALL.long$proc != "Rdm",]


#MNTD
SES.MNTD.ALL.long <- gather(SES.MNTD.ALL, "coarseness", "value", 1:6)
#only ef and comp
SES.MNTD.ALL.long.subset <- SES.MNTD.ALL.long[SES.MNTD.ALL.long$proc != "Rdm",]


#FD dendrogram
SES.FD.ALL.long <- gather(SES.FD.ALL, "coarseness", "value", 1:6)
#only ef and comp
SES.FD.ALL.long.subset <- SES.FD.ALL.long[SES.FD.ALL.long$proc != "Rdm",]

#MPD
SES.MPD.ALL.long <- gather(SES.MPD.ALL, "coarseness", "value", 1:6)
#only ef and comp
SES.MPD.ALL.long.subset <- SES.MPD.ALL.long[SES.MPD.ALL.long$proc != "Rdm",]

#F distinctness(dendrogram-based)
SES.FDI.ALL.long <- gather(SES.FDI.ALL, "coarseness", "value", 1:6)
#only ef and comp
SES.FDI.ALL.long.subset <- SES.FDI.ALL.long[SES.FDI.ALL.long$proc != "Rdm",]
#######################



#find the median of each group
med.MNTD <- group_by(SES.MNTD.ALL.long, proc, coarseness) %>% summarise(med=median(value))
med.MNTD
med.MNTD.subset <- med.MNTD[med.MNTD$proc != "Rdm",]
med.MNTD.subset

med.MPD <- group_by(SES.MPD.ALL.long, proc, coarseness) %>% summarise(med=median(value))
med.MPD
med.MPD.subset <- med.MPD[med.MPD$proc != "Rdm",]
med.MPD.subset

med.FD <- group_by(SES.FD.ALL.long, proc, coarseness) %>% summarise(med=median(value))
med.FD
med.FD.subset <- med.FD[med.FD$proc != "Rdm",]
med.FD.subset

med.FDis <- group_by(SES.FDis.ALL.long, proc, coarseness) %>% summarise(med=median(value))
med.FDis
med.FDis.subset <- med.FDis[med.FDis$proc != "Rdm",]
med.FDis.subset

med.FDI <- group_by(SES.FDI.ALL.long, proc, coarseness) %>% summarise(med=median(value))
med.FDI
med.FDI.subset <- med.FDI[med.FDI$proc != "Rdm",]
med.FDI.subset

#Combine ALL
medians.all <- cbind.data.frame(med.MNTD, med.MPD$med, med.FD$med, med.FDis$med, med.FDI$med, rep(tr, length(med.MNTD)), rep(p, length(med.MNTD)), rep(nc/p, length(med.MNTD)))
colnames(medians.all)<- c("proc", "coarseness", "MNND", "MPD", "FD", "FDis", "FDIavg", "tr", "p", "nc")
medians.all

#writeout and save
SimRunVar  #call and copy name as part of file name
write.csv(medians.all, file = "Medians_tr6p100nc1.4.csv", row.names = TRUE)


#######################


###Clean plots

# Overlay with transparent density plot
F1 <- ggplot(SES.MNTD.ALL.long.subset, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.MNTD.subset, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
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
  scale_x_continuous(limits = c(-10, 10)) +
  

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
F2<-ggplot(SES.FD.ALL.long.subset, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.FD.subset, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
  scale_shape_manual(values = c(24, 24))+
  
  ###AXES###
  #axis labels
  xlab("SES-FD") +
  ylab("Density") +
  
  #formatting axis titles: 
  theme(axis.title.x= element_text(face="bold", colour="black", size=22),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=0,size=14, colour="black")) + #adjusts size and color of y axis title
  theme(axis.title.y = element_text(face="bold", colour="black", size=22), #adjusts angles, size, etc of labels for y axis
        axis.text.y  = element_text(angle=0, vjust=0.5, size=14, colour="black")) +
  
  #axis scale adjustments to make axes always uniform
  scale_x_continuous(limits = c(-10, 10)) +
  
  
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
F3<-ggplot(SES.FDis.ALL.long.subset, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.FDis.subset, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
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
  scale_x_continuous(limits = c(-10, 10)) +
  
  
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
F4 <- ggplot(SES.FDI.ALL.long.subset, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.FDI.subset, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
  scale_shape_manual(values = c(24, 24))+
  
  ###AXES###
  #axis labels
  xlab("SES-FDIavg") +
  ylab("Density") +
  
  #formatting axis titles: 
  theme(axis.title.x= element_text(face="bold", colour="black", size=22),#adjusts size and color of x axis title
        axis.text.x  = element_text(angle=0,size=14, colour="black")) + #adjusts size and color of y axis title
  theme(axis.title.y = element_text(face="bold", colour="black", size=22), #adjusts angles, size, etc of labels for y axis
        axis.text.y  = element_text(angle=0, vjust=0.5, size=14, colour="black")) +
  
  #axis scale adjustments to make axes always uniform
  scale_x_continuous(limits = c(-10, 10)) +
  
  
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
F5 <- ggplot(SES.MPD.ALL.long.subset, aes(x=value, fill = coarseness, color = proc)) + 
  geom_density(size = 1) +
  
  scale_fill_manual(values = c('#f7f7f7','#bdbdbd','#969696','#737373','#525252', '#000000'), aesthetics = "fill") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "white")) +
  #or brewer(palette = "Reds", aesthetics = "color") +
  geom_point(data=med.MPD.subset, aes(x=med, y=-0.05, size = 2, shape=proc, fill = coarseness)) + #adds points at medians
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
  scale_x_continuous(limits = c(-10, 10)) +

  
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
require(cowplot)
#desired order: [dispersion = MNND, FDis], [richness/vol = MPD, FD], [redundancy = FDIavg] 
indexPlots <- plot_grid(F1, F3, F5, F2, F4, ncol = 2, nrow = 3, labels = "auto") + draw_figure_label(label = SimRunVar, position = "bottom.left")
indexPlots

# to export as graphic
ppi <- 300
png("indexPlots_tr6p100nc1.4.png", width = 8*ppi, height = 6*ppi, res = ppi, bg = "transparent")
plot(indexPlots)
dev.off()

###if want just the legend components
#legend <- get_legend(F1)
#plot(legend)
#png("indexPlots_legend.png", width = 1*ppi, height = 3*ppi, res = ppi)
#plot(legend)
#dev.off()
###
###
###

#######################################################################################



#######
#Wilcoxon Signed Rank Tests of differences between mean SES of distributions
####

#1.MNTD
##divergence
wt.d.rd <- wilcox.test(SES.MNTD.ALL$Rounded[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) #competition

wt.d.16 <-wilcox.test(SES.MNTD.ALL$Cat16[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.MNTD.ALL$Cat08[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.MNTD.ALL$Cat04[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.MNTD.ALL$Cat02[SES.MNTD.ALL$proc=="Div"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.MNTD.ALL$Rounded[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) #EF

wt.c.16 <- wilcox.test(SES.MNTD.ALL$Cat16[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.MNTD.ALL$Cat08[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.MNTD.ALL$Cat04[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.MNTD.ALL$Cat02[SES.MNTD.ALL$proc=="Conv"], SES.MNTD.ALL$TrContin[SES.MNTD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_MNND <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

t(WilcoxTests_MNND)

####################

#2.FD
##divergence
wt.d.rd <- wilcox.test(SES.FD.ALL$Rounded[SES.FD.ALL$proc=="Div"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) #competition

wt.d.16 <-wilcox.test(SES.FD.ALL$Cat16[SES.FD.ALL$proc=="Div"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.FD.ALL$Cat08[SES.FD.ALL$proc=="Div"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.FD.ALL$Cat04[SES.FD.ALL$proc=="Div"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.FD.ALL$Cat02[SES.FD.ALL$proc=="Div"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.FD.ALL$Rounded[SES.FD.ALL$proc=="Conv"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) #EF

wt.c.16 <- wilcox.test(SES.FD.ALL$Cat16[SES.FD.ALL$proc=="Conv"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.FD.ALL$Cat08[SES.FD.ALL$proc=="Conv"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.FD.ALL$Cat04[SES.FD.ALL$proc=="Conv"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.FD.ALL$Cat02[SES.FD.ALL$proc=="Conv"], SES.FD.ALL$TrContin[SES.FD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_FD <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

t(WilcoxTests_FD)


########################
#3.FDis
##divergence
wt.d.rd <- wilcox.test(SES.FDis.ALL$Rounded[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) #competition

wt.d.16 <-wilcox.test(SES.FDis.ALL$Cat16[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.FDis.ALL$Cat08[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.FDis.ALL$Cat04[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.FDis.ALL$Cat02[SES.FDis.ALL$proc=="Div"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.FDis.ALL$Rounded[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) #EF

wt.c.16 <- wilcox.test(SES.FDis.ALL$Cat16[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.FDis.ALL$Cat08[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.FDis.ALL$Cat04[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.FDis.ALL$Cat02[SES.FDis.ALL$proc=="Conv"], SES.FDis.ALL$TrContin[SES.FDis.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_FDis <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

t(WilcoxTests_FDis)

######################################################
####
#4.MPD
##divergence
wt.d.rd <- wilcox.test(SES.MPD.ALL$Rounded[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) #competition

wt.d.16 <-wilcox.test(SES.MPD.ALL$Cat16[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.MPD.ALL$Cat08[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.MPD.ALL$Cat04[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.MPD.ALL$Cat02[SES.MPD.ALL$proc=="Div"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.MPD.ALL$Rounded[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) #EF

wt.c.16 <- wilcox.test(SES.MPD.ALL$Cat16[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.MPD.ALL$Cat08[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.MPD.ALL$Cat04[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.MPD.ALL$Cat02[SES.MPD.ALL$proc=="Conv"], SES.MPD.ALL$TrContin[SES.MPD.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_MPD <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

t(WilcoxTests_MPD)

####################
####################

#5.FDI
##divergence
wt.d.rd <- wilcox.test(SES.FDI.ALL$Rounded[SES.FDI.ALL$proc=="Div"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 ) #competition

wt.d.16 <-wilcox.test(SES.FDI.ALL$Cat16[SES.FDI.ALL$proc=="Div"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.08 <-wilcox.test(SES.FDI.ALL$Cat08[SES.FDI.ALL$proc=="Div"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.04 <-wilcox.test(SES.FDI.ALL$Cat04[SES.FDI.ALL$proc=="Div"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

wt.d.02 <-wilcox.test(SES.FDI.ALL$Cat02[SES.FDI.ALL$proc=="Div"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Div"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95 )

##convergence
wt.c.rd <- wilcox.test(SES.FDI.ALL$Rounded[SES.FDI.ALL$proc=="Conv"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95) #EF

wt.c.16 <- wilcox.test(SES.FDI.ALL$Cat16[SES.FDI.ALL$proc=="Conv"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.08 <- wilcox.test(SES.FDI.ALL$Cat08[SES.FDI.ALL$proc=="Conv"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.04 <- wilcox.test(SES.FDI.ALL$Cat04[SES.FDI.ALL$proc=="Conv"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

wt.c.02 <- wilcox.test(SES.FDI.ALL$Cat02[SES.FDI.ALL$proc=="Conv"], SES.FDI.ALL$TrContin[SES.FDI.ALL$proc=="Conv"], alternative = "two.sided", paired = T, conf.int = T, conf.level = 0.95)

#combine
WilcoxTests_FDI <- cbind.data.frame(wt.d.rd$p.value, wt.d.rd$estimate, wt.d.16$p.value, wt.d.16$estimate, wt.d.08$p.value, wt.d.08$estimate, wt.d.04$p.value, wt.d.04$estimate, wt.d.02$p.value, wt.d.02$estimate, wt.c.rd$p.value, wt.c.rd$estimate, wt.c.16$p.value, wt.c.16$estimate, wt.c.08$p.value, wt.c.08$estimate, wt.c.04$p.value, wt.c.04$estimate, wt.c.02$p.value, wt.c.02$estimate)

t(WilcoxTests_FDI)
#####################


##########################

#Combine ALL
WilcoxTests.all <- cbind.data.frame(c(rep("Div", 10), rep("Conv", 10)), c(rep("Rounded", 2), rep("Cat16", 2), rep("Cat08",2), rep("Cat04",2), rep("Cat02",2),rep("Rounded", 2), rep("Cat16", 2), rep("Cat08",2), rep("Cat04",2), rep("Cat02",2)) ,t(WilcoxTests_MNND), t(WilcoxTests_FDis), t(WilcoxTests_MPD), t(WilcoxTests_FD), t(WilcoxTests_FDI), rep(tr, length(WilcoxTests_MNND)), rep(p, length(WilcoxTests_MNND)), rep(nc/p, length(WilcoxTests_MNND)))

colnames(WilcoxTests.all)<- c("proc", "coarseness", "MNND", "FDis", "MPD", "FD", "FDIavg", "tr", "p", "nc")

WilcoxTests.all

#writeout and save
SimRunVar  #call and copy name as part of file name
write.csv(WilcoxTests.all, file = "WilcoxSRTest_tr6p100nc1.4_estimates.csv", row.names = TRUE)

