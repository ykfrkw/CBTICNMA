## R code for 
## Components and delivery formats of cognitive-behavioral therapy for chronic insomnia: 
## a systematic review and component network meta-analysis
## 
## The codes are for Bayesian analysis with penalized regression. 
##

# libraries -------------------------------------------------------------
library(rjags)
library(stringr)
library(netmeta)
library(ggplot2)
library(xlsx)
library(readxl)
library(MCMCvis)
library(tidyr)
library(tidyverse)
library(igraph)
library(openxlsx)
options(digits=2)

remove(list=ls())

# load data  --------------------------------------------------------------------------------------
#setwd("G:/My Drive/PROJECT/APPLIED PROJECTS/cNMA insomnia/bayesian analysis")
Sys.setenv(LANGUAGE="en_US.UTF-8")
d1=read.csv("data_CBTICNMA.csv")
d1$r <- as.numeric(d1$r) 
d1$n <- as.numeric(d1$n)

data_cnma <- d1 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 
data_cnma

# perform cNMA frequentist  --------------------------------------------------------------------------------------
df_cnma <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma,
  studlab = study
)

dc1 <- discomb(TE, seTE, treat1, treat2, studlab, inactive=NULL,
               data = df_cnma, sm = "OR", fixed = FALSE, reference.group = "w")
dc1


# set up data for Bayesian cNMA --------------------------------------------------------------------------------------
data_cnma
data_cnma$r=floor(data_cnma$r)
data_cnma$n=floor(data_cnma$n)
components=c("se", "sd", "cr", "th", "cw",  "sr", "sc", "re", "pi", 
             "w",  "ns", "he", "tg", "ind", "gp", "ff", "ae")
studies=unique(data_cnma$study)
Ns=length(studies)
Nc=length(components)
dcom=as.data.frame(matrix(0, nrow=length(data_cnma$study), ncol=Nc))
colnames(dcom)=components  
d2=as.data.frame(cbind(data_cnma, dcom))

for (j in 1:length(d2$study)){
  for(i in 1:Nc){
    d2[j, which(colnames(d2)==components[i])] = grepl( components[i], d2$treatment_component[j], fixed = TRUE)
  }}

d2=d2[order(d2$study),]
d2$arm=0
studies=unique(d2$study)
Ns=length(studies)
na=c()
for (i in 1:Ns){j=length(d2$arm[d2$study==studies[i]])
na=c(na,j)}
for(i in 1:Ns){  d2$arm[d2$study==studies[i]]=1:na[i]}
d3=d2[,c(which(colnames(d2)=="study"),which(colnames(d2)=="arm"),which(colnames(d2)=="n"))]
n <- spread(d3, arm, n)[,2:6] 
d3=d2[,c(which(colnames(d2)=="study"),which(colnames(d2)=="arm"),which(colnames(d2)=="r"))]
r <- spread(d3, arm, r)[,2:6] 

for(i in 1:Nc){
  d6=d2[,c(which(colnames(d2)=="study"),which(colnames(d2)=="arm"),which(colnames(d2)==components[i]))]
  assign( paste("c",i, sep=""), spread(d6, arm,components[i])[,2:6] )}

#  NO INTERACTION MODEL  --------------------------------------------------------------------------------------
model1.string <-  "
model {
for(i in 1:Ns) { 
w[i,1]<- 0
theta[i,1]<- 0                                             
for (k in 1:na[i]) {r[i,k] ~ dbin(p[i,k],n[i,k])}                                                   
logit(p[i,1])<- u[i]		                    
for (k in 2:na[i]) {
logit(p[i,k])<- u[i] + theta[i,k]

##distribution of random effects
theta[i,k] ~ dnorm(md[i,k],precd[i,k])

## accounting for correlation between effect sizes estimated in multi-arm trials				             
md[i,k]<- mean[i,k]+ sw[i,k]                                   
w[i,k]<- (theta[i,k]  - mean[i,k])          
sw[i,k]<- sum(w[i,1:(k-1)])/(k-1)
precd[i,k]<- prec *2*(k-1)/k  

##consistency equations
mean[i,k] <-A1[i,k]-B1[i]
A1[i,k]<- d[1]*(1-equals(c1[i,k],0)) + d[2]*(1-equals(c2[i,k],0))+
d[3]*(1-equals(c3[i,k],0)) + d[4]*(1-equals(c4[i,k],0)) +
d[5]*(1-equals(c5[i,k],0)) + d[6]*(1-equals(c6[i,k],0)) + 
d[7]*(1-equals(c7[i,k],0)) + d[8]*(1-equals(c8[i,k],0)) + 
d[9]*(1-equals(c9[i,k],0)) + d[10]*(1-equals(c10[i,k],0))+
d[11]*(1-equals(c11[i,k],0)) + d[12]*(1-equals(c12[i,k],0))+
d[13]*(1-equals(c13[i,k],0)) + d[14]*(1-equals(c14[i,k],0))+
d[15]*(1-equals(c15[i,k],0)) + d[16]*(1-equals(c16[i,k],0)) +
d[17]*(1-equals(c17[i,k],0)) 
}

B1[i]<- 
d[1]*(1-equals(c1[i,1],0))  + d[2]*(1-equals(c2[i,1],0))+
d[3]*(1-equals(c3[i,1],0))  + d[4]*(1-equals(c4[i,1],0)) +
d[5]*(1-equals(c5[i,1],0))  + d[6]*(1-equals(c6[i,1],0)) + 
d[7]*(1-equals(c7[i,1],0))  + d[8]*(1-equals(c8[i,1],0)) + 
d[9]*(1-equals(c9[i,1],0))  + d[10]*(1-equals(c10[i,1],0)) + 
d[11]*(1-equals(c11[i,1],0))  + d[12]*(1-equals(c12[i,1],0)) +
d[13]*(1-equals(c13[i,1],0))  + d[14]*(1-equals(c14[i,1],0)) + 
d[15]*(1-equals(c15[i,1],0)) + d[16]*(1-equals(c16[i,1],0)) +
d[17]*(1-equals(c17[i,1],0))
}

##prior distribution for log-odds in baseline arm of study i
for (i in 1:Ns) {	u[i] ~ dnorm(0,.01)	}

## informative prior distribution for heterogeneity	
prec<-1/tau
tau~dlnorm(-1.67,inv.sd)
inv.sd<-1/(1.472*1.472)
	

## example for checking results - not needed for analysis
example1<-exp(d[2]+d[3]+d[4]+d[5]+d[8] - d[1])

## prior distribution for basic parameters		
for(k in 1:Nc) {d[k] ~ dnorm(0,.01)}
## calculate odd ratios
for(k in 1:Nc) {ORd[k]<- exp(d[k])}}
"

model1.spec<-textConnection(model1.string) 
data <- list(r=r,n=n, Ns=Ns, Nc=Nc, c1=c1, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6, c7=c7, c8=c8, 
             c9=c9, c10=c10, c11=c11, c12=c12, c13=c13, c14=c14, c15=c15, c16=c16, c17=c17, na=na)
jags.m=0
jags.m <- jags.model(model1.spec, data = data, n.chains =3, n.adapt = 2000)

params <- c("tau","d", "ORd") 
closeAllConnections()
samps<- coda.samples(jags.m, params, n.iter =20000)

MCMCtrace(samps,pdf = FALSE, params = "d") 

A1= MCMCsummary(samps)
rownames(A1)[1:Nc]=components
rownames(A1)[(Nc+1):(2*Nc)]=paste("logOR -", components)
round(A1[,c(4,3,5)], digits=3)

A1$results=paste(format(round(A1$`50%`,2), nsmall=2), "[", 
                 format(round(A1$`2.5%`,2), nsmall=2), ";", 
                 format(round(A1$`97.5%`,2), nsmall=2), "]", sep="")



cbind(components, A1$results[1:Nc])

samps.all.additive=as.data.frame(rbind(samps[[1]], samps[[2]], samps[[3]]))
options(digits=3)

## example for estimating results of combinations
# estimate OR between 2+3+4+5+8  vs  1
samps.all.additive$example<-with(samps.all.additive, exp(`d[2]`+`d[3]`+`d[4]`+`d[5]`+`d[8]` - `d[1]`))
quantile(samps.all.additive$example, c(.5, 0.025, 0.975))


# define useful functions  --------------------------------------------------------------------------------------

# this function provides the indicator of the interaction term between 2 components
place=function(i, j){
  if(i<j){pl=(2*Nc-i)*(i-1)/2+j-i}
  if(j<i){
    j1=j; i1=i;i=j1; j=i1 
    pl=(2*Nc-i)*(i-1)/2+j-i}
  return(pl)}
# e.g. place(2,3)

# this function provides the components comprising a specific interaction term
which.place=function(n){
  for (k in 1:(Nc-1)){
    for (l in (k+1):Nc){
      if(place(k,l)==n){pl=(paste(k,l))
      compo=paste(components[k],"-", components[l])
      }
      
    }}
  return(c(pl, compo))}
# e.g. which.place(23)

# this function provides all interactions needed for a combination of components
all.interactions=function(combination){
  nc=length(combination)
  combination=sort(combination)
  interactions=c()
  for ( i in 1:(nc-1)){
    for (j in (i+1):nc){
      interactions=c(interactions,place(combination[i],combination[j]))      }    }
  return(interactions)}

## eg all.interactions(c(2,4,5,11,12,15,16))


# define interactions -----------------------------------------------------
Ninter=Nc*(Nc-1)/2
interactions=array(NA, dim=c(Ns, 5, Ninter))
for (i in 1: Ns){
  for (j in 1:na[i]){
    place1=0
    for (k in 1:(Nc-1)){
      for (l in (k+1):Nc){
        place1=place1+1
        interactions[i, j, place1]=
          eval(parse(text=paste("c", k, sep="")))[i, j]*eval(parse(text=paste("c", l, sep="")))[i, j] }     }   }}



# SSVS INTERACTION MODEL, all interactions equiprobable -----------------------------------------------------
model1.SSVS <-  "
model {
for(i in 1:Ns) { 
w[i,1]<- 0
theta[i,1]<- 0                                             
for (k in 1:na[i]) {r[i,k] ~ dbin(p[i,k],n[i,k])}                                                   
logit(p[i,1])<- u[i]		                    
for (k in 2:na[i]) {
logit(p[i,k])<- u[i] + theta[i,k]

##distribution of random effects
theta[i,k] ~ dnorm(md[i,k],precd[i,k])

## accounting for correlation between effect sizes estimated in multi-arm trials				             
md[i,k]<- mean[i,k]+ sw[i,k]                                   
w[i,k]<- (theta[i,k]  - mean[i,k])          
sw[i,k]<- sum(w[i,1:(k-1)])/(k-1)
precd[i,k]<- prec *2*(k-1)/k  

##consistency equations
mean[i,k] <-A1[i,k]-B1[i]
A1[i,k]<-
d[1]*(1-equals(c1[i,k],0)) + d[2]*(1-equals(c2[i,k],0))+
d[3]*(1-equals(c3[i,k],0)) + d[4]*(1-equals(c4[i,k],0)) +
d[5]*(1-equals(c5[i,k],0)) + d[6]*(1-equals(c6[i,k],0)) + 
d[7]*(1-equals(c7[i,k],0)) + d[8]*(1-equals(c8[i,k],0)) + 
d[9]*(1-equals(c9[i,k],0)) + d[10]*(1-equals(c10[i,k],0))+
d[11]*(1-equals(c11[i,k],0)) + d[12]*(1-equals(c12[i,k],0)) +
d[13]*(1-equals(c13[i,k],0)) + d[14]*(1-equals(c14[i,k],0))+
d[15]*(1-equals(c15[i,k],0)) + d[16]*(1-equals(c16[i,k],0)) +
d[17]*(1-equals(c17[i,k],0))+
inprod(gamma[], interactions[i,k,]) }

B1[i]<-
d[1]*(1-equals(c1[i,1],0))  + d[2]*(1-equals(c2[i,1],0))+
d[3]*(1-equals(c3[i,1],0))  + d[4]*(1-equals(c4[i,1],0)) +
d[5]*(1-equals(c5[i,1],0))  + d[6]*(1-equals(c6[i,1],0)) + 
d[7]*(1-equals(c7[i,1],0))  + d[8]*(1-equals(c8[i,1],0)) + 
d[9]*(1-equals(c9[i,1],0))  + d[10]*(1-equals(c10[i,1],0)) + 
d[11]*(1-equals(c11[i,1],0))+ d[12]*(1-equals(c12[i,1],0))+ 
d[13]*(1-equals(c13[i,1],0))  + d[14]*(1-equals(c14[i,1],0)) + 
d[15]*(1-equals(c15[i,1],0)) + d[16]*(1-equals(c16[i,1],0)) +
d[17]*(1-equals(c17[i,1],0))+
inprod(gamma[], interactions[i,1,])  }

##prior distribution for log-odds in baseline arm of study i
for (i in 1:Ns) {	u[i] ~ dnorm(0,.01)	}

##prior distribution for heterogeneity	
prec<-1/tau
tau~dlnorm(-1.67,inv.sd)
inv.sd<-1/(1.472*1.472)


## SSVS
## to find the interactions for specific components use the place() function,
## e.g. run place(2,3) to find the interaction term between components 2 and 3
  for(k in 1:Ninter){
    IndA[k] ~ dcat(Pind[])
    Ind[k] <- IndA[k] - 1
    gamma[k] ~ dnorm(0, tauCov[IndA[k]])  }
  
zeta <- pow(eta, -2)
eta ~ dnorm(0,1000)I(0,)
tauCov[1] <- zeta
tauCov[2] <- zeta * 0.01  # g = 100

  Pind[1] <- 0.5 #P(I_j=1)= 0.5
  Pind[2] <- 0.5 
##prior distribution for basic parameters		
for(k in 1:Nc) {d[k] ~ dnorm(0,.01)}
for(k in 1:Nc) {ORd[k]<- exp(d[k])}

### example comparison for cheching results
### use the function all.interactions to find the which gamma to include
### e.g. run all.interactions(c(2,3,4))
####
#example1<-exp(d[2]+d[3]+d[4] +
#gamma[17]+gamma[18]+gamma[32]-d[1]) 
}
"
model1.spec<-textConnection(model1.SSVS) 
data <- list(r=r,n=n, Ns=Ns, Nc=Nc, c1=c1, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6, c7=c7, c8=c8, 
             c9=c9, c10=c10, c11=c11, c12=c12,c13=c13, c14=c14, c15=c15, c16=c16, c17=c17,  na=na, interactions=interactions, Ninter=Ninter  )
jags.m=0
jags.m <- jags.model(model1.spec, data = data, n.chains =4, n.adapt = 2000)

params <- c("tau", "ORd", "gamma", "Ind", "d", "eta") 
closeAllConnections()
samps2<- coda.samples(jags.m, params, n.iter =20000)

A2= MCMCsummary(samps2)
rownames(A2)[which(rownames(A2)=="d[1]"):which(rownames(A2)=="d[17]")]=paste("logOR -", components)
rownames(A2)[which(rownames(A2)=="ORd[1]"):which(rownames(A2)=="ORd[17]")]=components 

MCMCtrace(samps2,pdf = FALSE, params = "d") 
meanInd=c() 
for(i in 1:(Nc*(Nc-1)/2)){
  meanInd=c(meanInd, mean(rbind(samps2[[1]],samps2[[2]], samps2[[3]], samps2[[4]]) [, paste("Ind[", i,"]", sep="")]))}

medianG=c() 
for(i in 1:(Nc*(Nc-1)/2)){  medianG=c(medianG, median(rbind(samps2[[1]],samps2[[2]], samps2[[3]], samps2[[4]])[, paste("gamma[", i,"]", sep="")]))}
nmax=which(medianG==max(medianG))
nmin=which(medianG==min(medianG))
medianG[nmax]
medianG[nmin]
meanInd[nmax]
meanInd[nmin]

meanInd[order(medianG)[1]] 
meanInd[order(medianG)[2]]

which.place(nmax)
which.place(nmin)

dev.off()
ggplot(data.frame("gamma"=medianG, "Ind"=meanInd), aes(x=medianG, y=meanInd)) + 
  geom_point() + labs(x = "Estimated coefficent of interaction terms (SSVS)")+
  labs(y = "Variable selection frequency for interaction terms (SSVS)")


which.place(nmin) # "sd - th"
which.place(nmax) # "th - ind"

which.place(order(medianG)[1])
which.place(order(medianG)[4])

A2$results=paste(format(round(A2$`50%`,2), nsmall=2), "[", 
                 format(round(A2$`2.5%`,2), nsmall=2), ";", 
                 format(round(A2$`97.5%`,2), nsmall=2), "]", sep="")
options(max.print = 3000)
A2

### these are the strongest interactions found
# negative interactions
which.place(order(medianG)[1])   #"sd - th"
medianG[18]
exp(medianG[18])   #0.67
which.place(order(medianG)[2])   #"re - ae"
medianG[100]
exp(medianG[100])   #0.87
which.place(order(medianG)[3])   #"sr - sc"
medianG[32]
exp(medianG[32])   #0.87

# positive interactions
which.place(order(-medianG)[1])   #"th - ind"
medianG[55]   #"th - ind"
exp(medianG[55])   

## estimate treatment effects between combinations

samps.all.SSVS=as.data.frame(rbind(samps2[[1]], samps2[[2]], samps2[[3]],samps2[[4]]))
options(digits=3)


## example for estimating results of combinations
# estimate OR between 2+3+4  vs  1
all.interactions(c(2,3,4))
samps.all.SSVS$example<-with(samps.all.SSVS, exp(`d[2]`+`d[3]`+`d[4]`  #main effects
                                                 +`gamma[17]`+`gamma[18]` +`gamma[32]` #interactions
                                                 - `d[1]`))
quantile(samps.all.SSVS$example, c(.5, 0.025, 0.975))


# SSVS INTERACTION MODEL, with prior for I_k -----------------------------------------------------
active.components=c("se", "sd", "cr", "th", "sr", "sc", "re", "ff")
included=all.interactions(which(components %in% active.components))  # these interactions we want to have included
not.included=(1:(Nc*(Nc-1)/2))[!(1:(Nc*(Nc-1)/2) %in% included)]

model2.SSVS <-  "
################### AD part
model {
for(i in 1:Ns) { 
w[i,1]<- 0
theta[i,1]<- 0                                             
for (k in 1:na[i]) {r[i,k] ~ dbin(p[i,k],n[i,k])}                                                   
logit(p[i,1])<- u[i]		                    
for (k in 2:na[i]) {
logit(p[i,k])<- u[i] + theta[i,k]

##distribution of random effects
theta[i,k] ~ dnorm(md[i,k],precd[i,k])

## accounting for correlation between effect sizes estimated in multi-arm trials				             
md[i,k]<- mean[i,k]+ sw[i,k]                                   
w[i,k]<- (theta[i,k]  - mean[i,k])          
sw[i,k]<- sum(w[i,1:(k-1)])/(k-1)
precd[i,k]<- prec *2*(k-1)/k  

##consistency equations
mean[i,k] <-A1[i,k]-B1[i]
A1[i,k]<-
d[1]*(1-equals(c1[i,k],0)) + d[2]*(1-equals(c2[i,k],0))+
d[3]*(1-equals(c3[i,k],0)) + d[4]*(1-equals(c4[i,k],0)) +
d[5]*(1-equals(c5[i,k],0)) + d[6]*(1-equals(c6[i,k],0)) + 
d[7]*(1-equals(c7[i,k],0)) + d[8]*(1-equals(c8[i,k],0)) + 
d[9]*(1-equals(c9[i,k],0)) + d[10]*(1-equals(c10[i,k],0))+
d[11]*(1-equals(c11[i,k],0)) + d[12]*(1-equals(c12[i,k],0)) +
d[13]*(1-equals(c13[i,k],0)) + d[14]*(1-equals(c14[i,k],0))+
d[15]*(1-equals(c15[i,k],0)) + d[16]*(1-equals(c16[i,k],0)) +
d[17]*(1-equals(c17[i,k],0))+
inprod(gamma[], interactions[i,k,]) }

B1[i]<-
d[1]*(1-equals(c1[i,1],0))  + d[2]*(1-equals(c2[i,1],0))+
d[3]*(1-equals(c3[i,1],0))  + d[4]*(1-equals(c4[i,1],0)) +
d[5]*(1-equals(c5[i,1],0))  + d[6]*(1-equals(c6[i,1],0)) + 
d[7]*(1-equals(c7[i,1],0))  + d[8]*(1-equals(c8[i,1],0)) + 
d[9]*(1-equals(c9[i,1],0))  + d[10]*(1-equals(c10[i,1],0)) + 
d[11]*(1-equals(c11[i,1],0))+ d[12]*(1-equals(c12[i,1],0))+ 
d[13]*(1-equals(c13[i,1],0))  + d[14]*(1-equals(c14[i,1],0)) + 
d[15]*(1-equals(c15[i,1],0)) + d[16]*(1-equals(c16[i,1],0)) +
d[17]*(1-equals(c17[i,1],0))+
inprod(gamma[], interactions[i,1,])  }


##prior distribution for log-odds in baseline arm of study i
for (i in 1:Ns) {	u[i] ~ dnorm(0,.01)	}

##prior distribution for heterogeneity	
prec<-1/tau
tau~dlnorm(-1.67,inv.sd)
inv.sd<-1/(1.472*1.472)


## SSVS
## equiprobable interaction terms
  for(k in c(4, 8, 9, 10, 11, 12, 13, 14, 16, 19, 23, 24, 25, 
26, 27, 28, 29, 31, 33, 37, 38, 39, 40, 41, 42, 43, 
45, 46, 50, 51, 52, 53, 54, 55, 56, 58, 59, 60, 61, 
62, 63, 64, 65, 66, 67, 68, 69, 70, 73, 74, 75, 76, 
77, 78, 79, 81, 83, 84, 85, 86, 87, 88, 89, 91, 92, 
93, 94, 95, 96, 97, 98, 100, 101, 102, 103, 104, 105, 
106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 
117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 
128, 129, 130, 131, 132, 133, 134, 135, 136)){

    gamma[k] <-0 }
  
zeta <- pow(eta, -2)
eta ~ dnorm(0,1000)I(0,)
tauCov[1] <- zeta
tauCov[2] <- zeta * 0.01  # g = 100


## interaction terms with prior information
 for(k in c(1, 2, 3, 5, 6, 7, 15, 17, 18, 20, 21, 22, 30, 32, 34, 35, 36, 
44, 47, 48, 49, 57, 71, 72, 80, 82, 90, 99)){
    IndA[k] ~ dcat(Pind1[])
    Ind[k] <- IndA[k] - 1
    gamma[k] ~ dnorm(0, tauCov[IndA[k]])  }
  Pind1[1] <-0.2
  Pind1[2] <- 0.8 

##prior distribution for basic parameters		
for(k in 1:Nc) {d[k] ~ dnorm(0,.01)}
for(k in 1:Nc) {ORd[k]<- exp(d[k])}

}
"
model1.spec<-textConnection(model2.SSVS) 
data <- list(r=r,n=n, Ns=Ns, Nc=Nc, c1=c1, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6, c7=c7, c8=c8, 
             c9=c9, c10=c10, c11=c11, c12=c12,c13=c13, c14=c14, c15=c15, c16=c16, c17=c17,  na=na, interactions=interactions )
jags.m=0
jags.m <- jags.model(model1.spec, data = data, n.chains =4, n.adapt = 2000)

params <- c("tau", "ORd", "gamma", "Ind", "d") 
closeAllConnections()
samps3<- coda.samples(jags.m, params, n.iter =20000)

A3= MCMCsummary(samps3)
rownames(A3)[which(rownames(A3)=="d[1]"):which(rownames(A3)=="d[17]")]=paste("logOR -", components)
rownames(A3)[which(rownames(A3)=="ORd[1]"):which(rownames(A3)=="ORd[17]")]=components 
A3

MCMCtrace(samps3,pdf = FALSE, params = "d") 
meanInd2=c() 
for(i in included){
  meanInd2=c(meanInd2, mean(rbind(samps3[[1]],samps3[[2]], samps3[[3]], samps3[[4]]) [, paste("Ind[", i,"]", sep="")]))}
medianG2=c() 
for(i in included){
  medianG2=c(medianG2, median(rbind(samps3[[1]],samps3[[2]], samps3[[3]], samps3[[4]])[, paste("gamma[", i,"]", sep="")]))}


nmax=which(medianG2==max(medianG2))
nmin=which(medianG2==min(medianG2))
medianG2[nmax]
medianG2[nmin]

meanInd2[nmax]
meanInd2[nmin]

which.place(nmax) #se - ind
which.place(nmin) #se - w

dev.off()
ggplot(data.frame("gamma"=medianG2, "Ind"=meanInd2), aes(x=medianG2, y=meanInd2)) + 
  geom_point() + labs(x = "Estimated coefficent of interaction terms (SSVS)")+
  labs(y = "Variable selection frequency for interaction terms (SSVS)")


which.place(nmin) 
which.place(nmax) 
which.place(order(medianG2)[1])

A3$results=paste(format(round(A3$`50%`,2), nsmall=2), "[", 
                 format(round(A3$`2.5%`,2), nsmall=2), ";", 
                 format(round(A3$`97.5%`,2), nsmall=2), "]", sep="")




### these are the strongest interactions found
# negative interactions
which.place(included[order(medianG2)[1]])   #"sd - th"
medianG2[9]   #-0.437
exp(medianG2[9])    #0.646


which.place(included[order(medianG2)[2]])   #"sr - sc"
medianG2[23]   #-0.381
exp(medianG2[23])   #0.683

which.place(included[order(medianG2)[3]])   #"cr - th"
medianG2[14]   #-0.249
exp(medianG2[14])   #0.78

# positive interactions
which.place(included[order(-medianG2)[1]])
medianG2[13]

which.place(included[order(-medianG2)[2]])
medianG2[27]


## estimate treatment effects between combinations

samps.all.SSVS2=as.data.frame(rbind(samps3[[1]], samps3[[2]], samps3[[3]],samps3[[4]]))


## example for estimating results of combinations
# estimate OR between 2+3+4  vs  1
all.interactions(c(2,3,4))
samps.all.SSVS2$example<-with(samps.all.SSVS2, exp(`d[2]`+`d[3]`+`d[4]`  #main effects
                                                   +`gamma[17]`+`gamma[18]` +`gamma[32]` #interactions
                                                   - `d[1]`))
quantile(samps.all.SSVS2$example, c(.5, 0.025, 0.975))



# estimate OR between cr+tw  vs  w   (3+4-10)
all.interactions(c(3,4))
samps.all.SSVS2$example2<-with(samps.all.SSVS2, exp(`d[3]`+`d[4]`  #main effects
                                                    +`gamma[32]` #interactions
                                                    - `d[1]`))
quantile(samps.all.SSVS2$example2, c(.5, 0.025, 0.975))

# estimate OR between tw  vs  w   (4-10)
samps.all.SSVS2$example3<-with(samps.all.SSVS2, exp(`d[4]` - `d[1]`))
quantile(samps.all.SSVS2$example3, c(.5, 0.025, 0.975))



# summarize results -------------------------------------------------------
summ=cbind(
  A1[c(which(rownames(A1)=="se"):which(rownames(A1)=="ae"), which(rownames(A1)=="tau")),c("results")], 
  A2[c(which(rownames(A2)=="se"):which(rownames(A2)=="ae"), which(rownames(A2)=="tau")),c("results")],
  A3[c(which(rownames(A3)=="se"):which(rownames(A3)=="ae"), which(rownames(A3)=="tau")),c("results")]
)
rownames(summ)=c(components, "tau")
colnames(summ)=c("no interactions", "SSVS (equiprobable interactions)", "SSVS (informative priors about interactions")

write.xlsx2(summ, file="cNMA of CBT-I results.xlsx")



# estimating results of combinations ------

##### ind f2f full package vs ind f2f sleep hygiene
###ns+ff+ind+se+sd+cr+th+sr+sc vs ns+ff+ind+se (11+16+14+1+2+3+4+6+7  vs  11+16+14+1)
# the additive model
samps.all.additive$example<-with(samps.all.additive, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[6]`+`d[7]` - `d[1]`))
quantile(samps.all.additive$example, c(.5, 0.025, 0.975))
# the full interaction model
all.interactions(c(11,16,14,1,2,3,4,6,7))
all.interactions(c(11,16,14,1))
samps.all.SSVS$example<-with(samps.all.SSVS, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[6]`+`d[7]`  #main effects
                                                 +`gamma[1]`+`gamma[2]`+`gamma[3]`+`gamma[5]`+`gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[18]`+`gamma[20]`
                                                 +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                 +`gamma[32]`+`gamma[34]` +`gamma[35]`+`gamma[39]`+`gamma[42]`+`gamma[44]`+`gamma[47]`+`gamma[48]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                 +`gamma[71]`+`gamma[75]`+`gamma[78]`+`gamma[80]`+`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                 - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS$example, c(.5, 0.025, 0.975))
# the partially interactive model
samps.all.SSVS2$example<-with(samps.all.SSVS2, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[6]`+`d[7]`  #main effects
                                                   +`gamma[1]`+`gamma[2]`+`gamma[3]`+`gamma[5]`+`gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[18]`+`gamma[20]`
                                                   +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                   +`gamma[32]`+`gamma[34]` +`gamma[35]`+`gamma[39]`+`gamma[42]`+`gamma[44]`+`gamma[47]`+`gamma[48]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                   +`gamma[71]`+`gamma[75]`+`gamma[78]`+`gamma[80]`+`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                   - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS2$example, c(.5, 0.025, 0.975))



###scなし　ns+ff+ind+se+sd+cr+th+sr vs ns+ff+ind+se (11+16+14+1+2+3+4+6  vs  11+16+14+1)
# the additive model
samps.all.additive$example<-with(samps.all.additive, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[6]` - `d[1]`))
quantile(samps.all.additive$example, c(.5, 0.025, 0.975))
# the full interaction model
all.interactions(c(11,16,14,1,2,3,4,6))
all.interactions(c(11,16,14,1))
samps.all.SSVS$example<-with(samps.all.SSVS, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[6]`  #main effects
                                                 +`gamma[1]`+`gamma[2]`+`gamma[3]`+`gamma[5]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[18]`+`gamma[20]`
                                                 +`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                 +`gamma[32]`+`gamma[34]` +`gamma[39]`+`gamma[42]`+`gamma[44]`+`gamma[47]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                 +`gamma[75]`+`gamma[78]`+`gamma[80]`  #interactions
                                                 - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS$example, c(.5, 0.025, 0.975))
# the partially interactive model
samps.all.SSVS2$example<-with(samps.all.SSVS2, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[6]`  #main effects
                                                   +`gamma[1]`+`gamma[2]`+`gamma[3]`+`gamma[5]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[18]`+`gamma[20]`
                                                   +`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                   +`gamma[32]`+`gamma[34]` +`gamma[39]`+`gamma[42]`+`gamma[44]`+`gamma[47]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                   +`gamma[75]`+`gamma[78]`+`gamma[80]`  #interactions
                                                   - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS2$example, c(.5, 0.025, 0.975))


###srなし　ns+ff+ind+se+sd+cr+th+sc vs ns+ff+ind+se (11+16+14+1+2+3+4+7  vs  11+16+14+1)
# the additive model
samps.all.additive$example<-with(samps.all.additive, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[7]` - `d[1]`))
quantile(samps.all.additive$example, c(.5, 0.025, 0.975))
# the full interaction model
all.interactions(c(11,16,14,1,2,3,4,7))
all.interactions(c(11,16,14,1))
samps.all.SSVS$example<-with(samps.all.SSVS, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[7]`  #main effects
                                                 +`gamma[1]`+`gamma[2]`+`gamma[3]`+ `gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[18]`
                                                 +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                 +`gamma[32]` +`gamma[35]`+`gamma[39]`+`gamma[42]`+`gamma[44]`+`gamma[48]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                 +`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                 - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS$example, c(.5, 0.025, 0.975))
# the partially interactive model
samps.all.SSVS2$example<-with(samps.all.SSVS2, exp(`d[1]`+`d[2]`+`d[3]`+`d[4]`+`d[7]`  #main effects
                                                   +`gamma[1]`+`gamma[2]`+`gamma[3]`+ `gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[18]`
                                                   +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                   +`gamma[32]` +`gamma[35]`+`gamma[39]`+`gamma[42]`+`gamma[44]`+`gamma[48]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                   +`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                   - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS2$example, c(.5, 0.025, 0.975))

###thなし ns+ff+ind+se+sd+cr+sr+sc vs ns+ff+ind+se (11+16+14+1+2+3+4+6+7  vs  11+16+14+1)
# the additive model
samps.all.additive$example<-with(samps.all.additive, exp(`d[1]`+`d[2]`+`d[3]`+`d[6]`+`d[7]` - `d[1]`))
quantile(samps.all.additive$example, c(.5, 0.025, 0.975))
# the full interaction model
all.interactions(c(11,16,14,1,2,3,6,7)) # this shows the interactions of the components
all.interactions(c(11,16,14,1))
samps.all.SSVS$example<-with(samps.all.SSVS, exp(`d[1]`+`d[2]`+`d[3]`+`d[6]`+`d[7]`  #main effects
                                                 +`gamma[1]`+`gamma[2]`+`gamma[5]`+`gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[20]`
                                                 +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                 +`gamma[34]` +`gamma[35]`+`gamma[39]`+`gamma[42]`+`gamma[44]`
                                                 +`gamma[71]`+`gamma[75]`+`gamma[78]`+`gamma[80]`+`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                 - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS$example, c(.5, 0.025, 0.975))
# the partially interactive model
samps.all.SSVS2$example<-with(samps.all.SSVS2, exp(`d[1]`+`d[2]`+`d[3]`+`d[6]`+`d[7]`  #main effects
                                                   +`gamma[1]`+`gamma[2]`+`gamma[5]`+`gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[17]`+`gamma[20]`
                                                   +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                   +`gamma[34]` +`gamma[35]`+`gamma[39]`+`gamma[42]`+`gamma[44]`
                                                   +`gamma[71]`+`gamma[75]`+`gamma[78]`+`gamma[80]`+`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                   - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS2$example, c(.5, 0.025, 0.975))


###crなし ns+ff+ind+se+sd+th+sr+sc vs ns+ff+ind+se (11+16+14+1+2+3+4+6+7  vs  11+16+14+1)
# the additive model
samps.all.additive$example<-with(samps.all.additive, exp(`d[1]`+`d[2]`+`d[4]`+`d[6]`+`d[7]` - `d[1]`))
quantile(samps.all.additive$example, c(.5, 0.025, 0.975))
# the full interaction model
all.interactions(c(11,16,14,1,2,4,6,7))
all.interactions(c(11,16,14,1))
samps.all.SSVS$example<-with(samps.all.SSVS, exp(`d[1]`+`d[2]`+`d[4]`+`d[6]`+`d[7]`  #main effects
                                                 +`gamma[1]`+ `gamma[3]`+`gamma[5]`+`gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[18]`+`gamma[20]`
                                                 +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                 +`gamma[47]`+`gamma[48]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                 +`gamma[71]`+`gamma[75]`+`gamma[78]`+`gamma[80]`+`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                 - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS$example, c(.5, 0.025, 0.975))
# the partially interactive model
samps.all.SSVS2$example<-with(samps.all.SSVS2, exp(`d[1]`+`d[2]`+`d[4]`+`d[6]`+`d[7]`  #main effects
                                                   +`gamma[1]`+ `gamma[3]`+`gamma[5]`+`gamma[6]`+`gamma[10]`+`gamma[13]`+`gamma[15]`+`gamma[18]`+`gamma[20]`
                                                   +`gamma[21]`+`gamma[25]`+`gamma[28]`+`gamma[30]`
                                                   +`gamma[47]`+`gamma[48]`+`gamma[52]`+`gamma[55]`+`gamma[57]`
                                                   +`gamma[71]`+`gamma[75]`+`gamma[78]`+`gamma[80]`+`gamma[85]`+`gamma[88]`+`gamma[90]`  #interactions
                                                   - `d[1]`-`gamma[10]`-`gamma[13]`-`gamma[15]`))
quantile(samps.all.SSVS2$example, c(.5, 0.025, 0.975))
