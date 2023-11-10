# ABOUT -------------------------------------------------------------------
## R code for 
## Components and delivery formats of cognitive-behavioral therapy for chronic insomnia: 
## a systematic review and component network meta-analysis
##
## This is the code for primary analysis.
##

# PREPARE -----------------------------------------------------------------
## clear settings

Sys.setenv(LANGUAGE="en_US.UTF-8")
rm(list=ls())   ###clears memory
cat("\014")   ###clears console
if(!is.null(dev.list())) dev.off() ###clears plots

## library
library(dplyr)
#library(openxlsx)
library(meta)
library(netmeta)
library(tidyr)

# RAW ---------------------------------------------------------------------
## load data
data <- read.csv("data_CBTICNMA.csv")

## set data
data$r <- as.numeric(data$r) 
data$n <- as.numeric(data$n)

# NMA ---------------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
# NMA
## Format trial data

data_nma <- data %>% 
  dplyr::select(c(study, treatment_group,treatment_component,r,n))%>% # select columns
  dplyr::group_by(study,treatment_group) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) #combine the same treatment_group arms within a trial
head(data_nma)

df_nma <- netmeta::pairwise(
  treat = treatment_group,
  event = r,
  n = n,
  sm = "OR",
  data = data_nma,
  studlab = study
)
head(df_nma)

# +ANALYZE -----------------------------------------------------------------
# for NMA run:
net_nma <- netmeta(TE, seTE, treat1, treat2, studlab,
                   data = df_nma, ref = "PE",
                   sm = "OR", fixed = FALSE, small = "undesirable")
head(net_nma)

# p-score (SUCRA for frequentist NMA)
netrank(x=net_nma, 
        small.values="bad",
        sort = TRUE)

## Network graph
size_net_nma <- tapply(data_nma$n, data_nma$treatment_group, sum)

netgraph(net_nma,
         seq = "optimal",
         col = "black", plastic = FALSE,
         points = TRUE, pch = 21, cex.points = 2.5,
         col.points = "black",
         bg.points = "gray",
         # thickness = "se.fixed", #alternative choice
         thickness = "number.of.studies",
         multiarm = FALSE,
         number.of.studies = TRUE)

# global
decomp.design(net_nma)

# local
netsplit(net_nma) # here you can check which (and how many) loops show inconsistencies between direct and indirect TEs.

## Forest plot
forest(net_nma, xlim = c(0.2, 5),sortvar = -Pscore,
       smlab = paste("Remission \n (Random Effects Model)"),
       label.left = "Favours PE    ",
       label.right = "  Favors Intervention",)

## begin: prediction interval
net_nma_split<-netsplit(net_nma)
forest(net_nma_split,fontsize=6,spacing=0.5,addrow.subgroups=FALSE)

png("prediction_interval.png", width = 400, height = 1750)    # prepare device
forest(net_nma_split,fontsize=8,spacing=0.75,show="all", 
       only.reference=FALSE,prediction=TRUE, direct=FALSE,indirect=FALSE)
dev.off()   

## end: prediction interval

## League table
# Create a CSV file with league table for random effects model
league_nma <- netleague(net_nma, digits = 2, bracket = "(", separator = " to ",seq = netrank(net_nma))
write.table(league_nma$random, file = "league_nma_CBTI.csv",
            row.names = FALSE, col.names = FALSE, sep = ",")

# CNMA --------------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
## Format trial data

data_cnma <- data %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma,
  studlab = study
)
head(df_cnma)


# +ANALYZE ----------------------------------------------------------------
# CNMA
netconnection(df_cnma$treat1, df_cnma$treat2, df_cnma$treat1studlab) #disconnected

dc1 <- discomb(TE, seTE, treat1, treat2, studlab, 
               inactive="dt",
               data = df_cnma, sm = "OR", fixed = FALSE)
head(dc1)
dc1

#network graph
netgraph(dc1,plastic = FALSE, thickness = FALSE,
         points = FALSE, cex.points = 0.01)

#forest plot for cNMA
forest(netcomplex(dc1, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.2, 5),digits = 2)

# combinations -----
# you can calculate the iOR of a certain combination of components
## example: 
netcomplex(dc1,"cr+th+sr+sc",nchar.comps=4)
netcomplex(dc1,"cr+sr",nchar.comps=4)

# you can calculate the OR of a certain combination against another
## example: CBTI against PE
netcomparison(dc1, "se+sd+cr+th+sr+sc+ns+ind+ff", "se+ns+ind+ff")
## example: best package against the most common package
netcomparison(dc1, "se+sd+cr+th+sr+sc+ns+ind+ff", "se+sd+cr+sr+sc+re+ns+ind+ff")