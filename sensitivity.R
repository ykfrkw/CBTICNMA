# ABOUT -------------------------------------------------------------------
## R code for 
## Components and delivery formats of cognitive-behavioral therapy for chronic insomnia: 
## a systematic review and component network meta-analysis
## 
## The codes are for sensitivity analyses. 
## Please load the dataset and libarary using the primary analysis code first.
##

#1. Excluding studies without formal diagnosis of insomnia ------
## Format trial data
data_s1<- data[(data$Dx==1),]

data_cnma_s1 <- data_s1 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_s1 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_s1,
  studlab = study
)

# CNMA
dc1_s1 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
               data = df_cnma_s1, sm = "OR", fixed = FALSE)
dc1_s1

forest(netcomplex(dc1_s1, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#2. Excluding studies focusing on patients with comorbidities (both physical and psychological) ------
## Format trial data
data_s2 <- data[(data$Comorbidities<1.5),] 

data_cnma_s2 <- data_s2 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_s2 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_s2,
  studlab = study
)

# CNMA
dc1_s2 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                  data = df_cnma_s2, sm = "OR", fixed = FALSE)
dc1_s2

forest(netcomplex(dc1_s2, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#3. Excluding studies with overall high dropout rate (20% or more) -------
## Format trial data
data_s3<- data[(data$DropOutRate<20),] 

data_cnma_s3 <- data_s3 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_s3 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_s3,
  studlab = study
)

# CNMA
dc1_s3 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                  data = df_cnma_s3, sm = "OR", fixed = FALSE)
dc1_s3

forest(netcomplex(dc1_s3, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#4. Excluding studies at high overall risk of bias -------
## Format trial data
data_s4<- data[!(data$Overall=="H"),] 

data_cnma_s4 <- data_s4 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_s4 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_s4,
  studlab = study
)

# CNMA
dc1_s4 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                  data = df_cnma_s4, sm = "OR", fixed = FALSE)
dc1_s4

forest(netcomplex(dc1_s4, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#5. completer set analysis -----
# NA

#6. Using ISI only -------
## Format trial data
data_s6<- data[(data$Severity_scale=="ISI"),] 

data_cnma_s6 <- data_s6 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_s6 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_s6,
  studlab = study
)

# CNMA
dc1_s6 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                  data = df_cnma_s6, sm = "OR", fixed = FALSE)
dc1_s6

forest(netcomplex(dc1_s6, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#7. Excluding delivery format components (ind, gp, ff, he, ae, tg) -------
## Format trial data
data_s7 <- data

data_cnma_s7 <- data_s7 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

data_cnma_s7$treatment_component=
  gsub("\\+dt","", data_cnma_s7$treatment_component)
data_cnma_s7$treatment_component=
  gsub("\\+ind","", data_cnma_s7$treatment_component)
data_cnma_s7$treatment_component=
  gsub("\\+gp","", data_cnma_s7$treatment_component)
data_cnma_s7$treatment_component=
  gsub("\\+ff","", data_cnma_s7$treatment_component)
data_cnma_s7$treatment_component=
  gsub("\\+tg","", data_cnma_s7$treatment_component)
data_cnma_s7$treatment_component=
  gsub("\\+ae","", data_cnma_s7$treatment_component)
data_cnma_s7$treatment_component=
  gsub("\\+he","", data_cnma_s7$treatment_component)

data_cnma_s7 <- data_cnma_s7 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_s7 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_s7,
  studlab = study
)

# CNMA
dc1_s7 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                  data = df_cnma_s7, sm = "OR", fixed = FALSE)
dc1_s7

forest(netcomplex(dc1_s7, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w")),
       xlim = c(0.05, 20))

# 8. Continuous outcome ----------------------------------------------------------------
## Format trial data
data$Severity_ep_n <- as.numeric(data$Severity_ep_n) 
data$Severity_ep_mean <- as.numeric(data$Severity_ep_mean) 
data$Severity_ep_sd <- as.numeric(data$Severity_ep_sd) 

data_cnma_s8 <- data %>% 
  dplyr::select(c(study, treatment_component,Severity_ep_n,Severity_ep_mean,Severity_ep_sd)) %>% # select rows
  dplyr::group_by(study,treatment_component)

df_cnma_s8 <- netmeta::pairwise(
  treat = treatment_component,
  n = Severity_ep_n,
  mean = Severity_ep_mean,
  sd = Severity_ep_sd,
  sm = "SMD",
  data = data_cnma_s8,
  studlab = study,
)

# CNMA
dc1_cont_s8 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                       data = df_cnma_s8, sm = "SMD", fixed = FALSE)
dc1_cont_s8

forest(netcomplex(dc1_cont_s8, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(-2, 2))

#9. Excluding small trials -------
## Format trial data
data_s9<- data[!(data$n<10),] 

data_cnma_s9 <- data_s9 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_s9 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_s9,
  studlab = study
)

# CNMA
dc1_s9 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                  data = df_cnma_s9, sm = "OR", fixed = FALSE)
dc1_s9

forest(netcomplex(dc1_s9, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))