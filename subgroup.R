# ABOUT -------------------------------------------------------------------
## R code for 
## Components and delivery formats of cognitive-behavioral therapy for chronic insomnia: 
## a systematic review and component network meta-analysis
## 
## The codes are for subgroup analyses. 
## Please load the dataset and libarary using the primary analysis code first.
##

#r1. Age distribution --------
# age 18 to 35 -----
## Format trial data
data_r11<- data[(data$Age_mean<35),]

data_cnma_r11 <- data_r11 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r11 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r11,
  studlab = study
)

# CNMA
dc1_r11 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                  data = df_cnma_r11, sm = "OR", fixed = FALSE)
dc1_r11

forest(netcomplex(dc1_r11, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

# age 35 to 50 ----
data_r12<- data[(data$Age_mean>=35 & data$Age_mean<50),]

data_cnma_r12 <- data_r12 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r12 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r12,
  studlab = study
)

# CNMA
dc1_r12 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r12, sm = "OR", fixed = FALSE)
dc1_r12

forest(netcomplex(dc1_r12, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

# age 50 to 70 ------
data_r13<- data[(data$Age_mean>=50 & data$Age_mean<70),]

data_cnma_r13 <- data_r13 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r13 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r13,
  studlab = study
)

# CNMA
dc1_r13 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r13, sm = "OR", fixed = FALSE)
dc1_r13

forest(netcomplex(dc1_r13, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

# age 70+ ------
data_r14<- data[(data$Age_mean>=70),]

data_cnma_r14 <- data_r14 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r14 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r14,
  studlab = study
)

# CNMA
dc1_r14 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r14, sm = "OR", fixed = FALSE)
dc1_r14 # not estimable

forest(netcomplex(dc1_r14, 
                  c("se","sd","cr","th","cw","sr",
                    "sc",#"pi",
                    "re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#r2. Gender --------
data$F_n<-as.numeric(data$F_n)
data$Age_n<-as.numeric(data$Age_n)

#F>80% -----
data_r21<- data[(data$F_n/data$Age_n>0.8),]

data_cnma_r21 <- data_r21 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r21 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r21,
  studlab = study
)

# CNMA
dc1_r21 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r21, sm = "OR", fixed = FALSE)
dc1_r21

forest(netcomplex(dc1_r21, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

# F<20% -----
data_r22<- data[(data$F_n/data$Age_n<0.2),]

data_cnma_r22 <- data_r22 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r22 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r22,
  studlab = study
)

# CNMA
dc1_r22 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r22, sm = "OR", fixed = FALSE)
dc1_r22

forest(netcomplex(dc1_r22, 
                  c("se","sd","cr","th","sr", 
                    "sc","re","ns","w","ind","gp",
                    "ff","ae")),
       xlim = c(0.05, 20))

#r3. hypnotics  --------
# allowed -----
data_r31<- data[(data$dt==1),]

data_cnma_r31 <- data_r31 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r31 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r31,
  studlab = study
)

# CNMA

dc1_r31 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r31, sm = "OR", fixed = FALSE)
dc1_r31

forest(netcomplex(dc1_r31, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

# not allowed -----
data_r32<- data[(data$dt!=1),]

data_cnma_r32 <- data_r32 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r32 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r32,
  studlab = study
)

# CNMA
dc1_r32 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r32, sm = "OR", fixed = FALSE)
dc1_r32

forest(netcomplex(dc1_r32, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#r4. severity  --------
#ISI>=16 -----
data_r41<- data[(data$Severity_scale=="ISI"&data$Severity_bl_mean>=16),]

data_cnma_r41 <- data_r41 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r41 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r41,
  studlab = study
)

# CNMA
dc1_r41 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r41, sm = "OR", fixed = FALSE)
dc1_r41

forest(netcomplex(dc1_r41, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

#ISI<16 -----
data_r42<- data[(data$Severity_scale=="ISI"&data$Severity_bl_mean<16),]

data_cnma_r42 <- data_r42 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r42 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r42,
  studlab = study
)

# CNMA
dc1_r42 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r42, sm = "OR", fixed = FALSE)
dc1_r42

forest(netcomplex(dc1_r42, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))


#r5. comorbidities  --------
#Without -----
data_r51<- data[(substring(data$Comorbidities,first = 1,last = 1)==0),]

data_cnma_r51 <- data_r51 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r51 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r51,
  studlab = study
)

# CNMA
dc1_r51 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r51, sm = "OR", fixed = FALSE)
dc1_r51

forest(netcomplex(dc1_r51, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

# With -----
data_r52<- data[(substring(data$Comorbidities,first = 1,last = 1)>1),]

data_cnma_r52 <- data_r52 %>% 
  dplyr::select(c(study, treatment_component,r,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r = sum(r), n = sum(n)) 

df_cnma_r52 <- netmeta::pairwise(
  treat = treatment_component,
  event = r,
  n = n,
  sm = "OR",
  data = data_cnma_r52,
  studlab = study
)

# CNMA
dc1_r52 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                   data = df_cnma_r52, sm = "OR", fixed = FALSE)
dc1_r52

forest(netcomplex(dc1_r52, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

