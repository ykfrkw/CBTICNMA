# ABOUT -------------------------------------------------------------------
## R code for 
## Components and delivery formats of cognitive-behavioral therapy for chronic insomnia: 
## a systematic review and component network meta-analysis
## 
## The codes are for secondary outcomes. 
## Please load the dataset and libarary using the primary analysis code first.
##

# 2.Acceptability: dropouts for any reason (dichotomous) ---------------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
## Format trial data

data$dropout <- as.numeric(data$dropout) 
data$n <- as.numeric(data$n)

data_cnma_dropout <- data
data_cnma_dropout<- data_cnma_dropout[!(data_cnma_dropout$treatment_group=="WL")|
                                        (data_cnma_dropout$treatment_group=="NT")|
                                        (data_cnma_dropout$treatment_group=="TAU"),]
data_cnma_dropout<- data_cnma_dropout[!is.na(data_cnma_dropout$dropout),]

data_cnma_dropout <- data_cnma_dropout %>% 
  dplyr::select(c(study, treatment_component,dropout,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(dropout = sum(dropout), n = sum(n)) 

df_cnma_dropout <- netmeta::pairwise(
  treat = treatment_component,
  event = dropout,
  n = n,
  sm = "OR",
  data = data_cnma_dropout,
  studlab = study
)
#head(df_cnma_dropout)

# +ANALYZE ----------------------------------------------------------------
# CNMA
dc1_dropout <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
               data = df_cnma_dropout, sm = "OR", fixed = FALSE)
dc1_dropout

forest(netcomplex(dc1_dropout, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))

# 3. Sleep diary measures (continuous) ----------------------------------------------------------------
# 3.1. Sleep efficiency (SE, %) ----------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
## Format trial data

data$SE_n <- as.numeric(data$SE_n) 
data$SE_mean <- as.numeric(data$SE_mean) 
data$SE_sd <- as.numeric(data$SE_sd) 

data_cnma_31 <- data %>% 
  dplyr::select(c(study, treatment_component,SE_n,SE_mean,SE_sd)) %>% # select rows
  dplyr::group_by(study,treatment_component)

df_cnma_31 <- netmeta::pairwise(
  treat = treatment_component,
  n = SE_n,
  mean = SE_mean,
  sd = SE_sd,
  sm = "MD",
  data = data_cnma_31,
  studlab = study,
)

# +ANALYZE ----------------------------------------------------------------
# CNMA

dc1_cont_31 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                    data = df_cnma_31, sm = "MD", fixed = FALSE)
dc1_cont_31

forest(netcomplex(dc1_cont_31, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(-10, 10))

# 3.2. Total sleep time (TST, min) ----------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
## Format trial data

data$TST_n <- as.numeric(data$TST_n) 
data$TST_mean <- as.numeric(data$TST_mean) 
data$TST_sd <- as.numeric(data$TST_sd) 

data_cnma_32 <- data %>% 
  dplyr::select(c(study, treatment_component,TST_n,TST_mean,TST_sd)) %>% # select rows
  dplyr::group_by(study,treatment_component)

df_cnma_32 <- netmeta::pairwise(
  treat = treatment_component,
  n = TST_n,
  mean = TST_mean,
  sd = TST_sd,
  sm = "MD",
  data = data_cnma_32,
  studlab = study,
)

# +ANALYZE ----------------------------------------------------------------
# CNMA

dc1_cont_32 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                       data = df_cnma_32, sm = "MD", fixed = FALSE)
dc1_cont_32

forest(netcomplex(dc1_cont_32, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(-60, 60))

# 3.3. Sleep onset latency (SOL, min) ----------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
## Format trial data

data$SOL_n <- as.numeric(data$SOL_n) 
data$SOL_mean <- as.numeric(data$SOL_mean) 
data$SOL_sd <- as.numeric(data$SOL_sd) 

data_cnma_33 <- data %>% 
  dplyr::select(c(study, treatment_component,SOL_n,SOL_mean,SOL_sd)) %>% # select rows
  dplyr::group_by(study,treatment_component)

df_cnma_33 <- netmeta::pairwise(
  treat = treatment_component,
  n = SOL_n,
  mean = SOL_mean,
  sd = SOL_sd,
  sm = "MD",
  data = data_cnma_33,
  studlab = study,
)

# +ANALYZE ----------------------------------------------------------------
# CNMA

dc1_cont_33 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                       data = df_cnma_33, sm = "MD", fixed = FALSE)
dc1_cont_33

forest(netcomplex(dc1_cont_33, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(-30, 30))

#3.4. Wake after sleep onset (WASO, min) ----------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
## Format trial data

data$WASO_n <- as.numeric(data$WASO_n) 
data$WASO_mean <- as.numeric(data$WASO_mean) 
data$WASO_sd <- as.numeric(data$WASO_sd) 

data_cnma_34 <- data %>% 
  dplyr::select(c(study, treatment_component,WASO_n,WASO_mean,WASO_sd)) %>% # select rows
  dplyr::group_by(study,treatment_component)

df_cnma_34 <- netmeta::pairwise(
  treat = treatment_component,
  n = WASO_n,
  mean = WASO_mean,
  sd = WASO_sd,
  sm = "MD",
  data = data_cnma_34,
  studlab = study,
)

# +ANALYZE ----------------------------------------------------------------
# CNMA

dc1_cont_34 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
                       data = df_cnma_34, sm = "", fixed = FALSE)
dc1_cont_34

forest(netcomplex(dc1_cont_34, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(-60, 60))

# 4. Efficacy at long-term follow-up. (dichotomous, longest follow-up between 3 to 12 months)----------------------------------------------------------------
# +BUILD -------------------------------------------------------------------
## Format trial data
data$r_long <- as.numeric(data$r_long) 
data$n <- as.numeric(data$n)

data_cnma_4 <- data %>% 
  dplyr::select(c(study, treatment_component,r_long,n)) %>% # select rows
  dplyr::group_by(study,treatment_component) %>% # order
  dplyr::summarise(r_long = sum(r_long), n = sum(n)) 

df_cnma_4 <- netmeta::pairwise(
  treat = treatment_component,
  event = r_long,
  n = n,
  sm = "OR",
  data = data_cnma_4,
  studlab = study
)

# +ANALYZE ----------------------------------------------------------------
# CNMA
dc1_4 <- discomb(TE, seTE, treat1, treat2, studlab, inactive="dt",
               data = df_cnma_4, sm = "OR", fixed = FALSE)
dc1_4

forest(netcomplex(dc1_4, 
                  c("se","sd","cr","th","cw","sr",
                    "sc","pi","re","ns","w","ind","gp",
                    "ff","tg","he","ae")),
       xlim = c(0.05, 20))
