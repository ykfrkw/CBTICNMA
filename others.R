## R code for 
## Components and delivery formats of cognitive-behavioral therapy for chronic insomnia: 
## a systematic review and component network meta-analysis
## 
## The codes are for additional analyses. 
## Please load the dataset and library using the primary analysis code first.
##

# testing the imputation medhod -----

library(irr)
# remission + remission_long
data_remission   <- c(data$remission, data$remission_long)
data_remission <- as.numeric(data_remission)
# r_imputed + r_imputed_long
data_r_imputed   <- c(data$r_imputed, data$r_imputed_long)
data_r_imputed <- as.numeric(data_r_imputed)
# bind columns
correlation <- bind_cols(data_remission, data_r_imputed)

icc(correlation, model = "twoway",
    type = "consistency", unit = "single")
plot(correlation,
     xlab="raw",
     ylab="imputed",
     xlim=c(0,250),
     ylim=c(0,250))

# funnel plots  ------------------------------------------------------------------
# CBT vs PE -----
## Format trial data
df_binary <- df_nma %>%
  dplyr::filter(treat1=="CBT" & treat2=="PE")

head(df_binary)

# metabin
meta_binary <- metabin(event1,
                       n1,
                       event2,
                       n2,
                       data = df_binary,
                       studlab,
                       fixed = FALSE,
                       random = TRUE,
                       method.tau = "DL",
                       hakn = FALSE,
                       prediction = TRUE,
                       incr = 0.1,
                       sm = "OR",
                       digits.pval = 2) 

# Define fill colors for contour
col.contour = c("gray75", "gray85", "gray95")

# Generate funnel plot (we do not include study labels here)
funnel.meta(meta_binary, xlim = c(0.1, 1000),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour)

# Add a legend
legend("topright",#x = 1.6, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)

# Add a title
title("Contour-Enhanced Funnel Plot (CBT vs PE)")

# CBT vs WL -----
## Format trial data
df_binary <- df_nma %>%
  dplyr::filter(treat1=="CBT" & treat2=="WL")

head(df_binary)

# metabin
meta_binary <- metabin(event1,
                       n1,
                       event2,
                       n2,
                       data = df_binary,
                       studlab,
                       fixed = FALSE,
                       random = TRUE,
                       method.tau = "DL",
                       hakn = FALSE,
                       prediction = TRUE,
                       incr = 0.1,
                       sm = "OR",
                       digits.pval = 2) 

# Define fill colors for contour
col.contour = c("gray75", "gray85", "gray95")

# Generate funnel plot (we do not include study labels here)
funnel.meta(meta_binary, xlim = c(0.1, 1000),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour)

# Add a legend
legend("topright",#x = 1.6, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)

# Add a title
title("Contour-Enhanced Funnel Plot (CBT vs WL)")

# BT vs PE -----
## Format trial data
df_binary <- df_nma %>%
  dplyr::filter(treat1=="BT" & treat2=="PE")

head(df_binary)

# metabin
meta_binary <- metabin(event1,
                       n1,
                       event2,
                       n2,
                       data = df_binary,
                       studlab,
                       fixed = FALSE,
                       random = TRUE,
                       method.tau = "DL",
                       hakn = FALSE,
                       prediction = TRUE,
                       incr = 0.1,
                       sm = "OR",
                       digits.pval = 2) 

# Define fill colors for contour
col.contour = c("gray75", "gray85", "gray95")

# Generate funnel plot (we do not include study labels here)
funnel.meta(meta_binary, xlim = c(0.1, 1000),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour)

# Add a legend
legend("topright",#x = 1.6, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)

# Add a title
title("Contour-Enhanced Funnel Plot (BT vs PE)")

# BT vs WL -----
## Format trial data
df_binary <- df_nma %>%
  dplyr::filter(treat1=="BT" & treat2=="WL")

head(df_binary)

# metabin
meta_binary <- metabin(event1,
                       n1,
                       event2,
                       n2,
                       data = df_binary,
                       studlab,
                       fixed = FALSE,
                       random = TRUE,
                       method.tau = "DL",
                       hakn = FALSE,
                       prediction = TRUE,
                       incr = 0.1,
                       sm = "OR",
                       digits.pval = 2) 

# Define fill colors for contour
col.contour = c("gray75", "gray85", "gray95")

# Generate funnel plot (we do not include study labels here)
funnel.meta(meta_binary, xlim = c(0.1, 1000),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour)

# Add a legend
legend("topright",#x = 1.6, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)

# Add a title
title("Contour-Enhanced Funnel Plot (BT vs WL)")