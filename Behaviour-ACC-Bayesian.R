acc_allsessions_bayes_alltrials <- function(
  data_path,
  output_path,
  iter,
  cores,
  threads
) {
  # -------------------- Load libraries --------------------
  required_packages <- c("brms", "dplyr", "tidyr", "readxl", "posterior", "rstan")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed. Please install it before running on the cluster."))
    }
  }
  
  library(brms)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(posterior)
  
  # -------------------- Load data --------------------
  data_all <- readxl::read_xlsx(data_path)
  df_tus_eeg <- data_all%>%
    dplyr::filter(Stimulation < 5) 
  # -------------------- Preprocessing --------------------
  
  #Prepare TUS eeg data
  df_tus_eeg$Stimulation <- factor(df_tus_eeg$Stimulation)
  df_tus_eeg$Congruence <- ifelse(df_tus_eeg$Congruence == -1,-0.5, 0.5)
  df_tus_eeg$Feedback <- ifelse(df_tus_eeg$Feedback == 0,0.5, -0.5)
  df_tus_eeg$Subject <- as.factor(df_tus_eeg$Subject)
  #df_tus_eeg$Error <- ifelse(df_tus_eeg$Error == 0,-0.5,0.5)
  df_tus_eeg$BlockTrnr_z <- scale(df_tus_eeg$BlockTrnr, center = TRUE, scale = TRUE)
  df_tus_eeg$Trnr_z      <- scale(df_tus_eeg$Trnr, center = TRUE, scale = TRUE)
  
  
  df_tus_eeg <- df_tus_eeg %>%
    arrange(Subject, Trnr) %>%   # make sure trials are ordered
    group_by(Subject) %>%
    mutate(
      PrevError = lag(Error, default = -0.5),        # previous trial's Error
      PrevFeedback = lag(Feedback, default = -0.5),  # previous trial's Feedback
      PrevCongruence = lag(Congruence, default = -0.5) # previous trial's Congruence
    ) %>%
    ungroup()
  
      df_tus_eeg$Stimulation <- relevel(df_tus_eeg$Stimulation, "4")
      # -------------------- Priors --------------------
      priors_from_baseline <- c(
        prior(normal(3.14, 0.6), class = "Intercept"),
        prior(normal(1.55, 0.6), class = "b", coef = "Congruence"),
        prior(normal(0.0144, 0.6), class = "b", coef = "BlockTrnr_z"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation1"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation2"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation3"),
        prior(normal(0.037, 0.6), class = "b", coef = "PrevError"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation1:Congruence"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation2:Congruence"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation3:Congruence"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation1:PrevError"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation2:PrevError"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation3:PrevError"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation1:BlockTrnr_z"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation2:BlockTrnr_z"),
        prior(normal(0, 0.6), class = "b", coef = "Stimulation3:BlockTrnr_z")
      )
      # -------------------- Model --------------------
      model <- brm(
        Error ~ Stimulation*Congruence + Stimulation*PrevError + Stimulation*BlockTrnr_z +
          (1 + Stimulation*Congruence + Stimulation*PrevError + Stimulation*BlockTrnr_z || Subject),
        data   = df_tus_eeg,
        family = bernoulli(link = "logit"),
        prior  = priors_from_baseline,
        iter   = iter,
        cores  = cores,
        threads = threading(threads),
        control = list(adapt_delta = 0.99, max_treedepth = 12),
        save_pars = save_pars(all = TRUE)
      )
      # -------------------- Save --------------------
      save(model, file = output_path)
      return(model)
}