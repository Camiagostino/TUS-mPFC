acc_tus_bayes_all_trials <- function(
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
  df_tus_eeg <- data_all %>%
    dplyr::filter(Stimulation < 4)
  # -------------------- Preprocessing --------------------
df_tus_eeg <- df_tus_eeg %>%
  mutate(
    Stimulation = factor(Stimulation),
    Congruence = ifelse(Congruence == -1, -0.5, 0.5),
    Feedback   = ifelse(Feedback == 0, 0.5, -0.5),
    Subject    = as.factor(Subject),
    BlockTrnr_z = scale(BlockTrnr, center = TRUE, scale = TRUE)[,1],
    Trnr_z      = scale(Trnr, center = TRUE, scale = TRUE)[,1],
    # ---- BlockStage: early (1–4), middle (5–8), late (9–12) ----
    # Participants missing late blocks will have NA for those stages
    BlockStage = case_when(
      BlockTrnr >= 1  & BlockTrnr <= 4  ~ "Early",
      BlockTrnr >= 5  & BlockTrnr <= 8  ~ "Middle",
      BlockTrnr >= 9  & BlockTrnr <= 12 ~ "Late",
      TRUE ~ NA_character_              # safety net for unexpected values
    ),
    BlockStage = factor(BlockStage, levels = c("Early", "Middle", "Late"))
  ) %>%
  arrange(Subject, Trnr) %>%
  group_by(Subject) %>%
  mutate(
    PrevError      = lag(Error, default = -0.5),
    PrevFeedback   = lag(Feedback, default = -0.5),
    PrevCongruence = lag(Congruence, default = -0.5)
  ) %>%
  ungroup()
df_tus_eeg$Stimulation <- relevel(df_tus_eeg$Stimulation, "2")
df_tus_eeg$BlockStage <- relevel(df_tus_eeg$BlockStage, "Late")
# -------------------- Priors --------------------
priors_from_baseline <- c(
  prior(normal(3.14, 0.6), class = "Intercept"),
  prior(normal(1.55, 0.6), class = "b", coef = "Congruence"),
  prior(normal(0.0144, 0.6), class = "b", coef = "BlockStageMiddle"),
  prior(normal(0.0144, 0.6), class = "b", coef = "BlockStageEarly"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation1"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation3"),
  prior(normal(0.037, 0.6), class = "b", coef = "PrevError"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation1:Congruence"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation3:Congruence"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation1:PrevError"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation3:PrevError"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation1:BlockStageMiddle"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation3:BlockStageMiddle"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation1:BlockStageEarly"),
  prior(normal(0, 0.6), class = "b", coef = "Stimulation3:BlockStageEarly")
)
# -------------------- Model --------------------
model <- brm(
  Error ~ Stimulation*Congruence + Stimulation*PrevError + Stimulation*BlockStage +
    (1 + Stimulation*Congruence + Stimulation*PrevError + Stimulation*BlockStage || Subject),
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