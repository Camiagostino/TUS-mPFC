packages <- c("depmixS4", "dplyr", "ggplot2", "tidyr", "readxl", "tidyplots", "brms", "broom", 
              "lme4", "cmdstanr", "patchwork", "bayesplot", "reshape2", "truncnorm", "emmeans")

# Function to check if packages are installed, install them if not, and load them
load_packages <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      if (!require(package, character.only = TRUE)) {
        stop(paste("Package", package, "not found and installation failed. Please install it manually."))
      }
    }
    library(package, character.only = TRUE)
  }
}
load_packages(packages)


setwd("/Users/danieljanko/Desktop/Projects/TUS/data")

data <- read.csv("data_TUS.csv")

data$Error <- ifelse(data$Error == 0, -0.5, 0.5)
#data$Feedback <- ifelse(data$Feedback == 0, 0.5, -0.5)
data$Subject <- as.factor(data$Subject)
data$Stimulation <- as.factor(data$Stimulation)
#data$Trnr_z      <- scale(data$Trnr, center = TRUE, scale = TRUE)
data$Block <- scale(data$BlockTrnr, center = TRUE, scale = TRUE)

data <- data %>%
  arrange(Subject, Trnr) %>%   # make sure trials are ordered
  group_by(Subject) %>%
  mutate(
    PrevError = lag(ErrorMod, default = -0.5),        # previous trial's Error
    PrevFeedback = lag(Feedback, default = -0.5),   # previous trial's Feedback
    PrevCongruence = lag(Congruence, default = -0.5)
  ) %>%
  ungroup()
data$Stimulation <- relevel(data$Stimulation, "2")

priors_RT <- c(prior(normal(-0.145, 0.2), class = "b", coef = Congruence),
               prior(normal(0.2004,  0.2),  class = "b", coef = Error),
               prior(normal(-0.0015,  0.2),  class = "b", coef = PrevError),
               prior(normal(0.0069,  0.2),  class = "b", coef = Block),
               prior(normal(-1.119,  0.2),  class = Intercept),
               prior(normal(0, 0.2), class = b, coef = Stimulation1),
               prior(normal(0, 0.1), class = b, coef = Stimulation1:PrevFeedback),
               prior(normal(0, 0.1), class = b, coef = Stimulation1:PrevError),
               prior(normal(0, 0.1), class = b, coef = Stimulation1:Block),
               prior(normal(0, 0.1), class = b, coef = Stimulation1:Congruence),
               prior(normal(0, 0.1), class = b, coef = Stimulation1:Error),
               prior(normal(0, 0.2), class = b, coef = Stimulation3),
               prior(normal(0, 0.1), class = b, coef = Stimulation3:Block),
               prior(normal(0, 0.1), class = b, coef = Stimulation3:PrevError),
               prior(normal(0, 0.1), class = b, coef = Stimulation3:Congruence),
               prior(normal(0, 0.1), class = b, coef = Stimulation3:Error),
               prior(student_t(3, 0.4, 0.1), class = sd, coef = Intercept, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Block, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Congruence, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Error, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = PrevError, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation1, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation1:Block, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation1:PrevError, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation1:Congruence, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation1:Error, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation1:PrevFeedback, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation3, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation3:Block, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation3:PrevError, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation3:Congruence, group = Subject),
               prior(student_t(3, 0, 0.08), class = sd, coef = Stimulation3:Error, group = Subject)
               ) 

modelRT <- brm(RT ~ Stimulation*Congruence + Stimulation*Error + Stimulation*PrevError +Stimulation*Block + 
                 (1 + Stimulation*Congruence + Stimulation*Error +  Stimulation*PrevError + Stimulation*Block || Subject), 
               data = data, 
               prior = priors_RT,
               family = lognormal(),
               iter = 15000,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               threads = threading(4),
               cores = 4,
               save_pars = save_pars(all = TRUE)
)




