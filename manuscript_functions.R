# Packages required in the functions of this document
library(kableExtra)
library(tidyverse)
library(rms)
library(glmnet)
library(cowplot)
library(hrbrthemes)
library(gmodels)
library(pROC)


#----- Table functions-----

# Make show NA default in table and allow it to work with "tab" abbreviration
tab = function (..., useNA = 'ifany') base::table(..., useNA = useNA)
table = function (..., useNA = 'ifany') base::table(..., useNA = useNA)


# 2x2 table with row proportions 
row_tab = function(..., prop.c=FALSE, prop.t = FALSE, expected=FALSE, prop.chisq = FALSE) gmodels::CrossTable(..., prop.c=prop.c, prop.t=prop.t, expected=expected, prop.chisq = prop.chisq)


# 2x2 table with column proportions
col_tab = function(..., prop.r=FALSE, prop.t = FALSE, expected=FALSE, prop.chisq = FALSE) gmodels::CrossTable(..., prop.r=prop.r, prop.t=prop.t, expected=expected, prop.chisq = prop.chisq)


# ---- Formatting functions ----

# Re-code as factor using case_when
# PACAKGES REQUIRED: tidyverse
fct_case_when <- function(...) {
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

# Not in 
'%!in%' <- function(x,y)!('%in%'(x,y))

# Mode 
mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


#----- Validation functions -----

# --- Output brier score, calibration slope, and calibration intercept --

# INPUTS
#' validation_results: results from "lrm" validate function

save_val_results <- function(validation_results){
  brier <- validation_results["B", "index.orig"]
  slope <- validation_results["Slope", "index.orig"]
  intercept <- validation_results["Intercept", "index.orig"]
  
  save_results <- rbind("Brier score" = brier, "Calibration slope" = slope, "Calibration intercept" = intercept)
  
  return(save_results)
  
}

# --- Bootstrap CI for apparent c-statistic based only on predicted values ---
# PACAKGES REQUIRED: gmodels, pROC
# Automatically does 2000 replications

# INPUTS
# pred_vals_vector: vector of predicted values 
# response_vector: vector of response variables (outcome)

c_ci <- function(pred_vals_vector, response_vector){
  
  save_roc <- suppressMessages(roc(response_vector, pred_vals_vector))
  # Bootstrap CI with 2000 reps
  save_roc_ci <- suppressMessages(ci(save_roc, method="bootstrap"))
  c <- save_roc$auc
  roc_ci_lb <- save_roc_ci[1]
  roc_ci_ub <- save_roc_ci[3]
  
  c_stat_app <- cbind("c_stat" = c, "LB" = roc_ci_lb, "UB" = roc_ci_ub)
  rownames(c_stat_app) <- NULL
  
  return(c_stat_app)
}

# -- Boostrap CI for the c-statistic (apparent and optimism corrected) --
# PACAKGES REQUIRED: rms, proc
# The CI for the original c-statisic is estimated from the pROC package, which uses 2000 stratified bootstrap replicates
# The CI for the internal validation index corrected c-statistic is estimated with N samples (default is 2000) with replacement from the data and using Frank Harrell's validate function

# INPUTS
#' model: stored model fit from "lrm"
#' data: dataset
#' samples: number of re-samples 

c_stat_ci <- function(model, data, samples=2000){
  
  predvals <- predict(model)
  save_roc <- suppressMessages(roc(model$y, predvals))
  save_roc_ci <- suppressMessages(ci(save_roc, method="bootstrap"))
  roc_ci_lb <- save_roc_ci[1]
  roc_ci_ub <- save_roc_ci[3]
  
  c_stat_app <- cbind(roc_ci_lb, roc_ci_ub)
  
  # Bootstrap re-sampling for 95% CI of optimism corected c-statistic
  c_stat <- numeric(samples)  
  dxy <- numeric(samples)  
  c_stat_orig <- numeric(samples)  
  dxy_orig <- numeric(samples)  
  n <- nrow(data)
  
  for(i in 1 : samples) {
    v <- validate(model, B=1)
    dxy_orig[i] <- v['Dxy', 'index.orig']
    c_stat_orig[i] <- 0.5*(dxy_orig[i] + 1)
    dxy[i] <- v['Dxy', 'index.corrected']
    c_stat[i] <- 0.5*(dxy[i] + 1)
  }
  c_stat_mean <- mean(c_stat)
  c_stat_int <- quantile(c_stat, c(.025, .975))
  c_stat_correct <- c(c_stat_mean, c_stat_int)
  
  c_stat_mean_orig <- mean(c_stat_orig)
  c_stat_orig <- c(c_stat_mean_orig, c_stat_app)
  
  return(rbind(orig = c_stat_orig, corrected=c_stat_correct))
}




# -- External validation with LRM --
# PACKAGES REQUIRED: rms

# INPUTS
#' method: method for model updating, includes either "model recalibration", or "model revision" (unlike GLM below, this doesn't allow offset term, so we are limited to recalibration)
#' lp: linear predictor that is being validated
#' outcome: outcome of the model
#' formula: model formula when the model is being revised 
#' data: dataset 
#' samples: number of resamples to do for validation
#' boots: number of boostraps to do for estimating optimism
#' return: what to return, includes "model", "performance", "c_stat" or "cal_plot"


ev_lrm <- function(method, lp, outcome, formula=NULL, data, samples, boots, return){
  
  # Create forumula
  if (method=="model recalibration") formula <- as.formula(paste0(outcome, "~", lp))
  if (method=="model revision") formula <- formula
  
  # Run the logistic regression model
  model <- lrm(formula, data=data, x=TRUE, y=TRUE)
  
  # Validate the model
  validate <- validate(model, B=samples)
  
  # Output the c-statistic
  c_statistic <- 0.5 * (validate['Dxy', 'index.corrected'] + 1)
  
  # Bootstrap re-sampling for 95% CI of c-statistic
  c_stat <- numeric(boots)  
  dxy <- numeric(boots)  
  n <- nrow(data)
  
  for(i in 1 : boots) {
    g <- update(model, subset=sample(1:n, n, replace=TRUE))
    v <- validate(g, B=samples)
    dxy[i] <- v['Dxy', 'index.corrected']
    c_stat[i] <- 0.5*(dxy[i] + 1)
  }
  c_stat_median <- median(c_stat)
  dxy_ci <- quantile(dxy, c(.025, .975))
  c_stat_ci <- quantile(c_stat, c(.025, .975))
  
  # Calibration plot
  calibration <- calibrate(model, B=samples)
  
  # Setting which to return
  if (return=="model") return(print(model))
  if (return=="performance") return(validate)
  if (return=="c_stat") return(c(c_stat_median, c_stat_ci))
  if (return=="cal_plot") return(plot(calibration))
}


# -- External validation with GLM --
# PACKAGES REQUIRED: rms, tidyverse, ggplot2, hrbrthemes, cowplot, pROC

# INPUTS
#' method: method for model updating, includes either "original", "update intercept", "model recalibration", or "model revision"
#' lp: linear predictor that is being validated (relevant for all methods except model revision. should be predicted values not predicted risk)
#' outcome: outcome of the model (in quotes)
#' formula: model formula when the model is being revised 
#' data: dataset 
#' samples: number of bootstrap in bootstrap resampling (200 recommended)
#' return: what to return, includes "model", "performance", or "cal_plot"

ev_glm <- function(method, lp=NULL, outcome, formula=NULL, data, samples, return){
  
  # Create formula
  if (method=="original") formula <- as.formula(paste0(outcome, "~ -1 + offset(", lp, ")"))
  if (method=="update intercept") formula <- as.formula(paste0(outcome, "~ offset(", lp, ")"))
  if (method=="model recalibration") formula <- as.formula(paste0(outcome, "~", lp))
  if (method=="model revision") formula <- formula
  
  # Model 
  model <- glm(formula, data=data, family=binomial)
  model_outcome <- model$y
  data$out_var <- model$y
  
  # Predicted values
  data$pred_risk <- predict(model, type = "response")
  
  
  roc <- suppressMessages(roc(data, out_var, pred_risk, ci=TRUE))
  
  # Validate the model (apparent validation)
  performance <- bootstrap_validation(form=formula, df=data, boots=samples, predicted_risk=data$pred_risk, outcome=model_outcome) 
  
  calibration_plot <- suppressMessages(suppressWarnings(
    cal_plot(data=data, pred_risk=data$pred_risk, out_var=data$out_var)
  ))
  
  if (return=="model") return(print(model))
  if (return=="performance") return(performance)
  if (return=="cal_plot") return(calibration_plot)
  if (return=="roc") return(roc)
}

# -- Bootstrap internal validation with correction for optimism --
# Have to have run a model with lrm before running this code
# Runs inside ev_glm

# Inputs
#' form = formula for external validation
#' df = data
#' boots = number of bootstraps for internal validation
#' predicted_risks = stored vector of predicted risks 
#' outcome = model outcome as stored vector of outcome values

bootstrap_validation <- function(form=formula, df=data, boots=samples, predicted_risk=data$pred_risk, outcome=model_outcome) {
  
  # Set up matrix to store results
  row_names <- c("Dxy", "C-statistic", "R2", "Brier", "Intercept", "Slope", "E-max", "E-average")
  col_names <- c("apparent", "optimism", "corrected", "boots")
  performance <- matrix(nrow=8, ncol=4, dimnames = list(row_names, col_names))
  performance[,4] <- boots 
  
  validation <- val.prob(predicted_risk, outcome, pl=FALSE) %>% round(3)
  
  dxy <- validation["Dxy"]
  c_stat <- validation["C (ROC)"]
  r2 <- validation["R2"]
  brier <- validation["Brier"]
  intercept <- validation["Intercept"]
  slope <- validation["Slope"]
  e_max <- validation["Emax"]
  e_avg <- validation["Eavg"]
  
  apparent <- rbind(dxy, c_stat, r2,  brier, intercept, slope,  e_max, e_avg)
  
  # Add to first column of matrix
  performance[,1] <- apparent
  
  # Bootstrap validation and get confidence intervals
  
  # Place holders for bootstrap optimism 
  boot_sample_perf <- list()
  boot_orig_perf <- list()
  opp <- list()
  
  for(i in 1:boots){
    boot_sample <- sample_n(df, nrow(df), replace = TRUE)
    boot_model <- glm(form, data = boot_sample, family = binomial)
    boot_outcome <- boot_model$y 
    
    # Estimate performance measures using fitted model from bootstrap sample and bootstrap dataset
    # val.prob(predicted probability, vector of binary outcomes)
    boot_sample_perf[[i]] <- val.prob(predict(boot_model, data=boot_sample, type = "response"),
                                      boot_outcome, 
                                      pl = FALSE) %>% round(3)
    
    # Estimate performance measured with fitted model from bootstrap sample and original dataset
    boot_orig_perf[[i]] <- val.prob(predict(boot_model, newdata=df, type = "response"),
                                    outcome, 
                                    pl = FALSE) %>% round(3)
    
    # Estimate optimism
    opp[[i]] <- boot_sample_perf[[i]] - boot_orig_perf[[i]] 
  }
  
  dxy_opp <- mean(unlist(map(opp, "Dxy")))
  c_opp <- mean(unlist(map(opp, "C (ROC)")))
  r2_opp <- mean(unlist(map(opp, "R2")))
  brier_opp <- mean(unlist(map(opp, "Brier")))
  intercept_opp <- mean(unlist(map(opp, "Intercept")))
  slope_opp <- mean(unlist(map(opp, "Slope")))
  e_max_opp <- mean(unlist(map(opp, "Emax")))
  e_avg_opp <- mean(unlist(map(opp, "Eavg")))
  
  optimism <- rbind(dxy_opp, c_opp, r2_opp, brier_opp, intercept_opp, slope_opp, e_max_opp, e_avg_opp)
  
  # 2nd column of matrix
  performance[,2] <- optimism
  
  # 3rd column of matrix    
  performance[,3] <- performance[,1] - performance[,2]
  
  return(performance)
}


# -- Calibration plot --
# Runs inside ev_glm

# Inputs
#' data = data
#' pred_risk = vector of predicted risks
#' out_var = vector of outcome variable

cal_plot <- function(data=data, pred_risk=pred_risk, out_var=out_var){
  
  # Bin prediction into 10 groups
  g1 <- data %>% 
    mutate(bin = ntile(pred_risk, 10)) %>% 
    group_by(bin) %>%
    summarize(n = n(), # Get estimates and confidence intervals
           bin_pred = mean(pred_risk), 
           bin_prob = mean(out_var), 
           se = sqrt((bin_prob * (1 - bin_prob)) / n), 
           ul = bin_prob + 1.96 * se, 
           ll = bin_prob - 1.96 * se) %>%
    ggplot(aes(x=bin_pred, y=bin_prob, ymin = ll, ymax = ul)) +
    geom_pointrange(size =.75, fill="black", color = "grey", shape=23, fatten=5) +
    xlab("") +
    ylab("Observed probability") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    geom_abline(intercept = 0, slope = 1, size = 0.25, lty = 2) +
    geom_smooth(method = "loess", color="tomato1", se = FALSE, size=0.5, lty=1) +
    theme_ipsum(
      axis_title_size = 11,
      plot_margin = margin(10, 10, 0, 10)
    ) +
    ggtitle("Calibration plot")
  
  # The distribution plot        
  g2 <- ggplot(data, aes(x = pred_risk)) +
    geom_histogram(fill = "black", bins = 200) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    xlab("Predicted Probability") +
    ylab("") +
    theme_ipsum(
      axis_title_size = 11,
      plot_margin = margin(0, 10, 10, 10)
    ) +
    theme(panel.grid.minor = element_blank())
  
  # Combine them    
  graph <- plot_grid(g1, g2, nrow=2, rel_heights = c(5,2))
  
  return(graph)
}

# -- Calibration plot SOLO --
# Runs indepdently

# Inputs
#' data = data
#' pred_val_vec = vector of predicted values
#' out_var_vec = vector of outcome values
#' 
cal_plot_solo <- function(data=data, pred_val_vec=NULL, out_var_vec=NULL){

  graph_data <- data %>%  
    mutate(predval = pred_val_vec,
          pred_risk = 1/(1+exp(-predval)), 
           out_var = out_var_vec, 
           bin = ntile(pred_risk, 10)) %>% 
    select(pred_risk, out_var, bin)
  
  # Bin prediction into 10 groups
  g1 <- graph_data %>% 
    group_by(bin) %>%
    summarize(n = n(), # Get estimates and confidence intervals
              bin_pred = mean(pred_risk), 
              bin_prob = mean(out_var), 
              se = sqrt((bin_prob * (1 - bin_prob)) / n), 
              ul = bin_prob + 1.96 * se, 
              ll = bin_prob - 1.96 * se) %>%
    ggplot(aes(x=bin_pred, y=bin_prob, ymin = ll, ymax = ul)) +
    geom_pointrange(size =.75, fill="black", color = "grey", shape=23, fatten=5) +
    xlab("") +
    ylab("Observed probability") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    geom_abline(intercept = 0, slope = 1, size = 0.25, lty = 2) +
    geom_smooth(method = "loess", color="tomato1", se = FALSE, size=0.5, lty=1) +
    theme_ipsum(
      axis_title_size = 11,
      plot_margin = margin(10, 10, 0, 10)
    ) +
    ggtitle("Calibration plot")
  
  # The distribution plot        
  g2 <- ggplot(graph_data, aes(x = pred_risk)) +
    geom_histogram(fill = "black", bins = 200) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    xlab("Predicted Probability") +
    ylab("") +
    theme_ipsum(
      axis_title_size = 11,
      plot_margin = margin(0, 10, 10, 10)
    ) +
    theme(panel.grid.minor = element_blank())
  
  # Combine them    
  graph <- plot_grid(g1, g2, nrow=2, rel_heights = c(5,2))
  
  return(graph)
}

#----- Added value functions-----
# PACKAGES REQUIRED: tidyverse, kableExtra

# INPUTS
#' base_model: stored LRM results from the base model
#' add_model: stored LRM results from the added value model
#' Requires the two models have the same sample size

# Based on Frank's article 
added_value <- function(base_model, add_model, outcome_vector){
  
  # Outcome variable
  y = outcome_vector
  
  # Likelihood ratio test
  lrtest_results <- lrtest(base_model, add_model)
  
  # Save predicted risks (0-1)
  predrisk_base = predict(base_model, type="fitted.ind")
  predrisk_add = predict(add_model, type="fitted.ind")
  
  # Save predicted values
  predvals_base <- predict(base_model, type="fitted")
  predvals_add <- predict(add_model, type="fitted")
  
  # IDI and NRI 
  save <- improveProb(predrisk_base, predrisk_add, y)
  
  nri <- c(est = save$nri, lb = save$nri-1.96*save$se.nri, ub = save$nri+1.96*save$se.nri)
  idi <- c(est = save$idi, lb = save$idi-1.96*save$se.idi, ub = save$idi+1.96*save$se.idi)
  
  nri_idi <- rbind(nri=nri, idi=idi)
  
  # --- Table ---
  lra <- base_model$stats['Model L.R.']
  lrb <- add_model$stats['Model L.R.']
  base_adequacy <- lra/lrb
  frac_new_info_1 <- 1-base_adequacy
  
  ra  <- base_model$stats['R2']
  rb  <- add_model$stats['R2']
  
  var_base <- var(predvals_base)
  var_add <- var(predvals_add)
  rel_explain_var <- var_base/var_add
  frac_new_info_2 <- 1-rel_explain_var
  
  br2 <- function(p) var(p) / (var(p) + sum(p * (1 - p)) / length(p))
  
  frac_var_base <- br2(predvals_base)
  frac_var_add <- br2(predvals_add)
  rel_explain_var_frac <- frac_var_base/frac_var_add
  frac_new_info_3 = 1-rel_explain_var_frac 
  
  add_value_table <- as.data.frame(rbind("Base model LR X2"=lra, 
                                         "New model LR X2"=lrb, 
                                         "Adequacy of base model"=base_adequacy, 
                                         "Fraction of new information" = frac_new_info_1,
                                         "Base model Nagerlkerke pseudo R2"=ra, 
                                         "New model Nagerlkerke pseudo R2"=rb, 
                                         "Variance of base model risk"=var_base, 
                                         "Variance of new model risk"=var_add, 
                                         "Relative explained variation"=rel_explain_var, 
                                         "Fraction of new information"=frac_new_info_2, 
                                         "Base model fraction explained risk"=frac_var_base, 
                                         "New model fraction explained risk"=frac_var_add, 
                                         "Relative explained variation"=rel_explain_var_frac, 
                                         "Fraction of new information"=frac_new_info_3))  %>% 
    rownames_to_column("Index") %>% 
    rename(Value = "Model L.R.")
  
  formatted_table <- kable(add_value_table) %>% 
    kableExtra::kable_styling(font_size = 10) 
  
  return(list("table" = add_value_table,
                 "formatted_table" = formatted_table,
                 "lrtest" = lrtest_results,
                 "nri_idi" = nri_idi)
  )
}

#------ Analysis functions -----


# --- Model selection --
# PACKAGES REQUIRED: tidyverse, kableExtra


# INPUTS
#' model_formula: formula of full model
#' data: data frame
#' boots: number of bootstraps to repeat backward selection
#' return: includes "overview" or "frequency table" or "pairwise table" or "all"
#' pval_pair: pvalue for the pairwise comparison (recommend 0.01)

model_selection <- function(model_formula, data, boots=200, return, pval_pair=0.01, seed=5321) {

  # Estimate full model 
  formula <- model_formula
  full_mod <- glm(formula, data = data, family="binomial", x = T, y = T)
  full_est <- coef(full_mod)
  full_se <- coef(summary(full_mod))[, "Std. Error"]
  
  var_names <- names(full_est)
  outcome <- full_mod$y
  
  # Selected model 
  # From 1 round of backwards selection
  sel_mod <- stats::step(glm(formula, data = data, family="binomial", x=T, y=T), 
                  direction = "backward",
                  trace = 0, seed=seed)
  
  summary(sel_mod)
  
  sel_est <- stats::coef(sel_mod)[var_names]
  sel_est[is.na(sel_est)] <- 0
  names(sel_est) <- var_names
  sel_se <- coef(summary(sel_mod))[, "Std. Error"][var_names]
  sel_se[is.na(sel_se)] <- 0
  
  # Bootstrap selection
  bootnum <- boots
  boot_est <-  boot_se <- matrix(0, ncol = length(var_names), nrow = bootnum,
                                 dimnames = list(NULL, c(var_names)))
  
  set.seed(seed)
  for (i in 1:bootnum) {
    data_id <- sample(1:dim(data)[1], replace = T)
    boot_mod <- stats::step(glm(formula, data = data[data_id, ], family="binomial", x=T, y=T), 
                     direction = "backward", trace = 0)
   boot_est[i, names(stats::coef(boot_mod))] <- stats::coef(boot_mod)
   boot_se[i, names(stats::coef(boot_mod))] <- stats::coef(summary(boot_mod))[, "Std. Error"]
    
   boot_outcome <- boot_mod$y
   
  }

  
    
  # Bootstrap inclusion
  
  boot_01 <- (boot_est != 0) * 1
  boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)
  
  # TABLE
  sqe <- (t(boot_est) - full_est) ^ 2
  rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
  rmsdratio <- rmsd / full_se
  boot_mean <- apply(boot_est, 2, mean)
  boot_meanratio <- boot_mean / full_est
  boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
  boot_median <- apply(boot_est, 2, median)
  boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
  boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))
  
  overview <- round(cbind(full_est, full_se, sel_est, sel_se, 
                          rmsdratio, boot_inclusion, boot_relbias, boot_median, boot_025per, 
                          boot_975per), 4)
  overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
  
  # --- Overview table ---
  formatted_overview <- kable(overview, 
        col.names = c("Estimate","Standard error", "Estimate", "Standard error", "RMSD ratio", "Bootstrap inclusion frequency (%)", 
                      "Relative conditional bias (%)", "Median", "2.5th percentile", "97.5% percentile")) %>% 
    kableExtra::kable_styling(font_size = 10) %>% 
    add_header_above(c(" ", "Global Model"=2, "Selected Model"=2, " "=3, "Bootstrap"=3))
  
  
  # --- Model frequency table ---
  length_vars <- length(var_names)
  
  pred_ord <- names(boot_inclusion)[order(boot_inclusion, decreasing = T)]
  boot_01 <- boot_01[, pred_ord]
  boot_01 <- cbind(boot_01, count = rep(1, times = bootnum))
  boot_modfreq <- stats::aggregate(count ~ ., data = boot_01, sum)
  boot_modfreq[, "percent"] <- boot_modfreq$count / bootnum * 100
  boot_modfreq <- boot_modfreq[order(boot_modfreq[, "percent"], decreasing = T), ]
  boot_modfreq[, "cum_percent"] <- cumsum(boot_modfreq$percent)
  boot_modfreq <- boot_modfreq[boot_modfreq[, "cum_percent"] <= 80, ]
  
  if (dim(boot_modfreq)[1] > 20) boot_modfreq <- boot_modfreq[1:20, ]
  
  freq_table <- cbind("Predictors"= apply(boot_modfreq[,c(2:length_vars)], 1, 
                                          function(x) paste(names(x[x==1]), collapse=" ")),
                      boot_modfreq[,c("count", "percent", "cum_percent")])
  
  formatted_freq_table <- kable(freq_table) %>% 
    kableExtra::kable_styling(font_size = 10) 
  
  # --- Pairwise inclusion frequency in % ---
  # Table shows the variables that are often included together ("rope teams")or variables that may be substituted for one another ("competitors").
  # A > sign indicates that variables are selected together more often than would be expected by chance
  # A < sign indicates variables selected together less often than expected by chance
  
  pval <- pval_pair
  boot_pairfreq <- matrix(100, ncol = length(var_names)-1, nrow = length(var_names)-1,
                          dimnames = list(pred_ord[-1], 
                                          pred_ord[-1]))
  
  expect_pairfreq <- NULL
  combis <- combn(pred_ord[-1], 2)
  
  for (i in 1:dim(combis)[2]) {
    boot_pairfreq[combis[1, i], combis[2, i]] <-
      sum(apply(boot_01[, combis[, i]], 1, sum) == 2) / bootnum * 100
    expect_pairfreq[i] <-
      boot_inclusion[grepl(combis[1, i], pred_ord)][1] *
      boot_inclusion[grepl(combis[2, i], pred_ord)][1] / 100
    boot_pairfreq[combis[2, i], combis[1, i]] <-
      ifelse(is(suppressWarnings(try(chisq.test(boot_01[, combis[1, i]], 
                                                boot_01[, combis[2, i]]), 
                                     silent = T)), "try-error"),
             NA, ifelse(suppressWarnings(
               chisq.test(boot_01[, combis[1, i]],
                          boot_01[, combis[2, i]])$p.value) > pval,
               "", ifelse(as.numeric(boot_pairfreq[combis[1, i], combis[2, i]]) < 
                            expect_pairfreq[i], "<", ">")))
  }
  
  diag(boot_pairfreq) <- boot_inclusion[pred_ord][-1]
  
  pair_table <- kable(boot_pairfreq) %>% 
    kableExtra::kable_styling(font_size = 10) %>% 
    scroll_box(width="100%", height="600px")
  
  
  
  if (return=="overview") return(formatted_overview)
  if (return=="frequency table") return(formatted_freq_table)
  if (return=="pairwise table") return(pair_table)
  if (return=="all") 
    return(list("overview" = formatted_overview,
                "freq_table" = formatted_freq_table,
                "pair_table" = pair_table,
                "selected_model" = sel_mod,
                "raw_overview" = overview
           ))
  

}


# -- MODEL DEVELOPMENT --
# OTHER FUNCTIONS REQUIRED: save_val_results, c_ci, model_selection, cal_plot_solo

# INPUTS
#' data: data frame
#' outcome:  written as data$outcome
#' model_formula: formula of full model
#' approx_full_formula = formula with "predvals_full_model" as dependent variable
#' boostraps: number of boostraps for validation and bootstrap resampling
#' threshold: % of bootstraps that we require a variable to be selected in (default is 50)


model_development <- function(data, outcome, model_formula, approx_full_formula, bootstraps=200, threshold=50, seed=53421){
  
  # Outcome vector 
  outcome_vector <- outcome
  
  #--- Full model 
  full_model <- lrm(model_formula, data = data, x = T, y = T)
  
  # Save coefficients
  coef_full_model <- full_model$coefficients
  
  # Predicted values from full model
  predvals_full_model <- predict(full_model)
  
  # Performance
  full_val <- validate(full_model, B=bootstraps, seed=seed)
  perf_full_model <- save_val_results(full_val)
  c_full_model <- c_ci(predvals_full_model, outcome_vector)
  
  # C statistic - FULL MODEL (corrected)
  c_stat_full <- c_stat_ci(model=full_model, data=data)
  
  
  #---- Model approximation 
  approx_full_model <- ols(approx_full_formula, sigma=1, data = data)
  
  model_approximation_table <- fastbw(approx_full_model, aics=10000) 
  
  
  #------- Bootstrap resampling 
  test <- model_selection(model=model_formula, data = data, boot=bootstraps, return="all", pval_pair=0.01, seed=seed)
  
  raw_table <- as.data.frame(test$raw_overview) %>% rownames_to_column("variable")
  formatted_table <- test$overview
  boot_median_c_stat = test$c_stat_boot
  
  # Save coefficients from selected model and boot model
  coef_selected_model <- test$raw_overview[,"sel_est"]
  coef_boot_median <- test$raw_overview[,"boot_median"]
  
  #------- Selected model
  # Extract names of variables retained in selected model
  selected_vars_int <- raw_table %>% filter(coef_selected_model!=0) %>% pull(variable)
  #selected_vars_int <- names(coef_selected_model[which(coef_selected_model!=0)]) # old code
  selected_vars_int <- gsub("Former|Never|Current|Under|Normal|Over|25-35|35-45|45-55|55+", "", selected_vars_int)
  selected_vars_int <- gsub(")age|)age'", ")", selected_vars_int)
  selected_vars_int <- gsub("4)'", "4)", selected_vars_int)
  selected_vars_int <- unique(selected_vars_int)
  selected_vars <- selected_vars_int[-1]
  
  selected_form <- as.formula(paste0("unsuccessful ~", (paste(selected_vars, collapse = "+"))))
  
  selected_model <- lrm(selected_form, data = data, x = T, y = T)
  predvals_selected <- predict(selected_model)
  
  # Performance
  selected_val <- validate(selected_model, B=bootstraps, seed=seed)
  perf_selected_model <- save_val_results(selected_val)
  c_selected_model <- c_ci(predvals_selected, outcome_vector)
  
  # ------ Bootstrap model 
  
  # Extract names of variables retained in boot model
  boot_vars_int <- raw_table %>% filter(boot_inclusion >= threshold) %>% pull(variable)
  # boot_vars_int <- names(coef_boot_median[which(coef_boot_median!=0)])
  boot_vars_int <- gsub("Former|Never|Current|Under|Normal|Over|25-35|35-45|45-55|55+", "", boot_vars_int)
  boot_vars_int <- gsub(")age|)age'", ")", boot_vars_int)
  boot_vars_int <- gsub("4)'", "4)", boot_vars_int)
  boot_vars_int <- unique(boot_vars_int)
  boot_vars <- boot_vars_int[-1]
  
  boot_form <- as.formula(paste0("unsuccessful ~", (paste(boot_vars, collapse = "+"))))
  
  boot_model <- lrm(boot_form, data = data, x = T, y = T)
  predvals_boot_model <- predict(boot_model)
  
  # Coefficients from the model refit including only variables from the threshold of bootstrap %
  coef_boot_vars <- boot_model$coefficients
  
  # Performance
  boot_val <- validate(boot_model, B=bootstraps, seed=seed)
  perf_boot_model <- save_val_results(boot_val)
  c_boot_model <- c_ci(predvals_boot_model, outcome_vector)
  
  # C statistic - BOOT model (optimism corrected)
  c_stat_boot <- c_stat_ci(model=boot_model, data=data)
  
  # Calibration plot function
  data$predvals_boot_model <- predvals_boot_model
  calibration_plot_boot <- cal_plot_solo(data, predvals_boot_model, outcome_vector)
  
  
  # ------- Shrinkage 
  # Shrink coefficients from boot vars model (vars selected in threshold of boostraps) by the heuristic shrinkage factor
  
  # Heuristic shrinkage factor  
  shrinkage_factor <- (boot_model$stats["Model L.R."] - boot_model$stats["d.f."]) / (boot_model$stats["Model L.R."])
  
  # Predicted values from boot_model shrunk based on heuristic shrinkage factor
  predvals_shrink <- predict(boot_model)*shrinkage_factor
  predrisk_shrink <- 1/(1+exp(-predvals_shrink))
  
  # Shrink coefficients from threshold % of the bootstraps by shrinkage factor
  coef_shrink_model <- coef_boot_vars*shrinkage_factor
  
  data$predrisk_shrink = predrisk_shrink
  data$predvals_shrink = predvals_shrink
  data$outcome = outcome
  
  # Performance
  shrink_perf <- ev_glm(method="original", lp="predvals_shrink", outcome="outcome", data=data, samples=200, return="performance")
  
  brier <- shrink_perf["Brier", "apparent"]
  slope <- shrink_perf["Slope", "apparent"]
  intercept <- shrink_perf["Intercept", "apparent"]
  
  perf_shrink<- rbind("Brier score" = brier, "Calibration slope" = slope, "Calibration intercept" = intercept)
  
  c_shrink <- c_ci(predvals_shrink, outcome_vector)
  
  
  # ROC curve
  roc_shrink <- suppressMessages(roc(data, outcome, predrisk_shrink, ci=TRUE))
  
  
  # ----- Approximated model performance with boot vars 
  
  approx_boot_formula = as.formula(paste0("predvals_full_model ~ ",  (paste(boot_vars, collapse = "+"))))
  
  approx_boot_model <- ols(approx_boot_formula, data = data, x=TRUE, y=TRUE)
  
  # Save data from approximated model - this tells how well the boot variables predict the predicted values from the full model
  approx_boot_model_r2 <- approx_boot_model$stats["R2"]
  
  # Save predicted values
  predvals_approx_mod <- predict(approx_boot_model)
  data$predvals_approx_mod <- predvals_approx_mod
  
  # Performance
  approx_perf <- ev_glm(method="original", lp="predvals_approx_mod", outcome="outcome", data=data, samples=200, return="performance")
  
  brier <- approx_perf["Brier", "apparent"]
  slope <- approx_perf["Slope", "apparent"]
  intercept <- approx_perf["Intercept", "apparent"]
  
  perf_approx_model <- rbind("Brier score" = brier, "Calibration slope" = slope, "Calibration intercept" = intercept)
  c_approx_model <- c_ci(predvals_approx_mod, outcome_vector)
  
  
  # Coefficients from the approximated model
  coef_approx_model <- approx_boot_model$coefficients
  
  #------ Compare coefficients 
  
  coef_full <- data.frame(coef_full_model) %>% rownames_to_column("variable") %>% rename(full_model = coef_full_model) %>% 
    mutate(variable = gsub("[(]|[)]|=", "", variable))
  
  coef_selected <- data.frame(coef_selected_model)  %>% rownames_to_column("variable") %>% rename(selected_model = coef_selected_model) %>% 
    mutate(variable = gsub("[(]|[)]|=", "", variable),
           variable = gsub("rcsage, 4", "", variable))
  
  coef_boot_median <- data.frame(coef_boot_median)  %>% rownames_to_column("variable") %>% rename(coef_boot_median = coef_boot_median) %>% 
    mutate(variable = gsub("[(]|[)]|=", "", variable),
           variable = gsub("rcsage, 4", "", variable))
  
  coef_boot_model <- data.frame(coef_boot_vars)  %>% rownames_to_column("variable") %>% rename(coef_boot_model = coef_boot_vars) %>% 
    mutate(variable = gsub("[(]|[)]|=", "", variable))
  
  coef_shrink <- data.frame(coef_shrink_model)  %>% rownames_to_column("variable") %>% rename(coef_shrink = coef_shrink_model) %>% 
    mutate(variable = gsub("[(]|[)]|=", "", variable))
  
  coef_approx <- data.frame(coef_approx_model)  %>% rownames_to_column("variable") %>% rename(approx_model = coef_approx_model)  %>% 
    mutate(variable = gsub("[(]|[)]|=", "", variable))
  
  
  var_order <- coef_boot_median %>% pull(variable)
  
  coefficients_table <- Reduce(function(x,y) merge(x, y, by = "variable", all.x = TRUE, all.y = TRUE),
                               list(coef_full, coef_selected, coef_boot_median, coef_boot_model, coef_shrink, coef_approx)) %>% 
    mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
    arrange(factor(variable, var_order))
  
  
  #-------- Compare performance 
  
  # -- C statistics ---
  c_stat_list <- as.data.frame(rbind(c_full_model, c_selected_model, c_boot_model, c_shrink, c_approx_model))
  
  rownames(c_stat_list) <- c("full_model", "selected_model", "boot_model", "shrink", "approx")
  
  c_stats <- c_stat_list %>% 
    rownames_to_column("model") %>% 
    mutate_if(is.numeric, round, 2) %>% 
    mutate(c_ci = paste(c_stat, " (", LB, ", ", UB, ")", sep="")) %>% 
    select(model, c_ci)
  
  #--- Brier score, slope, and intercept ---
  performance_list <- as.data.frame(t(cbind(perf_full_model, perf_selected_model, perf_boot_model, perf_shrink, perf_approx_model)))
  
  rownames(performance_list) <- c("full_model", "selected_model",  "boot_model", "shrink", "approx")
  
  performance <- performance_list %>%
    rownames_to_column("model") %>% 
    rename(brier_score = "Brier score",
           cal_slope = "Calibration slope",
           cal_int = "Calibration intercept")
  
  
  #---- Combine data ---
  models_performance <- performance %>% 
    full_join(c_stats, by="model") 
  
  t_model_performance <- as.data.frame(t(models_performance)) %>% rownames_to_column("stat")
  
  
  # ----- Output 
  # C-statistics from:
  #   - Full model
  #   - boot model 
  return(list(
    "coefficients_table" = coefficients_table,
    "models_performance" = t_model_performance,
    "c_statistics" = list(full_model = c_stat_full, boot_model = c_stat_boot),
    "boot_calibration" = calibration_plot_boot,
    "roc_shrink" = roc_shrink,
    "raw_overview " = raw_table,
    "formatted_overview" = formatted_table,
    "r2_approx" = approx_boot_model_r2,
    "boot_model" = boot_model)
  )
}


