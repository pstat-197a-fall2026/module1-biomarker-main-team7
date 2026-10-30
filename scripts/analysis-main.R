# Questions 1 and 2 were completed by Stacy.
# Question 3 was completed by Lawrence.
# Question 4 was completed by Praveen.
# Each person set their own seeds in their analysis. 
# To match the results that each person observed, each section has its own seed.

library(tidyverse)
library(dplyr)
library(glue)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(rsample)
library(glmnet)

# Question 1:
set.seed(1)
getwd()
# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())

biomarker_clean_without_log_transform_outlier_trim <- read_csv('data/biomarker-raw.csv', 
                                                               skip = 2,
                                                               col_select = -2L,
                                                               col_names = c('group', 
                                                                             'empty',
                                                                             pull(var_names, abbreviation),
                                                                             'ados'),
                                                               na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  select(group, ados, everything())

random_protein_sample <- sample(setdiff(names(biomarker_clean_without_log_transform_outlier_trim), c("group", "ados")), size = 5)

random_protein_sample_tibble <- biomarker_clean_without_log_transform_outlier_trim |> 
  select(c("group", "ados"), random_protein_sample)

for (protein in random_protein_sample) {
  hist(random_protein_sample_tibble[[protein]], main = glue("{protein} distribution"), xlab = glue("{protein} values"))}

# Question 2: 
proteins_only <- biomarker_clean_without_log_transform_outlier_trim %>% select(-c(group, ados))
means <- colMeans(proteins_only , na.rm = TRUE)
sds   <- sapply(proteins_only , sd, na.rm = TRUE)

outliers_per_subject <- biomarker_clean_without_log_transform_outlier_trim  %>%
  rowwise() %>%
  mutate(outliers = sum(abs(c_across(names(proteins_only)) - means) > 3 * sds, na.rm = TRUE)) %>%
  ungroup()

barplot(outliers_per_subject %>% arrange(outliers %>% desc()) %>% pull(outliers), ylab = "Number of Outliers", xlab = "Participants", main = "Number of Outlier Protein Levels per Participant")

# Question 3
set.seed(101422)
load('data/biomarker-clean.RData')

# splitting data
split_biomarker_clean <- biomarker_clean %>%
  initial_split(strata = group, prop = 0.75)

train_biomarker_clean <- training(split_biomarker_clean)
test_biomarker_clean <- testing(split_biomarker_clean)

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- train_biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 15) %>%
  pull(protein)

# store predictors and response separately
predictors <- train_biomarker_clean %>%
  select(-c(group, ados))

response <- train_biomarker_clean %>% pull(group) %>% factor()

# fit RF
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 15) %>%
  pull(protein)

# select subset of interest
proteins_sstar <- unique(c(proteins_s1, proteins_s2))

# preprocessing function
prep_func <- function(.df){
  .df %>%
    select(group, any_of(proteins_sstar)) %>%
    mutate(class = (group == 'ASD')) %>%
    select(-group)
}

# prep_func applied to train and test
train_biomarker_sstar <- prep_func(train_biomarker_clean)
test_biomarker_sstar <- prep_func(test_biomarker_clean)

fit <- glm(class ~ ., data = train_biomarker_sstar, family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

test_biomarker_sstar %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = factor(pred > 0.5, levels = c(FALSE, TRUE)),
         tr = factor(class, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = tr,
                estimate = est,
                pred,
                event_level = 'second')

# Question 4
set.seed(67)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.8, strata = group)

train_data <- training(biomarker_split)
test_data <- testing(biomarker_split)

# function to compute tests
test_fn <- function(.df) {
  t_test(.df,
         formula = level ~ group,
         order = c("ASD", "TD"),
         alternative = "two-sided",
         var.equal = FALSE)
}

ttests_out <- train_data %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group,
               names_to = "protein",
               values_to = "level") %>%
  # nest by protein
  nest(data = c(level, group)) %>%
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1 / (2 * m) - digamma(1),
         rank = row_number(),
         p.adj = m * hm * p_value / rank)

# top 10 proteins by t-test
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

# store predictors and response separately
predictors <- train_data %>%
  select(-c(group, ados))

response <- train_data %>% pull(group) %>% factor()

# fit RF
set.seed(67)
rf_out <- randomForest(x = predictors,
                       y = response,
                       ntree = 1000,
                       importance = TRUE)

# top 10 proteins by RF importance
proteins_s2 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

# hard intersection like in class
proteins_sstar <- intersect(proteins_s1, proteins_s2)
proteins_sstar

# metrics to report
class_metrics <- metric_set(sensitivity,
                            specificity,
                            accuracy,
                            roc_auc)

# make the factor columns first, then compute metrics
eval_fn <- function(.df) {
  .df %>%
    mutate(truth = factor(class, levels = c(FALSE, TRUE)),
           estimate = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
    class_metrics(truth = truth,
                  estimate = estimate,
                  pred,
                  event_level = "second")
}

fit_panel <- function(proteins) {
  train_sub <- train_data %>%
    select(group, all_of(proteins)) %>%
    mutate(class = (group == "ASD")) %>%
    select(-group)
  
  test_sub <- test_data %>%
    select(group, all_of(proteins)) %>%
    mutate(class = (group == "ASD")) %>%
    select(-group)
  
  fit <- glm(class ~ ., data = train_sub, family = "binomial")
  
  test_sub %>%
    add_predictions(fit, type = "response") %>%
    eval_fn()
}

baseline_metrics <- fit_panel(proteins_sstar)
baseline_metrics

baseline_results <- baseline_metrics %>%
  select(.metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(method = "in-class (hard intersection)",
         n_proteins = length(proteins_sstar),
         proteins = paste(proteins_sstar, collapse = ", ")) %>%
  select(method, n_proteins, everything())

baseline_results

# order the baseline proteins by t-test p-value
proteins_ordered <- ttests_out %>%
  filter(protein %in% proteins_sstar) %>%
  arrange(p_value) %>%
  pull(protein)

# fit models using the top 1, 2, ..., k proteins
smaller_panels <- tibble(n_proteins = 1:length(proteins_ordered)) %>%
  mutate(metrics = map(n_proteins, ~ fit_panel(proteins_ordered[1:.x]))) %>%
  unnest(metrics) %>%
  select(n_proteins, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

smaller_panels

# 5-fold cross validation on the training data
train_folds <- vfold_cv(train_data, v = 5, strata = group)

# average cross validation AUC for the top k proteins
cv_auc <- function(k) {
  proteins <- proteins_ordered[1:k]
  
  fold_auc <- map_dbl(train_folds$splits, function(.split) {
    fold_train <- analysis(.split) %>%
      select(group, all_of(proteins)) %>%
      mutate(class = (group == "ASD")) %>%
      select(-group)
    
    fold_val <- assessment(.split) %>%
      select(group, all_of(proteins)) %>%
      mutate(class = (group == "ASD")) %>%
      select(-group)
    
    fit <- glm(class ~ ., data = fold_train, family = "binomial")
    
    fold_val %>%
      add_predictions(fit, type = "response") %>%
      eval_fn() %>%
      filter(.metric == "roc_auc") %>%
      pull(.estimate)
  })
  
  mean(fold_auc)
}

cv_results <- tibble(n_proteins = 1:length(proteins_ordered)) %>%
  mutate(cv_auc = map_dbl(n_proteins, cv_auc))

cv_results

n_simple <- cv_results %>%
  filter(cv_auc >= max(cv_auc) - 0.05) %>%
  slice_min(n_proteins, n = 1) %>%
  pull(n_proteins)

proteins_simple <- proteins_ordered[1:n_simple]
proteins_simple

simple_results <- smaller_panels %>%
  filter(n_proteins == n_simple) %>%
  mutate(method = "simpler panel",
         proteins = paste(proteins_simple, collapse = ", ")) %>%
  select(method, n_proteins, everything())

simple_results

x_train <- train_data %>% select(-c(group, ados)) %>% as.matrix()
y_train <- train_data %>% pull(group) == "ASD"

x_test <- test_data %>% select(-c(group, ados)) %>% as.matrix()

# pick lambda with 10-fold cross validation
set.seed(67)
lasso_cv <- cv.glmnet(x_train, y_train,
                      family = "binomial",
                      type.measure = "auc")

# proteins with nonzero coefficients
lasso_coef <- coef(lasso_cv, s = "lambda.1se")
proteins_lasso <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
proteins_lasso <- setdiff(proteins_lasso, "(Intercept)")
proteins_lasso

lasso_metrics <- test_data %>%
  mutate(class = (group == "ASD"),
         pred = predict(lasso_cv, newx = x_test,
                        s = "lambda.1se", type = "response")[, 1]) %>%
  eval_fn()

lasso_results <- lasso_metrics %>%
  select(.metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(method = "lasso",
         n_proteins = length(proteins_lasso),
         proteins = paste(proteins_lasso, collapse = ", ")) %>%
  select(method, n_proteins, everything())

lasso_results

comparison <- bind_rows(baseline_results, simple_results, lasso_results)

comparison %>% select(-proteins)

run_split <- function(seed) {
  set.seed(seed)
  split <- initial_split(biomarker_clean, prop = 0.8, strata = group)
  train <- training(split)
  test <- testing(split)
  
  # t-test selection
  tt <- train %>%
    select(-ados) %>%
    pivot_longer(-group, names_to = "protein", values_to = "level") %>%
    nest(data = c(level, group)) %>%
    mutate(ttest = map(data, test_fn)) %>%
    unnest(ttest) %>%
    arrange(p_value) %>%
    mutate(m = n(),
           hm = log(m) + 1 / (2 * m) - digamma(1),
           rank = row_number(),
           p.adj = m * hm * p_value / rank)
  s1 <- tt %>% slice_min(p.adj, n = 10) %>% pull(protein)
  
  # RF selection
  rf <- randomForest(x = train %>% select(-c(group, ados)),
                     y = factor(train$group),
                     ntree = 1000,
                     importance = TRUE)
  s2 <- rf$importance %>%
    as_tibble() %>%
    mutate(protein = rownames(rf$importance)) %>%
    slice_max(MeanDecreaseGini, n = 10) %>%
    pull(protein)
  
  sstar <- intersect(s1, s2)
  
  # simpler panel = top n_simple of the intersection by p-value
  simple <- tt %>% filter(protein %in% sstar) %>% pull(protein)
  simple <- simple[1:min(n_simple, length(simple))]
  
  # logistic regression on a panel
  glm_metrics <- function(proteins) {
    train_sub <- train %>%
      select(all_of(proteins)) %>%
      mutate(class = (train$group == "ASD"))
    test_sub <- test %>%
      select(all_of(proteins)) %>%
      mutate(class = (test$group == "ASD"))
    fit <- glm(class ~ ., data = train_sub, family = "binomial")
    test_sub %>%
      add_predictions(fit, type = "response") %>%
      eval_fn()
  }
  
  # lasso
  cv <- cv.glmnet(train %>% select(-c(group, ados)) %>% as.matrix(),
                  train$group == "ASD",
                  family = "binomial",
                  type.measure = "auc")
  lasso <- test %>%
    mutate(class = (group == "ASD"),
           pred = predict(cv,
                          newx = test %>% select(-c(group, ados)) %>% as.matrix(),
                          s = "lambda.1se", type = "response")[, 1]) %>%
    eval_fn()
  
  bind_rows(glm_metrics(sstar) %>% mutate(method = "in-class (hard intersection)"),
            glm_metrics(simple) %>% mutate(method = "simpler panel"),
            lasso %>% mutate(method = "lasso")) %>%
    mutate(seed = seed)
}

repeated_out <- map_dfr(67:86, run_split)

# average and sd of each metric across splits
repeated_summary <- repeated_out %>%
  group_by(method, .metric) %>%
  summarize(mean = mean(.estimate),
            sd = sd(.estimate),
            .groups = "drop")

repeated_summary %>%
  select(-sd) %>%
  pivot_wider(names_from = .metric, values_from = mean)

repeated_summary %>%
  select(-mean) %>%
  pivot_wider(names_from = .metric, values_from = sd)