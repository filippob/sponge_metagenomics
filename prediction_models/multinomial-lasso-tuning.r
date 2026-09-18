library("vip")
library("dplyr")
library("glmnet")
library("themis")
library("ggplot2")
library("parallel")
library("phyloseq")
library("tidyverse")
library("doParallel")
library("tidymodels")
library("data.table")


################################################################################
## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  if (tools::file_ext(args[1]) %in% c("r","R")) {
    
    source(args[1])
    # source("Analysis/config.R")
  } else {
    
    load(args[1])
  }
  
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    prjfolder = "Documents/metagenomics/sponge",
    repo = "Documents/metagenomics/sponge/sponge_metagenomics",
    output_folder = "analysis",
    input_data = "analysis/phyloseq.RData",
    target_variable = "island",
    nproc = 4,
    split_ratio = 0.80,
    method_cv = "repeatedcv",
    k_folds = 5,
    nrepeats_cv = 3,
    sampling_method = "up",
    force_overwrite = FALSE
  ))
}
################################################################################

HOME <- Sys.getenv("HOME")
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$output_folder)

if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

config_fname = file.path(outdir, "lasso.config.RData")
save(config, file = config_fname)

## loading data previously imported in phyloseq
writeLines(" - loading normalised phyloseq data")
fname = file.path(prjfolder, config$input_data)
load(fname)

otu <- as(otu_table(physeq_pruned), "matrix")
total_counts <- rowSums(otu)
prevalence <- rowSums(otu > 0)

keep <- total_counts >= 100 &
  prevalence >= 5

physeq_pruned <- prune_taxa(keep, physeq_pruned)
rm(otu)
gc()

## Import dataset
dataset <- otu_table(physeq_pruned) |> as.data.frame()
dataset <- t(as.matrix(dataset))
keep <- colMeans(dataset > 0) >= 0.10
dataset <- dataset[, keep]

dataset <- as.data.frame(dataset)
## subset random columns (thinning)
cols <- sample(names(dataset), 9000)
dataset <- dataset %>%
  select(all_of(cols))

dataset$sample_id = row.names(dataset)

metadata <- data.frame(sample_data(physeq_pruned))
metadata$sample_id = row.names(metadata)
metadata <- metadata |> select(sample_id,!!config$target_variable)

dataset <- metadata |> inner_join(dataset, by = "sample_id")

## Data preparation
dataset <- dataset %>% 
  mutate(
    !!config$target_variable := as.factor(.data[[config$target_variable]]), 
    sample_id = as.factor(sample_id)
    )

if(sum(is.na(dataset)) == 0) print("No missing data in the dataset: OK!")

## Setup parallel backend
writeLines(" - parallelising")
cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)

writeLines(" - splitting the data: training/validation/test")
dt <- select(dataset, -c(sample_id))
dt_split <- initial_split(dt, strata = !!config$target_variable, prop = config$split_ratio)
dt_train <- training(dt_split)
dt_test <- testing(dt_split)

#### Preprocessing

# We use tidymodels to build a recipe for data preprocessing:
#   
# -   remove correlated variables
# -   remove non informative variables (zero variance)
# -   standardize all variables
# -   impute missing data (Random Forest does not handle missing data)

lasso_recipe <- dt_train %>%
  recipe(reformulate(".", response = config$target_variable)) |>
  step_zv(all_predictors(), -all_outcomes()) %>%
  step_nzv(all_predictors()) |>
  step_normalize(all_numeric(), -all_outcomes()) |>
  step_upsample(all_outcomes(), over_ratio = 1.25)

temp <- recipes::prep(lasso_recipe) %>% juice()

#### Model building

# We now specify the structure of our model:
#   
# -   hyperparameters to tune: `mtry` (number of features to sample for each tree) and `min_n` (minimum number of data points in a node to allow further splitting)
# -   number of trees in the forest
# -   the problem at hand (classification)
# -   the engine (R package)
# 
# Then we put this in a workflow together with the preprocessing recipe

tune_spec <- multinom_reg(
  mode = "classification",
  penalty = tune(),
  mixture = 1      # α = 1 (Lasso)
) %>%
  set_engine("glmnet")

tune_wf <- workflow() %>%
  add_recipe(lasso_recipe) |>
  add_model(tune_spec)


#### Tuning the hyperparameters
writeLines(" - model tuning: hyperparameters")
# We use k-fold cross-validation to tune the hyperparameters in the training set
dt_folds <- vfold_cv(dt_train, v = config$k_folds, repeats = config$nrepeats_cv)

lasso_grid <- grid_regular(
  penalty(range = c(-10, -0.75)),
  levels = 14
)

regular_res <- tune_grid(
  tune_wf,
  metrics = metric_set(accuracy, kap, mcc),
  resamples = dt_folds,
  grid = lasso_grid
)

regular_res |>
  collect_metrics() |>
  filter(.metric == "mcc") |>
  print(n = 20)

library("repr")
p <- regular_res %>%
  collect_metrics() %>%
  # filter(.metric == "mcc") %>%
  ggplot(aes(penalty, mean, color = .metric)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "metric value")

fname = file.path(outdir, "figures/lasso_tuning.png")
ggsave(filename = fname, plot = p, device = "png")

best_auc <- select_best(x = regular_res, metric = "mcc")
show_best(regular_res, metric = "mcc") |> print()

# 2.  finalise the model:
final_lasso <- finalize_model(
  tune_spec,
  best_auc
)

print(final_lasso)

# 3.  finalise the workflow and fit it to the initial split (training and test data):
final_wf <- workflow() %>%
  add_recipe(lasso_recipe) %>%
  add_model(final_lasso)

final_res <- final_wf %>%
  last_fit(dt_split, metrics = metric_set(kap, accuracy, mcc, brier_class))

# 4.  evaluate the fine-tuned RF model:
writeLines(" - evaluating the final model")
print(final_res)
performance_metrics <- final_res %>%
  collect_metrics()

print(performance_metrics)
fname = file.path(outdir, "tables/performance_mertrics.csv")
fwrite(x = performance_metrics, file = fname, sep = ",")

# 5.  get variable importance:
final_res %>% 
  pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>% 
  vip()

coef_tbl <- final_res %>%
  pluck(".workflow", 1) %>%
  extract_fit_parsnip() %>%
  broom::tidy()

coef_tbl <- coef_tbl %>%
  filter(term != "(Intercept)") %>%
  mutate(importance = abs(estimate)) %>%
  arrange(desc(importance))

taxa <- data.frame(tax_table(physeq_pruned))
taxa$term <- row.names(taxa)
taxa <- select(taxa, c(term, Species))

coef_tbl <- coef_tbl %>%
  left_join(taxa, by = "term")

p <- coef_tbl %>%
  slice_max(importance, n = 20) %>%
  ggplot(aes(
    x = reorder(Species, importance),
    y = importance,
    fill = Species
  )) +
  xlab("species") + 
  geom_col() +
  coord_flip() + 
  guides(fill = "none") 

fname = file.path(outdir, "figures/variable_importance.png")
ggsave(filename = fname, plot = p, device = "png")

#### Predictions

# We collect the predictions on the test set: for each test observations we get the probabilities of belonging to each of the four classes.
preds <- final_res %>%
  collect_predictions()

fname = file.path(outdir, "tables/test_predictions.csv")
fwrite(x = preds, file = fname, sep = ",")

cm <- final_res %>%
  collect_predictions() %>%
  conf_mat(
    truth = !!sym(config$target_variable),
    estimate = .pred_class
  )

print(cm)

cm_df <- as.data.frame(cm$table)

p <- ggplot(cm_df, aes(x = Truth, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_discrete(limits = rev) +
  theme_bw()

fname <- file.path(outdir, "figures/confusion_matrix.png")
ggsave(filename = fname, plot = p, device = "png")

print("DONE!!")