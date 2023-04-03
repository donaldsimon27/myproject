#Load packages-------------
library(tidymodels)  
tidymodels::tidymodels_prefer()
options(tidymodels.dark = TRUE)
# Helper packages
library(rpart)                 #needed to use using rpart.plot
library(rpart.plot)            # for visualizing a decision tree
library(vip)                   # for variable importance plots
library(bestNormalize)         #To use step_orderNorm
library(future)
library(foreach)
library(doFuture)            
library(rngtools)             #requires rng tools to call up doRNG
library(doRNG)
registerDoFuture()
library(furrr)
library(parallel)
library(iterators)
library(doParallel)
library(kernlab)
library(themis)
library(ranger)
library(kknn)
library(Matrix)
library(glmnet)
library(C50)
library(rules)
library(pROC)
library(kableExtra)
library(patchwork)


#STEP1: LOAD FINAL BASELINE DATASET -----------
baselineX <- read.csv("~/Desktop/R.Projects/PredictLuminex/Modeling Datasets/baseline_nomsymoptoms.csv") |> 
  type.convert(as.is = FALSE) 


if (supportsMulticore()) {
  plan(multicore, workers=availableCores(omit=1))
} else {
  plan(multisession, workers=availableCores(omit=1))
}                              


#Step2: Remove unnecessary variables, create factors------------
baseline0 <- baselineX %>% 
  select(-1) 

data <-  baseline0 |>                #blx = baselineX
  select(-bodymass)      #remove bodymass since we already BMI in the dataset

#vars_to_factor <- c("CAVNUM_wk0")
#data[vars_to_factor] <- lapply(data[vars_to_factor], function(x) as.factor(x))
#rm(vars_to_factor)

data$Outcome <- data$Outcome |> relevel(ref='PoorOutcome')

#data <- data |> select(-c(1, 3:6, 70:84))      #removes Cavnum and symptoms

#Step3: Baseline recipe--------- 
normalized_recipe <- 
  recipe(Outcome ~ ., data = data) |> 
  update_role(PID, new_role = "ID") |> 
  step_zv(all_predictors()) |> 
  step_orderNorm(all_numeric_predictors()) |> 
  step_normalize(all_numeric_predictors()) |>  #newly_added_03Mar2023
  step_dummy(all_nominal_predictors()) |> 
  themis::step_downsample(Outcome)

#update_role(PID, new_role = "ID") |> 
#step_corr(all_predictors(), threshold = 0.9, method = "spearman")  #Can add as an additional step to recip

predictor_count <- sum(normalized_recipe$term_info$role == 'predictor')       #82 predictors
predictor_count

#Step4: Different Models--------------
C5.0_mod <-
  C5_rules(trees = tune(), min_n = tune()) |> 
  set_engine('C5.0') |> 
  set_mode("classification")
KNN_mod <-
  nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) |> 
  set_engine('kknn') |> 
  set_mode('classification')
Random_Forest_mod <-
  rand_forest(mtry = tune(), trees = tune(), min_n = tune()) |> 
  set_engine('ranger') |> 
  set_mode('classification')
Elastic_net_mod <-
  logistic_reg(penalty = tune(), mixture = tune()) |> 
  set_engine('glmnet') |> 
  set_mode("classification")


#Step5: CREATING THE WORKFLOW SET-----------
normalized_workflow <-
  workflowsets::workflow_set(
    preproc = list(normalised = normalized_recipe),
    models = list(C5.0 = C5.0_mod,
                  KNN = KNN_mod,
                  RF = Random_Forest_mod,
                  EN = Elastic_net_mod))

normalized_workflow <- normalized_workflow |> 
  mutate(wflow_id = gsub("(normalised_)|(numeric_)", "", wflow_id))           #to remove the normalised prefix


#Step6: Nested Crossvalidation---------------
#Split in training and and test set ??? D/W 
#??LOOCV

set.seed(14193)
folds <- nested_cv(data,
                   outside = vfold_cv(v = 10, repeats = 1, strata = "Outcome"),   #Fitting on this
                   inside = bootstraps(times = 50, strata = "Outcome"))
# inside = vfold_cv(v=5, repeats = 1, strata = TB))


#Step5: Model parameters and workflows------------
C5.0_params <- parameters(trees(range = c(1,100)), min_n())
RF_params <- parameters(mtry(range=c(1, predictor_count)), trees(range = c(1,100)), min_n()) 
KNN_params <- parameters(neighbors(), weight_func(), dist_power()) |> 
  recipes::update(weight_func = weight_func(c("triangular", "biweight", "triweight", "cos", "gaussian", "rank", "optimal")))
EN_params <- parameters(penalty(), mixture())



#Step6: Workflows---------
workflows <- normalized_workflow |>
  option_add(param_info = C5.0_params, id = "C5.0") |> 
  option_add(param_info = KNN_params, id = "KNN") |> 
  option_add(param_info = RF_params, id = "RF") |> 
  option_add(param_info = EN_params, id = "EN")


#Step7: Tuning-------------
bayes_ctrl <- control_bayes(no_improve = 15L, 
                            save_pred = TRUE, 
                            parallel_over = "resamples",    #"everything" -  success when parallel =  "resamples" instead of parallel =  "everything" is used 
                            save_workflow = TRUE, 
                            allow_par = TRUE, 
                            verbose = TRUE)


#tune_results <- foreach(i=1:length(folds$splits)) %dorng% {
tune_results <- foreach(i=1:length(folds$splits)) %do% {
  library(rules)
  workflows %>% workflow_map(seed = 15440,
                             fn = "tune_bayes",
                             resamples = folds$inner_resamples[[i]],
                             metrics = metric_set(roc_auc),
                             objective = exp_improve(),
                             iter = 50,
                             control = bayes_ctrl)
}


saveRDS(tune_results, "~/Desktop/R.Projects/PredictLuminex/Data/Exported Data/tune_results.rds")


#Step8: Tune results object----------
tune_results <- readRDS("~/Desktop/R.Projects/PredictLuminex/Data/Exported Data/tune_results.rds")


#Step9: Normalized recipe----------
normalised_training_recipe <- 
  recipe(Outcome ~ ., 
         data = data) |>
  update_role(PID, new_role = "ID") |>
  step_normalize(all_numeric_predictors()) |>
  step_dummy(all_nominal_predictors()) |> 
  step_smote(Outcome) 

#step_corr(all_predictors(), threshold = 0.9, method = "spearman")  #Can add as an additional step to recipe

#Step10: Fit models function-------
fit_models <- function(model, grid, data, type){
  
  best_result <-  grid |> 
    show_best(n=1)
  
  model_fit <-
    grid |> extract_workflow(model) |>
    finalize_workflow(best_result) |>
    update_recipe(normalised_training_recipe) |> 
    fit(data=analysis(data))
  
  pred <-
    predict(model_fit, assessment(data)) |> 
    bind_cols(predict(model_fit, assessment(data), type = "prob")) |> 
    bind_cols(assessment(data) |> select(Outcome))
  
  roc_auc <- pred |> roc_auc(truth = Outcome, .pred_PoorOutcome )
  sens <- pred|> sensitivity(Outcome, .pred_class)
  spec <- pred |> specificity(Outcome, .pred_class)
  
  perf <- tibble(Model = model, data ="train") |> 
    mutate(auc = roc_auc$.estimate, sens =  sens$.estimate, spec = spec$.estimate)
  
  
  if (type==1){
    return(perf)
  }
  else if(type==2){
    pred <- pred |> mutate(Resamples = data$id)
    return(pred)
  }
  else if(type==3){
    return(best_result)
  }
}


#Step9:Fit------------
training_results <- foreach(x=1:length(folds$splits)) %do% {
  model_results <- foreach(y=1:length(workflows$wflow_id)) %do% {
    fit_models(tune_results[[x]]$wflow_id[[y]], tune_results[[x]]$result[[y]], folds$splits[[x]], 1)
  }
  bind_rows(model_results)
}

training_average <- foreach(x=1:length(workflows$wflow_id), .combine=rbind) %dorng% {
  bind_rows(training_results) |> filter(Model==workflows$wflow_id[[x]]) |>
    group_by(Model, data) |> summarise_if(is.numeric, list(mean = mean, sd = sd), na.rm=TRUE)
}

training_average[3:8] <- training_average[3:8] |> round(3)

new_training_average <- foreach(x=1:length(workflows$wflow_id)) %do% {
  tibble(Model = workflows$wflow_id[[x]], data = "train") |>
    mutate(auc = paste0(training_average$auc_mean[[x]], "\u00B1", training_average$auc_sd[[x]]),
           sens = paste0(training_average$sens_mean[[x]], "\u00B1", training_average$sens_sd[[x]]),
           spec = paste0(training_average$spec_mean[[x]], "\u00B1", training_average$spec_sd[[x]]))
} |> bind_rows()

arrange(new_training_average, desc(auc))|> 
  kable(align=rep('c')) |> 
  kable_classic(full_width = F)

#Step10: Predictions----------
prediction_results <- foreach(x=1:length(folds$splits)) %do% {
  model_perf <- foreach(y=1:length(workflows$wflow_id)) %do% {
    fit_models(tune_results[[x]]$wflow_id[[y]], tune_results[[x]]$result[[y]], folds$splits[[x]], 2)
  }
}

roc_curves_folds <- foreach(x=1:length(workflows$wflow_id)) %do% {
  model_pred <- foreach(y=1:length(folds$splits)) %do% {
    prediction_results[[y]][[x]] |> 
      unnest(cols = c(Resamples))|> 
      rename(Resamples=id)
  }
  bind_rows(model_pred) |> 
    group_by(Resamples) |> 
    roc_curve(Outcome, .pred_PoorOutcome) |> 
    autoplot() + 
    labs(title = workflows$wflow_id[[x]]) + 
    theme_update(plot.title = element_text(hjust = 0.5)) +
    theme_bw() +
    theme(legend.position="left")        #alternative legend = "none
}


roc_curves_folds[[1]] + roc_curves_folds[[2]] + roc_curves_folds[[3]] + roc_curves_folds[[4]] + 
  plot_layout(nrow = 2, byrow = FALSE)


roc_curves <- foreach(x=1:length(workflows$wflow_id)) %do% {
  model_pred <- foreach(y=1:length(folds$splits), .combine = rbind) %do% {
    prediction_results[[y]][[x]] |> 
      unnest(cols = c(Resamples))|> 
      rename(Resamples=id)
  }
  pROC::roc(Outcome ~ .pred_PoorOutcome, data=model_pred, auc=T,levels=c('Cured','PoorOutcome'), ci=TRUE, of = "se", ci.type = "bars", ci.method = "bootstrap", boot.n = 5000, parallel = TRUE, plot=FALSE)
}
par(mfrow=c(2,2))
for (i in 1:4) {
  plot(roc_curves[[i]], main=workflows$wflow_id[i], print.auc=T)
}
par(mfrow=c(1,1))       #if you want individual plots