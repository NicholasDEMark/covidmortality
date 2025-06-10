library(svMisc)
library(dplyr)
library(CIPerm)

##### Prepare Data #####
# read data
df <- read.csv("/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_2006-2023.csv", 
               header=TRUE)

# rename
colnames(df) <-  c("cohort","cumu_month","cyear","cmonth","mortality")
# cumulative months
df <- df %>%
  arrange(cohort,cyear,cmonth) %>%
  group_by(cyear,cmonth) %>%
  mutate(cumu_month=cur_group_id())

##### Functions #########
# basic estimates with CI ##### 
# the data need to have 4 variables
# 1. cohort, 2. calendar year, 3. calendar month, 4. mx by month
# variable names should be following (cohort, cyear, cmonth,  mort)
# the function returns 1. cohort-specific birth loss (in absolute # and %)
# 2. the CI of the loss 
# 3. residuals from regression (a list object with length = # months * smooth cohort lag years)
# 4. mortality effect for each month for each cohort (a list object with length = # months)

# initialize
set.seed(978653421)
# years/months months during covid
cv_yr_beg <- 2020
cv_mn_beg <- 3
cv_yr_end <- 2020
cv_mn_end <- 12
t_delta   <- 12*cv_yr_end + cv_mn_end - (12*cv_yr_beg + cv_mn_beg) + 1
coh_beg   <- 1932
coh_end   <- 2013
n_coh     <- coh_end-coh_beg + 1
coh_lag1  <- 1

# set placebo quantities
plcb_n    <- 12
plcb_ybeg <- 2019
plcb_mbeg <- 1
ci_alpha  <- .99

cohort_discontinuity <-  function(data,
                                  shock_year_start=cv_yr_beg,
                                  shock_month_start=cv_mn_beg,
                                  shock_year_end=cv_yr_end,
                                  shock_month_end=cv_mn_end,
                                  cohort_start=coh_beg,
                                  cohort_end=coh_end,
                                  smooth_cohort_lag1=coh_lag1, ## c-2
                                  smooth_cohort_lag2=coh_lag2, ## c-6
                                  ci=ci_alpha,
                                  weighting=c("none", "triangular", "gaussian")) {
  
  # Match weighting argument
  weighting <- match.arg(weighting)
  
  # Cumulative months
  data <- data %>%
    arrange(cohort,cyear,cmonth) %>%
    group_by(cyear,cmonth) %>%
    mutate(cumu_month=cur_group_id())
  
  # Locate the starting and ending years
  cumu_start <- unique(data[data$cyear==shock_year_start&data$cmonth==shock_month_start,"cumu_month"])[[1]]
  cumu_end <- unique(data[data$cyear==shock_year_end&data$cmonth==shock_month_end,"cumu_month"])[[1]]
  
  # Kernel weighting function
  kernel_weight <- function(x, max_dist) {
    if (weighting == "none") {
      # No weighting - return constant weight of 1
      return(rep(1, length(x)))
    } else if (weighting == "triangular") {
      # Triangular kernel: weight decreases linearly with distance
      return(1 - abs(x) / max_dist)
    } else if (weighting == "gaussian") {
      # Gaussian kernel: weight decreases exponentially with distance
      return(exp(-0.5 * (x / max_dist)^2))
    }
  }
  
  # Loop over cohorts to measure counterfactuals and store residuals
  predict <- data.frame()
  for (i in 1:length(unique(data$cohort))) {
    # Locate target cohort
    c <- unique(data$cohort)[i]
    if (c >= cohort_start & c <= cohort_end) {
      # Targeted cohort and months whose counterfactuals are predicted
      target <- data[data$cohort==c & data$cumu_month>=cumu_start & data$cumu_month<=cumu_end,]
      target$fe <- seq(1, cumu_end-cumu_start+1)
      
      # Construct reference data (simplified for unweighted case)
      construct <- data.frame()
      for (t in smooth_cohort_lag1:smooth_cohort_lag2) {
        construct <- rbind(
          construct,
          data[data$cohort==unique(target$cohort)-t &
                 data$cyear %in% unique(target$cyear)-t &
                 data$cumu_month %in% (target$cumu_month-12*t),]
        )
      }
      
      # Only apply weighting calculations for weighted methods
      if (weighting != "none") {
        construct <- construct %>%
          arrange(cohort, cumu_month) %>%
          group_by(cohort) %>%
          mutate(
            fe = 1:n(),
            # Calculate distance from shock period
            dist_from_shock = abs(cumu_month - cumu_start),
            # Maximum distance for normalization
            max_dist = max(abs(construct$cumu_month - cumu_start))
          )
        
        # Apply kernel weighting
        construct$weight <- kernel_weight(construct$dist_from_shock, 
                                          construct$max_dist)
      } else {
        # For unweighted, just add a simple group identifier
        construct <- construct %>%
          arrange(cohort, cumu_month) %>%
          group_by(cohort) %>%
          mutate(fe = 1:n())
      }
      
      # Regression with optional weighting
      if (weighting == "none") {
        model <- lm(mortality ~ cumu_month + factor(fe), data = construct)
      } else {
        model <- lm(mortality ~ cumu_month + factor(fe), 
                    data = construct, 
                    weights = weight)
      }
      
      # Save prediction
      target$counterfactual <- predict.lm(model, target)
      
      # Save residuals
      target$residuals <- list(as.vector(model$residuals))
      predict <- rbind(predict, target)
    }
    
    # Monitor progress
    # progress(i, length(unique(data$cohort)))
  }
  
  
  # Rest of the function remains the same as in the original code
  # (ATE calculation, confidence intervals, etc.)
  cohort_ate <- predict %>%
    mutate(ate=mortality-counterfactual) %>%
    group_by(cohort) %>%
    summarize(delta_c=sum(ate),
              observed=sum(mortality),
              residuals_regression=list(unlist(unique(residuals))),
              delta_c_list=list(ate)) %>%
    mutate(percent_diff=delta_c/observed*100) %>%
    arrange(desc(cohort))
  
  # Confidence interval calculation (unchanged from original)
  cohort_ate$lower <- NA
  cohort_ate$higher <- NA
  for (i in 1:nrow(cohort_ate)) {
    CI <- cint(dset(cohort_ate[i,]$delta_c_list[[1]],
                    cohort_ate[i,]$residuals_regression[[1]]),
               conf.level = 1-(1-ci),#/(cohort_end-cohort_start+1), 
               tail = "Two")$conf.int*(cumu_end-cumu_start+1)
    
    cohort_ate[i,]$lower <- CI[1]
    cohort_ate[i,]$higher <- CI[2]
  }
  
  cohort_ate <- cohort_ate %>%
    mutate(percent_lower = lower/observed*100,
           percent_higher = higher/observed*100)
  
  cohort_ate <- cohort_ate %>%
    dplyr::select(cohort, delta_c, lower, higher, percent_diff,
                  percent_lower, percent_higher,
                  observed, residuals_regression, delta_c_list) %>%
    dplyr::rename(lower_c = lower,
                  higher_c = higher)
  
  return(cohort_ate)
}

# Loop over different coh_lag2 values
coh_lag2_values <- c(5, 6)

for (lag_value in coh_lag2_values) {
  # Set coh_lag2 for this iteration
  coh_lag2 <- lag_value
  
  # Print current value being processed
  cat("Processing coh_lag2 =", coh_lag2, "\n")
  
  # Use the function to derive results of cohort loss
  cohort_ate <- cohort_discontinuity(data=df,
                                     shock_year_start=cv_yr_beg,
                                     shock_month_start=cv_mn_beg,
                                     shock_year_end=cv_yr_end,
                                     shock_month_end=cv_mn_end,
                                     cohort_start=coh_beg,
                                     cohort_end=coh_end,
                                     smooth_cohort_lag1=coh_lag1, ## c-2
                                     smooth_cohort_lag2=coh_lag2, ## c-6
                                     ci=ci_alpha,
                                     weighting = "none")
  
  # Write estimate and percent change with CIs
  diff <- cohort_ate %>% select(-c("residuals_regression","delta_c_list"))
  # Create filename with the lag value
  filename <- paste0("/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/female_US", lag_value, ".csv")
  write.csv(diff, filename, row.names = FALSE)
  
  
  # Save the cohort_ate to a variable with the corresponding name
  assign(paste0("cohort_ate", lag_value), cohort_ate)
  
  cat("Completed coh_lag2 =", coh_lag2, "\n")
}

cat("All iterations completed.\n")


require(foreign)
ate5 <- subset(cohort_ate5, select = c(cohort, delta_c, lower_c, higher_c, percent_diff, percent_lower, percent_higher, observed))
ate6 <- subset(cohort_ate6, select = c(cohort, delta_c, lower_c, higher_c, percent_diff, percent_lower, percent_higher, observed))
write.dta(ate5, "/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/us_female_5.dta")
write.dta(ate6, "/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/us_female_6.dta")




################################

##### placebo diagnostics #####

## the function returns a list with 2 elements:
## 1. the main upper and lower-bound estimate of the bias for each cohort based on n placebo treatments
## 2. the full results of each placebo simulation, including percent difference of mortality its monthly variation (raw)
coh_lag1 <- 1
coh_lag2 <- 5

placebo_diagnostics <-
  function(data=df,
           placebo_total_months=plcb_n,
           placebo_year_start=plcb_ybeg,
           placebo_month_start=plcb_mbeg,
           treatment_delta=t_delta,
           cohort_start=coh_beg,
           cohort_end=coh_end,
           smooth_cohort_lag1=coh_lag1, ## c-2
           smooth_cohort_lag2=coh_lag2  ## c-6
  ) {
    
    ## loop over cohorts over placebo years
    placebo_results <- data.frame()
    
    for (p in 1:placebo_total_months){
      cumu_start <- seq(unique(data[data$cyear==placebo_year_start&
                                      data$cmonth==placebo_month_start,
                                    "cumu_month"])[[1]],
                        unique(data[data$cyear==placebo_year_start&
                                      data$cmonth==placebo_month_start,
                                    "cumu_month"])[[1]]+placebo_total_months-1,1)[p]
      cumu_end <- cumu_start + t_delta - 1
      
      ## loop over cohorts to measure counterfactuals and store residuals
      predict <- data.frame()
      for (i in 1:length(unique(data$cohort))){
        
        ## locate target cohort
        c <- unique(data$cohort)[i]
        if (c>=cohort_start&c<=cohort_end) {
          
          ## targeted cohort and months whose conterfactuals are predicted
          target <- data[data$cohort==c&data$cumu_month>=cumu_start&data$cumu_month<=cumu_end,]
          target$fe <- seq(1,cumu_end-cumu_start+1)
          construct <- data.frame()
          for (t in smooth_cohort_lag1:smooth_cohort_lag2){
            construct <- rbind(
              construct,
              data[data$cohort==unique(target$cohort)-t&
                     data$cyear %in% unique(target$cyear)-t&
                     data$cumu_month %in% (target$cumu_month-12*t),]
            )
          }
          
          ## fit the model
          construct <-
            construct %>%
            arrange(cohort,cumu_month) %>%
            group_by(cohort) %>%
            mutate(fe = 1:n())
          model <- lm(mortality~cumu_month+factor(fe),construct)
          
          ## save prediction
          target$counterfactual <- predict.lm(model, target)
          
          ## save residuals
          target$residuals <- list(as.vector(model$residuals))
          predict <- rbind(predict,target)
        } 
        
        ## monitor progress
        progress(i,length(unique(data$cohort)))
      }
      
      ## add ATE for each month
      cohort_ate <-
        predict %>%
        mutate(ate=mortality-counterfactual) %>%
        group_by(cohort) %>%
        summarize(delta_c=sum(ate),
                  observed=sum(mortality),
                  delta_c_list=list(ate)) %>%
        mutate(percent_diff=delta_c/observed*100) %>%
        arrange(desc(cohort))
      
      placebo_results <- rbind(placebo_results, cohort_ate)
      
      ## monitor progress
      ##    progress(p,placebo_total_months)
    }
    
    ## placebo results lower and upper bound
    placebo_results_main <- placebo_results %>%
      group_by(cohort) %>%
      summarize(min_percent_diff=min(percent_diff),
                max_percent_diff=max(percent_diff),
                median_percent_diff=median(percent_diff),
                min_delta_c=min(delta_c),
                max_delta_c=max(delta_c),
                median_delta_c=median(delta_c),
                list_percent_diff=list(percent_diff),
                list_delta_c=list(delta_c))
    
    return(list(placebo_results_main, placebo_results))
  }

placebo_results <-
  placebo_diagnostics(data=df,
                      placebo_total_months=plcb_n,
                      placebo_year_start=plcb_ybeg,
                      placebo_month_start=plcb_mbeg,
                      treatment_delta=t_delta,
                      cohort_start=coh_beg,
                      cohort_end=coh_end,
                      smooth_cohort_lag1=coh_lag1,
                      smooth_cohort_lag2=coh_lag2)


## save the placebo estimates - this includes min, max, and median for the placebo estimates for each cohort
placebo_estimates <- placebo_results[[1]] %>% dplyr::select(-c("list_percent_diff","list_delta_c"))
#write.csv(placebo_estimates,"placebo_est.csv",row.names = FALSE)

## kludge to write estimates using a format

tmp  <- matrix(0,nrow=n_coh,ncol=7)
tmp2 <- matrix(0,nrow=n_coh,ncol=4)

for (i in 1:7)
{
  tmp[,i] <- unlist((placebo_estimates[,i]),use.names=FALSE)
}

tmp2[,1] <- sprintf("%i",tmp[,1])

for (i in 2:4)
{
  tmp2[,i] <- sprintf("%10.3e",tmp[,i+3])
}

write.table(tmp2,file="placebo.sum",col.names=FALSE,row.names=FALSE,quote=FALSE)

## now check the multiple placebo values for each cohort
detailed_placebo_values <- placebo_results[[1]] %>% dplyr::select(c("cohort","list_percent_diff","list_delta_c"))

##### bias-corrected estimates with CI ######

## inputs of the function includes the second output of `placebo_diagnostics()`

bias_corrected_estimate <-
  function(placebo_results,
           placebo_total_months=plcb_n,
           treatment_delta=t_delta,
           ci=ci_alpha
  ) {
    
    ## pool the ate of each placebo simulation for each cohort
    ate_placebo <-
      placebo_results %>%
      group_by(cohort) %>%
      summarize(delta_c_placebo=list(unlist(delta_c_list)),
                percent_diff_placebo=list(unlist(delta_c_list)/observed))
    
    ## extract the median placebo ate of each month
    
    ## raw
    ate_placebo$median_delta_c_placebo <- NA
    for (r in 1:nrow(ate_placebo)){
      median_placebo <- c()
      for (i in 1:t_delta){
        median_placebo <- c(median_placebo,
                            median(ate_placebo[r,]$delta_c_placebo[[1]][seq(i,t_delta*placebo_total_months,t_delta)]) )
      } 
      ate_placebo[r,]$median_delta_c_placebo <- list(median_placebo)
    }
    
    ## percent
    ate_placebo$median_percent_placebo <- NA
    for (r in 1:nrow(ate_placebo)){
      median_placebo <- c()
      for (i in 1:t_delta){
        median_placebo <- c(median_placebo,
                            median(ate_placebo[r,]$percent_diff_placebo[[1]][seq(i,t_delta*placebo_total_months,t_delta)]) )
      } 
      ate_placebo[r,]$median_percent_placebo <- list(median_placebo)
    }
    
    ate_placebo <- ate_placebo %>%
      dplyr::select(cohort,median_delta_c_placebo,median_percent_placebo)
    
    ## merge with the original residuals from regression
    ate <- merge(cohort_ate %>%
                   group_by(cohort) %>%
                   mutate(percent_list=list(unlist(delta_c_list)/observed)) %>%
                   dplyr::select(cohort,delta_c_list,percent_list,
                                 residuals_regression,observed),
                 ate_placebo,
                 by = "cohort",
                 all.y = T)
    
    ## create the difference between all possible combinations of residuals (median bias vs original estimate)
    ## the number of possible combinations = t_delta^2
    ate$delta_c_adjust <- NA
    for (r in 1:nrow(ate)){
      ate[r,]$delta_c_adjust <- 
        list(as.vector(outer(ate[r,"delta_c_list"][[1]], ate[r,"median_delta_c_placebo"][[1]], "-")))
    }
    ate$percent_diff_adjust <- NA
    for (r in 1:nrow(ate)){
      ate[r,]$percent_diff_adjust <- 
        list(as.vector(outer(ate[r,"percent_list"][[1]], ate[r,"median_percent_placebo"][[1]], "-")))
    }
    
    ## add 99% CI for each cohort
    ate$lower <- NA
    ate$higher <- NA
    ate$residuals_regression_percent <- list(unlist(ate$residuals_regression)/ate$observed)
    for (i in 1:nrow(ate)){
      
      ## 99% CI from permutation
      CI <- cint(dset(unlist(ate[i,]$percent_diff_adjust),unlist(ate[i,]$residuals_regression_percent)), 
                 conf.level = 1-(1-ci)/(max(ate$cohort)-min(ate$cohort)+1), tail = "Two")$conf.int*(t_delta)
      
      ## add to the original data
      ate[i,]$lower <- CI[1]
      ate[i,]$higher <- CI[2]
    }
    
    ate <- ate %>%
      mutate(percent_lower = lower*100,
             percent_higher = higher*100) %>%
      select(-c("lower","higher"))
    
    ## add the original bias-corrected estimate
    ate$percent_diff <- NA
    for (r in 1:nrow(ate)){
      ate[r,]$percent_diff <- 
        mean(ate[r,]$percent_diff_adjust[[1]])*18
    }
    ate$percent_diff <- ate$percent_diff*100
    ## return the output
    return(ate)
  }

## estimate placebo adjusted ATEs for each cohort
bias_corrected_ate <- 
  bias_corrected_estimate(placebo_results[[2]],
                          placebo_total_months=plcb_n,
                          treatment_delta=t_delta,
                          ci=ci_alpha)

## add raw changes
bias_corrected_ate <- 
  bias_corrected_ate %>%
  mutate(delta_c_adjust = percent_diff*observed/100,
         delta_c_lower = percent_lower*observed/100,
         delta_c_higher = percent_higher*observed/100)

placebo_adj <- bias_corrected_ate %>% dplyr::select(
  cohort, percent_diff, percent_lower, percent_higher,
  delta_c_adjust, delta_c_lower, delta_c_higher)

require(foreign)
adj_ate5 <- subset(placebo_adj, select = c(cohort, delta_c_adjust, delta_c_lower, delta_c_higher, percent_diff, percent_lower, percent_higher ))
write.dta(adj_ate5, "/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/us_female_5_adj.dta")


