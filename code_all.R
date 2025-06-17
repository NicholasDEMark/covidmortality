rm(list = ls())
gc()
memory.limit(size=10^10)

library(svMisc)
library(dplyr)
library(CIPerm)

##### Prepare Data #####
# read data
df <- readRDS("~/NIH/Data/table_GOOD.rds")

df$mo <- as.numeric(df$mo)

df$c_mo <- (df$y - 2013)*12 + df$mo


##### Functions #########
# basic estimates with CI ##### 
# the data need to have six variables
# 1. cohort, 2. cumulated month, 
# 3.calendar year, 4. calendar month, 5. age. 6. mortality rate by month
# variable names should be following (cohort, cumu_month, cyear, cmonth, age, mortality)
# the function returns 1. cohort-specific mortality rate difference (in absolute # and %)
# 2. the CI of the difference 
# 3. residuals from regression (a list object with length = # months * smooth cohort lag years)
# 4. mortality rate difference for each month for each cohort (a list object with length = # months)

# initialize
set.seed(978653421)
# years/months are conceptions resulting in births during the GR
# assuming 9 months between a conception and birth
cv_yr_beg <- 2020
cv_mn_beg <- 1
cv_yr_end <- 2022
cv_mn_end <- 12
t_delta   <- 12*cv_yr_end + cv_mn_end - (12*cv_yr_beg + cv_mn_beg) + 1
coh_beg   <- 1924
coh_end   <- 2010
n_coh     <- coh_end-coh_beg + 1
coh_lag1  <- 1

# set placebo quantities
plcb_n    <- 12
plcb_ybeg <- 2018
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
                                  direction = c("horizontal", "oblique"),
                                  method = c("linear", "poisson"),
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
  predict_past <- data.frame()
  for (i in 1:length(unique(data$cohort))) {
    # Locate target cohort
    c <- unique(data$cohort)[i]
    if (c >= cohort_start & c <= cohort_end) {
      # Targeted cohort and months whose counterfactuals are predicted
      target <- data[data$cohort==c & data$cumu_month>=cumu_start & data$cumu_month<=cumu_end,]
      target$fe <- target$cmonth #seq(1, cumu_end-cumu_start+1)
      
      # Construct reference data (simplified for unweighted case)
      construct <- data.frame()
      for (t in smooth_cohort_lag1:smooth_cohort_lag2) {
        if (direction == "horizontal"){
          construct <- rbind(
            construct,
            data[data$cohort==unique(target$cohort)-t &
                   data$cyear %in% unique(target$cyear)-t &
                   data$cumu_month %in% (target$cumu_month[1:12]-12*t),]
          )
        }else{
          construct <- rbind(
            construct,
            data[data$cohort==unique(target$cohort) &
                   data$cyear == unique(target$cyear)[1]-t &
                   data$cumu_month %in% (target$cumu_month[1:12]-12*t),]
          )
        }
        
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
          group_by(age) %>%
          mutate(fe = cmonth)
      }
      
      # Regression with optional weighting
      if (weighting == "none") {
        if (method == "poisson"){
          model <- glm(deaths ~ offset(log(ep)) + cumu_month + factor(fe), data = construct, family=quasipoisson())
        }else{
          model <- lm(mortality ~ cumu_month + factor(fe), data = construct)
        }
      } else {
        if (method == "poisson"){
          model <- glm(deaths ~ offset(log(ep)) + cumu_month + factor(fe), data = construct, 
                       weights = weight, family=quasipoisson())
        }else{
          model <- lm(mortality ~ cumu_month + factor(fe), weights = weight, data = construct)
        }
      }
      
      if (method == "poisson"){
        # Save prediction
        target$counterfactual <- predict.glm(model, target, type = 'response')
        
        construct$counterfactual <- predict.glm(model, construct, type = 'response')
        
        # Save residuals
        target$residuals <- list(as.vector(model$residuals))
        predict <- rbind(predict, target)
        predict_past <- rbind(predict_past, construct)
      }else{
        # Save prediction
        target$counterfactual <- predict.lm(model, target, type = 'response')
        #print (c)
        
        # Save residuals
        target$residuals <- list(as.vector(model$residuals))
        predict <- rbind(predict, target)
      }
      
      
      
    }
    
    # Monitor progress
    #progress(i, length(unique(data$cohort)))
  }
  
  if (method == "poisson"){
    predict_past$ate <- predict_past$mortality - predict_past$counterfactual / predict_past$ep
  }
  
  
  # Rest of the function remains the same as in the original code
  # (ATE calculation, confidence intervals, etc.)
  
  if (method == "poisson"){
    if (direction == "horizontal"){
      cohort_ate <- predict %>%
        mutate(ate=mortality-counterfactual/ep) %>%
        group_by(cohort,age) %>%
        summarize(delta_c=sum(ate),
                  observed=sum(mortality),
                  #residuals_regression=list(unlist(unique(residuals))),
                  delta_c_list=list(ate)) %>%
        mutate(percent_diff=delta_c/observed*100) %>%
        arrange(desc(cohort))
    }else{
      cohort_ate <- predict %>%
        mutate(ate=mortality-counterfactual/ep) %>%
        group_by(cohort)%>%#,age) %>%
        summarize(delta_c=sum(ate),
                  observed=sum(mortality),
                  #residuals_regression=list(unlist(unique(residuals))),
                  delta_c_list=list(ate)) %>%
        mutate(percent_diff=delta_c/observed*100) %>%
        arrange(desc(cohort))
    }
  }else{
    if (direction == "horizontal"){
      cohort_ate <- predict %>%
        mutate(ate=mortality-counterfactual) %>%
        group_by(cohort, age) %>%
        summarize(delta_c=sum(ate),
                  observed=sum(mortality),
                  residuals_regression=list(unlist(unique(residuals))),
                  delta_c_list=list(ate)) %>%
        mutate(percent_diff=delta_c/observed*100) %>%
        arrange(desc(cohort))
    }else{
      cohort_ate <- predict %>%
        mutate(ate=mortality-counterfactual) %>%
        group_by(cohort)%>%#,age) %>%
        summarize(delta_c=sum(ate),
                  observed=sum(mortality),
                  residuals_regression=list(unlist(unique(residuals))),
                  delta_c_list=list(ate)) %>%
        mutate(percent_diff=delta_c/observed*100) %>%
        arrange(desc(cohort))
    }
  }

  
  if (method == "poisson"){
    cohort_ate$residuals_regression <- list(0)
    for (coh in unique(cohort_ate$cohort)){
      
      if (direction == "horizontal"){
        cohort_ate$residuals_regression [cohort_ate$cohort == coh] <-
          list(predict_past$ate [predict_past$age == cohort_ate$age [cohort_ate$cohort == coh]])
      }else{
        cohort_ate$residuals_regression [cohort_ate$cohort == coh] <-
          list(predict_past$ate [predict_past$cohort == coh])
      }
      
    }
  }
  
  
  
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

########################
#### COMPUTATION
#######################

coh_lag2_values <- c(5, 6)
sexes <- c("1", "2")
dirs <- c("horizontal", "oblique")
methods <- c("linear", "poisson")

for (method in methods){
  
  for (dir in dirs){
    
    for (sex in sexes){
      
      df2 <- subset(df, SEXE == sex)
      
      df2$ardy <- df2$y - df2$cohort
      
      df2 <- df2 [, c("cohort", "c_mo", "y", "mo", "ardy", "ep", "deaths", "mr")]
      
      df2 <- subset(df2, cohort > 1917 & cohort < 2020)
      
      # rename
      colnames(df2) <-  c("cohort","cumu_month","cyear","cmonth","age","ep", "deaths", "mortality")
      
      # cumulative months
      df2 <- df2 %>%
        arrange(cohort,cyear,cmonth) %>%
        group_by(cyear,cmonth) %>%
        mutate(cumu_month=cur_group_id())
      
      for (lag_value in coh_lag2_values) {
        # Set coh_lag2 for this iteration
        coh_lag2 <- lag_value
        
        # Print current value being processed
        cat("Processing coh_lag2 = ", coh_lag2, "\n", "Sex: ", sex,
            "\n", "Dir: ", dir, "\n", "Meth: ", method)
        
        # Use the function to derive results of cohort loss
        cohort_ate <- cohort_discontinuity(data=df2,
                                           shock_year_start=cv_yr_beg,
                                           shock_month_start=cv_mn_beg,
                                           shock_year_end=cv_yr_end,
                                           shock_month_end=cv_mn_end,
                                           cohort_start=coh_beg,
                                           cohort_end=coh_end,
                                           smooth_cohort_lag1=coh_lag1, ## c-2
                                           smooth_cohort_lag2=coh_lag2, ## c-6
                                           ci=ci_alpha,
                                           direction = dir,
                                           method = method,
                                           weighting = "none")
        
        # Write estimate and percent change with CIs
        loss <- cohort_ate %>% select(-c("residuals_regression","delta_c_list"))
        
        # Create filename with the lag value
        filename <- paste0("~/NIH/Results/FIN/dir_", dir,"_meth_", method,
                          "_sex_", sex,  "lag_", lag_value, ".csv")
        #write.csv(loss, filename, row.names = FALSE)
        
        # Save the cohort_ate to a variable with the corresponding name
        assign(paste0("cohort_ate", lag_value), cohort_ate)
        
        cat("Completed coh_lag2 =", coh_lag2, "\n")
      }
    }
  }
  
  
}


cat("All iterations completed.\n")
