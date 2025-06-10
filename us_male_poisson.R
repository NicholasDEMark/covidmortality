library(svMisc)
library(dplyr)
library(CIPerm)

##### Prepare Data #####
# read data
df <- read.csv("/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_2006-2023.csv", 
               header=TRUE)

# rename
colnames(df) <-  c("cohort","cumu_month","cyear","age","ep", "deaths","cmonth","mortality")
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
# years/months are conceptions resulting in births during the GR
# assuming 9 months between a conception and birth
cv_yr_beg <- 2020
cv_mn_beg <- 1
cv_yr_end <- 2020
cv_mn_end <- 12
t_delta   <- 12*cv_yr_end + cv_mn_end - (12*cv_yr_beg + cv_mn_beg) + 1
coh_beg   <- 1930
coh_end   <- 2019
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
        model <- glm(deaths ~ offset(log(ep)) + cumu_month + factor(fe), data = construct, family=quasipoisson())
      } else {
        model <- glm(mortality ~ offset(log(ep)) + cumu_month + factor(fe), 
                     data = construct, 
                     weights = weight, family=quasipoisson())
      }
      
      # Save prediction
      target$counterfactual <- predict.glm(model, target, type = 'response')
      
      construct$counterfactual <- predict.glm(model, construct, type = 'response')
      
      # Save residuals
      target$residuals <- list(as.vector(model$residuals))
      predict <- rbind(predict, target)
      predict_past <- rbind(predict_past, construct)
      
    }
    
    # Monitor progress
    #progress(i, length(unique(data$cohort)))
  }
  
  predict_past$ate <- predict_past$mortality - predict_past$counterfactual / predict_past$ep
  
  # Rest of the function remains the same as in the original code
  # (ATE calculation, confidence intervals, etc.)
  cohort_ate <- predict %>%
    mutate(ate=mortality-counterfactual/ep) %>%
    group_by(cohort,age) %>%
    summarize(delta_c=sum(ate),
              observed=sum(mortality),
              #residuals_regression=list(unlist(unique(residuals))),
              delta_c_list=list(ate)) %>%
    mutate(percent_diff=delta_c/observed*100) %>%
    arrange(desc(cohort))
  
  cohort_ate$residuals_regression <- list(0)
  
  for (coh in unique(cohort_ate$cohort)){
    cohort_ate$residuals_regression [cohort_ate$cohort == coh] <- 
      list(predict_past$ate [predict_past$age == cohort_ate$age [cohort_ate$cohort == coh]])
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


coh_lag2_values <- c(5, 6)#, 7)

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
  filename <- paste0("/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/male_US_poisson", lag_value, ".csv")
  write.csv(diff, filename, row.names = FALSE)
  
  # # Kludge to write estimates using a format
  # tmp <- matrix(0, nrow=n_coh, ncol=8)
  # for (i in 1:8) {
  #   tmp[,i] <- unlist(cohort_ate[,i], use.names=FALSE)
  # }
  # 
  # tmp2 <- tmp
  # tmp2[,1] <- sprintf("%i", tmp[,1])
  # tmp2[,8] <- sprintf("%10.3E", tmp[,8])
  # for (i in 2:4) {
  #   tmp2[,i] <- sprintf("%10.3E", tmp[,i])
  # }
  # for (i in 5:7) {
  #   tmp2[,i] <- sprintf("%6.2f", tmp[,i])
  # }
  # 
  # # Create formatted text filename with the lag value
  # txt_filename <- paste0("/Users/nicholasmark/Dropbox/gr1/runs2/births_estimates", lag_value, ".txt")
  # write.table(tmp2, file=txt_filename, col.names=FALSE, row.names=FALSE, quote=FALSE)
  # 
  # #Create formatted csv filename with the lag value
  # csv_filename <- paste0("/Users/nicholasmark/Dropbox/gr1/runs2/births_estimates", lag_value, ".csv")
  # write.table(tmp2, file=csv_filename, col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # Save the cohort_ate to a variable with the corresponding name
  assign(paste0("cohort_ate", lag_value), cohort_ate)

  # Create filename with the lag value
  filename <- paste0("/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/male_US", lag_value, ".csv")
  write.csv(diff, filename, row.names = FALSE)
  
  # Save the cohort_ate to a variable with the corresponding name
  assign(paste0("cohort_ate", lag_value), cohort_ate)
  
  cat("Completed coh_lag2 =", coh_lag2, "\n")
}

cat("All iterations completed.\n")


require(foreign)
ate5 <- subset(cohort_ate5, select = c(cohort, delta_c_adjust, lower_c, higher_c, percent_diff, percent_lower, percent_higher, observed))
ate6 <- subset(cohort_ate6, select = c(cohort, delta_c, lower_c, higher_c, percent_diff, percent_lower, percent_higher, observed))
write.dta(ate5, "/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/us_male_5.dta")
write.dta(ate6, "/Users/nicholasmark/Dropbox/covid/2025/runs/US/output/us_male_6.dta")


