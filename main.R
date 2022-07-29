library(tidyverse)
library(lubridate)
library(forecast)
inc_git <- read_csv("https://raw.githubusercontent.com/KITmetricslab/covid19-forecast-hub-de/master/data-truth/RKI/truth_RKI-Incident%20Cases_Germany.csv")


forecast_date <- max(inc_git$date)
print(forecast_date)
target_end_dates <- seq(forecast_date - weeks(2), forecast_date + weeks(3), by = "1 weeks") + days(5)


df <-  inc_git %>% group_by(location) %>% arrange(date) %>% 
  mutate(cum_cases = cumsum(value)
         ,acute_cases = cum_cases- dplyr::lag(cum_cases,7)
         , acute_cases = if_else(is.na(acute_cases), 0,acute_cases))

fc_level <- 95
h = as.numeric(difftime(max(target_end_dates),forecast_date, "days"))

cc = "Germany"
  
  # Estimate and predict
  mdl = df %>% 
    filter(location_name == cc) %>% 
    pull(acute_cases) %>% 
    tbats(seasonal.periods = c(7, 14), max.p = 7, max.q = 7) %>%
    forecast(level = fc_level, h = h)
  
  # If no variance in out of sample prediction (likely convergence issue),
  # estimate arima model instead
  if (sd(mdl$mean) < 0.5) {
    mdl = df %>% 
      filter(location_name == cc) %>% 
      pull(acute_cases) %>% 
      auto.arima(max.p = 14, max.q = 14) %>%
      forecast(level = fc_level, h = h)
  }
  
  quantiles <- c(0.01,0.025,seq(0.05,0.95,by =0.05),0.975,0.99)
  qf <- matrix(0, nrow=h, ncol=length(quantiles))
  m <- mdl$mean
  s <- (as.vector(mdl$upper)-as.vector(mdl$lower))/1.96/2
  for(h in 1:h){
    qf[h,] <- qnorm(quantiles, m[h], s[h])
  }
  qf[qf< 0] <- 0

  
  tempOS = tibble(
    Date = forecast_date + 1:h,
    yhat = mdl$mean,
    location = "DE"
  ) 

  
  tempCI <- tibble(Date = tempOS$Date,
                   as_tibble(qf)) 
  
  colnames(tempCI) <-  c("Date",quantiles)

fc = left_join(tempOS, tempCI)  %>% 
  rename(target_end_date = Date) %>% 
  mutate(yhat = if_else(yhat < 0, 0,yhat)
         ,yhat = round(yhat)
         , h = round(difftime(target_end_date,forecast_date, units = "weeks"))
  ) %>% 
  mutate(across(starts_with("0."), ~round(.) ) )%>% 
  filter(target_end_date %in% target_end_dates) %>% 
  pivot_longer(cols = c(-target_end_date,-location,-h,), names_to = "quantile") %>% 
  mutate(
    target = paste0(h," wk ahead inc case")
    , forecast_date = forecast_date
    , quantile = if_else(quantile == "yhat", NA_real_, as.double(quantile) ) 
    , type = if_else(is.na(quantile), "point", "quantile")
    , scenario_id = "forecast") %>% 
  filter(h > 0 | is.na(quantile)) %>% 
  select(forecast_date, target, target_end_date, location, type, quantile, value, scenario_id)

write_csv(fc,here::here("data-processed","CovidMetrics-epiBATS" , paste0(forecast_date,"-CovidMetrics-epiBATS.csv")))
