library(tidyverse)
library(lubridate)
library(forecast)

# Get data
df <- read_csv("https://raw.githubusercontent.com/covid19-forecast-hub-europe/covid19-forecast-hub-europe/main/data-truth/JHU/truth_JHU-Incident%20Cases.csv")%>% 
  filter(location == "DE") %>% 
  arrange(date) %>% 
  mutate(cum_cases = cumsum(value)
         , weekly_cases = cum_cases- dplyr::lag(cum_cases,7)
         , weekly_cases = if_else(is.na(weekly_cases), 0,weekly_cases))

# Prepare forecast
forecast_date <- max(df$date) + days(1)
target_end_dates <- seq(forecast_date , forecast_date + weeks(3), by = "1 weeks") + days(5)
fc_level <- 95
h = as.numeric(difftime(max(target_end_dates),forecast_date, "days"))

  
# Estimate and predict
mdl <-  df %>% 
    pull(weekly_cases) %>% 
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

# Compute quantile forecasts  
quantiles <- c(0.01,0.025,seq(0.05,0.95,by =0.05),0.975,0.99)
qf <- matrix(0, nrow=h, ncol=length(quantiles))
m <- mdl$mean
s <- (as.vector(mdl$upper)-as.vector(mdl$lower))/1.96/2
for(h in 1:h){
    qf[h,] <- qnorm(quantiles, m[h], s[h])
  }
qf[qf< 0] <- 0

  
fc = tibble(
    target_end_date = forecast_date + 1:h,
    yhat = as.vector(mdl$mean),
    location = "DE", 
    as_tibble(qf)
  ) 

colnames(fc)[-c(1:3)] <- quantiles


fc <- fc %>% 
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

dir.create(file.path("data-processed", "CovidMetrics-epiBATS"), recursive = TRUE)
write_csv(fc,here::here("data-processed","CovidMetrics-epiBATS" , paste0(forecast_date,"-CovidMetrics-epiBATS.csv")))
