res_010 <- readRDS(here::here("output", "res_010_50runs.rds"))
res_015 <- readRDS(here::here("output", "res_015_50runs.rds"))
res_020 <- readRDS(here::here("output", "res_020_50runs.rds"))
res_025 <- readRDS(here::here("output", "res_025_50runs.rds"))
res_030 <- readRDS(here::here("output", "res_030_50runs.rds"))
res_035 <- readRDS(here::here("output", "res_035_50runs.rds"))
source(here::here("R", "functions", "function_corEvalFish.R"))

res_correlation_010 <- lapply(res_010, function(x) correlation_eval(list_tree = x))
res_correlation_015 <- lapply(res_015, function(x) correlation_eval(list_tree = x))
res_correlation_020 <- lapply(res_020, function(x) correlation_eval(list_tree = x))
res_correlation_025 <- lapply(res_025, function(x) correlation_eval(list_tree = x))
res_correlation_035 <- lapply(res_035, function(x) correlation_eval(list_tree = x))
df_sensitivity <- data.frame(rbind(do.call(rbind, lapply(res_correlation_010, function(x) x$correlation)), 
                                  do.call(rbind, lapply(res_correlation_015, function(x) x$correlation)), 
                                  do.call(rbind, lapply(res_correlation_020, function(x) x$correlation)), 
                                  do.call(rbind, lapply(res_correlation_025, function(x) x$correlation))),
                            percent_insert = rep(c("ins10", "ins15", "ins20", "ins25"), each = 50))


# saving results ----------------------------------------------------------

saveRDS(df_sensitivity, file = here::here("output", "df_sensitivity_cor.rds"))
