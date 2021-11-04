
res_010 <- readRDS(here::here("output", "res_010_100runs.rds"))
res_015 <- readRDS(here::here("output", "res_015_100runs.rds"))
res_020 <- readRDS(here::here("output", "res_020_100runs.rds"))
res_025 <- readRDS(here::here("output", "res_025_100runs.rds"))
res_030 <- readRDS(here::here("output", "res_030_100runs.rds"))
res_035 <- readRDS(here::here("output", "res_035_100runs.rds"))
res_040 <- readRDS(here::here("output", "res_040_100runs.rds"))
res_045 <- readRDS(here::here("output", "res_045_100runs.rds"))
res_050 <- readRDS(here::here("output", "res_050_100runs.rds"))
res_055 <- readRDS(here::here("output", "res_055_100runs.rds"))
res_060 <- readRDS(here::here("output", "res_060_100runs.rds"))

source(here::here("R", "functions", "function_corEvalFish.R"))

df_sensitivity <- data.frame(rbind(do.call(rbind, lapply(lapply(res_010, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)), 
                                   do.call(rbind, lapply(lapply(res_015, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)), 
                                   do.call(rbind, lapply(lapply(res_020, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)), 
                                   do.call(rbind, lapply(lapply(res_025, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)),
                                   do.call(rbind, lapply(lapply(res_030, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)),
                                   do.call(rbind, lapply(lapply(res_035, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)),
                                   do.call(rbind, lapply(lapply(res_040, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)),
                                   do.call(rbind, lapply(lapply(res_045, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)),
                                   do.call(rbind, lapply(lapply(res_050, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)),
                                   do.call(rbind, lapply(lapply(res_055, function(x) correlation_eval(list_tree = x)), function(x) x$correlation)),
                                   do.call(rbind, lapply(lapply(res_060, function(x) correlation_eval(list_tree = x)), function(x) x$correlation))
                                   ),
                             percent_insert = rep(c("ins10", "ins15", "ins20", "ins25", "ins30", "ins35", "ins40", "ins45", "ins50",
                                                    "ins55", "ins60"), each = 100))



# saving results ----------------------------------------------------------

saveRDS(df_sensitivity, file = here::here("output", "df_sensitivity_cor.rds"))
