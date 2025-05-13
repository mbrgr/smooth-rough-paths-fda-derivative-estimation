#### Packages and Theme ####
library(biLocPol)
library(tidyverse)
library(lubridate)
library(hms)
library(interp)
library(reshape2)
library(plotly)
library(locpol)

deriv_est_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#### results ####
load("data/weather_covariance_derivative_estimation.RData") # for covariance estimation, Figure 13, Figure 16

#### Load and process the data ####
load("data/weather_in_nuremberg/weather_data_nuremberg.RData")
head(N, n = 7)
dim(N)

N0 = N[7:length(N$MESS_DATUM), ] 
which(N$MESS_DATUM == "2012-01-01 00:00:00")
which(N$MESS_DATUM == "2017-01-01 00:00:00")

N1 = N0 |> 
  mutate(time = as.character(time)) |> 
  dplyr::select(JAHR, MONAT, TAG, time, TT_10) |> 
  pivot_wider(names_from = time,
              values_from = TT_10)

# extend data to improve estimation at the boundaries
N2 = cbind(N1[, 1:3], 
           rbind(NA, N1[-length(N1$JAHR), -(1:3)]), 
           N1[, -(1:3)], 
           rbind(N1[-1, -(1:3)], NA))

# remove days that have wrong values by the previous procedure
N3 = N2[!(((N2[,3] %in% 28:29) & (N2[,2] == 2)) |
            (N2[,3] == 1) |
            ((N2[,3] == 31) & (N2[,2] %in% c(1,3,5,7,8,10,12))) |
            ((N2[,3] == 30) & (N2[,2] %in% c(4,6,9,11)))), ]
sum(is.na(N3))

# remove NA entries
to_rm = apply(N3, 1, function(v){any(is.na(v))}) %>% which()
N3 = N3[-to_rm, ]
sum(is.na(N3))

#### Mean Derivative Estimation ####

N_bar = matrix(0, 12, 3 * p)

for(m in 1:12){
  N_bar[m,] = N3[N3[,2] == m,-(1:3)] |> 
    apply(2, mean, na.rm = T)
}


x_temp = seq(0,1,length.out = 145)[-145]
x_test = seq(-1, 2, length.out = 3*p)

farben = c( "#a6cee3", "#03396c",  # Blau-Töne
            "#33a02c", "#66c21f", "#006400",  # Grün-Töne
            "#e31a1c", "#fb9a99", "#990000",  # Rot-Töne
            "#ff7f00", "#ffb300", "#b15928",  # Gelb-Orange-Töne
            "#1f78b4")

weights_mean_derivative = locPolWeights(x_test, x_temp, 3, 0.2, EpaK)$allWeig[,2,]
p = length(x_temp)
monthly_weather_deriv = matrix(0, p, 12)

for(m in 1:12){
  Y = N_bar[m,]
  monthly_weather_deriv[, m] = weights_mean_derivative %*% Y
}

dim(monthly_weather_deriv)
colnames(monthly_weather_deriv) = 1:12
deriv_tibble = monthly_weather_deriv |> 
  as_tibble() |> 
  cbind(x_temp) |> 
  pivot_longer(1:12, values_to = "deriv", names_to = "month") |> 
  mutate(month = as.numeric(month) ) 


time = N$UHRZEIT[7:150]
time 

##### Figure 11 #####
cbind(deriv_tibble, rep(time, each = 1728/144)) |> 
  mutate(month = as.factor(month)) |> 
  rename(time = "rep(time, each = 1728/144)") |> 
  ggplot(aes(x = time, y = deriv, col = month, linetype = month)) +
  geom_line() +
  labs(x = "time", y = NULL, title = "Derivatives of mean temperatur per month") + 
  scale_colour_manual(values = farben) + 
  scale_linetype_manual(values = c(2,3,4,4,4,1:3,5,5,5,1)) +
  deriv_est_theme

ggsave("grafics/derivatives_mean_temperature.png", device = "png",
       width = 5, height = 4, units = "in")

# test if the estimation produces reliable results 
rbind(deriv_tibble, deriv_tibble |> mutate(x_temp = x_temp +1)) |> 
  mutate(month = as.factor(month)) |> 
  ggplot(aes(x = x_temp, y = deriv, col = month, linetype = month)) +
  geom_line() +
  labs(x = "time", y = NULL, title = "Derivatives of mean temperatur per month") + 
  scale_colour_manual(values = farben) + 
  scale_linetype_manual(values = c(2,3,4,4,4,1:3,5,5,5,1)) +
  deriv_est_theme

##### Mean estimation #####
weights_mean = locPolWeights(x_test, x_temp, 3, 0.2, EpaK)$allWeig[,1,]
p = length(x_temp)
monthly_weather = matrix(0, p, 12)

for(m in 1:12){
  Y = N_bar[m,]
  monthly_weather[, m] = weights_mean %*% Y
}

dim(monthly_weather)
colnames(monthly_weather) = 1:12
mean_tibble = monthly_weather |> 
  as_tibble() |> 
  cbind(x_temp) |> 
  pivot_longer(1:12, values_to = "temp", names_to = "month") |> 
  mutate(month = as.numeric(month) ) 

##### Figure 12 #####
ggplot() +
  geom_line(data = cbind(mean_tibble, rep(time, each = 1728/144)) |> 
              mutate(month = as.factor(month)) |> 
              rename(time = "rep(time, each = 1728/144)"), 
            aes(x = time, y = temp, col = month, linetype = month)) +
  labs(x = "time", y = NULL, title = "Estimated daily mean temperature") + 
  scale_colour_manual(values = farben) + 
  scale_linetype_manual(values = c(2,3,4,4,4,1:3,5,5,5,1)) +
  deriv_est_theme + 
  geom_text(data = data.frame(x = c(hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0),
                                    hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -22, 0), hms(0, -22, 0), hms(0, -22, 0)), 
                              y = c(-.1, 0.3, 2.5, 6.6, 10.2, 14.25, 16, 15.5, 11.8, 7.6, 4.2, 1.5), 
                              month = gl(12,1)), aes(label = month, x = x, y = y), hjust = 0.4)

ggsave("grafics/mean_temperature.png", device = "png",
       width = 5, height = 4, units = "in")

#### (Co)-Variance Estimation ####
p.eval = 144
p = 144

#evaluation set with extended interval to next an previous day
x_test = seq(-1, 2, length.out = 3*p)

W = local_polynomial_weights(3*p, 0.3, p, T, m = 2, del = 1, eval.type = "diagonal", x.design.grid = x_test) # watch out for evaluation

# evaluation per month
g_hat = array(0, dim = c(p.eval, 3, 12))
for ( m in 1:12) {
  Y = N3[N3[,2] == m, -(1:3)] |>  
    observation_transformation(na.rm = T)
  g_hat[,,m] = eval_weights(W, Y)
}
est_per_month_G = g_hat[,1,]   # variance
est_per_month_G10 = g_hat[,2,] # del01 G
est_per_month_G01 = g_hat[,3,] # del10 G

colnames(est_per_month_G) = factor(1:12, labels = 1:12)
colnames(est_per_month_G10) = factor(1:12, labels = 1:12)
colnames(est_per_month_G01) = factor(1:12, labels = 1:12)
est_G = est_per_month_G %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12),
         g = rep("G", each =length(x))) 
est_G10 = est_per_month_G10 %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12), 
         g = rep("G10", each =length(x))) 
est_G01 = est_per_month_G01 %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12), 
         g = rep("G01", each = length(x))) 
est_per_month = rbind(est_G01, est_G10) %>% mutate(g = factor(g, levels = c("G10", "G01")))



est_per_month %>% 
  mutate(name = factor(name, levels = 1:12)) %>% 
  ggplot(aes( x = x, y = value, color = g)) + 
  geom_line() + 
  facet_wrap(name~.)


est_G %>% 
  mutate(name = factor(name, levels = 1:12)) %>% 
  ggplot(aes( x = x, y = value, color = g)) + 
  geom_line() + 
  facet_wrap(name~.)


time = N$time[7:150]
time 

final_est = est_per_month %>% 
  mutate(time = rep(rep(time, each = 12), 2)) %>% 
  mutate(name = factor(name, levels = 1:12)) 

##### Figure 16 #####
final_est %>% 
  ggplot(aes( x = time, y = value, color = g, linetype = g)) + 
  facet_wrap(name~.)  + 
  geom_line(linewidth = .9) +
  labs(y = NULL, x = "time", title = expression(Partial~derivatives~of~Gamma~on~the~diagonal) ) + 
  scale_discrete_manual(
    aesthetics = c("color", "linetype"),
    values = c("G10" = 2,"G01" = 4), 
    name = "Deriv.",
    labels = c("G10" = expression(italic(d)^{"(1,0)"}*Gamma), "G01" = expression(italic(d)^{"(0,1)"}*Gamma))
  ) +
  scale_x_continuous(breaks = c(time[25], time[121]))

ggsave("grafics/weather_part_deriv_gamma_diagonal_all_months.png", device = "png", width = 8, height =6, units = "in")

##### Figure 13 #####
final_est %>% 
  filter(name == 4 ) %>% 
  ggplot(aes( x = time, y = value, color = g, linetype = g)) + 
  geom_line(linewidth = .9) +
  deriv_est_theme + 
  labs(y = NULL, x = "time", title = expression(Partial~derivatives~of~Gamma~on~the~diagonal) ) + 
  scale_discrete_manual(
    aesthetics = c("color", "linetype"),
    values = c("G10" = 2,"G01" = 4), 
    name = "Deriv.",
    labels = c("G10" = expression(italic(d)^{"(1,0)"}*Gamma), "G01" = expression(italic(d)^{"(0,1)"}*Gamma))
  )  + 
  scale_x_continuous(breaks = c(time[25], time[73], time[121], time[169], time[217]))

ggsave("grafics/weather_part_deriv_gamma_diagonal.png", device = "png", width = 5, height = 4, units = "in")

#save.image("data/weather_covariance_derivative_estimation.RData")
