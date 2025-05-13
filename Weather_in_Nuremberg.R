#### Packages ####
library(biLocPol)
library(tidyverse)
library(lubridate)
library(hms)
library(locpol)
library(interp)
library(reshape2)
library(ffscb)
library(plotly)


deriv_est_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

#Z1 = read.table("data/weather_in_nuremberg_raw/produkt_zehn_min_tu_19950804_19991231_03668.txt", 
#                header = T, sep = ";")
Z2 = read.table("data/weather_in_nuremberg_raw/produkt_zehn_min_tu_20000101_20091231_03668.txt", 
                header = T, sep = ";")
Z3 = read.table("data/weather_in_nuremberg_raw/produkt_zehn_min_tu_20100101_20191231_03668.txt", 
                header = T, sep = ";")
Z4 = read.table("data/weather_in_nuremberg_raw/produkt_zehn_min_tu_20200101_20221231_03668.txt", 
                header = T, sep = ";")

#head(Z1)
#Z1 = Z1[-(1:48), ]
#Z1[Z1 == -999] = NA
Z2[Z2 == -999] = NA
Z3[Z3 == -999] = NA
Z4[Z4 == -999] = NA

#Z1$MESS_DATUM = ymd_hm(Z1$MESS_DATUM)
Z2$MESS_DATUM = ymd_hm(Z2$MESS_DATUM)
Z3$MESS_DATUM = ymd_hm(Z3$MESS_DATUM)
Z4$MESS_DATUM = ymd_hm(Z4$MESS_DATUM)
Z1$eor = NA

# Using only Z2, Z3 and Z4 to avoid to many NA
N =  rbind(Z2,Z3, Z4)
head(N)

N$JAHR = year(N$MESS_DATUM)
N$MONAT = month(N$MESS_DATUM)
N$TAG = day(N$MESS_DATUM)
N$UHRZEIT = as_hms(N$MESS_DATUM)

N =  N[!is.na(N$TT_10), ]
N

#load("R-Codes/data\\weather_data_nuremberg.RData")

N |> 
  filter(!is.na(TT_10)) |>
  summarise(.by = c(JAHR, MONAT, TAG), 
            n = n()) |> 
  filter(n == 144) |>  
  summarise(.by = MONAT, 
            m = n()) 

N$JAHRTAG = ymd(paste(year(N$MESS_DATUM), 1, day(N$MESS_DATUM)))
N_casted = acast(N, UHRZEIT  ~ JAHRTAG ~ MONAT, value.var = "TT_10" )
UHRZEIT = hms(hour = rep(seq(0,23,1), each = 6),
              minutes = rep(seq(0, 50, 10), 24), 
              seconds = rep(0, 144))
N_casted |>  dim() # Uhrzeit, Tag, Monat
N_casted[1:5, 1:100, 1]
shift = 48
N_add = array(NA, dim = c(144 + 2*shift, 682, 12))
N_add[shift + (1:144), , ] = N_casted
N_add[1:shift, -1, ] = N_casted[(145-shift):144, -682, ]
N_add[(145+shift):(shift*2 + 144), -682, ] = N_casted[1:shift, -1, ]
colnames(N_add) = colnames(N_casted)
head(N_add)
#pattern = "-(02|05|08|11|14|17|20|23|26)$"
pattern = "-(01|29|30|31)$"
N_final = N_add[,!grepl(pattern, colnames(N_add)),]
N_final

##### estimation procedure #####
N_Mittelwerte = N_final |> 
  apply(c(1,3), mean, na.rm = T)

N_Mittelwerte |> dim()

#### estimation of mean derivatives ####
x_temp = seq(0,1,length.out = 145)[-145]
x = c(-(1-x_temp[(145-shift):144]), x_temp, 1+x_temp[1:shift])

farben = c( "#a6cee3", "#03396c",  # Blau-Töne
            "#33a02c", "#66c21f", "#006400",  # Grün-Töne
            "#e31a1c", "#fb9a99", "#990000",  # Rot-Töne
            "#ff7f00", "#ffb300", "#b15928",  # Gelb-Orange-Töne
            "#1f78b4")

weights_mean_derivative = locPolWeights(x, x_temp, 3, 0.2, EpaK)$allWeig[,2,]
p = length(x_temp)
monthly_weather_deriv = matrix(0, p, 12)

for(m in 1:12){
  Y = N_Mittelwerte[, m]
  monthly_weather_deriv[, m] = weights_mean_derivative %*% Y
}

dim(monthly_weather_deriv)
colnames(monthly_weather_deriv) = 1:12
deriv_tibble = monthly_weather_deriv |> 
  as_tibble() |> 
  cbind(x_temp) |> 
  pivot_longer(1:12, values_to = "deriv", names_to = "month") |> 
  mutate(month = as.numeric(month) ) 

cbind(deriv_tibble, rep(UHRZEIT, each = 1728/144)) |> 
  mutate(month = as.factor(month)) |> 
  rename(UHRZEIT = "rep(UHRZEIT, each = 1728/144)") |> 
  ggplot(aes(x = UHRZEIT, y = deriv, col = month, linetype = month)) +
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

#### estimation of the covariance derivatives ####
p.eval = 72
p = 144 + 2*shift
x.design = (0:(144-1))/144
x.design.extended = c( x.design[ (144-shift+1):144 ] - 1, x.design, x.design[1:shift] + 1)

W = local_polynomial_weights(p, 0.3, p.eval, T, m = 2, del = 1, 
                             eval.type = "diagonal", x.design.grid = x.design.extended) # watch out for evaluation

g_hat = array(0, dim = c(p.eval, p.eval, 12))
for ( m in 1:12) {
  # has_na = apply(N_final[,,m], 2, function(x){sum(is.na(x)) > 0})
  Y = t(N_final[,,m]) |>  
    observation_transformation(na.rm = T)
  g_hat[,,m] = eval_weights(W, Y)[,2]
}


time = seq(0, 1, length.out = 72)

plot_ly() %>% 
  add_surface(x = ~ time, y = ~ time, z = g_hat[,,2]) 

deriv_cov_weather_diag = tibble(x = rep(time, 2), 
                                est = c(diag(g_hat[,,2]), diag(g_hat[,,3])), 
                                deriv = gl(2, 72, labels = c("G10", "G01")))
deriv_cov_weather_diag |> 
  ggplot(aes( x = x, y = est, color = deriv, linetype = deriv)) + 
  geom_line(linewidth = .9) +
  deriv_est_theme + 
  labs(y = NULL, x = "time", title = expression(Partial~derivatives~of~Gamma~on~the~diagonal) ) + 
  scale_discrete_manual(
    aesthetics = c("color", "linetype"),
    values = c(2,4), 
    name = "Deriv.",
    labels = c(expression(italic(d)^{"(1,0)"}*Gamma), expression(italic(d)^{"(0,1)"}*Gamma))
  ) 

#### estimation of mean itself #####
weights_mean = locPolWeights(x, x_temp, 3, 0.2, EpaK)$allWeig[,1,]

monthly_weather = matrix(0, p, 12)

for(m in 1:12){
  Y = N_Mittelwerte[,m]
  monthly_weather[, m] = weights_mean %*% Y
}

colnames(monthly_weather) = 1:12
mean_tibble = monthly_weather |> 
  as_tibble() |> 
  cbind(x_temp) |> 
  pivot_longer(1:12, values_to = "temp", names_to = "month") |> 
  mutate(month = as.numeric(month) ) 



ggplot() +
  geom_line(data = cbind(mean_tibble, rep(UHRZEIT, each = 1728/144)) |> 
              mutate(month = as.factor(month)) |> 
              rename(UHRZEIT = "rep(UHRZEIT, each = 1728/144)") , 
            aes(x = UHRZEIT, y = temp, col = month, lty = month)) +
  labs(x = "time", y = NULL, title = "Estimated daily mean temperatur") + 
  scale_colour_manual(values = farben) + 
  scale_linetype_manual(values = c(2,3,4,4,4,1:3,5,5,5,1)) + 
  geom_text(data = data.frame(x = c(hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0),
                                    hms(0, -15, 0), hms(0, -15, 0), hms(0, -15, 0), hms(0, -22, 0), hms(0, -22, 0), hms(0, -22, 0)), 
                              y = c(-.1, 0.3, 2.5, 6.6, 10.2, 14.25, 15.5, 16, 11.8, 7.6, 4.2, 1.5), 
                              month = gl(12,1)), aes(label = month, x = x, y = y), hjust = 0.4, size = 3) + 
  deriv_est_theme
ggsave("Grafics/mean_temperature.png", device = "png",
       width = 5, height = 4, units = "in")

