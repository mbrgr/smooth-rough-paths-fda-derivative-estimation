library(biLocPol)
library(tidyverse)
library(lubridate)
library(hms)
library(interp)
library(reshape2)
library(plotly)



deriv_est_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))


load("data/weather_data_nuremberg.RData")
head(N, n = 7)
dim(N)

#### Test mit reduziertem N ####
N0 = N[7:length(N$MESS_DATUM), ] # Jahre 2003 und 2004

which(N$MESS_DATUM == "2012-01-01 00:00:00")
which(N$MESS_DATUM == "2017-01-01 00:00:00")

#### deal with NAs ####
N0$TT_10 %>% is.na() %>% sum()

summary(N0)
N0$MONAT |> table()
tail(N0)
head(N0)

N1 = N0 |> 
  mutate(UHRZEIT = as.character(UHRZEIT)) |> 
  dplyr::select(JAHR, MONAT, TAG, UHRZEIT, TT_10) |> 
  pivot_wider(names_from = UHRZEIT,
              values_from = TT_10)
head(N0)


N2 = cbind(N1[, 1:3], 
           rbind(NA, N1[-length(N1$JAHR), -(1:3)]), 
           N1[, -(1:3)], 
           rbind(N1[-1, -(1:3)], NA))
#N2_alt = cbind(N1, N1[, -(1:3)], N1[, -(1:3)])

(N2- N2_alt)[1:5,280:300 ]
dim(N1)
dim(N2) # passt
head(N2)
sum(is.na(N2))
#view(t(N[1:3,]))
#view(t(N2[1:3, ]))
#view(t(N1[1:3, ]))

to_unit_interval = function(x){
  (x - min(x)) / (max(x) - min(x))
}


### entfernen der falschen Tage
N3 = N2[!(((N2[,3] %in% 28:29) & (N2[,2] == 2)) |
            (N2[,3] == 1) |
            ((N2[,3] == 31) & (N2[,2] %in% c(1,3,5,7,8,10,12))) |
            ((N2[,3] == 30) & (N2[,2] %in% c(4,6,9,11)))), ]
sum(is.na(N3))

to_rm = apply(N3, 1, function(v){any(is.na(v))}) %>% which()
N3 = N3[-to_rm, ]

sum(is.na(N3))

p.eval = 144
p = 144
#p_eval_ext = p.eval + 2*add
#x_design = seq(0, 1, length.out = p+1)[-(p+1)]
#x_design_extended = c(x_design[(p+1-add):p] - 1, x_design, x_design[1:add] + 1)
#x_design_extended_unit = to_unit_interval(x_design_extended)

x_test = seq(-1, 2, length.out = 3*p)

W = local_polynomial_weights(3*p, 0.3, p, T, m = 2, del = 1, eval.type = "diagonal", x.design.grid = x_test) # watch out for evaluation
#x_eval_orig = (1:p - 0.5)/p
#x_eval_extended = c(x_eval_orig[((p.eval + 1 )-(add/2)):p.eval] - 1, x_eval_orig, x_eval_orig[1:(add/2)] + 1)
#x_eval_extended %>% length()
#x_extended = c(x_orig[((p + 1 )-add):p] - 1, x_orig, x_orig[1:(add)] + 1)


#est_per_month_G10 = matrix(0, 116, 12)
#est_per_month_G01 = matrix(0, 116, 12)

g_hat = array(0, dim = c(p.eval, 3, 12))
for ( m in 1:12) {
  Y = N3[N3[,2] == m, -(1:3)] |>  
    observation_transformation(na.rm = T)
  g_hat[,,m] = eval_weights(W, Y)#[-c(1:(add/2), (p.eval + (add/2 + 1):add)),]
  #temp = eval_weights(W, Y)
  #est_per_month_G10[,m] = diag(temp[,,2])
  #est_per_month_G01[,m] = diag(temp[,,3])
}
est_per_month_G = g_hat[,1,]
est_per_month_G10 = g_hat[,2,]
est_per_month_G01 = g_hat[,3,]

colnames(est_per_month_G) = factor(1:12, labels = 1:12)
colnames(est_per_month_G10) = factor(1:12, labels = 1:12)
colnames(est_per_month_G01) = factor(1:12, labels = 1:12)
est_G = est_per_month_G %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12), #x = rep((1:(3*p)-0.5)/(3*p), each  = 12), 
         g = rep("G", each =length(x))) 
est_G10 = est_per_month_G10 %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12), #x = rep((1:(3*p)-0.5)/(3*p), each  = 12), 
         g = rep("G10", each =length(x))) 
est_G01 = est_per_month_G01 %>% 
  as_tibble() %>% 
  pivot_longer(cols = 1:12) %>% 
  mutate(x = rep(W$x.eval, each = 12), #x = rep((1:(3*p)-0.5)/(3*p), each  = 12), 
         g = rep("G01", each = length(x))) 
est_per_month = rbind(est_G01, est_G10) %>% mutate(g = factor(g, levels = c("G10", "G01")))

est_per_month %>% 
  # mutate(value = value/3) %>% 
  #  filter(x >= 1/3 & x <= 2/3) %>% 
  mutate(name = factor(name, levels = 1:12)) %>% 
  ggplot(aes( x = x, y = value, color = g)) + 
  geom_line() + 
  facet_wrap(name~.)


est_G %>% 
  mutate(name = factor(name, levels = 1:12)) %>% 
  ggplot(aes( x = x, y = value, color = g)) + 
  geom_line() + 
  facet_wrap(name~.)


time = N$UHRZEIT[7:150]
time 

final_est = est_per_month %>% 
  mutate(time = rep(rep(time, each = 12), 2)) %>% 
  mutate(name = factor(name, levels = 1:12)) 

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

ggsave("Grafics/weather_part_deriv_gamma_diagonal_all_months.png", device = "png", width = 8, height =6, units = "in")

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

ggsave("Grafics/weather_part_deriv_gamma_diagonal.png", device = "png", width = 5, height = 4, units = "in")


##### test for errors #####


N4 = colMeans(N3[N3[, 2] == 4, -(1:3)], na.rm = F)
plot(N4, type = "l")
length(N3)

mean_est = locpol::locPolSmootherC((1:(3*p)-0.5)/(3*p), N4, (1:(3*p)-0.5)/(3*p), 0.2, 2, kernel = locpol::EpaK)
plot(mean_est$beta0, type = "l")
plot(mean_est$beta1[(p+1):(2*p)], type = "l")


##### Estimation of full covariance of weather data #####
add = 48
N2 = cbind(N1[,1:3], 
           rbind(NA, N1[-length(N1$JAHR),(148-add):147]), 
           N1[, -(1:3)], 
           rbind(N1[-1, (1:add) + 3], NA))
N2[1:10, 1:100]

### entfernen der falschen Tage
N3 = N2[!(((N2[,3] %in% 28:29) & (N2[,2] == 2)) |
            (N2[,3] == 1) |
            ((N2[,3] == 31) & (N2[,2] %in% c(1,3,5,7,8,10,12))) |
            ((N2[,3] == 30) & (N2[,2] %in% c(4,6,9,11)))), ]
sum(is.na(N3))

to_rm = apply(N3, 1, function(v){any(is.na(v))}) %>% which()
N3 = N3[-to_rm, ]

rm(N, N0, N1, N2)

sum(is.na(N3))

p.eval = 100
p = 144 + 2*add
x.design = (0:(144-1))/144
x.design.extended = c( x.design[ (144-add+1):144 ] - 1, x.design, x.design[1:add] + 1)

W = local_polynomial_weights(p, 0.3, p.eval, T, m = 2, del = 1, 
                             eval.type = "full", x.design.grid = x.design.extended) # watch out for evaluation

g_hat = array(0, dim = c(p.eval, p.eval, 12))
for ( m in 1:12) {
  Y = N3[N3[,2] == m, -(1:3)] |>  
    observation_transformation()
  g_hat[,,m] = eval_weights(W, Y)[,,2]
}


time = (1:p.eval - 0.5)/p.eval * 24


cs2 = lisG10_UPcs2 = list(c(0, 1), c("lightblue", "darkred"))
UP = upper.tri(g_hat[,,3], diag = T)
estimate_03up = g_hat[,,3]
estimate_03dn = g_hat[,,3]
estimate_03dn[UP] = NA
estimate_03up[!UP] = NA
plot_ly() %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_03up, alpha = 0.9, colorscale = cs2, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_03dn, alpha = 0.9, colorscale = cs2, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  layout(
    scene = list(
      xaxis = list(title = "time"),
      yaxis = list(title = "time"),
      zaxis = list(title = "")
    )
  )%>%
  layout(
    scene = list(
      aspectratio = list(x = 1, y = 1, z = 1)  # Controls axis scaling
    )
  )
# save: cov01_march.png

UP = upper.tri(g_hat[,,11], diag = T)
estimate_11up = g_hat[,,11]
estimate_11dn = g_hat[,,11]
estimate_11dn[UP] = NA
estimate_11up[!UP] = NA
plot_ly() %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_11up, alpha = 0.9, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  add_surface(x = ~ time, y = ~ time, z = estimate_11dn, alpha = 0.9, lighting = list(
    ambient = 0.7, diffuse = 0.8, specular = 0.1, roughness = 0.9
  )) %>% 
  layout(
    scene = list(
      xaxis = list(title = "time"),
      yaxis = list(title = "time"),
      zaxis = list(title = "")
    )
  )  %>%
  layout(
    scene = list(
      aspectratio = list(x = 1, y = 1, z = 1)  # Controls axis scaling
    )
  )


#save: cov01_nov.png

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