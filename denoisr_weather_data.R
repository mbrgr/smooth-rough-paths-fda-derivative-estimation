#### installing the denoisr package ####
# install.packages("devtools")
# devtools::install_github("StevenGolovkine/denoisr")
library(denoisr)
library(somebm) # take care not available anymore on cran (only needed for the example)
library(ggplot2)

?estimate_H0_list()

deriv_est_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

# use colour theme from previous analysis
farben = c( "#a6cee3", "#03396c",  # Blau-Töne
            "#33a02c", "#66c21f", "#006400",  # Grün-Töne
            "#e31a1c", "#fb9a99", "#990000",  # Rot-Töne
            "#ff7f00", "#ffb300", "#b15928",  # Gelb-Orange-Töne
            "#1f78b4")

t0 = seq(0, 1, length.out = 24*6 + 1)[-145]
t0

#### Load and process the data ####
load("data/weather_in_nuremberg/weather_data_nuremberg.RData")

N0 = N[7:length(N$MESS_DATUM), ] 
N0$MESS_DATUM |> as.character()
class(N0$MESS_DATUM)
N0$time = N0$MESS_DATUM |> format(format = "%H:%M")
length(N0$time)

N1 = N0 |> 
  dplyr::select(JAHR, MONAT, TAG, time, TT_10) |> 
  tidyr::pivot_wider(names_from = time,
              values_from = TT_10)

#### for all days together ####

N1 |> dim() # 7542
days_with_NA = apply(N1[, -c(1,2,3)], 1, function(x){any(is.na(x))})
N1_woNA = N1[!days_with_NA, ]
N1_woNA |> dim() # 1054



# set M and k0 according to the paper
M = 6 * 24
k0 = M * exp( - log(log(M))^2 )
#t0 = seq(0,1, length.out = 12)
#k0 = 2

#### seperately for the months ####
N1_woNA |> dim() # 1054
N1_month = N1_woNA |> 
  dplyr::group_by(MONAT) |> 
  dplyr::group_map(.f = function(m, ...) {
    apply(m[,-(1:2)], 1, function(day) {
      tmp_list = list(t = t0, x = day)
    })
  })

# t0 =0.5
denoisr_estimate = lapply(1:12, function(i) {
  estimate_H0_list(N1_month[[i]], t0_list = t0, k0_list = floor(k0)) 
  }
)

denoisr_df = data.frame(est = unlist(denoisr_estimate), 
                        time = rep(t0, 12),
                        month = gl(12, 144, labels = 
              c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))
denoisr_df |> 
  ggplot(aes(x = time, y = est, col = month)) + 
  geom_line() + 
  labs(title = "Estimating the pointwise regularity with denoisr", 
       y = "estimated smoothness", x = "time" ) +
  facet_wrap(.~month) +
  scale_colour_manual(values = farben) + 
  scale_x_continuous(
    breaks = c(.25, .5, .75),
    labels = c("06:00", "12:00", "18:00")
  ) 
ggsave("grafics/denoisr_weater_all_months_k11.pdf", device = "pdf",
       width = 7, height = 5, units = "in")

#### example denoisr authors ####
par(mfrow = c(1,1))
X <- generate_fractional_brownian(N = 1000, M = 300, H = 0.5, sigma = 0.05)
M = 300
k0 = M * exp( - log(log(M))^2 )
H0 <- estimate_H0_list(X, t0_list = t0, k0_list = floor(k0))
plot(H0 ~ t0, type = "l")

X_frac <- generate_piecewise_fractional_brownian(N = 100, M = 300, 
                                            H = c(0.2, 0.5, 0.8), 
                                            sigma = 0.05)
H0_frac <- estimate_H0_list(X_frac, t0_list = t0, k0_list = rep(c(2,4,6), each = 144/3))
plot(H0_frac ~ t0, type = "l")
H0_frac_2 <- estimate_H0_list(X_frac, t0_list = t0, k0_list = k0)
plot(H0_frac_2 ~ t0, type = "l")
