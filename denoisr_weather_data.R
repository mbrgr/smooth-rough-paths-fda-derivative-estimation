#### installing the denoisr package ####
# install.packages("devtools")
# devtools::install_github("StevenGolovkine/denoisr")
library(denoisr)
library(somebm) # take care not available anymore on cran (only needed for the example)
library(ggplot2)

?estimate_H0_list()

t0 = seq(0, 1, length.out = 24*6 + 1)[-145]
t0

#### Load and process the data ####
load("data/weather_in_nuremberg/weather_data_nuremberg.RData")

N0 = N[7:length(N$MESS_DATUM), ] 
N0$MESS_DATUM |> as.character()
class(N0$MESS_DATUM)
N0$time = N0$MESS_DATUM |> format(format = "%H:%M")


N1 = N0 |> 
  dplyr::select(JAHR, MONAT, TAG, time, TT_10) |> 
  tidyr::pivot_wider(names_from = time,
              values_from = TT_10)

#### for all days together ####

N1 |> dim() # 7542
every_7 = seq(1, 7542, by = 10)
N1_7 = N1[every_7, ]
days_with_NA = apply(N1_7[, -c(1,2,3)], 1, function(x){any(is.na(x))})
N1_woNA = N1_7[!days_with_NA, ]
N1_woNA |> dim() # 1054

N_list = apply(N1_woNA, 1, function(day) {
  tmp_list = list(t = t0, x = day[-(1:3)])
})

names(N_list) <- apply(N1_woNA[, 1:3], 1, paste, collapse = "_")

M = 6 * 24
k0 = M * exp( - log(log(M))^2 )
#t0 = seq(0,1, length.out = 12)

denoisr_estimate_full = 
  estimate_H0_list(N_list, t0_list = t0, k0_list = floor(k0))

par(mfrow = c(1,1))
plot(denoisr_estimate_full ~ t0, type = "l", ylim = c(0,4))

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

plot(denoisr_estimate[[4]] ~ t0, type = "l", ylim =c(0,4))
par(mfrow = c(4,3))
for(i in 1:12){
  plot(denoisr_estimate[[i]] ~ t0, type = "l", ylim =c(0,4))
  abline(h = 1, lty = 2)
}

denoisr_df = data.frame(est = unlist(denoisr_estimate), 
                        time = rep(t0, 12),
                        month = gl(12, 144, labels = 
              c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))
denoisr_df |> 
  ggplot(aes(x = time, y = est, col = month)) + 
  geom_line()

#### example ####
par(mfrow = c(1,1))
X <- generate_fractional_brownian(N = 1000, M = 300, H = 0.5, sigma = 0.05)
H0 <- estimate_H0_list(X, t0_list = t0, k0_list = 14)
plot(H0, type = "l")

X_frac <- generate_piecewise_fractional_brownian(N = 100, M = 300, 
                                            H = c(0.2, 0.5, 0.8), 
                                            sigma = 0.05)
H0_frac <- estimate_H0_list(X_frac, t0_list = t0, k0_list = rep(c(2,4,6), each = 144/3))
plot(H0_frac, type = "l")
H0_frac_2 <- estimate_H0_list(X_frac, t0_list = t0, k0_list = 6)
plot(H0_frac_2, type = "l")
