library(denoisr)
library(ggplot2)
library(ffscb)

data("Biomechanics")
Biomechanics |> dim()
head(Biomechanics)
t0 = Biomechanics$Stance_Phase
bio_list_diff = lapply(2:19, function(i){
  tmp_list = list(t = Biomechanics[, 1], 
                  x = Biomechanics[, i] - Biomechanics[, 18 + i])
})
bio_list_diff


M = 201
k0 = M * exp( - log(log(M))^2 )
#t0 = seq(0,1, length.out = 12)

denoisr_estimate_full = 
  estimate_H0_list(bio_list_diff, t0_list = t0, k0_list = floor(k0))
plot(denoisr_estimate_full ~ t0, type = "l")

denoisr_estimate_full = 
  estimate_H0_list(bio_list_diff, t0_list = t0, k0_list = 2)
plot(denoisr_estimate_full ~ t0, type = "l")
