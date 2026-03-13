library(denoisr)
library(ggplot2)
library(ffscb)

deriv_est_theme = theme_grey(base_size = 15) + 
  theme(plot.title = element_text(size = 14))

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
data.frame(denoisr_estimate_full, t0*100) |> 
  ggplot(aes(x = t0, y = denoisr_estimate_full)) + 
  geom_line() + 
  labs(y = "estimated smoothness", x = "% of stance phase",  
       title = "estimating the pointwise smoothness with denoisr") + 
  deriv_est_theme
ggsave("grafics/denoisr_bio_k12.pdf", device = "pdf",
       width = 7, height = 5, units = "in")


plot(denoisr_estimate_full ~ t0, type = "l")

denoisr_estimate_full = 
  estimate_H0_list(bio_list_diff, t0_list = t0, k0_list = 2)
plot(denoisr_estimate_full ~ t0, type = "l")
