

x <- seq(min(S7), max(S7), length= 100)
f7 <- fitdist(S7[S7>0], "lnorm")
densi.f7 <- dlnorm(x, f7$estimate[1], f7$estimate[2])
lines(x, densi.f7, col ='green')