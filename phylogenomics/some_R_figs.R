# Rate heterogeneity graph
#gamma(x)

position <- seq(1:1000)
quant <- quantile (position)
#rgamma(nb observations, shape)
# 5 mutations
rg <- rgamma(1000, 1)


plot (position, rg)
abline(v=unlist(lowess(position, rg))) # better way ...
length(rg)
length(position)
length(rg) == length(position)
range(rg)
(
dgamma())


plot(pgamma(quant, shape = 1, rate = 5))
