m_m = displ$new(n)
m_m$getXmin()
m_m$getPars()
m_m$setXmin(2)
m_m$setPars(5)
(est = estimate_pars(m_m))
(est = estimate_xmin(m_m))
m_m$setXmin(est)
png("degree_distribution.png")
plot(m_m)
lines(m_m, col = 2)
dev.off()
dd = plot(m_m)
head(dd, 3)

# bs = bootstrap(m_m, no_of_sims = 10, threads = parallel::detectCores())
# hist(bs$bootstraps[, 2], breaks = "fd")
# hist(bs$bootstraps[, 3], breaks = "fd")
# plot(jitter(bs$bootstraps[, 2], factor = 1.2), bs$bootstraps[, 3])

bs_p = bootstrap_p(m_m, no_of_sims = 50, threads = parallel::detectCores())



data_alt <- dislnorm$new(n)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)


data_alt <- disexp$new(n)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)


data_alt <- dispois$new(n)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)


## continuous gives pval 0 for n
bw = df$bw[which(df$bw!=0)]
m_m = conpl$new(bw)
m_m$getXmin()
m_m$getPars()
m_m$setXmin(0.5)
m_m$setPars(5)
(est = estimate_pars(m_m))
(est = estimate_xmin(m_m))
m_m$setXmin(est)
png("para_betweenness_distribution.png")
plot(m_m)
lines(m_m, col = 2)
dev.off()
dd = plot(m_m)
head(dd, 3)

# bs = bootstrap(m_m, no_of_sims = 10, threads = parallel::detectCores())
# hist(bs$bootstraps[, 2], breaks = "fd")
# hist(bs$bootstraps[, 3], breaks = "fd")
# plot(jitter(bs$bootstraps[, 2], factor = 1.2), bs$bootstraps[, 3])

bs_p = bootstrap_p(m_m, no_of_sims = 50, threads = parallel::detectCores())

data_alt <- conlnorm$new(bw)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)


data_alt <- conexp$new(bw)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)


data_alt <- conpois$new(bw)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)

# check normality of closeness
shapiro.test(df$cl)
