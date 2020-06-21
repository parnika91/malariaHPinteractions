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

bs_p = bootstrap_p(m_m, no_of_sims = 10, threads = parallel::detectCores())
bs_p$p #0.6


data_alt <- dislnorm$new(n)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)
# test stat positive, p value not sig

data_alt <- disexp$new(n)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)
# test stat pos, pval sig

data_alt <- dispois$new(n)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)
# test stat pos, pval sig

## continuous gives pval 0 for n
bw = df$bw[which(df$bw!=0)]
m_m = conpl$new(bw)
m_m$getXmin()
m_m$getPars()
m_m$setXmin(0.1e-05)
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
bs_p$p # 0.06, means that it is probable that power law can describe the dist

data_alt <- conlnorm$new(bw)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)
# test statistic -2.4 means it is better described by lognorm and sig better (p twosided = 0.01)


data_alt <- conexp$new(bw)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)
# test stat: 4.3, better desc by power law, sig (P 2sided = 1.5e0.5)


data_alt <- conpois$new(bw)
data_alt$xmin <- est$xmin
data_alt$pars <- estimate_pars(data_alt)
comp <- compare_distributions(m_m, data_alt)
# same result as exp dist

# check normality of closeness
shapiro.test(ig_cl)
