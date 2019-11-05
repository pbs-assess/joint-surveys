library(sdmTMB)
library(ggplot2)
library(dplyr)

d <- tibble::tibble(X = runif(1000, 0, 1), Y = runif(1000, 0, 1))
d$year <- rep(1:10, each = 100)
d$type <- "a_obs"

nd <- expand.grid(X = seq(0, 1, length.out = 25), Y = seq(0, 1, length.out = 25), year = 1:10)
nd$type <- "b_grid"

.d <- rbind(d, nd)

# dat <- sim(x = d$X, y = d$Y, X = lapply(1:10, function(x) d$Y),
  # initial_betas = )

.rf_sim <- function(model, x, y) {
  out <- sdmTMB:::rf_sim(model, x, y)
  out - mean(out)
}

set.seed(1)
sigma_O <- 0.8
sigma_E <- 0.4
kappa <- 0.2
x = .d$X
y = .d$Y
rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1/kappa)
omega_s <- .rf_sim(model = rf_omega, x, y)
.d$omega_s <- omega_s

# ggplot(d, aes_string("X", "Y", colour = "omega_s")) +
#   geom_point() +
#   scale_color_gradient2()

rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1/kappa)

.d <- .d %>% group_by(year) %>%
  group_split() %>%
  purrr::map_dfr(
    ~tibble(X = .x$X, Y = .x$Y, type = .x$type, omega_s = .x$omega_s, year = .x$year,
      eps_st = .rf_sim(model = rf_epsilon, .x$X, .x$Y))) %>%
  ungroup() %>%
  arrange(type)

ggplot(.d, aes_string("X", "Y", colour = "eps_st")) +
  geom_point() +
  scale_color_gradient2() +
  facet_grid(type~year)

f <- ~ 0 + as.factor(year)
X_ij <- model.matrix(f, .d)

# b_j <- c(rlnorm(10, meanlog = 1, sdlog = 0.4), 0.2)
b_j <- c(rlnorm(10, meanlog = 1, sdlog = 0.4))
b_j

.d$eta <- as.numeric(.d$omega_s + .d$eps_st + X_ij %*% b_j)
.d$y <- exp(.d$eta)
.d$y[1:1000] <- rlnorm(1000, meanlog = .d$eta[1:1000], sdlog = 0.3) # observations

ggplot(.d, aes_string("X", "Y", colour = "y")) +
  geom_point() +
  scale_color_viridis_c() +
  facet_grid(type~year)

d <- filter(.d, type == "a_obs")
nd <- filter(.d, type != "a_obs")

sp <- make_spde(d$X, d$Y, n_knots = 100)
plot_spde(sp)
m <- sdmTMB(y ~ 0 + as.factor(year), data = d, spde = sp,
  family = sdmTMB::lognormal(link = "log"),
  silent = FALSE, time = "year")
m

p <- predict(m, newdata = nd, return_tmb_object = TRUE, xy_cols = c("X", "Y"))

ggplot(p$data, aes_string("X", "Y", fill = "est")) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~year)

.i <- get_index(p)

actual <- group_by(nd, year) %>%
  summarise(total = sum(y))

.f <- exp(mean(log(actual$total))) / exp(mean(log(.i$est)))

actual$total <- scale(actual$total)
.i$est <- scale(.i$est)

ggplot(.i, aes(year, est)) + geom_line() +
  geom_line(data = actual, mapping = aes(x = year, y = total),
    inherit.aes = FALSE, colour = "red")
