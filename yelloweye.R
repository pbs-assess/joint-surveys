library(dplyr)
library(ggplot2)
library(ggsidekick) # for fourth_root_power_trans
theme_set(ggsidekick::theme_sleek())
library(sdmTMB)
source("utils.R")
dir.create("data-generated", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

if (file.exists("../../gfs/report/data-cache/yelloweye-rockfish.rds")) {
  d <- readRDS("../../gfs/report/data-cache/yelloweye-rockfish.rds")
  d <- d$survey_sets
  d <- dplyr::filter(d, survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S"))
} else {
  d <- gfdata::get_survey_sets("yelloweye rockfish", ssid = c(22, 36))
}
d <- d %>% select(
  survey_abbrev, year, longitude, latitude, density_ppkm2,
  grouping_code, depth_m
) %>%
  rename(survey = survey_abbrev) %>%
  mutate(density_1000ppkm2 = density_ppkm2 / 1000)

# log and scale the depth predictor:
d$depth_log <- log(d$depth_m)
d$depth_centred <- d$depth_log - mean(d$depth_log)
d$depth_scaled <- d$depth_centred / sd(d$depth_centred)
d_utm <- convert2utm(d, coords = c("longitude", "latitude"))

ggplot(d_utm, aes(X, Y,
  size = density_1000ppkm2,
  colour = survey
)) +
  facet_wrap(~year) +
  geom_point(pch = 21) +
  scale_size_area()
ggsave("figs/hbll-joint-raw-data.pdf", width = 10, height = 10)

north_grid <- gfplot::hbll_n_grid$grid
south_grid <- gfplot::hbll_s_grid$grid
north_grid_utm <- convert2utm(north_grid)
south_grid_utm <- convert2utm(south_grid)

years <- sort(unique(d_utm$year))
north_grid_utm <- expand_prediction_grid(north_grid_utm, years = years) %>%
  mutate(survey = "HBLL OUT N") %>%
  mutate(depth_centred = log(depth) - mean(d$depth_log)) %>%
  mutate(depth_scaled = depth_centred / sd(d$depth_centred))
south_grid_utm <- expand_prediction_grid(south_grid_utm, years = years) %>%
  mutate(survey = "HBLL OUT S") %>%
  mutate(depth_centred = log(depth) - mean(d$depth_log)) %>%
  mutate(depth_scaled = depth_centred / sd(d$depth_centred))
joint_grid_utm <- dplyr::bind_rows(north_grid_utm, south_grid_utm)

sp <- make_spde(d_utm$X, d_utm$Y, n_knots = 175)
plot_spde(sp)

model_file <- "data-generated/hbll-out-joint.rds"
if (!file.exists(model_file)) {
  tictoc::tic()
  m <- sdmTMB(
    formula = density_1000ppkm2 ~ 0 +
      as.factor(year) + depth_scaled + I(depth_scaled^2),
    data = d_utm,
    spde = sp,
    time = "year",
    silent = FALSE,
    anisotropy = FALSE,
    ar1_fields = FALSE,
    include_spatial = TRUE,
    family = tweedie(link = "log")
  )
  tictoc::toc()
  saveRDS(m, file = model_file)
} else {
  m <- readRDS(model_file)
}
m

d_utm$resids <- residuals(m) # randomized quantile residuals
hist(d_utm$resids)
pdf("figs/hbll-joint-residuals-qq.pdf", width = 5, height = 5)
par(cex = 0.75)
qqnorm(d_utm$resids)
abline(a = 0, b = 1)
dev.off()

ggplot(d_utm, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave("figs/hbll-joint-residual-map.pdf", width = 10, height = 10)

predictions <- predict(m,
  newdata = joint_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y")
)

plot_map <- function(dat, column) {
  ggplot() +
    geom_point(
      data = dat,
      aes_string("X", "Y", colour = column), size = 0.5
    ) +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(predictions$data, "exp(est)") +
  scale_colour_viridis_c(option = "C") +
  ggtitle("Prediction (fixed effects + all random effects)")
ggsave("figs/hbll-joint-prediction.pdf", width = 10, height = 10)

plot_map(predictions$data, "exp(est)") +
  scale_colour_viridis_c(trans = "fourth_root_power", option = "C") +
  ggtitle("Prediction (fourth root transformed colour; fixed effects + all random effects)")
ggsave("figs/hbll-joint-prediction-sqrt.pdf", width = 10, height = 10)

plot_map(predictions$data, "exp(est_non_rf)") +
  ggtitle("Fixed effects only") +
  scale_colour_viridis_c(trans = "fourth_root_power", option = "C")
ggsave("figs/hbll-joint-non-rf.pdf", width = 10, height = 10)

plot_map(predictions$data, "est_rf") +
  ggtitle("All spatial and spatiotemporal random effects") +
  scale_colour_gradient2()
ggsave("figs/hbll-joint-rf.pdf", width = 10, height = 10)

plot_map(filter(predictions$data, year == 2018), "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_colour_gradient2()
ggsave("figs/hbll-joint-omega.pdf", width = 5, height = 5)

plot_map(predictions$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_colour_gradient2()
ggsave("figs/hbll-joint-epsilon.pdf", width = 10, height = 10)

# not bias correcting for speed for now:
ind <- get_index(predictions, bias_correct = FALSE)

scale <- 2 * 2 # 2 x 2 km grid
ggplot(ind, aes(year, est * scale)) + geom_line() +
  geom_ribbon(aes(ymin = lwr * scale, ymax = upr * scale), alpha = 0.4) +
  xlab("Year") + ylab("Estimated density (1000s of fish)") +
  geom_vline(xintercept = seq(2006, 2018, 2), lty = 2, alpha = 0.2)
ggsave("figs/hbll-index.pdf", width = 8, height = 5)

# what about the individual surveys? -----------------------

pred_north <- predict(m,
  newdata = north_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y")
)
ind_north <- get_index(pred_north, bias_correct = FALSE)

pred_south <- predict(m,
  newdata = south_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y")
)
ind_south <- get_index(pred_south, bias_correct = FALSE)

ind_north$type <- "HBLL N"
ind_south$type <- "HBLL S"
ind$type <- "HBLL all"

bind_rows(ind_north, ind_south) %>%
  bind_rows(ind) %>%
  ggplot(aes(year, est * scale)) + geom_line(aes(colour = type)) +
  geom_ribbon(aes(ymin = lwr * scale, ymax = upr * scale, fill = type), alpha = 0.4, colour = NA) +
  facet_wrap(~type, ncol = 1) +
  xlab("Year") + ylab("Estimated density (1000s of fish)") +
  geom_vline(xintercept = seq(2006, 2018, 2), lty = 2, alpha = 0.2) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")
ggsave("figs/hbll-index-components.pdf", width = 5.5, height = 8)


# # -----------------------
#
#
# d_utm_north <- dplyr::filter(d_utm, survey == "HBLL OUT N")
# sp_north <- make_spde(d_utm_north$X, d_utm_north$Y, n_knots = 125)
# plot_spde(sp_north);axis(1);axis(2)
#
# ggplot(d_utm_north, aes(X, Y, size = density_1000ppkm2)) +
#   facet_wrap(~year) +
#   geom_point(pch = 21) +
#   scale_size_area()
#
# ggplot(d_utm_north, aes(X, Y, size = depth_scaled)) +
#   facet_wrap(~year) +
#   geom_point(pch = 21) +
#   scale_size_area()
#
# m_north <- sdmTMB(
#   formula = density_1000ppkm2 ~ 0 +
#     as.factor(year) + depth_scaled + I(depth_scaled^2),
#   data = d_utm_north,
#   spde = sp_north,
#   time = "year",
#   silent = FALSE,
#   anisotropy = FALSE,
#   ar1_fields = FALSE,
#   family = tweedie(link = "log")
# )
# m_north
#
# north_grid_utm <- gfplot::hbll_n_grid$grid %>%
#   convert2utm() %>%
#   expand_prediction_grid(years = unique(d_utm_north$year)) %>%
#   mutate(survey = "HBLL OUT N") %>%
#   mutate(depth_centred = log(depth) - mean(d$depth_log)) %>%
#   mutate(depth_scaled = depth_centred / sd(d$depth_centred))
#
# predictions_north <- predict(m_north, newdata = north_grid_utm,
#   return_tmb_object = TRUE, xy_cols = c("X", "Y"))
#
# plot_map <- function(dat, column) {
#   ggplot() +
#     geom_point(data = filter(dat, survey == "HBLL OUT N"),
#       aes_string("X", "Y", colour = column), size = 0.5) +
#     facet_wrap(~year) +
#     coord_fixed()
# }
#
# # predictions_north$data %>%
# #   group_by(year) %>%
# #   mutate(est_centred = est - mean(est)) %>%
# #   plot_map("est_centred") +
# #   scale_colour_gradient2() +
# #   ggtitle("Prediction (fixed effects + all random effects)")
#
# library(gfsynopsis) # for fourth_root_power_trans
# trans <- "fourth_root_power"
#
# predictions_north$data %>%
#   plot_map("exp(est)") +
#   scale_colour_viridis_c(trans = trans) +
#   ggtitle("Prediction (fixed effects + all random effects)")
#
# plot_map(predictions_north$data, "exp(est_non_rf)") +
#   ggtitle("Fixed effects only") +
#   scale_colour_viridis_c(trans = trans)
#
# plot_map(filter(predictions_north$data, year %in% 2017), "omega_s") +
#   ggtitle("Spatial random effects only") +
#   scale_colour_gradient2()
#
# plot_map(predictions_north$data, "epsilon_st") +
#   ggtitle("Spatiotemporal random effects only") +
#   scale_colour_gradient2()
