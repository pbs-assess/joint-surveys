library(dplyr)
library(ggplot2)
# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick) # for fourth_root_power_trans
theme_set(ggsidekick::theme_sleek())
library(sdmTMB)
source("utils.R")
dir.create("data-generated", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

f <- "data-generated/yelloweye-rockfish-inside.rds"
if (file.exists(f)) {
  d <- readRDS(f)
} else {
  d <- gfdata::get_survey_sets("yelloweye rockfish",
    ssid = c(22, 36, 39, 40))
  saveRDS(d, file = f)
}
d <- d %>%
  filter(survey_series_id %in% c(39, 40)) %>% # inside only
  select(
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

joint_grid <- readRDS("data-generated/hbll-inside-grid.rds")
joint_grid_utm <- convert2utm(joint_grid)
years <- sort(unique(d_utm$year))
joint_grid_utm <- expand_prediction_grid(joint_grid_utm, years = years) %>%
  mutate(depth_centred = log(depth) - mean(d$depth_log)) %>%
  mutate(depth_scaled = depth_centred / sd(d$depth_centred))
joint_grid_utm <- mutate(joint_grid_utm, year_fake = ifelse(year == 2003, 2004, year))
joint_grid_utm <- mutate(joint_grid_utm, Y_cent = Y - mean(d_utm$Y))

# n_years <- filter(d_utm, survey %in% "HBLL INS N") %>% pull(year) %>% unique()
# s_years <- filter(d_utm, survey %in% "HBLL INS S") %>% pull(year) %>% unique()
# north_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS N", year %in% n_years)
# south_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS S", year %in% s_years)
north_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS N")
south_grid_utm <- filter(joint_grid_utm, survey %in% "HBLL INS S")

sp <- make_spde(d_utm$X, d_utm$Y, n_knots = 175)
plot_spde(sp)

# d_utm <- mutate(d_utm, year_fake = ifelse(year %in% c(2003, 2004), 2005, year))

d_utm$Y_cent <- d_utm$Y - mean(d_utm$Y)
model_file <- "data-generated/hbll-inside-joint.rds"
if (!file.exists(model_file)) {
  tictoc::tic()
  m <- sdmTMB(
    formula = density_1000ppkm2 ~ 0 +
      Y_cent + I(Y_cent^2) +
      as.factor(year) + depth_scaled + I(depth_scaled^2),
    data = d_utm,
    spde = sp,
    time = "year",
    silent = FALSE,
    anisotropy = TRUE,
    ar1_fields = TRUE,
    include_spatial = FALSE,
    control = sdmTMBcontrol(step.min = 0.5),
    family = tweedie(link = "log")
  )
  tictoc::toc()
  saveRDS(m, file = model_file)
} else {
  m <- readRDS(model_file)
}
m

predictions <- predict(m,
  newdata = joint_grid_utm,
  return_tmb_object = TRUE, xy_cols = c("X", "Y")
)
# not bias correcting for speed for now:
ind <- get_index(predictions, bias_correct = FALSE)

d_utm$resids <- residuals(m) # randomized quantile residuals
hist(d_utm$resids)
pdf("figs/hbll-joint-residuals-qq.pdf", width = 5, height = 5)
par(cex = 0.75)
qqnorm(d_utm$resids)
abline(a = 0, b = 1)
dev.off()

ggplot(d_utm, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point(size = 0.5) + facet_wrap(~year) + coord_fixed()
ggsave("figs/hbll-joint-residual-map.pdf", width = 10, height = 10)

plot_map <- function(dat, column) {
  ggplot() +
    geom_point(data = dat, aes_string("X", "Y", colour = column), size = 0.5) +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(predictions$data, "exp(est)") +
  scale_colour_viridis_c(option = "C") +
  ggtitle("Prediction (fixed effects + all random effects)")
ggsave("figs/hbll-joint-prediction.pdf", width = 10, height = 10)

plot_map(predictions$data, "exp(est)") +
  scale_colour_viridis_c(trans = "fourth_root_power", option = "C") +
  ggtitle(paste0("Prediction (fourth root transformed colour; ",
  "fixed effects + all random effects)"))
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

ind_north$type <- "HBLL INS N"
ind_south$type <- "HBLL INS S"
ind$type <- "HBLL INS all"

n_years <- filter(d_utm, survey %in% "HBLL INS N") %>% pull(year) %>% unique()
s_years <- filter(d_utm, survey %in% "HBLL INS S") %>% pull(year) %>% unique()

ind_north_plot <- filter(ind_north, year %in% n_years)
ind_south_plot <- filter(ind_south, year %in% s_years)
all_plot <- bind_rows(ind_north_plot, ind_south_plot) %>%
  bind_rows(ind)

bind_rows(ind_north, ind_south) %>%
  bind_rows(ind) %>%
  ggplot(aes(year, est * scale)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(colour = type), pch = 21) +
  geom_point(data = all_plot, aes(colour = type)) +
  geom_ribbon(aes(ymin = lwr * scale, ymax = upr * scale, fill = type),
    alpha = 0.2, colour = NA) +
  facet_wrap(~type, ncol = 1) +
  xlab("Year") + ylab("Estimated density (1000s of fish)") +
  geom_vline(xintercept = seq(2004, 2018, 2), lty = 2, alpha = 0.2, lwd = 0.2) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(breaks = seq(2004, 2018, 2))
ggsave("figs/hbll-index-components.pdf", width = 5.5, height = 8.5)
