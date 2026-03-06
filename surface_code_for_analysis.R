# ------------------------------
# Packages
# ------------------------------
library(terra)
library(tmap)
library(tidyverse)
library(tidyr)
library(purrr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(lattice)
library(MuMIn)
library(flextable)
library(broom.mixed)


# ------------------------------
# Function 1: compute mean slope time-series for one run
# ------------------------------
get_slope_ts <- function(dem_dir, run_id) {
  
  keep_names <- c("original.tif","20m.tif","40m.tif","60m.tif",
                  "80m.tif","100m.tif","120m.tif")
  dem_files  <- file.path(dem_dir, keep_names)
  stopifnot(all(file.exists(dem_files)))
  
  # parse time from filenames
  base <- tools::file_path_sans_ext(basename(dem_files))
  time_min <- ifelse(base == "original", 0,
                     as.numeric(sub("m$", "", base)))
  stopifnot(!any(is.na(time_min)))
  
  # order by time
  ord <- order(time_min)
  dem_files <- dem_files[ord]
  time_min  <- time_min[ord]
  
  # read DEMs and compute slope
  dem_list <- lapply(dem_files, rast)
  dem_list <- lapply(dem_list, function(r) {
    names(r) <- "elevation"; r
  })
  
  slope_list <- lapply(dem_list, function(dem_r) {
    terrain(dem_r, v = "slope",
            unit = "degrees", neighbors = 8)
  })
  
  # mean slope per time step
  slope_mean <- vapply(
    slope_list,
    function(s) global(s, "mean", na.rm = TRUE)[1, 1],
    numeric(1)
  )
  
  data.frame(
    run        = run_id,
    time_min   = time_min,
    slope_mean = slope_mean
  )
}

# ------------------------------
# Function 2: visualisation at 0 and 120 minutes
# ------------------------------
plot_slope_0_120 <- function(dem_dir, run_id,
                             hs_angle = 45, hs_dir = 315,
                             slope_alpha = 0.65) {
  
  tmap_options(max.raster = 1e8)
  
  keep_names <- c("original.tif","20m.tif","40m.tif","60m.tif",
                  "80m.tif","100m.tif","120m.tif")
  dem_files  <- file.path(dem_dir, keep_names)
  stopifnot(all(file.exists(dem_files)))
  
  base <- tools::file_path_sans_ext(basename(dem_files))
  time_min <- ifelse(base == "original", 0, as.numeric(sub("m.*$", "", base)))
  stopifnot(!any(is.na(time_min)))
  
  f0   <- dem_files[time_min == 0][1]
  f120 <- dem_files[time_min == 120][1]
  
  dem0   <- rast(f0)
  dem120 <- rast(f120)
  
  # ---- SLOPE for colouring (degrees) ----
  slope0_deg   <- terrain(dem0,   v = "slope", unit = "degrees", neighbors = 8)
  slope120_deg <- terrain(dem120, v = "slope", unit = "degrees", neighbors = 8)
  
  # ---- SLOPE/ASPECT for hillshade MUST be radians ----
  slope0_rad   <- terrain(dem0,   v = "slope",  unit = "radians", neighbors = 8)
  aspect0_rad  <- terrain(dem0,   v = "aspect", unit = "radians", neighbors = 8)
  slope120_rad <- terrain(dem120, v = "slope",  unit = "radians", neighbors = 8)
  aspect120_rad<- terrain(dem120, v = "aspect", unit = "radians", neighbors = 8)
  
  hs0   <- shade(slope0_rad,   aspect0_rad,   angle = hs_angle, direction = hs_dir)
  hs120 <- shade(slope120_rad, aspect120_rad, angle = hs_angle, direction = hs_dir)
  
  # Optional: boost hillshade contrast a bit (often helps a lot)
  hs0   <- clamp((hs0   - 0.4) * 1.4 + 0.4, 0, 1)
  hs120 <- clamp((hs120 - 0.4) * 1.4 + 0.4, 0, 1)
  
  # ---- Paper-like breaks (more resolution at low slopes) ----
  breaks <- c(0, 2, 5, 8, 12, 16, 20, 25, 30, 35, 45, 60, 90)
  
  # ---- Paper-like palette (green -> yellow -> orange -> brown -> dark) ----
  slope_cols <- c(
    "#1a9850", "#66bd63", "#a6d96a", "#d9ef8b",
    "#fee08b", "#fdae61", "#f46d43", "#d73027",
    "#8c510a", "#3b1f0b"
  )
  
  make_map <- function(hs, slope_deg, title_txt) {
    tm_shape(hs) +
      tm_raster(
        col.scale  = tm_scale_continuous(values = "grey10"),
        col.legend = tm_legend(show = FALSE)
      ) +
      tm_shape(slope_deg) +
      tm_raster(
        col.scale  = tm_scale_intervals(values = slope_cols, breaks = breaks, style = "fixed"),
        col_alpha  = slope_alpha,
        col.legend = tm_legend(title = "Slope (deg)")
      ) +
      tm_layout(frame = FALSE, 
                legend.text.size = 0.5,
                legend.title.size = 0.6,
                legend.position = c("right", "bottom"),
                legend.bg.color = "white",
                legend.bg.alpha = 0.7,
                legend.frame = F,
                legend.just = "right",
                inner.margins = c(0.02, 0.02, 0.01, 0.2)) +
      tm_title(title_txt)
  }
  
  map_slope0   <- make_map(hs0,   slope0_deg,   paste0("Slope - ", run_id, " (0 min)"))
  map_slope120 <- make_map(hs120, slope120_deg, paste0("Slope - ", run_id, " (120 min)"))
  
  print(tmap_arrange(map_slope0, map_slope120))
}


# ------------------------------
# PHASE 1: load all runs
# ------------------------------
ts_all <- bind_rows(
  get_slope_ts("./Surface_analysis/Experiment_1/DEM/GRA-1-LO/", "GRA-1-LO"),
  get_slope_ts("./Surface_analysis/Experiment_2/DEM/GRA-2-LO/", "GRA-2-LO"),
  get_slope_ts("./Surface_analysis/Experiment_3/DEM/GRA-3-LO/", "GRA-3-LO"),
  get_slope_ts("./Surface_analysis/Experiment_4/DEM/GRA-4-LO/", "GRA-4-LO"),
  get_slope_ts("./Surface_analysis/Experiment_5/DEM/GRA-5-LO/", "GRA-5-LO"),
  get_slope_ts("./Surface_analysis/Experiment_6/DEM/GRA-6-LO/", "GRA-6-LO"),
  get_slope_ts("./Surface_analysis/Experiment_7/DEM/GRA-7-LO/", "GRA-7-LO"),
  get_slope_ts("./Surface_analysis/Experiment_8/DEM/GRA-8-LO/", "GRA-8-LO")
)

# plot figures using tmap to visualise the change in the surface
plot_slope_0_120("./Surface_analysis/Experiment_1/DEM/GRA-1-LO/", "GRA-1-LO")
plot_slope_0_120("./Surface_analysis/Experiment_2/DEM/GRA-2-LO/", "GRA-2-LO")
plot_slope_0_120("./Surface_analysis/Experiment_3/DEM/GRA-3-LO/", "GRA-3-LO")
plot_slope_0_120("./Surface_analysis/Experiment_4/DEM/GRA-4-LO/", "GRA-4-LO")
plot_slope_0_120("./Surface_analysis/Experiment_5/DEM/GRA-5-LO/", "GRA-5-LO")
plot_slope_0_120("./Surface_analysis/Experiment_6/DEM/GRA-6-LO/", "GRA-6-LO")
plot_slope_0_120("./Surface_analysis/Experiment_7/DEM/GRA-7-LO/", "GRA-7-LO")
plot_slope_0_120("./Surface_analysis/Experiment_8/DEM/GRA-8-LO/", "GRA-8-LO")


# ------------------------------
# PHASE 2: add categorical variables
# ------------------------------
ts_all_df <- ts_all %>%
  mutate(
    contact_area = case_when(
      run %in% c("GRA-1-LO","GRA-2-LO","GRA-5-LO","GRA-6-LO") ~ "Large",
      TRUE ~ "Small"
    ),
    material = case_when(
      run %in% c("GRA-1-LO","GRA-3-LO","GRA-5-LO","GRA-7-LO") ~ "Dolerite",
      TRUE ~ "Quartzite"
    )
  )

ts_all_df$contact_area <- factor(ts_all_df$contact_area,
                                 levels = c("Large","Small"))
ts_all_df$material <- factor(ts_all_df$material,
                             levels = c("Dolerite","Quartzite"))



# ------------------------------
# PHASE 3: per-run linear models and extract stats
# ------------------------------

pearson_table <- ts_all_df %>%
  group_by(run) %>%
  summarise(
    cor_test = list(cor.test(time_min, slope_mean, method = "pearson")),
    lm_fit   = list(lm(slope_mean ~ time_min)),
    .groups  = "drop"
  ) %>%
  mutate(
    pearson_r = map_dbl(cor_test, ~ unname(.x$estimate)),     # correlation coefficient r
    df        = map_dbl(cor_test, ~ unname(.x$parameter)),    # df = n - 2
    t_value   = map_dbl(cor_test, ~ unname(.x$statistic)),    # t statistic
    p_value   = map_dbl(cor_test, ~ .x$p.value),              # p value
    r2        = map_dbl(lm_fit,   ~ summary(.x)$r.squared)    # regression R^2
  ) %>%
  select(run, pearson_r, df, t_value, r2, p_value) %>%
  arrange(run)

pearson_table

## print the table to word
pearson_table_word <- pearson_table %>%
  mutate(
    pearson_r = round(pearson_r, 3),
    t_value   = round(t_value, 2),
    r2        = round(r2, 3),
    p_value   = ifelse(p_value < 0.001, "< 0.001",
                       sprintf("%.3f", p_value))
  ) %>%
  rename(
    'Granite sample' = run,
    `Pearson r` = pearson_r,
    df = df,
    t = t_value,
    `R²` = r2,
    p = p_value
  )

ft <- flextable(pearson_table_word)
ft <- ft %>% 
  theme_booktabs()  %>% 
  autofit() %>% 
  bold(part = "header") %>% 
  align(align = "center", part = "all") %>% 
  set_caption(
    caption = "Table X. Per-run Pearson correlation between grinding intensity and MSI."
  )

doc <- read_docx()
doc <- body_add_docx(doc, ft)
# save_as_docx(ft, path = "/Folder_for_Honour_experimentation (UWA 2025-2026)/tabe1.docx")

ann_run <- ts_alpearson_table_wordann_run <- ts_all_df %>%
  group_by(run) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(slope_mean ~ time_min, data = .x)),
    intercept = map_dbl(fit, ~ coef(.x)[1]),
    slope     = map_dbl(fit, ~ coef(.x)[2]),
    r2        = map_dbl(fit, ~ summary(.x)$r.squared),
    p_value   = map_dbl(fit, ~ summary(.x)$coefficients["time_min", "Pr(>|t|)"])
  ) %>%
  mutate(
    label = sprintf("Intercept = %.3f\nSlope = %.5f\nR² = %.3f\np = %.3g",
                    intercept, slope, r2, p_value)
  ) %>%
  ungroup()

ann_run <- ann_run %>%
  mutate(
    x = case_when(
      run %in% c("GRA-6-LO", "GRA-8-LO") ~ 105,  # move further right
      TRUE ~ 5                                # left group stays near left
    ),
    y = case_when(
      run %in% c("GRA-6-LO", "GRA-8-LO") ~ Inf,  # top for those two
      TRUE ~ -Inf                          # bottom for the rest
    ),
    hjust = case_when(
      run %in% c("GRA-6-LO", "GRA-8-LO") ~ 1.05, # push text slightly inward from right
      TRUE ~ 0
    ),
    vjust = case_when(
      run %in% c("GRA-6-LO", "GRA-8-LO") ~ 1.3,  # push downward from top edge
      TRUE ~ -0.2                              # push further downward at bottom
    )
  )


# Plot all runs using scatter plot.
ggplot(ts_all_df, aes(time_min, slope_mean, colour = material, shape = contact_area)) +
  geom_point(size = 2.5) +
  theme(legend.position = "none") +
  facet_wrap(~ run, nrow = 3, scales = "free_y") +
  geom_smooth(method = "lm", colour = "dodgerblue3", se = FALSE, linetype = "dashed",
              size = 0.8) +
  scale_colour_manual(values = c(
    "Dolerite"  = "#f97b4f",
    "Quartzite" = "#fcbc66"
  )) +
  labs(x = "Time (min)",
       y = "Mean slope value for surface characteristics",
       colour = "Material of active-stone") +
  scale_shape_discrete(name = "Contact Area of active-stone") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(vjust = 3)) + 
  scale_x_continuous(breaks = seq(0,120, 20)) + 
  geom_text(
    data = ann_run,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    size = 2.5
  )

# ------------------------------
# PHASE 4: correlation of intercept and slope
# ------------------------------
coef_df <- ts_all_df %>%
  split(.$run) %>%
  map(~ lm(slope_mean ~ time_min, data = .x)) %>%
  map(coef) %>%
  do.call(rbind, .) %>%
  as.data.frame()

names(coef_df) <- c("Intercept","Slope")

ggplot(coef_df, aes(Intercept, Slope, colour = row.names(coefficients))) +
  geom_point() +
  scale_color_discrete(name = "Granite sample") +
  geom_smooth(method = "lm", colour = "black", se = FALSE) +
  theme(axis.title.y = element_text(vjust = 3),
        legend.title = element_text(size = 8, face = 2),
        legend.text  = element_text(size = 6)) +
  theme_bw() +
  labs(title = "Relationship between intercept and slope (x = strokes)",
       fill = "run")

cor.test(coef_df$Intercept, coef_df$Slope)



# ------------------------------
# PHASE 5: Fitting a LMM 
# ------------------------------

# Mixed linear model: random intercept + random slope 
## using the maximum to minimum process  
# First we fitted the most complex model  
options(contrasts = c("contr.sum", "contr.poly"))
m_ri <- lmer(slope_mean ~ time_min*contact_area*material + (1 | run)+ (0 + time_min| run), data = ts_all_df) 
summary(m_ri)

# Make a table for reporting the model

fe <- tidy(m_ri, effects = "fixed", conf.int = T, conf.method="Wald") %>% 
  transmute(
    Predictor = term,
    Estimate  = estimate,
    SE        = std.error,
    df        = df,
    t         = statistic,
    p         = p.value,
    CI_low    = conf.low,
    CI_high   = conf.high
  ) %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE       = round(SE, 3),
    df       = round(df, 2),
    t        = round(t, 2),
    CI_low   = round(CI_low, 3),
    CI_high  = round(CI_high, 3),
    CI       = paste0(CI_low, " – ", CI_high),
    p        = ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
  ) %>%
  select(Predictor, Estimate, SE, df, t, CI, p)

fe <- fe %>%
  mutate(Predictor = recode(Predictor,
                            "(Intercept)" = "Intercept",
                            "time_min" = "Grinding intensity (min)",
                            "contact_area1" = "Contact area",
                            "material1" = "Material",
                            "time_min:contact_area1" = "Grinding intensity (min) × Contact area",
                            "time_min:material1" = "Grinding intensity (min) × Material",
                            "contact_area1:material1" = "Contact area × Material",
                            "time_min:contact_area1:material1" = "Grinding intensity (min) × Contact area × Material"
  ))

# R2
r2 <- MuMIn::r.squaredGLMM(m_ri)
R2m <- round(r2[1, "R2m"], 3)
R2c <- round(r2[1, "R2c"], 3)

# Random effects SDs
vc <- as.data.frame(VarCorr(m_ri))
re_sd <- vc %>%
  transmute(
    Group = grp,
    Effect = var1,
    SD = round(sdcor, 3)
  )

# sample sizes
n_obs <- nobs(m_ri)
n_runs <- length(unique(ts_all_df$run))

# Make a flextable object
ft <- flextable(fe) %>%
  theme_booktabs() %>%
  autofit() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 11, part = "all") %>%
  set_caption("Table X. Linear mixed-effects model predicting mean surface slope from time, contact area, and material.")

# Add a short model summary under the table
summary_lines <- c(
  paste0("Observations: ", n_obs),
  paste0("Runs (random effect): ", n_runs),
  paste0("Random effects SD (Intercept | run): ", re_sd$SD[re_sd$Group == "run" & re_sd$Effect == "(Intercept)"]),
  paste0("Random effects SD (Time | run): ", re_sd$SD[re_sd$Group == "run.1" & re_sd$Effect == "time_min"]),
  paste0("Residual SD: ", re_sd$SD[re_sd$Group == "Residual"]),
  paste0("Marginal R²: ", R2m, " ; Conditional R²: ", R2c)
)
doc <- read_docx()
doc <- body_add_par(doc, "Mixed-model results", style = "heading 1")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "Model summary:", style = "Normal")
for (ln in summary_lines) {
  doc <- body_add_par(doc, ln, style = "Normal")
}

# print(doc, target = "m_ri_model_table.docx")

# now we can see how this model predict 
pred_grid <- ts_all_df %>%
  distinct(run, material, contact_area) %>%
  tidyr::crossing(
    time_min = seq(min(ts_all_df$time_min),
                   max(ts_all_df$time_min),
                   by = 1)
  )

# Get mixed-model predictions
pred_grid$fit <- predict(m_ri, newdata = pred_grid, re.form = NULL)

lm_pred <- ts_all_df %>%
  group_by(run) %>%
  summarise(
    model = list(lm(slope_mean ~ time_min, data = cur_data())),
    .groups = "drop"
  ) %>%
  tidyr::crossing(
    time_min = seq(min(ts_all_df$time_min),
                   max(ts_all_df$time_min),
                   by = 1)
  ) %>%
  mutate(
    fit_lm = map2_dbl(model, time_min,
                      ~ predict(.x, newdata = data.frame(time_min = .y)))
  )

ggplot(ts_all_df, aes(time_min, slope_mean,
                      colour = material,
                      shape  = contact_area)) +
  
  geom_point(size = 2.5) +
  
  facet_wrap(~ run, nrow = 3, scales = "free_y") +
  
  # Mixed model conditional fits
  geom_line(data = pred_grid,
            aes(y = fit, colour = material),
            linewidth = 0.9) +
  
  # ---- NEW: simple LM per run ----
geom_line(data = lm_pred,
          aes(x = time_min, y = fit_lm),
          colour = "blue",
          linetype = "dashed",
          linewidth = 0.8,
          inherit.aes = FALSE) +
  
  scale_colour_manual(values = c(
    "Dolerite"  = "#f97b4f",
    "Quartzite" = "#fcbc66"
  )) +
  
  labs(x = "Time (min)",
       y = "Mean slope value for surface characteristics",
       colour = "Material of active-stone") +
  
  scale_shape_discrete(name = "Contact Area of active-stone") +
  
  theme_bw() +
  
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(vjust = 3)
  ) +
  
  scale_x_continuous(breaks = seq(0,120,20))


## Now we reduce the model 
m1 <- lmer(slope_mean ~ time_min*as.factor(contact_area) + (1 | run) + (0 + time_min | run),
           data = ts_all_df)
anova(m1)
summary(m1)

# Make a table for this model 

# ---- fixed effects table ----
fe_m1 <- broom.mixed::tidy(
  m1,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
) %>%
  transmute(
    Predictor = term,
    Estimate  = estimate,
    SE        = std.error,
    df        = df,
    t         = statistic,
    p         = p.value,
    CI_low    = conf.low,
    CI_high   = conf.high
  ) %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE       = round(SE, 3),
    df       = round(df, 2),
    t        = round(t, 2),
    CI_low   = round(CI_low, 3),
    CI_high  = round(CI_high, 3),
    CI       = paste0(CI_low, " – ", CI_high),
    p        = ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
  ) %>%
  select(Predictor, Estimate, SE, df, t, CI, p)

# rename
fe_m1 <- fe_m1 %>%
  mutate(Predictor = recode(Predictor,
                            "(Intercept)" = "Intercept",
                            "time_min" = "Grinding intensity (min)",
                            "as.factor(contact_area)1" = "Contact area",
                            "time_min:as.factor(contact_area)1" = "Grinding intensity × Contact area"
  ))

# ---- R2 (marginal/conditional) ----
r2_m1 <- MuMIn::r.squaredGLMM(m1)
R2m <- round(r2_m1[1, "R2m"], 3)
R2c <- round(r2_m1[1, "R2c"], 3)

# ---- random effects SDs + model sizes ----
vc <- as.data.frame(VarCorr(m1))
re_sd <- vc %>%
  transmute(
    Group  = grp,
    Effect = var1,
    SD     = round(sdcor, 3)
  )

n_obs  <- nobs(m1)
n_runs <- length(unique(ts_all_df$run))

# ---- build flextable ----
ft_m1 <- flextable(fe_m1) %>%
  theme_booktabs() %>%
  autofit() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 11, part = "all") %>%
  set_caption("Table X. Linear mixed-effects model (m1) predicting mean surface slope from time and contact area.")

# ---- model summary lines ----
# Note: group names may appear as 'run' and 'run.1' depending on VarCorr() output
sd_int <- re_sd$SD[re_sd$Group == "run"   & re_sd$Effect == "(Intercept)"]
sd_slo <- re_sd$SD[re_sd$Effect == "time_min"][1]        # works even if group label differs
sd_res <- re_sd$SD[re_sd$Group == "Residual"]

summary_lines <- c(
  paste0("Observations: ", n_obs),
  paste0("Runs (random effect): ", n_runs),
  paste0("Random intercept SD (run): ", sd_int),
  paste0("Random slope SD (time): ", sd_slo),
  paste0("Residual SD: ", sd_res),
  paste0("Marginal R²: ", R2m, " ; Conditional R²: ", R2c)
)

# ---- export to Word ----
doc <- read_docx()
doc <- body_add_par(doc, "Mixed-model results (m1)", style = "heading 1")
doc <- body_add_flextable(doc, ft_m1)
doc <- body_add_par(doc, "Model summary:", style = "Normal")
for (ln in summary_lines) {
  doc <- body_add_par(doc, ln, style = "Normal")
}

print(doc, target = "m1_model_table.docx")


# ------------------------------
# PHASE 4: model diagnostics
# ------------------------------

# Based on the warning from the lmer() function there might such a small deviation in the random effect of slope so we need to test of a simpler model can fit the data just as good 
m2 <- lmer(slope_mean ~ time_min*as.factor(contact_area) + (1 | run),
           data = ts_all_df) 

anova(m1,m2)

anova(m_ri, m2)

anova(m1,m_ri)


# m_ri is the best model



### Diagnostics of the best models. 
# Random effects dotplot
dotplot(ranef(m_ri, condVar = TRUE),
        strip = strip.custom(bg = "grey90"),
        par.settings = list(
          plot.symbol = list(col = "black", cex = 0.9),
          plot.line   = list(col = "darkblue", lwd = 2)
        ))


# Diagnostics for m1
fm4 <- fortify.merMod(m_ri)

ggplot(fm4, aes(.fitted, .resid)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  xlab("Fitted values") + ylab("Residuals")

ggplot(fm4, aes(run, .resid)) +
  geom_boxplot() +
  coord_flip() +
  ylab("Residuals")

ggplot(fm4, aes(.fitted, slope_mean)) +
  geom_point(colour = "blue") +
  facet_wrap(~ run, nrow = 3, scales = "free_y") +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Fitted values") + ylab("Observed values")

ggplot(fm4, aes(sample = .resid)) +
  geom_qq() + geom_qq_line()

ggplot(ranef(m1)$run, aes(sample = `(Intercept)`)) +
  geom_qq() + geom_qq_line()

ggplot(ranef(m1)$run, aes(sample = time_min)) +
  geom_qq() + geom_qq_line()


rs_m1 <- data.frame(resid_m1 = fm4$.resid)
ggplot(rs_m1, aes(x = resid_m1)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15,
                 fill = "grey70",
                 color = "black") +
  stat_function(fun = dnorm,
                args = list(mean = mean(rs_m1$resid_m1),
                            sd   = sd(rs_m1$resid_m1)),
                color = "red", linewidth = 1) +
  labs(title = "Histogram of model residuals",
       x = "Residuals", y = "Density") +
  theme_minimal() 


