library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(caret)
library(randomForest)
library(broom)
library(ggplot2)
library(cluster)
library(factoextra)
library(dplyr)
library(readr)   # for parse_number()

run_metric_regression <- function(data, metric_col) {
  formula <- as.formula(
    paste0(metric_col, " ~ angle + v0 + spin_rpm + hb + vb")
  )
  
  model <- lm(formula, data = data)
  s <- summary(model)
  
  tibble(
    term      = rownames(s$coefficients),
    estimate  = s$coefficients[,1],
    std_error = s$coefficients[,2],
    t_value   = s$coefficients[,3],
    p_value   = s$coefficients[,4],
    r_squared = s$r.squared
  )
}

simulate_pitch_nathan_R <- function(
    v0        = 95,   # mph
    theta     = -3,   # deg
    phi       = 0,    # deg
    spin      = c(back = 2300, side = 0, gyro = 0),  # rpm
    position0 = c(x = 0, y = 55, z = 6),             # ft
    dt        = 0.002,
    t_max     = 1.0
) {
  ## -------- helpers ----------
  DEG2RAD    <- function(d) d * pi / 180
  RPM_TO_RAD <- function(rpm) rpm * 2 * pi / 60
  MPH_TO_FTS <- function(mph) mph * 5280 / 3600
  
  ## -------- sanitize inputs (ignore names internally) ----------
  position0 <- as.numeric(position0)
  if (length(position0) != 3) {
    stop("position0 must have 3 elements (x, y, z)")
  }
  # pos indices: 1=x, 2=y, 3=z
  
  spin <- as.numeric(spin)
  if (length(spin) != 3) {
    stop("spin must have 3 elements (back, side, gyro) in that order")
  }
  back_rpm <- spin[1]
  side_rpm <- spin[2]
  gyro_rpm <- spin[3]
  
  ## -------- constants (from JS) ----------
  cd0       <- 0.3008
  cdspin    <- 0.0292         # Cd = cd0 + cdspin * (spinRPM/1000)
  a0        <- 0.6
  a1        <- 1.8
  rho       <- 0.0023769
  g         <- 32.174
  mass_oz   <- 5.125
  radius_in <- 1.45
  
  lb_per_slug <- 32.174
  mass_lb     <- mass_oz / 16
  m           <- mass_lb / lb_per_slug     # slugs
  R           <- radius_in / 12            # ft
  A           <- pi * R * R                # ft^2
  
  ## -------- initial velocity (world coordinates) ----------
  v  <- MPH_TO_FTS(v0)
  th <- DEG2RAD(theta)
  ph <- DEG2RAD(phi)
  
  vx0 <-  v * cos(th) * sin(ph)
  vy0 <-  v * cos(th) * cos(ph)   # matches JS (note: + not -)
  vz0 <-  v * sin(th)
  
  # vel indices: 1=x, 2=y, 3=z
  vel <- c(vx0, vy0, vz0)
  
  ## -------- internal position (y = -y_world) ----------
  pos <- c(
    position0[1],   # x
    -position0[2],  # internal y
    position0[3]    # z
  )
  
  ## -------- spin basis from vhat ----------
  vmag0 <- sqrt(sum(vel^2))
  vhat  <- vel / vmag0
  
  # zhat = world up = (0,0,1)
  zhat <- c(0, 0, 1)
  
  # back axis = vhat × zhat
  back <- c(
    vhat[2] * zhat[3] - vhat[3] * zhat[2],
    vhat[3] * zhat[1] - vhat[1] * zhat[3],
    vhat[1] * zhat[2] - vhat[2] * zhat[1]
  )
  back_mag <- sqrt(sum(back^2))
  if (back_mag < 1e-12) back <- c(1, 0, 0) else back <- back / back_mag
  
  # side axis = projection of zhat onto plane ⟂ vhat
  dot_z_v <- sum(zhat * vhat)
  side <- zhat - dot_z_v * vhat
  side_mag <- sqrt(sum(side^2))
  if (side_mag < 1e-12) side <- c(0, 0, 1) else side <- side / side_mag
  
  gyro <- vhat
  
  # omega_rpm = back*back_axis + side*side_axis + gyro*gyro_axis
  omega_rpm <- back_rpm * back + side_rpm * side + gyro_rpm * gyro
  omega     <- RPM_TO_RAD(omega_rpm)
  totalSpinRPM <- sqrt(sum(omega_rpm^2))
  
  ## -------- storage ----------
  t <- 0
  
  # world y is -pos[2]
  out <- data.frame(
    time = 0,
    x    = pos[1],
    y    = -pos[2],
    z    = pos[3]
  )
  
  y_world_at <- function(pos_vec) {
    -pos_vec[2] + 1.5   # same as JS: -pos.y + 1.5
  }
  
  ## -------- integrate until plate or t_max ----------
  while (t < t_max && y_world_at(pos) > 17/12) {  # 17/12 ft ≈ front of plate
    vmag <- sqrt(sum(vel^2))
    if (vmag < 1e-6) break
    
    vhat_i <- vel / vmag
    
    # omega_perp magnitude
    dot_ov <- sum(omega * vhat_i)
    omega_perp <- sqrt(max(sum(omega^2) - dot_ov^2, 0))
    
    # spin parameter S
    S <- (R * omega_perp) / vmag
    
    # Cd, Cl
    Cd <- cd0 + cdspin * (totalSpinRPM / 1000)
    Cl <- (a0 * S) / (1 + a1 * S)
    
    q  <- 0.5 * rho * vmag^2
    Fd <- q * Cd * A
    Fl <- q * Cl * A
    
    # lift direction ∥ omega × vhat
    cross <- c(
      omega[2] * vhat_i[3] - omega[3] * vhat_i[2],
      omega[3] * vhat_i[1] - omega[1] * vhat_i[3],
      omega[1] * vhat_i[2] - omega[2] * vhat_i[1]
    )
    cmag <- sqrt(sum(cross^2))
    if (cmag > 1e-12) {
      lift_dir <- cross / cmag
    } else {
      lift_dir <- c(0, 0, 0)
    }
    
    # accelerations
    a_drag <- (-Fd * vhat_i) / m
    a_lift <- ( Fl * lift_dir) / m
    a_grav <- c(0, 0, -g)
    
    a   <- a_drag + a_lift + a_grav
    vel <- vel + a * dt
    pos <- pos + vel * dt
    t   <- t + dt
    
    out <- rbind(
      out,
      data.frame(
        time = t,
        x    = pos[1],
        y    = -pos[2],
        z    = pos[3]
      )
    )
  }
  
  rownames(out) <- NULL
  out
}

compute_break_metrics_R <- function(traj_spin, traj_nospin = NULL) {
  if (is.null(traj_spin) || nrow(traj_spin) < 2) {
    return(list(horzBreak = 0, vertBreak = 0))
  }
  
  last  <- tail(traj_spin, 1)
  first <- head(traj_spin, 1)
  
  if (!is.null(traj_nospin) && nrow(traj_nospin) > 1) {
    last_ns <- tail(traj_nospin, 1)
    return(list(
      horzBreak = (last$x - last_ns$x) * 12,  # inches
      vertBreak = (last$z - last_ns$z) * 12   # inches
    ))
  }
  
  # fallback: displacement from straight line
  list(
    horzBreak = (last$x - first$x) * 12,
    vertBreak = (last$z - first$z) * 12
  )
}


calc_metrics = function(v0=95,theta=-1,phi=0,spin=c(back = 3000, side = 0, gyro = 0)){
  
  traj_spin <- simulate_pitch_nathan_R(
    v0        = v0,
    theta     = theta,
    phi       = 0,
    spin      = spin,
    position0 = c(x = 0, y = 55, z = 6),
    dt        = 0.002,
    t_max     = 1.0
  )
  
  traj_nospin <- simulate_pitch_nathan_R(
    v0        = v0,
    theta     = theta,
    phi       = phi,
    spin      = c(back = 0, side = 0, gyro = 0),
    position0 = c(x = 0, y = 55, z = 6),
    dt        = 0.002,
    t_max     = 1.0
  )
  
  breaks <- compute_break_metrics_R(traj_spin, traj_nospin)
  
  return(c(breaks$horzBreak,breaks$vertBreak))
  
}

# Helper: JS tag rule (p3, m3, p500, m500)
fmtTag <- function(x) {
  if (x >= 0) paste0("p", x) else paste0("m", abs(x))
}

# -----------------------------
# Base pitch constants (JS-matching)
# -----------------------------

# -----------------------------
# Base parameters (as you had)
# -----------------------------
baseV0 <- c(
  fb = 96,
  cb = 86,
  sl = 88,
  ch = 87
)

# base spin vectors (JS values: x=gyro, y=back, z=side)
baseSpin <- list(
  fb = c(0, 3000, 0),
  cb = c(0, -3000, 0),
  sl = c(600, -700, -2800),
  ch = c(400, 1800, 700)
)

spinMag <- function(v) sqrt(sum(v * v))
baseSpinMag <- sapply(baseSpin, spinMag)

# helper: format +/- tags like JS ("p3", "m500", etc.)
fmtTag <- function(d) ifelse(d >= 0, paste0("p", d), paste0("m", abs(d)))

# -----------------------------
# Variation grids (JS-matching)
# -----------------------------
pitchTypes  <- c("fb", "cb", "sl", "ch")
angles      <- c(0, 45, 90)
speedDelta  <- c(-3, 0, 3)
spinDelta   <- c(-500, 0, 500)

# -----------------------------
# Put base spin into a tibble
# -----------------------------
spin_df <- tibble(
  pitchType    = names(baseSpin),
  baseSpinX    = sapply(baseSpin, function(v) v[1]),  # gyro
  baseSpinY    = sapply(baseSpin, function(v) v[2]),  # back
  baseSpinZ    = sapply(baseSpin, function(v) v[3]),  # side
  baseSpinMag  = baseSpinMag[names(baseSpin)]
)

# -----------------------------
# Build presets with full spin breakdown
# -----------------------------
df_presets <- expand_grid(
  pitchType = pitchTypes,
  angleDeg  = angles,
  dSpeed    = speedDelta,
  dSpin     = spinDelta
) %>%
  left_join(spin_df, by = "pitchType") %>%
  mutate(
    # actual velocity
    v0 = baseV0[pitchType] + dSpeed,
    
    # target total spin magnitude
    spin_rpm = baseSpinMag + dSpin,
    
    # scale factor for [gyro, back, side] vector
    spin_scale = spin_rpm / baseSpinMag,
    
    gyro_rpm = baseSpinX * spin_scale,
    back_rpm = baseSpinY * spin_scale,
    side_rpm = baseSpinZ * spin_scale,
    
    speedTag = fmtTag(dSpeed),
    spinTag  = fmtTag(dSpin),
    
    # JS ID format: `${pitchType}_${angleDeg}_${speedTag}_${spinTag}`
    id         = paste(pitchType, angleDeg, speedTag, spinTag, sep = "_"),
    labelShort = id,
    label = paste0(
      toupper(pitchType), " – ",
      angleDeg, "° – dSpeed=", dSpeed, " – dSpin=", dSpin,
      " (v0=", v0, "mph, spin=", round(spin_rpm), "rpm, ",
      "gyro=", round(gyro_rpm), ", back=", round(back_rpm),
      ", side=", round(side_rpm), ")"
    )
  ) %>%
  select(
    pitchType, angleDeg, dSpeed, dSpin,
    v0, spin_rpm,
    gyro_rpm, back_rpm, side_rpm,
    speedTag, spinTag,
    id, labelShort, label
  )

hb = list()
vb = list()
for(i in 1:nrow(df_presets)){
  
  gyro = df_presets$gyro_rpm[i]
  back = df_presets$back_rpm[i]
  side = df_presets$side_rpm[i]
  v0 = df_presets$v0[i]
  
  breaks = calc_metrics(v0,theta=-1,phi = 0, spin=c(back=back,side=side,gyro=gyro))
  hb[[i]] = breaks[1]
  vb[[i]] = breaks[2]
}

hbs = unlist(hb)
vbs = unlist(vb)

df_presets$hb = hbs
df_presets$vb = vbs

# Starting point: df and df_presets already in memory
df = read.csv('./Documents/github/allPitchCones.csv') 
df = df %>%
  group_by(pitch,spin,angle,speed,frame) %>%
  arrange(run) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  filter(row == 1)


df_mapped <- df %>%
  mutate(
    # match JS pitchType (fb/cb/sl/ch)
    pitchType = as.character(pitch),
    
    # "0deg" -> 0, "45deg" -> 45, etc.
    angleDeg  = parse_number(angle),
    
    # "speed+3" -> 3, "speed-3" -> -3, "speed0" -> 0
    dSpeed    = parse_number(speed),
    
    # "spin+500" -> 500, "spin-500" -> -500, "spin0" -> 0
    dSpin     = parse_number(spin)
  ) %>%
  # join on the 4 defining knobs
  left_join(
    df_presets,
    by = c("pitchType", "angleDeg", "dSpeed", "dSpin")
  ) %>%
  filter(pitchType != 'custom')

N_FRAMES <- 80  

df_mapped <- df_mapped %>%
  group_by(id) %>%  # id uniquely defines a preset (pitchType/angle/dSpeed/dSpin)
  mutate(
    # normalized time from 0 (release) to 1 (plate)
    t_norm = (frame - min(frame)) / (max(frame) - min(frame)),
    
    # distance from home plate in feet (linear in normalized time)
    dist_from_plate_ft   = 55 * (1 - t_norm),
    dist_from_release_ft = 55 - dist_from_plate_ft
  ) %>%
  ungroup()

#starting analysis
dfm <- df_mapped %>%
  group_by(id) %>%
  mutate(
    t_norm = (frame - min(frame)) / (max(frame) - min(frame))
  ) %>%
  ungroup()


#PCA BULLSHIT
df_pca <- df_mapped %>% 
  ungroup() %>% 
  select(fovea_total, mid_total, far_total,
         total, spread, meanSignal, delta_total, corr_prev) %>%
  mutate_all(as.numeric)

# remove any missing
df_clean <- df_pca %>% drop_na()

df_nz <- df_clean[, sapply(df_clean, function(x) sd(x, na.rm=TRUE) > 1e-12)]

# Run global PCA on ALL frames
global_pca <- prcomp(df_nz, center = TRUE, scale. = TRUE)
global_scores <- global_pca$x


df_master <- df_mapped %>%
  drop_na(fovea_total, mid_total, far_total,
          total, spread, meanSignal, delta_total, corr_prev) %>%
  mutate(
    PC1 = global_scores[,1],
    PC2 = global_scores[,2],
    PC3 = global_scores[,3],
    PC4 = global_scores[,4],
    PC5 = global_scores[,5]
  )

#Measuring seam orientation effects
pitch_types <- unique(df_master$pitch)

run_pc_regression <- function(data, pc_col) {
  
  formula <- as.formula(
    paste0(pc_col, " ~ v0 + spin_rpm + angle + hb + vb")
  )
  
  model <- lm(formula, data = data)
  s <- summary(model)
  
  # s$coefficients is a matrix with columns:
  # Estimate, Std. Error, t value, Pr(>|t|)
  tibble(
    term     = rownames(s$coefficients),
    estimate = s$coefficients[, 1],
    std_error = s$coefficients[, 2],
    t_value   = s$coefficients[, 3],
    p_value   = s$coefficients[, 4],
    r_squared = s$r.squared
  )
}


df_master <- df_master %>%
  mutate(
    dist_bin = cut(
      dist_from_release_ft,
      breaks = seq(60, 0, by = -5),
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  # convert bin labels → actual numeric bin midpoints
  mutate(
    dist_from_plate = {
      # turn factor levels like "[0,5)" into a numeric matrix
      mm <- read.table(
        text = gsub("\\[|\\)|\\]", "", as.character(dist_bin)),
        sep = ",",
        col.names = c("lo", "hi")
      )
      (mm$lo + mm$hi) / 2   # midpoint
    }
  ) %>%
  mutate(dist_from_plate = dist_from_plate_ft)


#metric regressions
metrics <- c("fovea_total", "mid_total", "far_total",
             "meanSignal")

results_metric_seam <- df_master %>%
  group_by(pitch, dist_bin) %>%
  group_modify(~ {
    data_bin <- .x
    map_dfr(metrics, function(mname) {
      run_metric_regression(data_bin, metric_col = mname) %>%
        mutate(metric = mname)
    })
  }) %>%
  ungroup()

seam_effects <- results_metric_seam %>%
  filter(term %in% c("angle45deg", "angle90deg"))

seam_effects <- seam_effects %>%
  mutate(
    dist_low  = as.numeric(stringr::str_extract(dist_bin, "(?<=\\[)\\d+")),
    dist_high = as.numeric(stringr::str_extract(dist_bin, "(?<=,)\\d+")),
    dist_mid  = (dist_low + dist_high)/2
  )

seam_effects %>%
  filter(metric == "fovea_total") %>%
  ggplot(aes(dist_mid, estimate, color = term)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(aes(shape = p_value < 0.05)) +
  facet_wrap(~ pitch) +
  scale_x_continuous(breaks = seq(0, 60, 5)) +  # <-- forward direction
  labs(
    x = "Distance from Release (ft)",
    y = "Effect of Seam Angle on Foveal Activity",
    color = "Seam Angle",
    shape = "p < 0.05"
  ) +
  theme_bw()


#Regression per pitch type
pcs <- paste0("PC", 1:5)

results_pitch_seam <- df_master %>%
  group_by(pitch, dist_bin) %>%
  group_modify(~ {
    data_bin <- .x
    
    map_dfr(pcs, function(pcname) {
      run_pc_regression(data_bin, pc_col = pcname) %>%
        mutate(PC = pcname)
    })
  }) %>%
  ungroup()

seam_effects <- seam_effects %>%
  mutate(
    dist_low  = as.numeric(stringr::str_extract(dist_bin, "(?<=\\[)\\d+")),
    dist_high = as.numeric(stringr::str_extract(dist_bin, "(?<=,)\\d+")),
    dist_mid  = (dist_low + dist_high)/2,
    dist_from_plate = 60 - dist_mid     # <-- NEW METRIC
  )

seam_effects %>%
  filter(metric == "fovea_total") %>%
  ggplot(aes(dist_from_plate, estimate, color = term)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(aes(shape = p_value < 0.05)) +
  facet_wrap(~ pitch) +
  scale_x_reverse(breaks = seq(0, 60, by = 5)) +  # <-- CLEAN, CORRECT
  labs(
    x = "Distance From Home Plate (ft)",
    y = "Effect of Seam Angle on Foveal Activity",
    color = "Seam Angle",
    shape = "p < 0.05"
  ) +
  theme_bw()


library(gt)

summary_table_clean <- seam_effects %>%
  filter(metric == "fovea_total") %>%
  group_by(pitch, term) %>%
  summarize(
    peak_effect   = max(abs(estimate), na.rm = TRUE),
    significant   = if_else(any(p_value < 0.05), "Yes", "No"),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = term,
    values_from = c(peak_effect, significant)
  )

# --- Professional formatting ---
gt(summary_table_clean) %>%
  tab_header(
    title = md("**Foveal Effects of Seam Orientation**"),
    subtitle = "Peak magnitude and statistical significance"
  ) %>%
  fmt_number(
    columns = contains("peak_effect"),
    decimals = 0
  ) %>%
  tab_style(
    style = cell_borders(color = "gray85", sides = "all", weight = px(1)),
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "gray95"),
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_options(
    table.border.top.width = px(2),
    table.border.bottom.width = px(2),
    table.font.size = 12,
    data_row.padding = px(4)
  )


df_master <- df_master %>%
  mutate(
    angle = factor(angle, levels = c("0deg", "45deg", "90deg"))
  )

run_interaction_regression <- function(df, metric_col) {
  fml <- as.formula(
    paste0(metric_col, " ~ angle * spin_rpm * v0")
  )
  broom::tidy(lm(fml, data = df))
}

interaction_results <- df_master %>%
  group_by(pitch, dist_bin) %>%
  group_modify(~ run_interaction_regression(.x, "fovea_total")) %>%
  ungroup()


library(stringr)

# prepare distance metric
interaction_plot_df <- interaction_results %>%
  filter(str_detect(term, "angle")) %>%   # keep angle, angle:spin, angle:v0, angle:spin:v0
  mutate(
    dist_low  = as.numeric(str_extract(dist_bin, "(?<=\\[)\\d+")),
    dist_high = as.numeric(str_extract(dist_bin, "(?<=,)\\d+")),
    dist_mid  = (dist_low + dist_high)/2,
    dist_from_plate = 60 - dist_mid,
    effect_mag = abs(estimate)
  )

ggplot(interaction_plot_df,
       aes(x = dist_from_plate,
           y = term,
           fill = effect_mag)) +
  geom_tile(color = "gray80") +
  scale_x_reverse(breaks = seq(0, 60, by = 5)) +
  scale_fill_viridis_c(option = "plasma") +
  facet_wrap(~ pitch, ncol = 2) +
  labs(
    title = "Seam Orientation × Spin Rate × Velocity Interaction Strength",
    x = "Distance From Home Plate (ft)",
    y = "Interaction Term",
    fill = "| β |"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8)
  )


#Pitch Classifier based on retinal activity
# choose retinal features to feed into the classifier
retina_feats <- c(
  "fovea_total", 
  "mid_total", 
  "far_total",
  # "spread",
  # "delta_total",
  "corr_prev"
)

df_master <- df_master %>%
  mutate(
    condition = sub("^[a-z]{2}_", "", id)   # removes cb_, fb_, sl_, ch_
  )

valid_conditions <- df_master %>%
  distinct(pitch, condition) %>%
  count(condition) %>%
  filter(n == 4) %>%        # must appear in fb, cb, sl, ch
  pull(condition)

df_balanced <- df_master %>%
  filter(condition %in% valid_conditions)

valid_triplets <- df_master %>%
  distinct(pitch, condition, dist_bin) %>%
  count(condition, dist_bin) %>%
  filter(n == 4) %>%     # must appear for fb, cb, sl, ch
  dplyr::select(condition, dist_bin)

df_balanced <- df_master %>%
  inner_join(valid_triplets, by = c("condition", "dist_bin"))

retina_feats <- c("fovea_total", "mid_total", "far_total", "corr_prev","dist_from_plate")

df_runlevel <- df_balanced %>%
  ungroup() %>%
  filter(angle == '0deg') %>%
  # mutate(dist_from_plate = dist_from_plate_ft) %>%
  group_by(pitch, condition, dist_bin, dist_from_plate) %>%
  dplyr::summarise(
    across(all_of(c("fovea_total", "mid_total", "far_total", "corr_prev", "hb","vb","v0","spin_rpm")), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pitch = factor(pitch),
    dist_bin = factor(dist_bin)
  )



# Sort by condition + decreasing distance
ggplot(df, aes(x = dist_from_plate, y = corr_prev, color = pitch)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 0.7, alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Distance From Plate (ft)",
    y = "Fovea Total",
    color = "Pitch",
    title = "Fovea Total vs Distance From Plate"
  )

df <- df_runlevel %>%
  arrange(condition, pitch, desc(dist_from_plate))

df <- df %>%
  group_by(condition, pitch) %>%   # <-- KEY FIX
  mutate(
    trial_break = dist_from_plate > lag(dist_from_plate, default = first(dist_from_plate)),
    trial_id = cumsum(trial_break)
  ) %>%
  ungroup()

# Create global trial UID
df$trial_uid <- paste(df$condition, df$pitch, df$trial_id, sep = "_")

set.seed(42)

trials <- unique(df$trial_uid)

test_trials  <- sample(trials, size = floor(length(trials)*0.30))
train_trials <- setdiff(trials, test_trials)

train_df <- df %>% filter(trial_uid %in% train_trials)
test_df  <- df %>% filter(trial_uid %in% test_trials)

trials_df <- df %>%
  distinct(trial_uid, pitch)

set.seed(42)

test_trials <- trials_df %>%
  group_by(pitch) %>%
  sample_frac(0.30) %>%
  pull(trial_uid)

train_trials <- setdiff(trials_df$trial_uid, test_trials)

df2 <- df %>%
  mutate(
    dist_bin = cut(
      dist_from_plate,
      breaks = seq(60, 0, by = -5),
      include.lowest = TRUE,
      right = FALSE
    )
  )

# Train/test split by true trial_id (from previous code)
train_df2 <- df2 %>% filter(trial_uid %in% train_trials)
test_df2  <- df2 %>% filter(trial_uid %in% test_trials)

# RF (bucketed distance instead of numeric)
rf_bin <- randomForest(
  x = train_df2[, c("dist_bin", "fovea_total","corr_prev")],
  y = as.factor(train_df2$pitch),
  ntree = 500
)

test_df2$pred <- predict(rf_bin, test_df2[, c("dist_bin", "fovea_total","corr_prev")])
test_df2$prob <- predict(rf_bin, test_df2[, c("dist_bin", "fovea_total","corr_prev")], type = "prob")
# Accuracy by distance bin
acc_by_bin2 <- test_df2 %>%
  group_by(dist_bin) %>%
  summarize(
    n = n(),
    acc = mean(pred == pitch)
  )

print(acc_by_bin2)
plot_df <- test_df2 %>%
  group_by(dist_bin, pitch, pred) %>%
  summarize(n = n(), .groups = "drop")

ggplot(plot_df, aes(x = pred, y = pitch, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  facet_wrap(~ dist_bin, ncol = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "RF: Predicted vs Actual Pitch Types Across Distance",
    x = "Predicted",
    y = "Actual",
    fill = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Reverse the factor order (so [55,60) is left, [0,5) is right)
acc_by_bin2 <- acc_by_bin2 %>%
  mutate(dist_bin = factor(dist_bin, levels = rev(levels(dist_bin))))

ggplot(acc_by_bin2, aes(x = dist_bin, y = acc, group = 1)) +
  geom_line(size = 1.3, color = "firebrick") +
  geom_point(size = 3, color = "firebrick") +
  ylim(0.25, 1.0) +
  labs(
    title = "Classification Accuracy vs. Distance",
    x = "Distance Bin from Home Plate (ft)",
    y = "Accuracy"
  ) +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Example: density of predicted probabilities by distance
library(reshape2)

prob_df <- cbind(
  test_df2[, c("dist_bin", "pitch")],
  as.data.frame(test_df2$prob)
)

# Long format: one row per (trial, class)
prob_long <- prob_df %>%
  pivot_longer(
    cols = c(cb, ch, fb, sl),   # adjust if your column names differ
    names_to = "class",
    values_to = "prob"
  )

# Keep ONLY the probability of the *true* class
prob_true <- prob_long %>%
  filter(class == pitch)

# Now plot: confidence in the correct class vs distance
ggplot(prob_true, aes(prob, color = pitch)) +
  geom_density() +
  facet_wrap(~ dist_bin) +
  labs(
    title = "RF Confidence in the TRUE Pitch Class by Distance",
    x = "Predicted Probability of True Class",
    y = "Density"
  )

# prob_true has: dist_bin, pitch, prob  (prob = P(true class))

summary_table <- prob_true %>%
  group_by(dist_bin, pitch) %>%
  summarise(mean_prob = mean(prob), .groups = "drop") %>%
  pivot_wider(
    names_from = pitch,
    values_from = mean_prob
  ) %>%
  arrange((dist_bin))   # optional: put 60ft → 0ft downward

summary_table

library(knitr)

summary_table_pretty <- prob_true %>%
  group_by(dist_bin, pitch) %>%
  summarise(mean_prob = mean(prob), .groups = "drop") %>%
  pivot_wider(
    names_from = pitch,
    values_from = mean_prob
  ) %>%
  arrange((dist_bin)) %>%   # 60ft → 0ft
  mutate(across(cb:sl, ~ scales::percent(.x, accuracy = 0.1)))

num_mat <- summary_table_pretty[, 2:5] %>%
  mutate(across(everything(),
                ~ as.numeric(gsub("%", "", .x)) / 100))

blue_red <- colorRampPalette(c("blue", "red"))

tab <- summary_table_pretty %>%
  kable(
    caption = "Mean RF Confidence (True-Class Probability)",
    align   = "c"
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    font_size  = 14
  ) %>%
  row_spec(0, bold = TRUE, background = "#f0f0f0") %>%
  column_spec(1, bold = TRUE, border_right = TRUE)

# for each pitch column (cb, ch, fb, sl), add its own color vector
pitch_cols <- names(num_mat)

for (j in seq_along(pitch_cols)) {
  vals <- num_mat[[j]]              # numeric 0–1
  # map 0–1 → 100-step blue→red
  idx  <- as.numeric(cut(vals,
                         breaks = seq(0, 1, length.out = 101),
                         include.lowest = TRUE))
  col_vec <- blue_red(100)[idx]
  
  tab <- tab %>%
    column_spec(
      column    = j + 1,            # +1 because col 1 is dist_bin
      background = col_vec,
      color      = "white"
    )
}

tab


#Correlation analysis
cor_df <- prob_true %>%
  select(prob, spin_rpm, v0, hb, vb) %>%
  cor(use = "complete.obs")

cor_df


df_pt <- prob_true %>%
  select(pitch, prob, spin_rpm, v0, hb, vb)

# Function to compute a correlation row for one pitch
cor_for_pitch <- function(df_pitch) {
  cor_mat <- cor(df_pitch[, c("prob","spin_rpm","v0","hb","vb")], use="pairwise")
  # extract correlations with "prob"
  tibble(
    spin_rpm = cor_mat["prob","spin_rpm"],
    v0       = cor_mat["prob","v0"],
    hb       = cor_mat["prob","hb"],
    vb       = cor_mat["prob","vb"]
  )
}

# Run for each pitch type
cor_by_pitch <- df_pt %>%
  group_by(pitch) %>%
  group_modify(~ cor_for_pitch(.x)) %>%
  ungroup()

cor_by_pitch


library(dplyr)
library(kableExtra)

# Your correlation table
# Blue → white → red color map
col_fun <- colorRampPalette(c("blue", "white", "red"))

# ===========================
# 1. CREATE PERCENT TEXT
# ===========================
cor_percent <- cor_by_pitch %>%
  mutate(across(-pitch, ~ sprintf("%.1f%%", .x * 100)))

# ===========================
# 2. COLOR SCALE (numeric)
# ===========================
col_fun <- colorRampPalette(c("blue", "white", "red"))

num_vals <- cor_by_pitch %>%
  select(-pitch) %>%
  as.matrix()

# Map numeric values to colors
color_map <- function(x) {
  if (is.na(x)) return("transparent")     # keep NA clean
  idx <- round(((x + 1) / 2) * 199) + 1   # scale -1→1 into 1→200
  col_fun(200)[idx]
}

color_matrix <- apply(num_vals, c(1,2), color_map)

# ===========================
# 3. BUILD THE TABLE
# ===========================
k <- cor_percent %>%
  kable("html",
        caption = "Correlation Between TRUE-Class RF Probability and Pitch Traits (Percent Scale)",
        align = "c") %>%
  kable_styling(full_width = FALSE, font_size = 14) %>%
  row_spec(0, bold = TRUE, background = "#f0f0f0") %>%
  column_spec(1, bold = TRUE, border_right = TRUE)

# apply colors column by column
for (j in 2:ncol(cor_percent)) {
  k <- k %>% column_spec(j, background = color_matrix[, j-1])
}

k

ggplot(prob_true, aes(spin_rpm, prob, color = pitch)) +
  geom_point(alpha = 0.15) +
  geom_smooth(se=FALSE) +
  facet_wrap(~ pitch) +
  labs(title = "RF Confidence vs Spin Rate", y="True-Class Probability")


