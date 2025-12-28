suppressPackageStartupMessages({
  library(tidyverse)
  library(haven)
  library(survey)
  library(Hmisc)
  library(broom)
})

options(survey.lonely.psu = "adjust")

# ==========================================================
# 0) Choose base directory ONCE
# Pick any *.xpt in your NHANES folder (e.g., DEMO_G.xpt)
# ==========================================================
xpt_dir <- dirname(file.choose())
if (!dir.exists(xpt_dir)) stop("目录不存在：", xpt_dir)
message("Using directory: ", xpt_dir)

# ==========================================================
# Common helpers
# ==========================================================
read_xpt_safe <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  haven::read_xpt(path) %>% as_tibble()
}

read_GH <- function(stem) {
  g_path <- file.path(xpt_dir, paste0(stem, "_G.xpt"))
  h_path <- file.path(xpt_dir, paste0(stem, "_H.xpt"))
  bind_rows(
    read_xpt_safe(g_path) %>% mutate(cycle = "2011-2012"),
    read_xpt_safe(h_path) %>% mutate(cycle = "2013-2014")
  )
}

most_frequent <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  tab <- table(x)
  names(tab)[which.max(tab)]
}

collapse_to_unique <- function(df, df_name = "component") {
  if (!all(c("SEQN", "cycle") %in% names(df))) {
    stop(sprintf("[%s] missing SEQN/cycle columns.", df_name))
  }
  dup_groups <- df %>% count(SEQN, cycle) %>% filter(n > 1)
  if (nrow(dup_groups) == 0) return(df)
  message(sprintf("[%s] Found %d duplicated SEQN+cycle IDs -> collapsing.", df_name, nrow(dup_groups)))
  same_rows <- df %>%
    group_by(SEQN, cycle) %>%
    summarise(.same = n_distinct(across(everything())) == 1, .groups = "drop")
  if (all(same_rows$.same)) {
    message(sprintf("[%s] Duplicates are identical rows -> using distinct().", df_name))
    return(df %>% distinct(SEQN, cycle, .keep_all = TRUE))
  }
  df %>%
    group_by(SEQN, cycle) %>%
    summarise(
      across(
        .cols = -all_of(c("SEQN", "cycle")),
        .fns = ~{
          if (is.numeric(.x)) {
            if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
          } else {
            most_frequent(.x)
          }
        }
      ),
      .groups = "drop"
    )
}

# ==========================================================
# STEP1) Read + merge NHANES components (2011–2014), make WTMEC4YR
# ==========================================================
message("\n===== STEP1 =====")
DEMO <- read_GH("DEMO") %>% collapse_to_unique("DEMO")
CFQ <- read_GH("CFQ") %>% collapse_to_unique("CFQ")
DR1 <- read_GH("DR1TOT") %>% collapse_to_unique("DR1TOT")
DR2 <- read_GH("DR2TOT") %>% collapse_to_unique("DR2TOT")
PAQ <- read_GH("PAQ") %>% collapse_to_unique("PAQ")
BMX <- read_GH("BMX") %>% collapse_to_unique("BMX")
ALQ <- read_GH("ALQ") %>% collapse_to_unique("ALQ")
MCQ <- read_GH("MCQ") %>% collapse_to_unique("MCQ")
DPQ <- read_GH("DPQ") %>% collapse_to_unique("DPQ")
BPX <- read_GH("BPX") %>% collapse_to_unique("BPX")
DIQ <- read_GH("DIQ") %>% collapse_to_unique("DIQ")
GLU <- read_GH("GLU") %>% collapse_to_unique("GLU")
GHB <- read_GH("GHB") %>% collapse_to_unique("GHB")
COT <- read_GH("COT") %>% collapse_to_unique("COT")

master <- DEMO %>%
  left_join(CFQ, by = c("SEQN","cycle")) %>%
  left_join(DR1, by = c("SEQN","cycle")) %>%
  left_join(DR2, by = c("SEQN","cycle")) %>%
  left_join(PAQ, by = c("SEQN","cycle")) %>%
  left_join(BMX, by = c("SEQN","cycle")) %>%
  left_join(ALQ, by = c("SEQN","cycle")) %>%
  left_join(MCQ, by = c("SEQN","cycle")) %>%
  left_join(DPQ, by = c("SEQN","cycle")) %>%
  left_join(BPX, by = c("SEQN","cycle")) %>%
  left_join(DIQ, by = c("SEQN","cycle")) %>%
  left_join(GLU, by = c("SEQN","cycle")) %>%
  left_join(GHB, by = c("SEQN","cycle")) %>%
  left_join(COT, by = c("SEQN","cycle"))

dup_master <- master %>% count(SEQN, cycle) %>% filter(n > 1)
if (nrow(dup_master) > 0) {
  print(dup_master)
  stop("master still has duplicated SEQN+cycle after collapsing components.")
} else {
  message("OK: master has unique SEQN+cycle (one row per participant).")
}

master <- master %>% mutate(WTMEC4YR = WTMEC2YR * 0.5)
saveRDS(master, file = file.path(xpt_dir, "master_NHANES_2011_2014_withCOT.rds"))
message("Saved: master_NHANES_2011_2014_withCOT.rds")

# ==========================================================
# STEP2) Build analytic sample + flow
# ==========================================================
message("\n===== STEP2 =====")
df0 <- readRDS(file.path(xpt_dir, "master_NHANES_2011_2014_withCOT.rds"))

flow <- tibble(step = character(), n = integer())
add_flow <- function(label, data) {
  flow <<- bind_rows(flow, tibble(step = label, n = nrow(data)))
  data
}

DR1_RELIABLE_VAR <- "DR1DRSTZ"
if (!DR1_RELIABLE_VAR %in% names(df0)) {
  stop("Required dietary reliability variable not found: ", DR1_RELIABLE_VAR,
       ". Please confirm DR1TOT file was merged.")
}

kcal_candidates_day1 <- c("DR1TKCAL", "DR1IKCAL")
kcal_candidates_day2 <- c("DR2TKCAL", "DR2IKCAL")
kcal1_found <- kcal_candidates_day1[kcal_candidates_day1 %in% names(df0)]
kcal2_found <- kcal_candidates_day2[kcal_candidates_day2 %in% names(df0)]

if (length(kcal1_found) == 0) {
  stop("No Day1 energy intake variable found (expected DR1TKCAL or DR1IKCAL).")
}

df0 <- df0 %>%
  mutate(
    kcal_day1 = .data[[kcal1_found[1]]],
    kcal_day2 = if (length(kcal2_found) > 0) .data[[kcal2_found[1]]] else NA_real_,
    kcal_excl = if_else(!is.na(kcal_day2), (kcal_day1 + kcal_day2) / 2, kcal_day1)
  )

KCAL1_VAR <- "kcal_excl"

df <- df0 %>% add_flow("0) NHANES 2011–2014 merged (STEP1)", .)

df <- df %>%
  filter(!is.na(RIDAGEYR), RIDAGEYR >= 60) %>%
  add_flow("1) Age ≥ 60 years", .)

cog_vars <- c("CFDCST1","CFDCST2","CFDCST3","CFDCSR","CFDAST","CFDDS")
missing_cog <- setdiff(cog_vars, names(df))
if (length(missing_cog) > 0) stop("Cognitive test variables missing: ", paste(missing_cog, collapse = ", "))

df <- df %>%
  mutate(cerad_immediate = CFDCST1 + CFDCST2 + CFDCST3) %>%
  filter(across(all_of(cog_vars), ~ !is.na(.x))) %>%
  add_flow("2) Completed all cognitive tests", .)

df <- df %>%
  filter(.data[[DR1_RELIABLE_VAR]] == 1) %>%
  add_flow("3) Day1 dietary recall reliable (DR1DRSTZ == 1)", .)

df <- df %>%
  filter(!is.na(.data[[KCAL1_VAR]])) %>%
  filter(
    (RIAGENDR == 2 & .data[[KCAL1_VAR]] >= 500 & .data[[KCAL1_VAR]] <= 5000) |
      (RIAGENDR == 1 & .data[[KCAL1_VAR]] >= 500 & .data[[KCAL1_VAR]] <= 8000)
  ) %>%
  add_flow("4) Excluded extreme energy intake", .)

if (!"MCQ160F" %in% names(df)) {
  warning("MCQ160F not found. Stroke exclusion step skipped.")
} else {
  df <- df %>%
    filter(is.na(MCQ160F) | MCQ160F != 1) %>%
    add_flow("5) Excluded self-reported stroke (MCQ160F != 1)", .)
}

key_covars <- c("RIDAGEYR","RIAGENDR","RIDRETH1","DMDEDUC2",
                "INDFMPIR","BMXBMI","LBXCOT",
                "SDMVSTRA","SDMVPSU","WTMEC4YR")
missing_cov <- setdiff(key_covars, names(df))
if (length(missing_cov) > 0) stop("Key covariates missing: ", paste(missing_cov, collapse = ", "))

df_analytic <- df %>%
  filter(across(all_of(key_covars), ~ !is.na(.x))) %>%
  add_flow("6) Complete key covariates", .)

print(flow)
saveRDS(df_analytic, file.path(xpt_dir, "df_analytic_step2.rds"))
write.csv(flow, file.path(xpt_dir, "flow_step2.csv"), row.names = FALSE)
message("Saved: df_analytic_step2.rds, flow_step2.csv")

# ==========================================================
# STEP3.1) Dietary OBS
# ==========================================================
message("\n===== STEP3.1 =====")
df <- readRDS(file.path(xpt_dir, "df_analytic_step2.rds"))

mean_day12 <- function(d1, d2) ifelse(!is.na(d2), (d1 + d2) / 2, d1)

score_tertile <- function(x, reverse = FALSE) {
  q <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE, type = 2)
  s <- case_when(
    is.na(x) ~ NA_real_,
    x <= q[1] ~ 0,
    x <= q[2] ~ 1,
    x > q[2] ~ 2
  )
  if (reverse) s <- 2 - s
  s
}

anti_base <- c("DR1TVC","DR1TATOC","DR1TFIBE","DR1TSELE","DR1TZINC","DR1TMAGN",
               "DR1TCALC","DR1TCOPP","DR1TVB6","DR1TVB12","DR1TFOLA","DR1TVB2","DR1TNIAC")
carotenoid_vars <- c("DR1TACAR","DR1TBCAR","DR1TCRYP","DR1TLYCO","DR1TLZ")
pro_fat1 <- "DR1TTFAT"; pro_fat2 <- "DR2TTFAT"
pro_iron1 <- "DR1TIRON"; pro_iron2 <- "DR2TIRON"

for (v1 in anti_base) {
  v2 <- sub("^DR1", "DR2", v1)
  df[[paste0(v1, "_mean")]] <- if (v2 %in% names(df)) mean_day12(df[[v1]], df[[v2]]) else df[[v1]]
}

caro1 <- rowSums(df[, carotenoid_vars], na.rm = TRUE)
caro2_vars <- sub("^DR1", "DR2", carotenoid_vars)
caro2 <- if (all(caro2_vars %in% names(df))) rowSums(df[, caro2_vars], na.rm = TRUE) else NA_real_
df$DR1TCARO_mean <- mean_day12(caro1, caro2)

df$fat_mean <- mean_day12(df[[pro_fat1]], df[[pro_fat2]])
df$iron_mean <- mean_day12(df[[pro_iron1]], df[[pro_iron2]])

df <- df %>%
  group_by(RIAGENDR) %>%
  mutate(
    across(
      ends_with("_mean") & !any_of(c("fat_mean", "iron_mean")),
      ~ score_tertile(.x, reverse = FALSE),
      .names = "{.col}_score"
    ),
    fat_mean_score = score_tertile(fat_mean, reverse = TRUE),
    iron_mean_score = score_tertile(iron_mean, reverse = TRUE)
  ) %>%
  ungroup()

diet_score_vars <- c(
  paste0(paste0(anti_base, "_mean"), "_score"),
  "DR1TCARO_mean_score","fat_mean_score","iron_mean_score"
)

df <- df %>% mutate(OBS_dietary = rowSums(across(all_of(diet_score_vars)), na.rm = FALSE))
saveRDS(df, file.path(xpt_dir, "df_step3_dietOBS.rds"))
message("Saved: df_step3_dietOBS.rds")

# ==========================================================
# STEP3.2) Lifestyle OBS + Total OBS
# ==========================================================
message("\n===== STEP3.2 =====")
df <- readRDS(file.path(xpt_dir, "df_step3_dietOBS.rds"))

df <- df %>%
  mutate(
    smoke_score = case_when(
      LBXCOT < 0.015 ~ 2,
      LBXCOT < 10 ~ 1,
      LBXCOT >= 10 ~ 0,
      TRUE ~ NA_real_
    )
  )

score_tertile <- function(x, reverse = FALSE) {
  q <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE, type = 2)
  s <- case_when(
    is.na(x) ~ NA_real_,
    x <= q[1] ~ 0,
    x <= q[2] ~ 1,
    x > q[2] ~ 2
  )
  if (reverse) s <- 2 - s
  s
}

df <- df %>%
  group_by(RIAGENDR) %>%
  mutate(bmi_score = score_tertile(BMXBMI, reverse = TRUE)) %>%
  ungroup()

df <- df %>%
  mutate(
    alcohol_score = case_when(
      ALQ130 == 0 ~ 2,
      RIAGENDR == 1 & ALQ130 < 30 ~ 1,
      RIAGENDR == 2 & ALQ130 < 15 ~ 1,
      TRUE ~ 0
    )
  )

to_na_days <- function(x) dplyr::case_when(x %in% c(77, 99) ~ NA_real_, TRUE ~ as.numeric(x))
to_na_mins <- function(x) dplyr::case_when(x %in% c(7777, 9999, 777, 999, 77, 99) ~ NA_real_, TRUE ~ as.numeric(x))

df <- df %>%
  mutate(
    vig_yes = case_when(PAQ650 == 1 ~ 1, PAQ650 == 2 ~ 0, TRUE ~ NA_real_),
    vig_days = to_na_days(PAQ655),
    vig_mins = to_na_mins(PAD660),
    vig_min_week = case_when(
      vig_yes == 0 ~ 0,
      vig_yes == 1 ~ vig_days * vig_mins,
      TRUE ~ NA_real_
    ),
    mod_yes = case_when(PAQ665 == 1 ~ 1, PAQ665 == 2 ~ 0, TRUE ~ NA_real_),
    mod_days = to_na_days(PAQ670),
    mod_mins = to_na_mins(PAD675),
    mod_min_week = case_when(
      mod_yes == 0 ~ 0,
      mod_yes == 1 ~ mod_days * mod_mins,
      TRUE ~ NA_real_
    ),
    PA_MET = mod_min_week * 4.0 + vig_min_week * 8.0
  )

df <- df %>%
  group_by(RIAGENDR) %>%
  mutate(pa_score = score_tertile(PA_MET, reverse = FALSE)) %>%
  ungroup()

df <- df %>%
  mutate(
    OBS_lifestyle = smoke_score + bmi_score + alcohol_score + pa_score,
    OBS_total = OBS_dietary + OBS_lifestyle,
    OBS_z = as.numeric(scale(OBS_total))
  )

saveRDS(df, file.path(xpt_dir, "df_step3_OBS_final.rds"))
df_final <- df %>% filter(!is.na(OBS_total))
saveRDS(df_final, file.path(xpt_dir, "df_analysis_OBS_final.rds"))
message("Saved: df_step3_OBS_final.rds, df_analysis_OBS_final.rds")

# ==========================================================
# STEP4) Survey-weighted models + RCS + tables/figures
# ==========================================================
message("\n===== STEP4 =====")
df <- readRDS(file.path(xpt_dir, "df_step3_OBS_final.rds")) %>% filter(!is.na(OBS_total))

na_if_any <- function(x, miss = c(7, 9, 77, 99, 777, 999, 7777, 9999)) {
  if (is.numeric(x)) x[x %in% miss] <- NA
  x
}

df <- df %>%
  mutate(
    RIAGENDR = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
    RIDRETH1 = factor(na_if_any(RIDRETH1)),
    DMDEDUC2 = factor(na_if_any(DMDEDUC2)),
    DMDMARTL = if ("DMDMARTL" %in% names(.)) factor(na_if_any(DMDMARTL)) else DMDMARTL
  )

df <- df %>%
  mutate(
    z_cerad_immediate = as.numeric(scale(cerad_immediate)),
    z_cerad_delayed = as.numeric(scale(CFDCSR)),
    z_aft = as.numeric(scale(CFDAST)),
    z_dsst = as.numeric(scale(CFDDS)),
    cognition_global = rowMeans(across(c(z_cerad_immediate, z_cerad_delayed, z_aft, z_dsst)), na.rm = FALSE)
  )

nhanes_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights= ~WTMEC4YR,
  nest = TRUE,
  data = df
)

cov_m1 <- c("RIDAGEYR","RIAGENDR","RIDRETH1")
cov_m2 <- c(cov_m1, "DMDEDUC2","INDFMPIR")
cov_m3 <- c(cov_m2, "BMXBMI","LBXCOT","ALQ130")

keep_in_design <- function(vars, design) vars[vars %in% names(design$variables)]
cov_m1 <- keep_in_design(cov_m1, nhanes_design)
cov_m2 <- keep_in_design(cov_m2, nhanes_design)
cov_m3 <- keep_in_design(cov_m3, nhanes_design)

make_rcs_cols <- function(x, knots, xvar) {
  mat <- Hmisc::rcspline.eval(x, knots = knots, inclx = TRUE)
  colnames(mat) <- c(paste0(xvar, "_lin"),
                     paste0(xvar, "_nl", seq_len(ncol(mat) - 1)))
  as_tibble(mat)
}

fit_models_with_rcs <- function(design, outcome, xvar, covars,
                                knot_probs = c(0.05,0.35,0.65,0.95)) {
  dat <- design$variables
  x <- dat[[xvar]]
  knots <- as.numeric(quantile(x, probs = knot_probs, na.rm = TRUE, type = 2))
  knots <- unique(knots)
  if (length(knots) < 4) {
    knots <- unique(as.numeric(quantile(x, probs = c(0.10,0.50,0.90), na.rm = TRUE, type = 2)))
  }
  rcs_df <- make_rcs_cols(x, knots, xvar)
  dat2 <- bind_cols(dat, rcs_df)
  design2 <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights= ~WTMEC4YR,
    nest = TRUE,
    data = dat2
  )
  cov_part <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
  nl_cols <- names(rcs_df)[grepl(paste0("^", xvar, "_nl"), names(rcs_df))]
  f_lin <- as.formula(paste(outcome, "~", xvar, "+", cov_part))
  f_rcs <- as.formula(paste(outcome, "~",
                            paste(c(paste0(xvar, "_lin"), nl_cols), collapse = " + "),
                            "+", cov_part))
  m_lin <- svyglm(f_lin, design = design2)
  m_rcs <- svyglm(f_rcs, design = design2)
  p_overall <- as.numeric(regTermTest(m_rcs,
                                      as.formula(paste("~", paste(c(paste0(xvar, "_lin"), nl_cols), collapse = " + ")))
  )$p)
  p_nonlinear <- as.numeric(regTermTest(m_rcs,
                                        as.formula(paste("~", paste(nl_cols, collapse = " + ")))
  )$p)
  list(knots = knots, model_linear = m_lin, model_rcs = m_rcs,
       p_overall = p_overall, p_nonlinear = p_nonlinear)
}

exposure <- "OBS_total"
outcomes <- c("cerad_immediate", "CFDCSR", "CFDAST", "CFDDS", "cognition_global")

run_all <- function(outcome) {
  r1 <- fit_models_with_rcs(nhanes_design, outcome, exposure, cov_m1)
  r2 <- fit_models_with_rcs(nhanes_design, outcome, exposure, cov_m2)
  r3 <- fit_models_with_rcs(nhanes_design, outcome, exposure, cov_m3)
  tibble(
    outcome = outcome,
    model = c("Model1","Model2","Model3"),
    p_overall = c(r1$p_overall, r2$p_overall, r3$p_overall),
    p_nonlinear = c(r1$p_nonlinear, r2$p_nonlinear, r3$p_nonlinear)
  )
}

p_table <- purrr::map_dfr(outcomes, run_all)
print(p_table)
write.csv(p_table, file.path(xpt_dir, "step4_p_overall_p_nonlinear.csv"), row.names = FALSE)

if (!("OBS_z" %in% names(df))) df$OBS_z <- as.numeric(scale(df$OBS_total))

fit_linear_table <- function(outcome, covars, model_name) {
  cov_part <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
  f <- as.formula(paste(outcome, "~ OBS_z +", cov_part))
  m <- svyglm(f, design = nhanes_design)
  tt <- broom::tidy(m) %>% filter(term == "OBS_z")
  tt %>% transmute(
    outcome = outcome, model = model_name,
    beta = estimate, se = std.error, p = p.value,
    ci_low = estimate - 1.96*std.error,
    ci_high= estimate + 1.96*std.error
  )
}

lin_table <- bind_rows(
  purrr::map_dfr(outcomes, ~fit_linear_table(.x, cov_m1, "Model1")),
  purrr::map_dfr(outcomes, ~fit_linear_table(.x, cov_m2, "Model2")),
  purrr::map_dfr(outcomes, ~fit_linear_table(.x, cov_m3, "Model3"))
)

print(lin_table)
write.csv(lin_table, file.path(xpt_dir, "step4_linear_beta_perSD.csv"), row.names = FALSE)

plot_rcs <- function(outcome, covars, xvar = "OBS_total",
                     knot_probs = c(0.05,0.35,0.65,0.95),
                     ref = c("median","mean","p10"),
                     file = NULL) {
  fit <- fit_models_with_rcs(nhanes_design, outcome, xvar, covars, knot_probs = knot_probs)
  m <- fit$model_rcs
  dat <- nhanes_design$variables
  x <- dat[[xvar]]
  x_grid <- seq(quantile(x, 0.01, na.rm=TRUE), quantile(x, 0.99, na.rm=TRUE), length.out = 200)
  rcs_grid <- as_tibble(Hmisc::rcspline.eval(x_grid, knots = fit$knots, inclx = TRUE))
  colnames(rcs_grid) <- c(paste0(xvar, "_lin"),
                          paste0(xvar, "_nl", seq_len(ncol(rcs_grid) - 1)))
  mode_val <- function(v) {
    u <- unique(v[!is.na(v)])
    if (length(u) == 0) return(NA)
    u[which.max(tabulate(match(v, u)))]
  }
  pred_base <- tibble(!!xvar := x_grid) %>% bind_cols(rcs_grid)
  for (cv in covars) {
    if (!cv %in% names(dat)) next
    if (is.numeric(dat[[cv]])) pred_base[[cv]] <- median(dat[[cv]], na.rm = TRUE)
    else pred_base[[cv]] <- mode_val(dat[[cv]])
  }
  pr <- predict(m, newdata = pred_base, se.fit = TRUE, type = "link")
  if (is.list(pr) && !is.null(pr$fit)) {
    pred_base$fit <- as.numeric(pr$fit)
    pred_base$se <- as.numeric(pr$se.fit)
  } else {
    pred_base$fit <- as.numeric(pr)
    X <- model.matrix(delete.response(terms(m)), data = pred_base)
    V <- vcov(m)
    pred_base$se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  }
  pred_base$ci_low <- pred_base$fit - 1.96 * pred_base$se
  pred_base$ci_high <- pred_base$fit + 1.96 * pred_base$se
  ref_type <- ref[1]
  x_ref <- switch(ref_type,
                  median = median(x, na.rm=TRUE),
                  mean = mean(x, na.rm=TRUE),
                  p10 = as.numeric(quantile(x, 0.10, na.rm=TRUE, type=2)),
                  median(x, na.rm=TRUE))
  ref_rcs <- as_tibble(Hmisc::rcspline.eval(x_ref, knots = fit$knots, inclx = TRUE))
  colnames(ref_rcs) <- colnames(rcs_grid)
  ref_row <- tibble(!!xvar := x_ref) %>% bind_cols(ref_rcs)
  for (cv in covars) {
    if (!cv %in% names(dat)) next
    if (is.numeric(dat[[cv]])) ref_row[[cv]] <- median(dat[[cv]], na.rm = TRUE)
    else ref_row[[cv]] <- mode_val(dat[[cv]])
  }
  ref_pred <- as.numeric(predict(m, newdata = ref_row, type = "link"))
  pred_base <- pred_base %>%
    mutate(
      fit_centered = fit - ref_pred,
      ci_low_centered = ci_low - ref_pred,
      ci_high_centered = ci_high - ref_pred
    )
  g <- ggplot(pred_base, aes(x = .data[[xvar]], y = fit_centered)) +
    geom_ribbon(aes(ymin = ci_low_centered, ymax = ci_high_centered), alpha = 0.2) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      x = xvar,
      y = paste0(outcome, " (difference vs ", ref_type, " ", xvar, ")"),
      title = paste0("Survey-weighted RCS (Model 3): ", outcome),
      subtitle = paste0("p_overall=", signif(fit$p_overall,3), ", p_nonlinear=", signif(fit$p_nonlinear,3))
    ) +
    theme_classic()
  if (!is.null(file)) ggsave(file, g, width = 7, height = 5, dpi = 300)
  g
}

g1 <- plot_rcs("cognition_global", cov_m3, xvar = "OBS_total",
               file = file.path(xpt_dir, "rcs_Model3_cognition_global.png"))
print(g1)

fit_linear_table_per10 <- function(outcome, covars, model_name) {
  cov_part <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
  f <- as.formula(paste(outcome, "~ I(OBS_total/10) +", cov_part))
  m <- svyglm(f, design = nhanes_design)
  tt <- broom::tidy(m) %>% filter(term == "I(OBS_total/10)")
  tt %>% transmute(
    outcome = outcome, model = model_name,
    beta_per10 = estimate, se = std.error, p = p.value,
    ci_low = estimate - 1.96*std.error,
    ci_high= estimate + 1.96*std.error
  )
}

lin_table_per10 <- bind_rows(
  purrr::map_dfr(outcomes, ~fit_linear_table_per10(.x, cov_m1, "Model1")),
  purrr::map_dfr(outcomes, ~fit_linear_table_per10(.x, cov_m2, "Model2")),
  purrr::map_dfr(outcomes, ~fit_linear_table_per10(.x, cov_m3, "Model3"))
)

print(lin_table_per10)
write.csv(lin_table_per10, file.path(xpt_dir, "step4_linear_beta_per10pts.csv"), row.names = FALSE)

# ==========================================================
# STEP5) Quartiles vs global cognition trend plot
# IMPORTANT: use nhanes_design in memory (no re-read)
# ==========================================================
message("\n===== STEP5 =====")
df5 <- nhanes_design$variables %>%
  mutate(
    OBS_q = ntile(OBS_total, 4),
    OBS_q = factor(OBS_q, levels = 1:4,
                   labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"))
  )

design_q <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC4YR,
  nest = TRUE,
  data = df5
)

mean_q <- svyby(~cognition_global, ~OBS_q, design = design_q, svymean, vartype = "se")

p_trend <- ggplot(mean_q, aes(x = OBS_q, y = cognition_global)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = cognition_global - se, ymax = cognition_global + se), width = 0.15) +
  geom_line(aes(group = 1), linewidth = 0.8) +
  theme_classic()

print(p_trend)
ggsave(file.path(xpt_dir, "OBS_quartiles_vs_global_cognition.png"),
       p_trend, width = 6, height = 4.5, dpi = 300)

trend_model <- svyglm(cognition_global ~ as.numeric(OBS_q), design = design_q)
print(summary(trend_model))

# ==========================================================
# STEP6) Forest plot (Q2–Q4 vs Q1) for global cognition
# IMPORTANT: use design_q from STEP5 (no re-read)
# ==========================================================
message("\n===== STEP6 =====")
forest_model <- svyglm(
  cognition_global ~ OBS_q + RIDAGEYR + RIAGENDR + RIDRETH1 +
    DMDEDUC2 + INDFMPIR + BMXBMI + LBXCOT + ALQ130,
  design = design_q
)

forest_df <- tidy(forest_model) %>%
  filter(grepl("OBS_q", term)) %>%
  mutate(
    term = recode(term,
                  "OBS_qQ2" = "Q2 vs Q1",
                  "OBS_qQ3" = "Q3 vs Q1",
                  "OBS_qQ4 (Highest)" = "Q4 vs Q1")
  )

p_forest <- ggplot(forest_df, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = estimate - 1.96*std.error,
                     xmax = estimate + 1.96*std.error),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "β (95% CI) for global cognition", y = "",
       title = "Association between OBS quartiles and global cognition (Model 3)") +
  theme_classic()

print(p_forest)
ggsave(file.path(xpt_dir, "OBS_quartiles_forest_global_cognition.png"),
       p_forest, width = 6, height = 4, dpi = 300)

print(summary(forest_model))

message("\nALL DONE. Outputs saved to: ", xpt_dir)

# ==========================================================
# STEP 3.5) Baseline Characteristics Table (Table 1)
# 这一步生成按 OBS 四分位数分组的基线特征表
# ==========================================================
message("\n===== STEP 3.5: Generating Table 1 =====")

if (!requireNamespace("tableone", quietly = TRUE)) {
  warning("需要 'tableone' 包来生成基线表。请运行: install.packages('tableone')")
} else {
  library(tableone)
}

df_tbl <- readRDS(file.path(xpt_dir, "df_step3_OBS_final.rds")) %>%
  filter(!is.na(OBS_total))

df_tbl <- df_tbl %>%
  mutate(
    OBS_Quartile = ntile(OBS_total, 4),
    OBS_Quartile = factor(OBS_Quartile, levels = 1:4,
                          labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)")),
   
    Gender = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
   
    Race = factor(RIDRETH1, levels = c(1, 2, 3, 4, 5),
                  labels = c("Mexican American", "Other Hispanic",
                             "Non-Hispanic White", "Non-Hispanic Black", "Other Race")),
   
    Education = factor(DMDEDUC2, levels = c(1, 2, 3, 4, 5),
                       labels = c("< 9th Grade", "9-11th Grade", "High School Grad",
                                  "Some College", "College Grad")),
   
    Smoking_Status = case_when(
      LBXCOT < 0.015 ~ "Non-smoker",
      LBXCOT >= 0.015 & LBXCOT < 10 ~ "Passive/Light",
      LBXCOT >= 10 ~ "Active Smoker",
      TRUE ~ NA_character_
    ),
    Smoking_Status = factor(Smoking_Status, levels = c("Non-smoker", "Passive/Light", "Active Smoker"))
  )

cont_vars <- c("RIDAGEYR", "BMXBMI", "INDFMPIR", "OBS_total",
               "cognition_global", "cerad_immediate", "CFDCSR", "CFDAST", "CFDDS")

cat_vars <- c("Gender", "Race", "Education", "Smoking_Status")

my_vars <- c(cont_vars, cat_vars)

design_tbl <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC4YR,
  nest = TRUE,
  data = df_tbl
)

message("正在计算加权基线表 (可能需要几秒钟)...")

tab1 <- svyCreateTableOne(
  vars = my_vars,
  strata = "OBS_Quartile",
  data = design_tbl,
  test = TRUE,
  addOverall = TRUE
)

tab1_mat <- print(tab1, printToggle = FALSE, showAllLevels = TRUE, smd = FALSE)

write.csv(tab1_mat, file.path(xpt_dir, "Table1_Baseline_Characteristics.csv"))
message("Table 1 已保存至: Table1_Baseline_Characteristics.csv")
