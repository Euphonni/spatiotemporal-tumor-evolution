suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(FNN)
  library(stringr)
  library(grid)
  library(Ckmeans.1d.dp)
})

## ==================================================
## parameters
## ==================================================
input_dir <- "data/example_input"
file_pattern <- "^Integrated_Data_.*\\.csv$"

REGION_COL  <- "region"
BARCODE_COL <- "barcode"

alpha <- 0.85
N_CP_BANDS <- 7
KNN_K <- 10
FRONT_Q <- 0.85

FLOW_SHOW_Q <- 0.70
FLOW_MAX_ARROWS <- 800
FLOW_SCALE <- 800

out_root <- "results"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

## ==================================================
## unified theme
## ==================================================
theme_efm <- theme_void() +
  theme(
    plot.title = element_blank()
  )

## ==================================================
## input files
## ==================================================
input_files <- list.files(
  input_dir,
  pattern = file_pattern,
  full.names = TRUE
)

if (length(input_files) == 0) {
  stop(
    "No input files found in: ", input_dir,
    "\nExpected files like: Integrated_Data_region5.csv"
  )
}

## ==================================================
## main loop
## ==================================================
for (input_file in input_files) {
  
  message("Processing: ", basename(input_file))
  
  REGION_NAME <- str_remove(basename(input_file), "^Integrated_Data_|\\.csv$")
  prefix <- REGION_NAME
  
  out_dir <- file.path(out_root, REGION_NAME)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  cells <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  ## ----------------------------
  ## required columns check
  ## ----------------------------
  need_cols <- c(
    BARCODE_COL, REGION_COL, "pseudotime",
    "pxl_row_in_fullres", "pxl_col_in_fullres"
  )
  miss <- setdiff(need_cols, colnames(cells))
  if (length(miss) > 0) {
    message("SKIP: missing columns: ", paste(miss, collapse = ", "))
    next
  }
  
  ## ----------------------------
  ## basic cleanup
  ## ----------------------------
  cells <- cells %>%
    filter(
      !is.na(.data[[BARCODE_COL]]),
      !is.na(.data[[REGION_COL]]),
      !is.na(pseudotime),
      !is.na(pxl_row_in_fullres),
      !is.na(pxl_col_in_fullres)
    )
  
  if (nrow(cells) == 0) {
    message("SKIP: no valid rows after cleanup")
    next
  }
  
  ## ----------------------------
  ## pseudotime / coordinates
  ## ----------------------------
  cells$pseudotime <- as.numeric(cells[["pseudotime"]])
  cells$x <- as.numeric(cells[["pxl_row_in_fullres"]])
  cells$y <- as.numeric(cells[["pxl_col_in_fullres"]])
  
  ## ----------------------------
  ## pseudotime normalization
  ## ----------------------------
  pt_min <- min(cells$pseudotime, na.rm = TRUE)
  pt_max <- max(cells$pseudotime, na.rm = TRUE)
  
  if (!is.finite(pt_min) || !is.finite(pt_max) || pt_min == pt_max) {
    message("SKIP: pseudotime invalid or range = 0")
    next
  }
  
  cells$t_global <- (cells$pseudotime - pt_min) / (pt_max - pt_min)
  
  ## ----------------------------
  ## r_dir by microzone
  ## ----------------------------
  res_list <- list()
  
  for (mz in unique(cells[[REGION_COL]])) {
    
    df <- cells %>% filter(.data[[REGION_COL]] == mz)
    
    if (nrow(df) < 10) next
    
    cx <- median(df$x, na.rm = TRUE)
    cy <- median(df$y, na.rm = TRUE)
    
    dx <- df$x - cx
    dy <- df$y - cy
    
    df$d <- sqrt(dx^2 + dy^2)
    df$theta <- atan2(dy, dx)
    
    df$theta_bin <- cut(
      df$theta,
      breaks = seq(-pi, pi, length.out = 73),
      include.lowest = TRUE
    )
    
    R_theta <- df %>%
      group_by(theta_bin) %>%
      summarise(R = max(d, na.rm = TRUE), .groups = "drop")
    
    df <- df %>% left_join(R_theta, by = "theta_bin")
    df$r_dir <- ifelse(df$R > 0, df$d / df$R, 0)
    
    res_list[[as.character(mz)]] <- df
  }
  
  if (length(res_list) == 0) {
    message("SKIP: all microzones have <10 cells for r_dir")
    next
  }
  
  cells <- bind_rows(res_list)
  
  ## ----------------------------
  ## CP
  ## ----------------------------
  cells$CP <- alpha * cells$t_global + (1 - alpha) * cells$r_dir
  
  ## ----------------------------
  ## KNN
  ## ----------------------------
  coords <- as.matrix(cells[, c("x", "y")])
  
  if (nrow(coords) <= KNN_K) {
    message("SKIP: not enough cells for KNN (n <= KNN_K)")
    next
  }
  
  knn <- get.knn(coords, k = KNN_K)
  
  ## ----------------------------
  ## CP slope
  ## ----------------------------
  cells$CP_slope <- sapply(seq_len(nrow(cells)), function(i) {
    nb <- knn$nn.index[i, ]
    mean(abs(cells$CP[i] - cells$CP[nb]), na.rm = TRUE)
  })
  
  ## ----------------------------
  ## Evolutionary flow
  ## ----------------------------
  cells$flow_vx <- NA_real_
  cells$flow_vy <- NA_real_
  
  for (i in seq_len(nrow(cells))) {
    nb <- knn$nn.index[i, ]
    
    dx <- cells$x[nb] - cells$x[i]
    dy <- cells$y[nb] - cells$y[i]
    dist <- sqrt(dx^2 + dy^2) + 1e-6
    
    dCP <- cells$CP[nb] - cells$CP[i]
    
    ux <- dx / dist
    uy <- dy / dist
    w  <- dCP / dist
    
    cells$flow_vx[i] <- mean(w * ux, na.rm = TRUE)
    cells$flow_vy[i] <- mean(w * uy, na.rm = TRUE)
  }
  
  cells$flow_mag <- sqrt(cells$flow_vx^2 + cells$flow_vy^2)
  
  ## ----------------------------
  ## evolutionary front
  ## ----------------------------
  cp_thr <- quantile(cells$CP, FRONT_Q, na.rm = TRUE)
  slope_thr <- quantile(cells$CP_slope, FRONT_Q, na.rm = TRUE)
  
  cells$is_front_CP <- cells$CP >= cp_thr
  cells$is_front_slope <- cells$CP_slope >= slope_thr
  cells$is_front <- cells$is_front_CP | cells$is_front_slope
  
  ## ==================================================
  ## CP bands by ckmeans
  ## ==================================================
  x_cp <- cells$CP
  x_valid <- x_cp[is.finite(x_cp)]
  
  if (length(x_valid) < 5) {
    message("SKIP: too few valid CP values for ckmeans")
    next
  }
  
  max_k <- min(N_CP_BANDS, length(unique(x_valid)))
  if (max_k < 2) {
    message("SKIP: not enough unique CP values for clustering")
    next
  }
  
  sse <- numeric(max_k)
  
  for (k in seq_len(max_k)) {
    ck_tmp <- Ckmeans.1d.dp(x_valid, k = k)
    groups <- split(x_valid, ck_tmp$cluster)
    sse[k] <- sum(vapply(groups, function(g) sum((g - mean(g))^2), numeric(1)))
  }
  
  if (length(sse) < 3) {
    k_opt <- max_k
  } else {
    d2 <- diff(diff(sse))
    if (length(d2) == 0) {
      k_opt <- max_k
    } else {
      k_opt <- which.max(-d2) + 1
      k_opt <- max(2, min(k_opt, max_k))
    }
  }
  
  ck <- Ckmeans.1d.dp(x_valid, k = k_opt)
  
  cells$CP_band <- NA_integer_
  cells$CP_band[is.finite(cells$CP)] <- ck$cluster
  
  ## ==================================================
  ## plots
  ## ==================================================
  p1 <- ggplot(cells, aes(x = x, y = y, color = CP)) +
    geom_point(size = 0.4) +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_efm
  
  p2 <- ggplot(cells, aes(x = x, y = y, color = CP_slope)) +
    geom_point(size = 0.4) +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_efm
  
  p3 <- ggplot(cells, aes(x = x, y = y)) +
    geom_point(color = "grey85", size = 0.3) +
    geom_point(
      data = subset(cells, is_front),
      color = "#F8766D",
      size = 0.5
    ) +
    coord_fixed() +
    theme_efm
  
  p4 <- ggplot(cells, aes(x = x, y = y, color = factor(CP_band))) +
    geom_point(size = 0.4) +
    coord_fixed() +
    theme_efm +
    labs(color = "EFM Subclone\n")
  
  set.seed(123)
  flow_thr <- quantile(cells$flow_mag, FLOW_SHOW_Q, na.rm = TRUE)
  
  flow_df <- cells %>%
    filter(is.finite(flow_mag), flow_mag >= flow_thr)
  
  if (nrow(flow_df) > 0) {
    flow_df <- flow_df %>%
      slice_sample(n = min(FLOW_MAX_ARROWS, nrow(flow_df)))
  }
  
  p5 <- ggplot(cells, aes(x = x, y = y)) +
    geom_point(color = "grey80", size = 0.3, alpha = 0.6) +
    geom_segment(
      data = flow_df,
      aes(
        xend = x + flow_vx * FLOW_SCALE,
        yend = y + flow_vy * FLOW_SCALE
      ),
      arrow = arrow(type = "open", angle = 25, length = unit(0.12, "cm")),
      linewidth = 0.3,
      color = "#D55E00",
      alpha = 0.8
    ) +
    coord_fixed() +
    theme_efm
  
  ## ==================================================
  ## save outputs
  ## ==================================================
  ggsave(
    filename = paste0(prefix, "_EFM_height_CP.png"),
    plot = p1,
    path = out_dir,
    width = 6, height = 6, dpi = 300
  )
  
  ggsave(
    filename = paste0(prefix, "_EFM_slope_CP.png"),
    plot = p2,
    path = out_dir,
    width = 6, height = 6, dpi = 300
  )
  
  ggsave(
    filename = paste0(prefix, "_EFM_front.png"),
    plot = p3,
    path = out_dir,
    width = 6, height = 6, dpi = 300
  )
  
  ggsave(
    filename = paste0(prefix, "_EFM_subclones.png"),
    plot = p4,
    path = out_dir,
    width = 6, height = 6, dpi = 300
  )
  
  ggsave(
    filename = paste0(prefix, "_EFM_flow.png"),
    plot = p5,
    path = out_dir,
    width = 6, height = 6, dpi = 300
  )
  
  write.csv(
    cells,
    file.path(out_dir, paste0(prefix, "_EFM_cells.csv")),
    row.names = FALSE
  )
  
  write.csv(
    cells %>% select(all_of(BARCODE_COL), CP, CP_band),
    file.path(out_dir, paste0(prefix, "_EFM_Loupe_CP_subclone.csv")),
    row.names = FALSE
  )
  
  summary_df <- data.frame(
    region_file = prefix,
    n_cells = nrow(cells),
    n_microzones = dplyr::n_distinct(cells[[REGION_COL]]),
    alpha = alpha,
    knn_k = KNN_K,
    front_quantile = FRONT_Q,
    cp_band_k = k_opt,
    cp_threshold = cp_thr,
    slope_threshold = slope_thr
  )
  
  write.csv(
    summary_df,
    file.path(out_dir, paste0(prefix, "_summary.csv")),
    row.names = FALSE
  )
  
  graphics.off()
}

message("✅ All regions processed successfully.")