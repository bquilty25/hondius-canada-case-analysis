args <- commandArgs(trailingOnly = TRUE)

date_input_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  "output/canada_case_probability/baseline_date_distribution.csv"
}
output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  dirname(date_input_path)
}

probability_input_path <- file.path(output_dir, "baseline_draws.csv")
source_input_path <- file.path(output_dir, "baseline_source_distribution.csv")
linelist_input_path <- "external/Hondius_hantavirus_h2026/data/linelist/2026_hantavirus.csv"
plot_end_date <- as.Date("2026-05-17")

reference_lines <- data.frame(
  line_label = c("Peak date", "Main disembarkation"),
  exposure_date = as.Date(c(NA, "2026-05-10"))
)

if (!file.exists(date_input_path)) {
  stop(sprintf("Input file not found: %s", date_input_path))
}

suppressPackageStartupMessages(library(ggplot2))

weighted_quantile <- function(values, weights, probability = 0.5) {
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord] / sum(weights)
  cumulative <- cumsum(weights)
  values[which(cumulative >= probability)[1]]
}

date_summary <- read.csv(date_input_path, check.names = FALSE)
date_summary$exposure_date <- as.Date(date_summary$exposure_date)

required_cols <- c(
  "hypothesis",
  "hypothesis_label",
  "exposure_date",
  "probability_mean",
  "probability_lower",
  "probability_upper"
)
missing_cols <- setdiff(required_cols, names(date_summary))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Missing required columns in %s: %s",
    date_input_path,
    paste(missing_cols, collapse = ", ")
  ))
}

date_summary$hypothesis_plot <- factor(
  ifelse(date_summary$hypothesis == "secondary", "Secondary", "Tertiary"),
  levels = c("Secondary", "Tertiary")
)

summary_data <- do.call(
  rbind,
  lapply(split(date_summary, date_summary$hypothesis_plot), function(df) {
    df <- df[order(df$exposure_date), ]
    out <- df[which.max(df$probability_mean), c("hypothesis_plot", "exposure_date")]
    out$line_label <- "Peak date"
    out
  })
)
palette <- c(Secondary = "#0B6E4F", Tertiary = "#C84C09")
line_palette <- c("Peak date" = "#666666", "Main disembarkation" = "#3B5BA5")
line_types <- c("Peak date" = "dashed", "Main disembarkation" = "dotdash")
line_breaks <- c("Peak date", "Main disembarkation")
plot_start_date <- min(date_summary$exposure_date, na.rm = TRUE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

date_density_plot <- ggplot(
  date_summary,
  aes(x = exposure_date, colour = hypothesis_plot, fill = hypothesis_plot)
) +
  geom_ribbon(
    aes(ymin = probability_lower, ymax = probability_upper),
    alpha = 0.16,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_area(aes(y = probability_mean), alpha = 0.18, position = "identity", show.legend = FALSE) +
  geom_line(aes(y = probability_mean), linewidth = 1.1) +
  geom_vline(
    data = summary_data,
    aes(xintercept = exposure_date, linetype = line_label, colour = line_label),
    linewidth = 0.8,
    show.legend = TRUE
  ) +
  geom_vline(
    data = data.frame(exposure_date = as.Date("2026-05-10"), line_label = "Main disembarkation"),
    aes(xintercept = exposure_date, linetype = line_label, colour = line_label),
    linewidth = 0.8,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = palette) +
  scale_linetype_manual(values = line_types, breaks = line_breaks, name = NULL) +
  scale_x_date(
    date_labels = "%d %b",
    date_breaks = "3 days",
    limits = c(plot_start_date, plot_end_date),
    expand = c(0.01, 0.01)
  ) +
  scale_colour_manual(values = c(palette, line_palette), breaks = c("Secondary", "Tertiary"), name = NULL) +
  labs(
    title = "Joint posterior support across exposure dates",
    subtitle = "Bayesian support by exposure date from the fitted transmission-timing and incubation posteriors",
    x = "Exposure date",
    y = "Joint posterior mean probability mass",
    colour = NULL,
    fill = NULL
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(linetype = "solid", shape = 16, linewidth = 1.1)),
    linetype = guide_legend(order = 2, override.aes = list(colour = unname(line_palette[line_breaks]), linewidth = 0.9)),
    fill = "none"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  filename = file.path(output_dir, "date_overlap_density.png"),
  plot = date_density_plot,
  width = 8,
  height = 5,
  dpi = 320
)

if (file.exists(linelist_input_path)) {
  onset_raw <- read.csv(linelist_input_path, check.names = FALSE)

  parse_date_col <- function(x) as.Date(ifelse(x %in% c("", "NA", NA), NA, x))
  onset_raw$symptom_onset <- parse_date_col(onset_raw$symptom_onset)
  onset_raw$confirmation_date <- parse_date_col(onset_raw$confirmation_date)
  onset_raw$outcome_date <- parse_date_col(onset_raw$outcome_date)

  # Keep only cases that have at least one dateable event
  has_date <- !is.na(onset_raw$symptom_onset) |
    !is.na(onset_raw$confirmation_date) |
    !is.na(onset_raw$outcome_date)
  onset_raw <- onset_raw[has_date, ]

  # Add Canadian case as a synthetic row
  canada_row <- data.frame(
    Gh_ID = "CA",
    status = "confirmed",
    symptom_onset = as.Date("2026-05-14"),
    confirmation_date = as.Date("2026-05-15"),
    outcome_date = as.Date(NA),
    stringsAsFactors = FALSE
  )
  onset_raw <- rbind(
    onset_raw[, intersect(names(onset_raw), c("Gh_ID", "status", "symptom_onset", "confirmation_date", "outcome_date"))],
    canada_row
  )

  # Order cases by symptom onset (earliest at bottom); Canadian case at bottom
  case_order <- c(
    "CA",
    rev(onset_raw$Gh_ID[onset_raw$Gh_ID != "CA"][
      order(
        onset_raw$symptom_onset[onset_raw$Gh_ID != "CA"],
        onset_raw$Gh_ID[onset_raw$Gh_ID != "CA"]
      )
    ])
  )
  onset_raw$case_label <- ifelse(
    onset_raw$Gh_ID == "CA", "Canadian case",
    paste("Hondius case", onset_raw$Gh_ID)
  )
  label_order <- ifelse(case_order == "CA", "Canadian case", paste("Hondius case", case_order))
  onset_raw$case_label <- factor(onset_raw$case_label, levels = label_order)

  # Pivot to long format
  onset_long <- rbind(
    data.frame(
      case_label = onset_raw$case_label, status = onset_raw$status,
      date = onset_raw$symptom_onset, date_type = "Symptom onset",
      stringsAsFactors = FALSE
    ),
    data.frame(
      case_label = onset_raw$case_label, status = onset_raw$status,
      date = onset_raw$confirmation_date, date_type = "Confirmation / positive test",
      stringsAsFactors = FALSE
    ),
    data.frame(
      case_label = onset_raw$case_label, status = onset_raw$status,
      date = onset_raw$outcome_date, date_type = "Outcome",
      stringsAsFactors = FALSE
    )
  )
  onset_long <- onset_long[!is.na(onset_long$date), ]
  onset_long$date <- as.Date(onset_long$date)
  onset_long$date_type <- factor(
    onset_long$date_type,
    levels = c("Symptom onset", "Confirmation / positive test", "Outcome")
  )
  onset_long$status_group <- ifelse(
    onset_long$status %in% c("confirmed", "probable"), onset_long$status, "other"
  )

  onset_plot_start_date <- min(onset_long$date, na.rm = TRUE) - 2

  date_type_shapes <- c(
    "Symptom onset"                   = 16,
    "Confirmation / positive test"    = 15,
    "Outcome"                         = 17
  )
  onset_status_palette <- c(confirmed = "#0077BB", probable = "#EE7733", other = "#BBBBBB")

  onset_plot <- ggplot(
    onset_long,
    aes(x = date, y = case_label, colour = status_group, shape = date_type)
  ) +
    geom_point(size = 3) +
    geom_vline(
      data = data.frame(date = as.Date("2026-05-10"), line_label = "Main disembarkation"),
      aes(xintercept = date, linetype = line_label),
      colour = line_palette[["Main disembarkation"]],
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    scale_colour_manual(values = onset_status_palette, name = "Status") +
    scale_shape_manual(values = date_type_shapes, name = "Date type") +
    scale_linetype_manual(values = c("Main disembarkation" = line_types[["Main disembarkation"]]), name = NULL) +
    scale_x_date(
      date_labels = "%d %b",
      date_breaks = "3 days",
      limits = c(onset_plot_start_date, plot_end_date),
      expand = c(0.01, 0.01)
    ) +
    labs(
      title = "Case event timeline",
      subtitle = "Symptom onset, confirmation, and outcome dates through 17 May 2026",
      x = "Date",
      y = NULL
    ) +
    guides(
      colour   = guide_legend(order = 1),
      shape    = guide_legend(order = 2),
      linetype = guide_legend(order = 3, override.aes = list(colour = line_palette[["Main disembarkation"]], linewidth = 0.9))
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(
    filename = file.path(output_dir, "onset_timings.png"),
    plot = onset_plot,
    width = 9,
    height = 6,
    dpi = 320
  )
}

if (file.exists(source_input_path)) {
  source_summary <- read.csv(source_input_path, check.names = FALSE)
  source_summary$hypothesis_plot <- factor(
    ifelse(source_summary$hypothesis == "secondary", "Generation 2", "Generation 3"),
    levels = c("Generation 2", "Generation 3")
  )
  source_summary$source_label <- factor(
    paste("Case", source_summary$source_id),
    levels = paste("Case", sort(unique(source_summary$source_id)))
  )

  source_plot <- ggplot(
    source_summary,
    aes(x = source_label, y = probability_mean, fill = hypothesis_plot)
  ) +
    geom_col(position = "dodge", width = 0.72) +
    geom_errorbar(
      aes(ymin = probability_lower, ymax = probability_upper),
      position = position_dodge(width = 0.72),
      width = 0.2,
      linewidth = 0.6
    ) +
    scale_fill_manual(values = c("Generation 2" = palette[["Secondary"]], "Generation 3" = palette[["Tertiary"]])) +
    labs(
      title = "Posterior support for candidate source cases",
      subtitle = "Joint posterior probability that each Hondius case seeded the Canadian case",
      x = "Candidate source case",
      y = "Joint posterior mean probability",
      fill = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot"
    )

  ggsave(
    filename = file.path(output_dir, "source_posterior_bars.png"),
    plot = source_plot,
    width = 8,
    height = 5,
    dpi = 320
  )
}

if (file.exists(probability_input_path)) {
  draws <- read.csv(probability_input_path, check.names = FALSE)
  prob_required_cols <- c("draw", "prob_secondary", "prob_tertiary")
  prob_missing_cols <- setdiff(prob_required_cols, names(draws))
  if (length(prob_missing_cols) == 0) {
    if (!"posterior_draw_weight" %in% names(draws)) {
      draws$posterior_draw_weight <- 1 / nrow(draws)
    }

    probability_plot_data <- rbind(
      data.frame(
        draw = draws$draw,
        hypothesis = "Secondary",
        probability = draws$prob_secondary,
        weight = draws$posterior_draw_weight
      ),
      data.frame(
        draw = draws$draw,
        hypothesis = "Tertiary",
        probability = draws$prob_tertiary,
        weight = draws$posterior_draw_weight
      )
    )

    probability_plot_data$hypothesis <- factor(
      probability_plot_data$hypothesis,
      levels = c("Secondary", "Tertiary")
    )

    probability_medians <- do.call(
      rbind,
      lapply(split(probability_plot_data, probability_plot_data$hypothesis), function(df) {
        data.frame(
          hypothesis = df$hypothesis[[1]],
          probability = weighted_quantile(df$probability, df$weight, 0.5)
        )
      })
    )

    legacy_probability_plot <- ggplot(
      probability_plot_data,
      aes(x = probability, colour = hypothesis, fill = hypothesis, weight = weight)
    ) +
      geom_density(alpha = 0.25, linewidth = 1.1, adjust = 1.1) +
      geom_vline(
        data = probability_medians,
        aes(xintercept = probability, linetype = "Peak value"),
        linewidth = 0.8,
        colour = "#666666",
        show.legend = TRUE
      ) +
      scale_colour_manual(values = palette) +
      scale_fill_manual(values = palette) +
      scale_linetype_manual(values = c("Peak value" = "dashed"), name = NULL) +
      scale_x_continuous(expand = c(0.01, 0.01)) +
      coord_cartesian(xlim = c(0, 1)) +
      labs(
        title = "Posterior probability distributions for the Canada case",
        subtitle = "Weighted by the posterior draw support under the fully Bayesian Canada-case update",
        x = "Posterior probability",
        y = "Density",
        colour = NULL,
        fill = NULL
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        panel.grid.minor = element_blank(),
        plot.title.position = "plot"
      )

    ggsave(
      filename = file.path(output_dir, "probability_overlap_density.png"),
      plot = legacy_probability_plot,
      width = 8,
      height = 5,
      dpi = 320
    )
  }
}

message(sprintf("Wrote ggplot outputs to %s", output_dir))
