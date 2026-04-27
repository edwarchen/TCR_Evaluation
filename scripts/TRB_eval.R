library(openPrimeR)
library(Biostrings)

hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")

cfg <- list(
  TRBV = list(
    ref_path = "data/TRBV_clean.fasta",
    primer_path = "data/primer_set/TRBV_primers.fasta",
    direction = "fw",
    primer_id_suffix = "_fw"
  ),
  TRBJ = list(
    ref_path = "data/TRBJ_clean.fasta",
    primer_path = "data/primer_set/TRBJ_primers.fasta",
    direction = "rev",
    primer_id_suffix = "_rev"
  )
)

load_ref <- function(ref_path, direction) {
  ref <- read_templates(ref_path, hdr.structure, delim = "|", id.column = "GROUP")
  ref$Group <- gsub("[-*].*", "", ref$Group)

  # 为 forward 引物显式设置允许结合区间，避免不同模板格式导致的不一致。
  if (direction == "fw") {
    ref$Allowed_Start_fw <- 1
    ref$Allowed_End_fw <- nchar(ref$Sequence)
    ref$Allowed_fw <- ref$Sequence
  }
  ref
}

load_primers_by_direction <- function(primer_path, direction, primer_id_suffix) {
  if (direction == "fw") {
    read_primers(primer_path, fw.id = primer_id_suffix)
  } else {
    read_primers(primer_path, rev.id = primer_id_suffix)
  }
}

run_coverage <- function(one_cfg, eval_settings) {
  ref <- load_ref(one_cfg$ref_path, one_cfg$direction)
  primers <- load_primers_by_direction(
    one_cfg$primer_path,
    one_cfg$direction,
    one_cfg$primer_id_suffix
  )
  constraint.df <- check_constraints(
    primers,
    ref,
    eval_settings,
    active.constraints = "primer_coverage"
  )
  template.df <- update_template_cvg(ref, constraint.df)
  list(
    constraint.df = constraint.df,
    template.df = template.df
  )
}

summarize_and_plot <- function(
  name,
  one_result,
  plot_dir = "results/coverage_plots"
) {
  constraint.df <- one_result$constraint.df
  template.df <- one_result$template.df

  ratio <- as.numeric(get_cvg_ratio(constraint.df, template.df))
  cvg.stats <- get_cvg_stats(constraint.df, template.df, for.viewing = TRUE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_file <- file.path(plot_dir, paste0(name, "_cov.png"))

  grDevices::png(plot_file, width = 1400, height = 900, res = 120)
  print(plot_template_cvg(constraint.df, template.df))
  grDevices::dev.off()

  list(
    name = name,
    ratio = ratio,
    cvg.stats = cvg.stats,
    plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE)
  )
}

save_cov_stats <- function(name, one_summary, output_dir = "results/coverage_tables") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stats_file <- file.path(output_dir, paste0(name, "_cov.csv"))
  write.csv(one_summary$cvg.stats, stats_file, row.names = FALSE)
  normalizePath(stats_file, winslash = "/", mustWork = FALSE)
}

xml_file <- system.file("extdata", "settings", "A_Taq_PCR_design.xml", package = "openPrimeR")
eval_settings <- read_settings(xml_file)
conOptions(eval_settings)$allowed_other_binding_ratio <- c(max = 1.0)

results <- lapply(cfg, run_coverage, eval_settings = eval_settings)
if (is.null(names(results))) {
  names(results) <- paste0("target_", seq_along(results))
}
coverage_summary <- mapply(
  summarize_and_plot,
  names(results),
  results,
  SIMPLIFY = FALSE
)
saved_stats <- mapply(
  save_cov_stats,
  names(coverage_summary),
  coverage_summary,
  SIMPLIFY = FALSE
)

for (target_name in names(results)) {
  cat("\n===", target_name, "primer coverage ===\n")
  cat(target_name, "coverage ratio:", coverage_summary[[target_name]]$ratio, "\n")
  print(coverage_summary[[target_name]]$cvg.stats)
  cat(target_name, "coverage plot saved:", coverage_summary[[target_name]]$plot_file, "\n")
  cat(target_name, "coverage table saved:", saved_stats[[target_name]], "\n")
}
