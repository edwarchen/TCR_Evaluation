library(openPrimeR)
library(Biostrings)

hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")

cfg <- list(
  TRBV = list(
    ref_path = "data/IGHV.fasta",
    primer_path = "data/primer_set/IGHV_primers.fasta",
    direction = "fw",
    primer_id_suffix = "_fw"
  ),
  TRBJ = list(
    ref_path = "data/IGHJ.fasta",
    primer_path = "data/primer_set/IGHJ_primers.fasta",
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
  check_constraints(
    primers,
    ref,
    eval_settings,
    active.constraints = "primer_coverage"
  )
}

xml_file <- system.file("extdata", "settings", "A_Taq_PCR_design.xml", package = "openPrimeR")
eval_settings <- read_settings(xml_file)
conOptions(eval_settings)$allowed_other_binding_ratio <- c(max = 1.0)

results <- lapply(cfg, run_coverage, eval_settings = eval_settings)

cat("=== TRBV primer coverage ===\n")
print(results$TRBV$primer_coverage)
cat("\n=== TRBJ primer coverage ===\n")
print(results$TRBJ$primer_coverage)

