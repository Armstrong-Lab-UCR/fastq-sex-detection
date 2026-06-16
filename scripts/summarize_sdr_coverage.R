library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

coverage_dir <- args[1]
samplesheet_file <- args[2]
output_csv <- args[3]
output_plot <- args[4]

samplesheet <- read.csv(samplesheet_file, header = TRUE)
names(samplesheet) <- tolower(names(samplesheet))

if (!"sample" %in% names(samplesheet)) {
  stop("samplesheet must contain a 'sample' column")
}

# Optional sex column
sex_col <- intersect(names(samplesheet), c("sample_sex", "sex", "metadata_sex"))

if (length(sex_col) > 0) {
  sex_col <- sex_col[1]
  metadata <- samplesheet %>%
    dplyr::select(sample, Sample_Sex = all_of(sex_col))
} else {
  metadata <- samplesheet %>%
    dplyr::select(sample) %>%
    dplyr::mutate(Sample_Sex = "Unknown")
}

cov_files <- list.files(
  coverage_dir,
  pattern = "\\.sdr\\.cov$",
  full.names = TRUE
)

if (length(cov_files) == 0) {
  stop("No .sdr.cov files found in coverage directory")
}

sdr_cov <- purrr::map_dfr(cov_files, function(file) {

  sample_name <- basename(file) %>%
    stringr::str_remove("\\.sdr\\.cov$")

  x <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    comment.char = "",
    check.names = FALSE
  )

  names(x) <- gsub("^#", "", names(x))

  x %>%
    dplyr::mutate(sample = sample_name)
})

per_ind <- sdr_cov %>%
  dplyr::mutate(
    SDR_Coverage = as.numeric(coverage),
    SDR_Depth = as.numeric(meandepth)
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    mean_SDR_Coverage = mean(SDR_Coverage, na.rm = TRUE),
    mean_SDR_Depth = mean(SDR_Depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(metadata, by = "sample")

write.csv(per_ind, output_csv, row.names = FALSE)

p <- ggplot(
  per_ind,
  aes(
    x = mean_SDR_Coverage,
    y = mean_SDR_Depth,
    color = Sample_Sex
  )
) +
  geom_point(size = 1.5, alpha = 0.85) +
  labs(
    x = "Mean SDR Coverage",
    y = "Mean SDR Depth",
    color = "Sample Sex"
  ) +
  theme_bw()

ggsave(
  filename = output_plot,
  plot = p,
  width = 6,
  height = 5
)
