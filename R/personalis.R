# Contains functions to read data from Personalis folders

GUESS_MAX <- 21474836 # largest value without warning - this is anyway more than an excel sheet can have rows

#' Read a set of personalis results folders into a MultiAssayExperiment object
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @param pathlist List with paths to Personalis results folders
#' @param report_type Type of small variant report that should be read (default: preferred)
#' @param read_vcf_data Boolean to specify if VCF file should be read in as well (default: FALSE)
#' @return MultiAssayExperiment
#' @export
read_personalis <- function(pathlist, report_type = "preferred", read_vcf_data = FALSE) {
  stopifnot("`report_type` must be one of 'preferred' or 'all'." = report_type %in% c("preferred", "all"))

  experiments <- list()
  message(">>> Reading gene expression...")
  experiments[["gene_expression"]] <- read_personalis_gene_expression(pathlist)

  message(">>> Reading DNA small variant data... ")
  message("Reading DNA small variant data (somatic)... ")
  experiments[["dna_small_variants_somatic"]] <- read_personalis_small_variant_reports(pathlist, "DNA", "somatic", report_type)
  message("Reading DNA small variant data (tumor)... ")
  experiments[["dna_small_variants_tumor"]] <- read_personalis_small_variant_reports(pathlist, "DNA", "tumor", report_type)
  message("Reading DNA small variant data (normal)... ")
  experiments[["dna_small_variants_normal"]] <- read_personalis_small_variant_reports(pathlist, "DNA", "normal", report_type)

  if (read_vcf_data) {
    message("Reading DNA small variant data (somatic) from VCF files... ")
    experiments[["dna_small_variants_somatic_vcf"]] <- read_personalis_vcf_files(pathlist, "DNA", "somatic")
    message("Reading DNA small variant data (tumor) from VCF files... ")
    experiments[["dna_small_variants_tumor_vcf"]] <- read_personalis_vcf_files(pathlist, "DNA", "tumor")
    message("Reading DNA small variant data (normal) from VCF files... ")
    experiments[["dna_small_variants_normal_vcf"]] <- read_personalis_vcf_files(pathlist, "DNA", "normal")
  }

  message(">>> Reading RNA small variant data... ")
  message("Reading RNA small variant data (somatic)... ")
  experiments[["rna_small_variants_somatic"]] <- read_personalis_small_variant_reports(pathlist, "RNA", "somatic", report_type)
  message("Reading RNA small variant data (tumor)... ")
  experiments[["rna_small_variants_tumor"]] <- read_personalis_small_variant_reports(pathlist, "RNA", "tumor", report_type)
  message("Reading RNA small variant data (normal)... ")
  experiments[["rna_small_variants_normal"]] <- read_personalis_small_variant_reports(pathlist, "RNA", "normal", report_type)

  # for RNA we only have tumor or somatic variant VCF data
  if (read_vcf_data) {
    message("Reading RNA small variant data (somatic) from VCF files... ")
    experiments[["rna_small_variants_somatic_vcf"]] <- read_personalis_vcf_files(pathlist, "RNA", "somatic")
    message("Reading RNA small variant data (tumor) from VCF files... ")
    experiments[["rna_small_variants_tumor_vcf"]] <- read_personalis_vcf_files(pathlist, "RNA", "tumor")
  }

  message(">>> Reading microsatellite stability (MSI) data... ")
  experiments[["msi"]] <- read_personalis_msi_reports(pathlist)

  message(">>> Reading copy number variation (CNV) data... ")
  experiments[["cnv"]] <- read_personalis_cnv_reports(pathlist)

  message(">>> Reading HLA types... ")
  experiments[["hla"]] <- read_personalis_hla_reports(pathlist)

  message(">>> Reading TCR clone reports... ")
  experiments[["tcr"]] <- read_personalis_tcr_reports(pathlist)

  message(">>> Reading BCR clone reports... ")
  experiments[["bcr"]] <- read_personalis_bcr_reports(pathlist)

  mae <- MultiAssayExperiment(
    # remove experiments with no samples (=NULL)
    experiments = experiments[!is.null(experiments)]
  )
  return(mae)
}


#
# -------------- GENE EXPRESSION -------------------
#

#' Read in gene expression data from personalis folders
#' @param sample_paths A vector of paths to personalis folders
#' @return SummarizedExperiment
#' @export
read_personalis_gene_expression <- function(sample_paths) {
  gex_list <- read_samples(sample_paths, read_personalis_gene_expression_sample, "gene_expression")
  do.call(SummarizedExperiment::cbind, gex_list)
}

#' Read a single gene expression sample from a personalis folder
#'
#' @param sample_folder path to personalis folder
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom readxl read_excel
#' @keywords internal
read_personalis_gene_expression_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  gex_path <- file.path(sample_folder, "RNA_Pipeline", "Gene_Expression_Reports", sprintf("RNA_%s_tumor_rna_expression_report.xlsx", sample_name))
  if (!file.exists(gex_path)) {
    # try alternative gene expression folder
    gex_path <- file.path(sample_folder, "RNA_Pipeline", "Expression_Reports", sprintf("RNA_%s_tumor_rna_gene_expression_report.xlsx", sample_name))
  }
  gex_df <- read_excel(gex_path, progress = FALSE) |> tibble::column_to_rownames("NCBI Gene ID")
  tpm_mat <- as.matrix(gex_df[, "TPM", drop = FALSE])
  cpm_mat <- as.matrix(gex_df[, "CPM", drop = FALSE])
  count_mat <- as.matrix(gex_df[, "RNA-Seq Raw Counts", drop = FALSE])
  colnames(tpm_mat) <- sample_name
  colnames(cpm_mat) <- sample_name
  colnames(count_mat) <- sample_name
  row_data <- gex_df[, "Gene Symbol", drop = FALSE]
  se <- SummarizedExperiment(
    assays = list(
      counts = count_mat,
      tpm = tpm_mat,
      normCounts = cpm_mat
    ),
    rowData = row_data
  )
}

#
# -------------------- SMALL VARIANTS --------------------
#

#' Read in small variant data from personalis folder
#' We only read in the "Preferred Transcript" report here:
#'
#' Preferred Transcript: The RefSeq accession.version for the transcript used for variant analysis.
#' Personalis uses a curated list of transcripts, which is based on the number of
#' times a transcript (accession.version) is referred to in COSMIC. If not present in
#' COSMIC, the default transcript would be the one corresponding to the longest CDS.
#' (source: Personalis Analysis_Pipeline_Documentation)
#'
#' We also do not read the `cancer_clinical_research`, the `cancer_research` and the `lowpupulationfreq` reports,
#' because they are subsets of the full report.
#'
#' In addition, Personalis also provides raw VCF files with all (unfiltered) variants. We currently don't
#' read them in because without additional filtering they are not very useful. If you need this level of information,
#' feel free to start from the raw vcf files or even run your own variant calling based on the FASTQ files.
#'
#' @param sample_paths A vector of paths to personalis folders
#' @param sample_type Can be one or multiple of of 'tumor', 'normal', or 'somatic'.
#'   'tumor' refers to tumor sample vs. genome reference (i.e. somatic+germline mutations),
#'   'normal' refers to normal sample vs. genome reference (i.e. germline mutations) and
#'   'somatic' refers to tumor vs. normal (i.e. somatic mutations only).
#' @param modality modality from which the variants were called. Can be either 'DNA' or 'RNA'
#' @param report_type report_type which should be read. Can be either 'preferred' or 'all'
#' @return SummarizedExperiment
#' @importFrom dplyr select distinct bind_rows
#' @importFrom BumpyMatrix splitAsBumpyMatrix
#' @importFrom purrr map
#' @export
read_personalis_small_variant_reports <- function(sample_paths, modality, sample_type, report_type) {
  read_all <- function(sample, ...) {
    list(
      variant_report = read_personalis_small_variant_report_sample(sample, ...),
      summary_stats = read_personalis_somatic_variants_summary_statistics(sample, ...)
    )
  }
  variant_list <- read_samples(
    sample_paths,
    read_all,
    sprintf("%s small variant (%s)", modality, sample_type),
    modality = modality,
    sample_type = sample_type,
    report_type = report_type
  )

  if (!length(variant_list)) {
    return(NULL)
  }

  col_data <- map(variant_list, "summary_stats") |>
    bind_rows() |>
    tibble::column_to_rownames("sample")

  all_variants <- map(variant_list, "variant_report") |> bind_rows()
  row_data <- all_variants |>
    select(
      mut_id,
      Chromosome,
      POS,
      `Variant Type`,
      `Genomic Variant`
    ) |>
    distinct()
  stopifnot("mut_id is not a unique identifier" = !anyDuplicated(row_data$mut_id))

  # ignore columns that are already in `row_data` except for the `mut_id` identifier
  row_data_cols <- colnames(row_data)[colnames(row_data) != "mut_id"]
  all_variants <- all_variants[, !(colnames(all_variants) %in% row_data_cols)]
  all_variants <- add_dummy_entry(all_variants, col_data)
  variants_bm <- splitAsBumpyMatrix(
    all_variants,
    row = all_variants$mut_id,
    column = all_variants$sample,
    sparse = TRUE
  )
  # sort row and col data to be consistent with BM
  row_data <- tibble::column_to_rownames(row_data, "mut_id")[rownames(variants_bm), ]
  col_data <- col_data[colnames(variants_bm), ]

  se <- SummarizedExperiment(
    assays = list(
      variants = variants_bm
    ),
    rowData = row_data,
    colData = col_data
  )
}


#' Read personalis Small Variant reports for a single sample.
#'
#' We only read in the "Preferred Transcript" report here:
#'
#' Preferred Transcript: The RefSeq accession.version for the transcript used for variant analysis.
#' Personalis uses a curated list of transcripts, which is based on the number of
#' times a transcript (accession.version) is referred to in COSMIC. If not present in
#' COSMIC, the default transcript would be the one corresponding to the longest CDS.
#' (source: Personalis Analysis_Pipeline_Documentation)
#'
#' We also do not read the `cancer_clinical_research`, the `cancer_research` and the `lowpupulationfreq` reports,
#' because they are subsets of the full report.
#'
#' In addition, Personalis also provides raw VCF files with all (unfiltered) variants. If you need this level of information,
#' please specify read_vcf_data=TRUE in your read_personalis() call.
#'
#' @param sample_folder path to personalis folder
#' @param sample_types Can be one or multiple of of 'tumor', 'normal', or 'somatic'.
#'   'tumor' refers to tumor sample vs. genome reference (i.e. somatic+germline mutations),
#'   'normal' refers to normal sample vs. genome reference (i.e. germline mutations) and
#'   'somatic' refers to tumor vs. normal (i.e. somatic mutations only).
#' @param modality modality from which the variants were called. Can be either 'DNA' or 'RNA'
#' @param report_type report_type which should be read. Can be either 'preferred' or 'all'
#' @return data frame with all variants
#' @importFrom dplyr mutate rename_with
#' @keywords internal
read_personalis_small_variant_report_sample <- function(sample_folder, modality, sample_type, report_type) {
  stopifnot("`modality` must be one of 'DNA' or 'RNA'." = modality %in% c("DNA", "RNA"))
  stopifnot("`sample_types` must be in 'tumor', 'normal', or 'somatic'." = sample_type %in% c("tumor", "normal", "somatic"))
  stopifnot("`report_type` must be one of 'preferred' or 'full'." = report_type %in% c("preferred", "all"))
  sample_name <- basename(sample_folder)
  # determine report folder based on parameter
  report_folder <- ifelse(report_type == "preferred", "Preferred_Transcripts", "All_Transcripts")

  # the prefix of the normal filename ends with N rather than T, e.g. K13T -> tumor, K13N -> normal
  # even though they are all in the K13T folder.
  tmp_sample_name <- if_else(sample_type == "normal", sub("T$", "N", sample_name), sample_name)
  variant_table <- read_excel(
    file.path(
      sample_folder,
      sprintf("%s_Pipeline", modality),
      "Annotated_SmallVariant_Reports",
      report_folder,
      sprintf("%s_%s_%s_%s_small_variant_report_%s.xlsx", modality, tmp_sample_name, sample_type, tolower(modality), report_type)
    ),
    guess_max = GUESS_MAX
  ) |>
    mutate(sample = sample_name) |>
    # in older versions, the "Chromosome" column is called "Sequence"
    rename_with(\(x) if_else(x == "Sequence", "Chromosome", x)) |>
    mutate(mut_id = sprintf("%s_%s_%s", Chromosome, `Genomic Variant`, `Variant Type`)) |>
    mutate(`Variant ID` = as.character(`Variant ID`)) |>
    mutate(`dbSNP Build` = as.character(`dbSNP Build`))


  variant_table
}

#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr bind_rows
#' @importFrom purrr keep
#' @keywords internal
read_personalis_somatic_variants_summary_statistics <- function(sample_folder, modality, sample_type, report_type) {
  stopifnot("`modality` must be one of 'DNA' or 'RNA'." = modality %in% c("DNA", "RNA"))
  sample_name <- basename(sample_folder)

  # Return empty tibble for sample types other than somatic. For now we only have metrics for somatic.
  if (sample_type != "somatic") {
    return(tibble(sample = sample_name))
  }

  html_file <- file.path(
    sample_folder,
    "QC_REPORT",
    sprintf("%s_%s_%s_statistics.html", modality, sample_name, tolower(modality))
  )
  tables <- read_html(html_file) |>
    html_elements("#somatic_variant_annotation") |>
    html_elements("table") |>
    html_table(na.strings = "N/A")
  tables |>
    # some reports contain two such tables, some only one
    keep(\(x) "SNVs" %in% colnames(x)) |>
    lapply(function(df) {
      colnames(df) <- make.names(colnames(df))
      colnames(df)[1] <- "metric"
      df |>
        mutate(across(c(SNVs, Indels, Total), fix_thousands_separator))
    }) |>
    bind_rows() |>
    pivot_longer(-metric, names_to = "mut_type", values_to = "value") |>
    mutate(sample = sample_name) |>
    mutate(var_name = sprintf("%s (%s)", metric, mut_type)) |>
    select(sample, var_name, value) |>
    pivot_wider(id_cols = sample, names_from = "var_name", values_from = "value")
}

#' Read HTML report on variant calling summary statistics and store information
#' based on given modality and sample_type:
#' DNA -> Normal DNA or Tumor DNA "Summary Small Variants" table
#' RNA -> Tumor RNA "Summary Small Variants" table
#' somatic -> "Tumor/Normal Concordance" table
#'
#' @param sample_folder path to personalis folder
#' @param sample_type Can be one or multiple of of 'tumor', 'normal', or 'somatic'.
#'   'tumor' refers to tumor sample vs. genome reference (i.e. somatic+germline mutations),
#'   'normal' refers to normal sample vs. genome reference (i.e. germline mutations) and
#'   'somatic' refers to tumor vs. normal (i.e. somatic mutations only).
#' @param modality modality from which the variants were called. Can be either 'DNA' or 'RNA'
#' @importFrom rvest read_html html_elements html_nodes html_text
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_to_title
#' @importFrom dplyr all_of across contains
read_personalis_variant_calling_summary_statistics <- function(sample_folder, modality, sample_type) {
  stopifnot("`modality` must be one of 'DNA' or 'RNA'." = modality %in% c("DNA", "RNA"))
  stopifnot("`sample_type` must be one of 'tumor', 'normal', or 'somatic'." = sample_type %in% c("tumor", "normal", "somatic"))
  sample_name <- basename(sample_folder)

  html_file <- file.path(
    sample_folder,
    "QC_REPORT",
    sprintf("%s_%s_%s_statistics.html", modality, sample_name, tolower(modality))
  )

  # decide based on modality and sample_type which content to read from html
  html_section <- if_else(sample_type == "somatic", "#concordance", sprintf("#%s_%s", str_to_title(sample_type), modality))
  table_number <- if_else(sample_type == "somatic", 1, 2)
  columns_to_fix <- if (sample_type == "somatic") c() else c("SNVs", "Indels", "Total")

  tables <- read_html(html_file) |>
    html_elements(html_section) |>
    html_elements("table") |>
    html_table(na.strings = "N/A")

  if (!length(tables)) {
    return(tibble())
  } else {
    tes <- tables[table_number] |>
      lapply(function(df) {
        colnames(df) <- make.names(colnames(df))
        colnames(df)[1] <- "metric"
        df |>
          mutate(across(all_of(columns_to_fix), fix_thousands_separator))
      }) |>
      bind_rows() |>
      pivot_longer(-metric, names_to = "mut_type", values_to = "value") |>
      mutate(sample = sample_name) |>
      mutate(var_name = sprintf("%s (%s)", metric, mut_type)) |>
      select(sample, var_name, value) |>
      pivot_wider(id_cols = sample, names_from = "var_name", values_from = "value") |>
      mutate(across(contains("Number"), fix_thousands_separator))
  }
}

#' Read in (unfiltered) small variant data from VCF files from personalis folders
#' @param sample_paths A vector of paths to personalis folders
#' @param sample_type Can be one or multiple of of 'tumor', 'normal', or 'somatic'.
#'   'tumor' refers to tumor sample vs. genome reference (i.e. somatic+germline mutations),
#'   'normal' refers to normal sample vs. genome reference (i.e. germline mutations) and
#'   'somatic' refers to tumor vs. normal (i.e. somatic mutations only).
#' @param modality modality from which the variants were called. Can be either 'DNA' or 'RNA'
#' @return SummarizedExperiment
#' @importFrom BumpyMatrix splitAsBumpyMatrix
#' @export
read_personalis_vcf_files <- function(sample_paths, modality, sample_type) {
  stopifnot("`modality` must be one of 'DNA' or 'RNA'." = modality %in% c("DNA", "RNA"))

  read_all <- function(sample, ...) {
    list(
      vcf_data = read_personalis_vcf_files_sample(sample, ...),
      summary_stats = read_personalis_variant_calling_summary_statistics(sample, ...)
    )
  }

  variant_list <- read_samples(
    sample_paths,
    read_all,
    sprintf("%s small variant VCF data (%s)", modality, sample_type),
    modality = modality,
    sample_type = sample_type
  )
  if (!length(variant_list)) {
    return(NULL)
  }

  col_data <- bind_rows(map(variant_list, "summary_stats"))
  if (nrow(col_data)) {
    col_data <- col_data |>
      tibble::column_to_rownames("sample")
  }

  all_variants <- map(variant_list, "vcf_data") |> bind_rows()
  row_data <- all_variants |>
    select(
      mut_id,
      CHROM,
      POS,
      REF,
      ALT
    ) |>
    distinct()
  stopifnot("mut_id is not a unique identifier" = !anyDuplicated(row_data$mut_id))

  # ignore columns that are already in `row_data` except for the `mut_id` identifier
  row_data_cols <- colnames(row_data)[colnames(row_data) != "mut_id"]
  all_variants <- all_variants[, !(colnames(all_variants) %in% row_data_cols)]
  all_variants <- add_dummy_entry(all_variants, col_data)
  variants_bm <- splitAsBumpyMatrix(
    all_variants,
    row = all_variants$mut_id,
    column = all_variants$sample,
    sparse = TRUE
  )
  # sort row and col data to be consistent with BM
  row_data <- tibble::column_to_rownames(row_data, "mut_id")[rownames(variants_bm), ]
  col_data <- col_data[colnames(variants_bm), ]

  se <- SummarizedExperiment(
    assays = list(
      variants = variants_bm
    ),
    rowData = row_data,
    colData = col_data
  )
}

#' Read Personalis (raw) VCF files for a provided folders/multiple samples and merge to one table.
#'
#' @param sample_folder base path to personalis folders
#' @param sample_type Can be one or multiple of of 'tumor', 'normal', or 'somatic'.
#'   'tumor' refers to tumor sample vs. genome reference (i.e. somatic+germline mutations),
#'   'normal' refers to normal sample vs. genome reference (i.e. germline mutations) and
#'   'somatic' refers to tumor vs. normal (i.e. somatic mutations only).
#' @param modality modality from which the variants were called. Can be either 'DNA' or 'RNA'
#' @return data frame with all variants
read_personalis_vcf_files_sample <- function(sample_folder, modality, sample_type) {
  stopifnot("`modality` must be one of 'DNA'/'dna' or 'RNA'/'rna'." = modality %in% c("DNA", "RNA", "dna", "rna"))
  stopifnot("`sample_type` must be in 'tumor', 'normal', or 'somatic'." = sample_type %in% c("tumor", "normal", "somatic"))
  sample_name <- basename(sample_folder)
  # the prefix of the normal filename ends with N rather than T, e.g. K13T -> tumor, K13N -> normal
  # even though they are all in the K13T folder.
  tmp_sample_name <- if_else(sample_type == "normal", sub("T$", "N", sample_name), sample_name)

  variant_table <- parse_vcf_to_df(
    file.path(
      sample_folder,
      sprintf("%s_Pipeline", modality),
      "Variants",
      sprintf("%s_%s_%s_%s.%s", modality, tmp_sample_name, sample_type, tolower(modality), "vcf.gz")
    )
  )

  if (nrow(variant_table)) {
    variant_table <- variant_table |>
      mutate(sample = sample_name) |>
      mutate(mut_id = sprintf("%s_%s_%s_%s", CHROM, POS, REF, ALT))
  }

  variant_table
}


#
# -------------------- CNV --------------------------
#

#' Read Personalis CNV data for a list of samples
#' @return SummarizedExperiment
#' @param sample_paths List of directories with Personalis samples
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom BumpyMatrix splitAsBumpyMatrix
#' @export
read_personalis_cnv_reports <- function(sample_paths) {
  read_all <- function(path) {
    list(
      cnv_report = read_personalis_cnv_report_sample(path),
      summary_stats = read_personalis_cnv_summary_statistics(path)
    )
  }
  cnv_list <- read_samples(sample_paths, read_all, "CNV")

  if (!length(cnv_list)) {
    return(NULL)
  }

  col_data <- bind_rows(map(cnv_list, "summary_stats"))
  if (nrow(col_data)) {
    col_data <- col_data |> tibble::column_to_rownames("sample")
  }

  all_cnv <- bind_rows(map(cnv_list, "cnv_report"))
  row_data <- all_cnv |>
    select(cnv_id, `Gene Symbol`, `Chromosome`, `Segment Start`, `Segment End`) |>
    distinct()
  stopifnot("cnv_id is not a unique identifier" = !any(duplicated(row_data$cnv_id)))

  # ignore columns that are already in `row_data` except for the `mut_id` identifier
  row_data_cols <- colnames(row_data)[colnames(row_data) != "cnv_id"]
  all_cnv <- all_cnv[, !(colnames(all_cnv) %in% row_data_cols)]

  all_cnv <- add_dummy_entry(all_cnv, col_data)
  cnv_bm <- splitAsBumpyMatrix(
    all_cnv,
    row = all_cnv$cnv_id,
    column = all_cnv$sample,
    sparse = TRUE
  )
  col_data <- col_data[colnames(cnv_bm), ]
  row_data <- tibble::column_to_rownames(row_data, "cnv_id")[rownames((cnv_bm)), ]

  se <- SummarizedExperiment(
    assays = list(cnv = cnv_bm),
    rowData = row_data,
    colData = col_data
  )
}

#' Read personalis CNV data for a single sample
#' @return data frame with all CNV
#' @importFrom dplyr any_of across cur_column bind_rows
#' @importFrom purrr keep
#' @keywords internal
read_personalis_cnv_report_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  cnv_file <- file.path(
    sample_folder,
    "DNA_Pipeline",
    "Annotated_CopyNumber_Reports",
    sprintf("DNA_%s_somatic_dna_gene_cna_report.xlsx", sample_name)
  )
  COL_TYPES <- list(
    "NCBI Gene ID" = as.numeric,
    "Gene Symbol" = as.character,
    "CNA Type" = as.character,
    "AbsoluteCN" = as.numeric,
    "Chromosome" = as.character,
    "Segment Start" = as.numeric,
    "Segment End" = as.numeric,
    # "Estimated Sample purity" = as.numeric,
    # "Estimated Sample Ploidy" = as.numeric,
    "Percent of Gene in Event" = \(x) as.numeric(sub("%", "", x))
  )
  suppressWarnings({
    cnv_table_list <- list(
      # `skip` raises an error when the table is empty. Workaround:
      # `guess` and then remove the unnecessary columns
      # we also can't specify the columns at import time, because in some personalis versions, some columns
      # are omitted.
      amp = read_excel(cnv_file, sheet = "AMP", col_types = NULL) |>
        # In older reports the "Chromosome" column is called sequence
        rename_with(\(x) if_else(x == "Sequence", "Chromosome", x)) |>
        select(-any_of(c("log posterior probability", "B-allele Frequency", "Allelotype", "Mean_log2Ratio"))) |>
        mutate(across(names(COL_TYPES), \(x) COL_TYPES[[cur_column()]](x))),
      del = read_excel(cnv_file, sheet = "DEL", col_types = NULL) |>
        # In older reports the "Chromosome" column is called sequence
        rename_with(\(x) if_else(x == "Sequence", "Chromosome", x)) |>
        select(-any_of(c("Wilcoxon pvalue", "KS pvalue"))) |>
        mutate(across(names(COL_TYPES), \(x) COL_TYPES[[cur_column()]](x)))
    )
  })

  # it seems some types get inferred incorrectly when a table is empty, therefore
  # only concatenate tables that have at least one row
  cnv_table <- bind_rows(keep(cnv_table_list, \(x) nrow(x) > 0))
  if (nrow(cnv_table)) {
    cnv_table <- cnv_table |>
      mutate(sample = sample_name) |>
      # if a segment spans multiple genes, there will be multiple rows per gene
      mutate(cnv_id = sprintf("%s_%i_%i_%s", Chromosome, `Segment Start`, `Segment End`, `Gene Symbol`))
  }

  cnv_table
}

#' @importFrom rvest read_html html_elements html_nodes html_text
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @keywords internal
read_personalis_cnv_summary_statistics <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  html_file <- file.path(
    sample_folder,
    "QC_REPORT",
    sprintf("DNA_%s_dna_statistics.html", sample_name)
  )
  # unfortunately, this is not a table, but a div of divs that looks like a table
  table <- (read_html(html_file) |> html_elements("#copy_number"))
  # unfortunately, it seems missing in newer versions of the report
  if (!length(table)) {
    return(tibble())
  } else {
    table <- table[[1]]
    titles <- table |>
      html_nodes(".title") |>
      html_text()
    values <- table |>
      html_nodes(".value") |>
      html_text()

    cnv_metrics <- tibble(metric = titles[1:5], value = values[1:5]) |>
      mutate(sample = sample_name) |>
      pivot_wider(id_cols = sample, names_from = metric, values_from = value)
    cnv_metrics
  }
}


#
# -------------------- HLA types --------------------------
#

#' Read Personalis HLA types for a list of samples
#' @param sample_paths List of paths to read
#' @return SummarizedExperiment
#' @export
read_personalis_hla_reports <- function(sample_paths) {
  hla_list <- read_samples(sample_paths, read_personalis_hla_sample, "HLA")
  do.call(SummarizedExperiment::cbind, hla_list)
}

#' Read HLA types for a single sample
#' @return data frame with one row per HLA gene per Allele.
#' @keywords internal
read_personalis_hla_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  # The HLAs are inferred on the "normal" sample, therefore the final T needs to be replaced with N
  tmp_sample_name <- sub("T$", "N", sample_name)
  hla_file <- file.path(
    sample_folder,
    "HLA",
    sprintf("DNA_%s_hla.xlsx", tmp_sample_name)
  )
  hla_tab <- read_excel(hla_file)
  # Put the data in one column (treating Allele 1 and 2 as separate variables)
  hla_tab <- dplyr::bind_rows(
    hla_tab |> select(hla_gene = Gene, hla_type = Allele1) |> mutate(locus = sprintf("%s.1", hla_gene), allele = 1),
    hla_tab |> select(hla_gene = Gene, hla_type = Allele2) |> mutate(locus = sprintf("%s.2", hla_gene), allele = 2)
  )
  row_data <- hla_tab |>
    select(locus, hla_gene, allele) |>
    tibble::column_to_rownames("locus")
  hla_type_mat <- as.matrix(hla_tab[, "hla_type", drop = FALSE])
  colnames(hla_type_mat) <- sample_name
  se <- SummarizedExperiment(
    assays = list(hla_types = hla_type_mat), rowData = row_data
  )

  se
}


#
# -------------------- TCR --------------------------
#

#' Read Personalis TCR clone report for a list of samples
#' @param sample_paths List of paths to read
#' @return SummarizedExperiment
#' @export
#' @importFrom purrr map
#' @importFrom BumpyMatrix splitAsBumpyMatrix
read_personalis_tcr_reports <- function(sample_paths) {
  read_all <- function(path) {
    list(
      clone_report = read_personalis_tcr_clone_report_sample(path),
      summary_stats = read_personalis_tcr_summary_statistics(path)
    )
  }
  tcr_list <- read_samples(sample_paths, read_all, "TCR")

  col_data <- dplyr::bind_rows(map(tcr_list, "summary_stats"))

  if (nrow(col_data) == 0) {
    return(NULL)
  }

  col_data <- col_data |> tibble::column_to_rownames("sample")

  all_tcr <- dplyr::bind_rows(map(tcr_list, "clone_report"))
  tcr_bm <- splitAsBumpyMatrix(
    all_tcr,
    row = all_tcr$locus,
    column = all_tcr$sample,
    sparse = FALSE
  )
  col_data <- col_data[colnames(tcr_bm), ]
  se <- SummarizedExperiment(
    assays = list(tcr = tcr_bm),
    colData = col_data
  )
}


#' Read TCRs for a single sample
#' @return data frame with one row per TCR clone
#' @importFrom dplyr if_else bind_rows
#' @keywords internal
read_personalis_tcr_clone_report_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  tcr_df <- lapply(c("TRA", "TRB"), function(locus) {
    locus_file_name <- if_else(locus == "TRA", "tcr_alpha", "tcr_beta")

    # Try RNA_TCR_BCR folder first
    tcr_file <- file.path(
      sample_folder,
      "RNA_TCR_BCR",
      "Immune_Repertoire",
      sprintf("RNA_%s_%s_clone_report.xlsx", sample_name, locus_file_name)
    )

    if (!file.exists(tcr_file)) {
      # Try RNA_TCR folder
      tcr_file <- file.path(
        sample_folder,
        "RNA_TCR",
        "Immune_Repertoire",
        sprintf("RNA_%s_rna_%s_clone_report.xlsx", sample_name, locus_file_name)
      )
    }

    if (!file.exists(tcr_file)) {
      # try alternative path in RNA_TCR
      tcr_file <- file.path(
        sample_folder,
        "RNA_TCR",
        "Immune_Repertoire",
        sprintf("RNA_%s_%s_clone_report.xlsx", sample_name, locus_file_name)
      )
    }

    read_excel(tcr_file, guess_max = GUESS_MAX) |>
      mutate(locus = locus, sample = sample_name)
  }) |> bind_rows()

  tcr_df
}

#' Extract summary statistics from the TCR HTML report
#' @return data frame
#' @importFrom rvest read_html html_elements html_table
#' @keywords internal
read_personalis_tcr_summary_statistics <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  html_file <- file.path(
    sample_folder,
    "QC_REPORT",
    sprintf("RNA_%s_tcr_statistics.html", sample_name)
  )
  html <- read_html(html_file)
  clonality_alpha <- (html |> html_elements("#tcr_alpha_overview") |> html_elements("table") |> html_table())[[1]]
  clonality_beta <- (html |> html_elements("#tcr_beta_overview") |> html_elements("table") |> html_table())[[1]]

  tibble(
    sample = sample_name,
    clonality_alpha = clonality_alpha$X2[1],
    clonality_beta = clonality_beta$X2[1],
  )
}


#
# -------------------- BCR --------------------------
#

#' Read BCR clone report for a list of samples
#' @param sample_paths List of paths to read
#' @return SummarizedExperiment
#' @export
#' @importFrom purrr map
#' @importFrom BumpyMatrix splitAsBumpyMatrix
read_personalis_bcr_reports <- function(sample_paths) {
  read_all <- function(path) {
    list(
      clone_report = read_personalis_bcr_clone_report_sample(path),
      isotype_composition = read_personalis_bcr_isotype_composition_sample(path),
      summary_stats = read_personalis_bcr_summary_statistics(path)
    )
  }
  bcr_list <- read_samples(sample_paths, read_all, "BCR")

  # Merge summary stats and isotype composition
  isotype_composition <- dplyr::bind_rows(map(bcr_list, "isotype_composition"))
  summary_stats <- dplyr::bind_rows(map(bcr_list, "summary_stats"))
  col_data <- dplyr::left_join(summary_stats, isotype_composition, by = "sample")

  if (nrow(col_data) == 0) {
    return(NULL)
  }

  col_data <- col_data |> tibble::column_to_rownames("sample")

  all_bcr <- dplyr::bind_rows(map(bcr_list, "clone_report"))
  bcr_bm <- splitAsBumpyMatrix(
    all_bcr,
    row = all_bcr$chain,
    column = all_bcr$sample,
    sparse = FALSE
  )
  col_data <- col_data[colnames(bcr_bm), ]
  se <- SummarizedExperiment(
    assays = list(bcr = bcr_bm),
    colData = col_data
  )
}


#' Read BCRs for a single sample.
#' @return data frame with one row per BCR clone
#' @importFrom dplyr bind_rows
#' @keywords internal
read_personalis_bcr_clone_report_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)

  # Try RNA_TCR_BCR folder first
  bcr_file <- file.path(
    sample_folder,
    "RNA_TCR_BCR",
    "Immune_Repertoire",
    sprintf("RNA_%s_bcr_heavy_chain_report.xlsx", sample_name)
  )

  if (!file.exists(bcr_file)) {
    # Try RNA_BCR folder
    bcr_file <- file.path(
      sample_folder,
      "RNA_BCR",
      "Immune_Repertoire",
      sprintf("RNA_%s_bcr_heavy_chain_report.xlsx", sample_name)
    )
  }

  if (!file.exists(bcr_file)) {
    # Try alternative naming pattern in RNA_TCR_BCR
    bcr_file <- file.path(
      sample_folder,
      "RNA_TCR_BCR",
      "Immune_Repertoire",
      sprintf("RNA_%s_rna_bcr_heavy_chain_report.xlsx", sample_name)
    )
  }

  if (!file.exists(bcr_file)) {
    return(tibble(chain = character(), sample = character()))
  }

  bcr_df <- read_excel(bcr_file, guess_max = GUESS_MAX) |>
    mutate(chain = "IGH", sample = sample_name)

  bcr_df
}

#' Read BCR isotype composition for a single sample
#' @return data frame with isotype composition statistics
#' @importFrom readxl read_excel
#' @importFrom tidyr pivot_wider
#' @keywords internal
read_personalis_bcr_isotype_composition_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)

  # Try RNA_TCR_BCR folder first
  isotype_file <- file.path(
    sample_folder,
    "RNA_TCR_BCR",
    "Immune_Repertoire",
    sprintf("RNA_%s_bcr_isotype_composition.xlsx", sample_name)
  )

  if (!file.exists(isotype_file)) {
    # Try RNA_BCR folder
    isotype_file <- file.path(
      sample_folder,
      "RNA_BCR",
      "Immune_Repertoire",
      sprintf("RNA_%s_bcr_isotype_composition.xlsx", sample_name)
    )
  }

  if (!file.exists(isotype_file)) {
    # Try alternative naming pattern in RNA_TCR_BCR
    isotype_file <- file.path(
      sample_folder,
      "RNA_TCR_BCR",
      "Immune_Repertoire",
      sprintf("RNA_%s_rna_bcr_isotype_composition.xlsx", sample_name)
    )
  }

  if (!file.exists(isotype_file)) {
    return(tibble(sample = sample_name))
  }

  isotype_df <- read_excel(isotype_file, guess_max = GUESS_MAX)

  # Convert isotype composition to wide format for colData
  isotype_wide <- isotype_df |>
    pivot_wider(
      names_from = names(isotype_df)[2], # Second column is isotype name
      values_from = names(isotype_df)[1] # First column is the fraction value
    ) |>
    mutate(sample = sample_name)

  isotype_wide
}

#' Extract summary statistics from the BCR HTML report
#' @return data frame
#' @importFrom rvest read_html html_elements html_table
#' @keywords internal
read_personalis_bcr_summary_statistics <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  html_file <- file.path(
    sample_folder,
    "QC_REPORT",
    sprintf("RNA_%s_bcr_statistics.html", sample_name)
  )
  html <- read_html(html_file)
  clonality_bcrh <- (html |> html_elements("#bcr_overview") |> html_elements("table") |> html_table())[[1]]
  clonality_iga <- (html |> html_elements("#iga_overview") |> html_elements("table") |> html_table())[[1]]
  clonality_igd <- (html |> html_elements("#igd_overview") |> html_elements("table") |> html_table())[[1]]
  clonality_igg <- (html |> html_elements("#igg_overview") |> html_elements("table") |> html_table())[[1]]
  clonality_igm <- (html |> html_elements("#igm_overview") |> html_elements("table") |> html_table())[[1]]

  tibble(
    sample = sample_name,
    clonality_bcrh = clonality_bcrh$X2[1],
    clonality_iga = clonality_iga$X2[1],
    clonality_igd = clonality_igd$X2[1],
    clonality_igg = clonality_igg$X2[1],
    clonality_igm = clonality_igm$X2[1]
  )
}


#
# -------------------- MSI --------------------------
#

#' Read MSI information for a list of samples
#'
#' MSI information is exclusively available form the HTML reports.
#' The information for individual loci will be stored in an assay of a SummarizedExperiment.
#' Summary statistics will be available from colData
#'
#' @param sample_paths Path to personalis samples
#' @return SummarizedExperiment
#' @export
read_personalis_msi_reports <- function(sample_paths) {
  msi_list <- read_samples(sample_paths, read_personalis_msi_sample, "MSI")
  do.call(SummarizedExperiment::cbind, msi_list)
}

#' Extract MSI locus information and summary statistics from a single MSI HTML report
#' @return SummarizedExperiment
#' @keywords internal
read_personalis_msi_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  html_file <- file.path(
    sample_folder,
    "QC_REPORT",
    sprintf("DNA_%s_msi_statistics.html", sample_name)
  )
  html <- read_html(html_file)
  locus_table <- (html |>
    html_elements("#msi_table_by_loc") |>
    html_table())[[1]]
  stats_table <- (html |> html_elements("#msi_table") |> html_table())[[1]]

  locus_mat <- locus_table |>
    mutate(sample = sample_name) |>
    pivot_wider(id_cols = sample, names_from = "Microsatellite locus", values_from = "Stability status") |>
    tibble::column_to_rownames("sample") |>
    as.matrix() |>
    t()

  col_data <- stats_table |>
    mutate(sample = sample_name) |>
    tibble::column_to_rownames("sample")

  se <- SummarizedExperiment(
    assays = list(
      bethesda_consensus_panel = locus_mat
    ),
    colData = col_data
  )
  se
}
