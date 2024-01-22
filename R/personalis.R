# Contains functions to read data from Personalis folders

GUESS_MAX <- 21474836 # largest value without warning - this is anyway more than an excel sheet can have rows

#' Read a set of personalis results folders into a MultiAssayExperiment object
#' @param pathlist List with paths to Personalis results folders
#' @return MultiAssayExperiment
#' @export
read_personalis <- function(pathlist) {
  experiments <- list()
  message(">>> Reading gene expression...")
  experiments[["gene_expression"]] <- read_personalis_gene_expression(pathlist)

  message(">>> Reading DNA small variant data... ")
  message("Reading DNA small variant data (somatic)... ")
  experiments[["dna_small_variants_somatic"]] <- read_personalis_small_variant_reports(pathlist, "DNA", "somatic")
  message("Reading DNA small variant data (tumor)... ")
  experiments[["dna_small_variants_tumor"]] <- read_personalis_small_variant_reports(pathlist, "DNA", "tumor")
  message("Reading DNA small variant data (normal)... ")
  experiments[["dna_small_variants_normalr"]] <- read_personalis_small_variant_reports(pathlist, "DNA", "normal")

  message(">>> Reading RNA small variant data... ")
  message("Reading RNA small variant data (somatic)... ")
  experiments[["rna_small_variants_somatic"]] <- read_personalis_small_variant_reports(pathlist, "RNA", "somatic")
  message("Reading RNA small variant data (tumor)... ")
  experiments[["rna_small_variants_tumor"]] <- read_personalis_small_variant_reports(pathlist, "RNA", "tumor")
  message("Reading RNA small variant data (normal)... ")
  experiments[["rna_small_variants_normal"]] <- read_personalis_small_variant_reports(pathlist, "RNA", "normal")

  message(">>> Reading microsatellite stability (MSI) data... ")
  experiments[["msi"]] <- read_personalis_msi_reports(pathlist)

  message(">>> Reading copy number variation (CNV) data... ")
  experiments[["cnv"]] <- read_personalis_cnv_reports(pathlist)

  message(">>> Reading HLA types... ")
  experiments[["hla"]] <- read_personalis_hla_reports(pathlist)

  message(">>> Reading TCR clone reports... ")
  experiments[["tcr"]] <- read_personalis_tcr_reports(pathlist)

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
#' @return SummarizedExperiment
#' @importFrom dplyr select
#' @importFrom purrr map
#' @export
read_personalis_small_variant_reports <- function(sample_paths, modality, sample_type) {
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
    sample_type = sample_type
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
      Sequence,
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
#' In addition, Personalis also provides raw VCF files with all (unfiltered) variants. We currently don't
#' read them in because without additional filtering they are not very useful. If you need this level of information,
#' feel free to start from the raw vcf files or even run your own variant calling based on the FASTQ files.
#'
#' @param sample_folder path to personalis folder
#' @param sample_types Can be ome or multiple of of 'tumor', 'normal', or 'somatic'.
#'   'tumor' refers to tumor sample vs. genome reference (i.e. somatic+germline mutations),
#'   'normal' refers to normal sample vs. genome reference (i.e. germline mutations) and
#'   'somatic' refers to tumor vs. normal (i.e. somatic mutations only).
#' @param modality modality from which the variants were called. Can be either 'DNA' or 'RNA'
#' @return data frame with all variants
#' @importFrom dplyr mutate
#' @keywords internal
read_personalis_small_variant_report_sample <- function(sample_folder, modality, sample_type) {
  stopifnot("`modality` must be one of 'DNA' or 'RNA'." = modality %in% c("DNA", "RNA"))
  stopifnot("`sample_types` must be in 'tumor', 'normal', or 'somatic'." = sample_type %in% c("tumor", "normal", "somatic"))
  sample_name <- basename(sample_folder)
  # the prefix of the normal filename ends with N rather than T, e.g. K13T -> tumor, K13N -> normal
  # even though they are all in the K13T folder.
  tmp_sample_name <- if_else(sample_type == "normal", sub("T$", "N", sample_name), sample_name)
  variant_table <- read_excel(
    file.path(
      sample_folder,
      sprintf("%s_Pipeline", modality),
      "Annotated_SmallVariant_Reports",
      "Preferred_Transcripts",
      sprintf("%s_%s_%s_%s_small_variant_report_preferred.xlsx", modality, tmp_sample_name, sample_type, tolower(modality))
    ),
    guess_max = GUESS_MAX
  ) |>
    mutate(sample = sample_name) |>
    mutate(mut_id = sprintf("%s_%s_%s", Sequence, `Genomic Variant`, `Variant Type`))

  variant_table
}

#'
#' @importFrom tidyr pivot_longer
#' @keywords internal
read_personalis_somatic_variants_summary_statistics <- function(sample_folder, modality, sample_type) {
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
  tables[2:3] |>
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


#
# -------------------- CNV --------------------------
#

#' Read Personalis CNV data for a list of samples
#' @return SummarizedExperiment
#' @param sample_paths List of directories with Personalis samples
#' @importFrom purrr map
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

  col_data <- bind_rows(map(cnv_list, "summary_stats")) |>
    tibble::column_to_rownames("sample")

  all_cnv <- bind_rows(map(cnv_list, "cnv_report"))
  row_data <- all_cnv |>
    select(cnv_id, `Gene Symbol`, `Sequence`, `Segment Start`, `Segment End`) |>
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
    "Sequence" = as.character,
    "Segment Start" = as.numeric,
    "Segment End" = as.numeric,
    "Estimated Sample purity" = as.numeric,
    "Estimated Sample Ploidy" = as.numeric,
    "Percent of Gene in Event" = \(x) as.numeric(sub("%", "", x))
  )
  suppressWarnings({
    cnv_table_list <- list(
      # `skip` raises an error when the table is empty. Workaround:
      # `guess` and then remove the unnecessary columns
      # we also can't specify the columns at import time, because in some personalis versions, some columns
      # are omitted.
      amp = read_excel(cnv_file, sheet = "AMP", col_types = NULL) |>
        select(-any_of(c("log posterior probability", "B-allele Frequency", "Allelotype", "Mean_log2Ratio"))) |>
        mutate(across(names(COL_TYPES), \(x) COL_TYPES[[cur_column()]](x))),
      del = read_excel(cnv_file, sheet = "DEL", col_types = NULL) |>
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
      mutate(cnv_id = sprintf("%s_%i_%i_%s", Sequence, `Segment Start`, `Segment End`, `Gene Symbol`))
  }

  cnv_table
}

#' @importFrom rvest read_html html_elements html_nodes html_text
#' @keywords internal
read_personalis_cnv_summary_statistics <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  html_file <- file.path(
    sample_folder,
    "QC_REPORT",
    sprintf("DNA_%s_dna_statistics.html", sample_name)
  )
  # unfortunately, this is not a table, but a div of divs that looks like a table
  table <- (read_html(html_file) |> html_elements("#copy_number"))[[1]]
  titles <- table %>%
    html_nodes(".title") %>%
    html_text()
  values <- table %>%
    html_nodes(".value") %>%
    html_text()

  cnv_metrics <- tibble(metric = titles[1:5], value = values[1:5]) |>
    mutate(sample = sample_name) |>
    pivot_wider(id_cols = sample, names_from = metric, values_from = value)
  cnv_metrics
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
  hla_tab <- bind_rows(
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
read_personalis_tcr_reports <- function(sample_paths) {
  read_all <- function(path) {
    list(
      clone_report = read_personalis_tcr_clone_report_sample(path),
      summary_stats = read_personalis_tcr_summary_statistics(path)
    )
  }
  tcr_list <- read_samples(sample_paths, read_all, "TCR")

  col_data <- bind_rows(map(tcr_list, "summary_stats"))

  if (nrow(col_data) == 0) {
    return(NULL)
  }

  col_data <- col_data |> tibble::column_to_rownames("sample")

  all_tcr <- bind_rows(map(tcr_list, "clone_report"))
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


#' Read TCRs for a single sample.
#' @return data frame with one row per TCR clone
#' @keywords internal
read_personalis_tcr_clone_report_sample <- function(sample_folder) {
  sample_name <- basename(sample_folder)
  tcr_df <- lapply(c("TRA", "TRB"), function(locus) {
    locus_file_name <- if_else(locus == "TRA", "tcr_alpha", "tcr_beta")
    tcr_file <- file.path(
      sample_folder,
      "RNA_TCR",
      "Immune_Repertoire",
      sprintf("RNA_%s_rna_%s_clone_report.xlsx", sample_name, locus_file_name)
    )
    if (!file.exists(tcr_file)) {
      # try alternative path
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
