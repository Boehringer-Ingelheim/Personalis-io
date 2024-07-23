#' Remove thousands separator (',') in numbers
#' and convert to numeric.
#' @param col a character vector
#' @return numeric vector
#' @keywords internal
fix_thousands_separator <- function(col) {
  as.numeric(sub(",", "", col))
}

#' Catch file not found errors
#'
#' If an error occurs:
#'   - return NULL if it is due to a "path does not exist" error
#'   - stop otherwise
#'
#' @param path path to read
#' @param callback IO function to use
#' @param ... additional aruments passed to callback
#' @return Returns the value of `callback(path, ...)` or NULL if a file not found error occured.
#' @keywords internal
catch_file_not_found <- function(path, callback, ...) {
  tryCatch(
    {
      callback(path, ...)
    },
    error = function(e) {
      if (
        # List of erorr messages to catch
        grepl("`path` does not exist: .*", conditionMessage(e)) || # from read_excel
          grepl("'.*?' does not exist.", conditionMessage(e)) || # from read_html
          grepl("cannot open the connection", conditionMessage(e)) # from read.vcfR
      ) {
        return(NULL)
      } else {
        stop(e)
      }
    }
  )
}

#' Print a warning message for each missing sample
#'
#' @param sample_list List of samples. Missing samples are represented by NULL
#' @param sample_path List of filenames from which the samples were read in (same order as sample_list)
#' @param mod_name Name of the data modality to be used for the error message. e.g. "gene expression"
#' @return a boolean vector indicating which elements in `sample_list` are NULL
#' @keywords internal
warn_missing_samples <- function(sample_list, sample_paths, mod_name) {
  failed_sample_mask <- vapply(sample_list, is.null, logical(1))
  lapply(sample_paths[failed_sample_mask], function(path) {
    message(sprintf("\u26A0\uFE0F No gene %s found for %s. Sample ignored.", mod_name, basename(path)))
  })
  if (sum(failed_sample_mask)) {
    warning(sprintf("%i %s sample(s) were not found. Please check the log for details.", sum(failed_sample_mask), mod_name))
  }
  return(failed_sample_mask)
}

#' Apply IO function for single sample to multiple samples
#'
#' Ignores samples that cannot be found and prints a warning message for missing samples.
#'
#' @param sample_paths {chracter} List of paths to samples
#' @param io_func {function} Callback function to apply to a single directory
#' @param description {character} Name of the modality that is read in (e.g. DNA). Used in the warning message
#'   for missing samples.
#'
#' @return {list} List with one element per successfully read sample. Each element is a return value of `io_func`.
#' @keywords internal
read_samples <- function(sample_paths, io_func, description, ...) {
  sample_list <- pbapply::pblapply(
    sample_paths,
    catch_file_not_found,
    io_func,
    ...
  )
  failed_sample_mask <- warn_missing_samples(sample_list, sample_paths, description)
  sample_list[!failed_sample_mask]
}


#' Ensure that every col in col data has an entry in the bumpy matrix
#'
#' When splitting data frames into a bumpy matrix, samples that don't have a row in the
#' data frame, won't result in a column of the bumpy matrix. This is problematic,
#' as it a sample might have been successfully measured, but simply not resulted in a call of a
#' CNV, mutation, etc.
#'
#' Loosing the sample would result in loss of information:
#'  * we loose the information that there is "no mutation"
#'  * there might be metadata in colData that would be lost
#'
#' While this is not an ideal solution, this function provides a workaround by adding
#' rows that are all NAs for the missing samples. When converting the bumpy matrix
#' back to a data frame, the dummy entries are lost, so we end up with the original
#' data frame, which is nice!.
#'
#' @param df {data.frame} data frame that will be converted into a BumpyMatrix
#' @param col_data {data.frame} data frame that is used as colData (must have rownames that are sample identifiers!)
#' @param sample_col {character} column in `df` that contains the sample identifier
#' @return {tibble} new data frame with dummy entries added
#' @importFrom tibble as_tibble
#' @keywords internal
add_dummy_entry <- function(df, col_data, sample_col = "sample") {
  missing_samples <- setdiff(rownames(col_data), unique(df[[sample_col]]))
  dummy_entries <- lapply(
    missing_samples, function(sample) {
      # List of NAs for each column in `df`
      dummy_list <- sapply(colnames(df), \(x) NA, USE.NAMES = TRUE, simplify = FALSE)
      dummy_list[[sample_col]] <- sample
      as_tibble(dummy_list)
    }
  )
  dplyr::bind_rows(
    df,
    dummy_entries
  )
}

#' Parse VCF files for a provided path and construct data frame.
#'
#' @param path path to VCF file in `*.vcf` or `*.vcf.gz` format
#' @return {tibble} new data frame with all variants (fixed field and genotype information)
#' @importFrom dplyr mutate left_join
#' @importFrom vcfR read.vcfR vcfR2tidy
#' @importFrom stringr str_split_i str_replace
#' @importFrom tibble as_tibble
parse_vcf_to_df <- function(path) {
  # parse VCF file
  vcf_content <- tryCatch({
    read.vcfR(path) 
  }, error = function(e) {
      read.vcfR(str_replace(path, "vcf.gz", "vcf"))
    }
  )
  
  # fixed field content to data frame
  fixed_df <- vcfR2tidy(vcf_content)$fix

  # GT content to data frame
  gt_df <- vcfR2tidy(vcf_content)$gt
  
  # create addition column with observed nucleotides in order to avoid collisions when we do the left_join
  gt_df <- gt_df |>
    dplyr::mutate(ALT = str_split_i(gt_GT_alleles, "/", 2))
  
  # next use ChromKey, POS and ALT for joining vcf content data frames
  joined_vcf_df <- fixed_df |>
    dplyr::left_join(gt_df, by = c("ChromKey", "POS", "ALT"))

  as_tibble(joined_vcf_df)
}
