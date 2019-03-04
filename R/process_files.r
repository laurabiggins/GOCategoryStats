#' Import and process GMT file.
#'
#' Import GMT file containing gene sets.
#'
#' The GMT format is described \href{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{GMT format}.
#' The first column should contain gene set names, the 2nd column contains a brief description and the 3rd column contains a whitespace separated list of genes.
#'
#' @param file file path for gmt file containing gene sets
#' @param min_genes minimum number of genes required in the functional category for the category to be loaded.
#' @param max_genes maximum number of genes required in the functional category for the category to be loaded.
#' @return A named list containing all the gene sets imported from \code{gmt_file_path}
#' @examples
#' file <- "http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/Mouse_GO_AllPathways_no_GO_iea_December_01_2018_symbol.gmt"
#' x <- process_GMT(file, 5, 5000)
process_GMT <- function(file, min_genes = 3, max_genes = 5000) {

  gmt_file <- data.table::fread(
    file,
    sep = "\n",
    header = FALSE,
    data.table = FALSE
  )[, 1]

  categories <- strsplit(gmt_file, "\t")

  names(categories) <- sapply(categories, `[[`, 1)

  genes_in_categories <- lapply(categories, `[`, -(1:2))

  genes_in_categories <- sapply(genes_in_categories, toupper)

  category_lengths <- sapply(categories, length)

  genes_in_categories[category_lengths >= min_genes & category_lengths <= max_genes]
}


#' Import GTF file.
#'
#' @param file file path for gtf file
#' @param genes_only if true, only return the entries annotated as "gene"
#' @return A data frame containing all the genes from \code{gtf_file}
#' @examples
#' import_GTF("~/file.gtf")
import_GTF <- function(file, genes_only = TRUE) {

  gtf <- data.table::fread(
    gtf_file,
    header = FALSE,
    select = c(1,3:5,9),
    data.table = FALSE
  )

  colnames(gtf) <- c("chr","type","start","end","info")

  if (!(
    is.character(gtf$type) &
    is.integer(gtf$start) &
    is.integer(gtf$end) &
    is.character(gtf$info)
  )) {
    stop("gtf file format is not as expected")
  }

  if (genes_only) return(gtf[gtf$type == "gene", ])

  gtf
}


#' Process GTF file.
#'
#' parses the info column of a GTF file to extract the gene name and length.
#'
#' @param gtf gtf file that has been imported using \code{import_GTF}.
#' The gtf file must include a column named info that contains data in the following format: a character vector in the
#' format "gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"
#' @param gene_names if true, add a column containing gene name
#' @param gene_length if true, add a column containing gene length
#' @param remove_duplicates remove rows if gene names are duplicated. The row that is kept
#' is the first instance of that gene. To remove duplicates based on something other that
#' gene name, use the separate \code{remove_duplicates} function.
#' @return A data frame containing all the genes from \code{gtf_file}
#' @examples
#' parse_GTF_info(gtf_file)
parse_GTF_info <- function(gtf, gene_names = TRUE, gene_length = TRUE, remove_duplicates = TRUE) {

  if (gene_names) {
    if (is.null(gtf$info)) stop("The gtf file must include a column named 'info'")
    if (sum(grepl("gene_name", x = gtf$info)) == 0) {
      stop("The info column in the gtf file must include gene names when the gene_names argument is set to TRUE,
           in the format specified in the documentation for the function process_GTF")
    }
    if (sum(grepl("gene_name", x = gtf$info)) < nrow(gtf)) {
      n_names <- sum(grepl("gene_name", x = gtf$info))
      warning("a gene name could not be found for every gene: n_names found from nrow(gtf) entries")
    }
    gtf$gene_name <- extract_gene_names(gtf$info)
  }
  if (gene_length) gtf$length <- gtf$end - gtf$start
  ifelse(remove_duplicates, return(remove_duplicates(gtf, column = "gene_name")), return(gtf))
}


#' Extract gene names from GTF
#'
#' Parse info column from GTF file to extract gene names
#'
#' @param info_column character vector from the info column of a gtf file, probably the output from \code{\link{process_GTF}}, in the
#' format "gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"
#' @return A data frame containing the gene sets of interest along with information on
#' gene lengths etc. \url{http://rstudio.com} some more text \link{process_GMT}
#' @examples
#' x <- extract_gene_names(m_musculus$info)
extract_gene_names <- function(info_column) {
  gene_names <- sapply(info_column, function(info) {
    sapply(strsplit(info, split = ";"), function(x) {
      x[grep(pattern = "gene_name", x)]
    })
  })
  gene_names <- gsub(gene_names, pattern = "gene_name", replacement = "")
  # trim white space, apostrophes etc
  clean_text(gene_names, remove_empty = FALSE, remove_dup = FALSE)
}


#' get_GO_sets.
#'
#' gets sets of genes
#'
#' @param selected_go_names names of gene sets of interest
#' @param all_go_categories gene sets that have been imported from \code{\link{process_GMT}}
#' @param gtf dataframe containing genes and information, probably the output from \code{\link{process_GTF}}

get_GO_sets <- function(selected_go_names, all_go_categories) {
  all_go_categories[match(selected_go_names, names(all_go_categories))]
}


#' Get information on selected gene ontology categories.
#'
#' Uses information from gtf file to get lengths of genes etc.
#'
#' @param selected_go_names names of gene sets of interest
#' @param all_go_categories gene sets that have been imported from \code{\link{process_GMT}}
#' @param gtf dataframe containing genes and information, probably the output from \code{\link{process_GTF}}
#' @return A data frame containing the gene sets of interest along with information on
#' gene lengths etc. \url{http://rstudio.com} some more text \link{process_GMT}
#' @examples
#' get_GO_info(selected_go_names, all_go_categories, gtf)
get_GO_info <- function(selected_go_names, all_go_categories, gtf) {
  if (is.null(gtf$gene_name)) gtf <- process_gene_names(gtf)

  go_sets <- get_GO_sets(selected_go_names, all_go_categories)
  go_category_lengths <- sapply(go_sets, length)

  go_info <- data.frame(selected_go_names)

  # add column showing total number of genes in GO category
  go_info$total_genes <- go_category_lengths[(match(selected_go_names, as.vector(names(go_category_lengths))))]

  go_info$genes_found_in_gtf <- lapply(lapply(go_sets, fastmatch::fmatch, gtf$gene_name), function(x) sum(!is.na(x)))
  # It's about 40 times faster using fmatch than normal match, so definitely worth doing
  # I tried this matching a few different ways and this was the fastest way that I tried.
  if (is.null(gtf$length) & (is.null(gtf$start) | is.null(gtf$end))) {
    stop("No length information could be found in the gtf file")
  } else if (is.null(gtf$length) & (!is.null(gtf$start) | !is.null(gtf$end))) {
    gtf$length <- gtf$end - gtf$start
  }

  gene_lengths_for_categories <- lapply(lapply(go_sets, fastmatch::fmatch, gtf$gene_name), function(x) gtf$length[x])

  # it's quicker to do sum(x)/length(x) than mean() but we need to remove the NAs so it's easier to use mean
  go_info$mean_gene_length <- lapply(gene_lengths_for_categories, mean, na.rm = TRUE)
  go_info$median_gene_length <- lapply(gene_lengths_for_categories, median, na.rm = TRUE)

  go_info
}

