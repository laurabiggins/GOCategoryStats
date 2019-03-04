# https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html

#' Gene set functional categories.
#'
#' A dataset containing genes grouped in to functional sets. Each functional
#' group is named. This dataset is processed from a gmt file.
#' Filtering for reasonable sized categories has been turned off. This set of
#' categories should mainly be used for testing for biases.
#'
#' @format A named list with x number of components
#' @source http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/Mouse_GO_AllPathways_no_GO_iea_February_01_2019_symbol.gmt
#' @details all_go_categories <- process_GMT(file_location, min_genes = 3, max_genes = 500000)
#> [1] "all_go_categories"


#' Gene set functional categories.
#'
#' A dataset containing genes grouped in to functional sets. Each functional
#' group is named. This dataset is processed from a gmt file.
#' The
#'
#'
#' @format A named list with x number of components
#' @source http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/Mouse_GO_AllPathways_no_GO_iea_February_01_2019_symbol.gmt
#' @details all_go_categories <- process_GMT(file_location)
#> [1] "all_go_categories"


#' Set of suspect gene ontology categories
#'
#' These categories are derived from an analysis of public datasets ....
#' Essentially these GO categories appeared more than expected when running GO
#' overrepresentation tests
#'
#' @format Tabled data. The names are the GO categories. The numbers are the
#' proportion of times the category was identified as significant from a set of
#' 373 analyses
"suspects1"


#' Set of chromosome data
#'
#' A dataset derived from 2 sets of genes, containing the proportion of genes that
#' were located on each chromosome. The query set are \code{genes_1} and the
#' background set \code{bg_genes}.
#'
#' @format A data frame containing 2 columns. The row names are the chromosome names
"chromosomes"


#' Set of mouse genes
#'
#' A dataset containing a set of 100 mouse genes.
#'
#' @format A character vector containing 100 gene names
#'
#' @examples
#' x <- head(genes_1)
#' x
"genes_1"


#' Set of background genes
#'
#' A dataset containing a set of background mouse genes.
#' All protein coding genes in a certain release?
#'
#' @format A character vector containing 29125 gene names
#' @examples
#' head(bg_genes)
"bg_genes"


#' Mouse genome annotations
#'
#' A dataset containing gene annotations from Ensembl for the mouse genome.
#'
#' This is a gtf file from Ensembl that has been imported into R and saved as a dataframe.
#' Each line is an annotation for a region in the genome.
#' The data has been filtered to only include the gene annotations. The original
#' file also included transcripts, coding sequences, 3' and 5' UTRs etc.
#'
#' @format A data frame containing information on the Mus musculus genome.
#' @examples
#' head(m_musculus)
#' @source <ftp://ftp.ensembl.org/pub/release-94/gtf/mus_musculus/Mus_musculus.GRCm38.94.gtf.gz>
"m_musculus"




