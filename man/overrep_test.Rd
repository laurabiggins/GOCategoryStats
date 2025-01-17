% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overrepresentation_test.r
\name{overrep_test}
\alias{overrep_test}
\title{overrep_test}
\usage{
overrep_test(categories, query_genes, bg_genes, min_query = 3,
  pval_threshold = 0.05, ease = TRUE, sig_digits = 4,
  mult_test = TRUE, super_strict = FALSE)
}
\arguments{
\item{categories}{list of named character vectors containing the functional groups.
Each vector should contain gene names or IDs. The name of each vector should
be the functional category.}

\item{query_genes}{character vector of gene names}

\item{bg_genes}{character vector of gene names to use as the background set
containing genes and associated information including location of genes}

\item{min_query}{minimum number of query genes in the category for it be tested,
no point having a category with one query gene in it}

\item{pval_threshold}{p value threshold. Only results with p-value/corrected p-value
less than this thrreshold will be returned.}

\item{mult_test}{apply multiple testing correction (Benjamini-Hochberg FDR is used)
and use the corrected value to filter for significant results.
This should usually be set to TRUE (default). If set to false, the correction is
still applied but the uncorrected pvalue is used to filter by.}

\item{super_strict}{stricter pvalue correction where it takes the number of
tests as being the total number of functional categories. By default the number
of tests corrected for is only the number of functional categories that contain

min_query genes}
}
\value{
results of functional overrepresentation test. If no categories have a
p-value <= pval_threshold a NULL object will be returned.
}
\description{
Overrepresentation test - Functional (usually gene ontology) analysis
}
\examples{
go_results <- overrep_test(all_go_categories, genes_1, bg_genes)
head(go_results)
}
