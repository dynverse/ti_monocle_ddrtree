#!/usr/local/bin/Rscript

task <- dyncli::main()
task <- dyncli::main(
  c("--dataset", "/code/example.h5", "--output", "./output.h5"),
  "/code/definition.yml"
)

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(monocle)

#   ____________________________________________________________________________
#   Load data                                                               ####

params <- task$params
counts <- as.matrix(task$counts)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# just in case
if (is.factor(params$norm_method)) {
  params$norm_method <- as.character(params$norm_method)
}

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# load in the new dataset
pd <- Biobase::AnnotatedDataFrame(data.frame(row.names = rownames(counts)))
fd <- Biobase::AnnotatedDataFrame(data.frame(row.names = colnames(counts), gene_short_name = colnames(counts)))
cds <- monocle::newCellDataSet(t(counts), pd, fd)

# estimate size factors and dispersions
cds <- BiocGenerics::estimateSizeFactors(cds)
cds <- BiocGenerics::estimateDispersions(cds)

# filter features if requested
if (params$filter_features) {
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= params$filter_features_mean_expression)
  cds <- setOrderingFilter(cds, ordering_genes)

  print(nrow(ordering_genes))
}

# if low # cells or features, do not select the number of cluster centers automatically -> https://github.com/cole-trapnell-lab/monocle-release/issues/26
# this avoids the error "initial centers are not distinct."
if (ncol(counts) < 500 || nrow(counts) < 500) {
  params$auto_param_selection <- FALSE
}

# reduce the dimensionality
cds <- monocle::reduceDimension(
  cds,
  max_components = params$max_components,
  reduction_method = params$reduction_method,
  norm_method = params$norm_method,
  auto_param_selection = params$auto_param_selection
)

# order the cells
cds <- monocle::orderCells(cds)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# extract the igraph and which cells are on the trajectory
gr <- cds@auxOrderingData[[params$reduction_method]]$cell_ordering_tree
to_keep <- setNames(rep(TRUE, nrow(counts)), rownames(counts))

# convert to milestone representation
cell_graph <- igraph::as_data_frame(gr, "edges") %>% mutate(directed = FALSE)

if ("weight" %in% colnames(cell_graph)) {
  cell_graph <- cell_graph %>% rename(length = weight)
} else {
  cell_graph <- cell_graph %>% mutate(length = 1)
}

cell_graph <- cell_graph %>% select(from, to, length, directed)

dimred <- t(cds@reducedDimS)
colnames(dimred) <- paste0("Comp", seq_len(ncol(dimred)))

#   ____________________________________________________________________________
#   Save output                                                             ####

output <- dynwrap::wrap_data(cell_ids = rownames(dimred)) %>%
  dynwrap::add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep
  ) %>%
  dynwrap::add_dimred(dimred) %>%
  dynwrap::add_timings(checkpoints)

dyncli::write_output(output, task$output)
