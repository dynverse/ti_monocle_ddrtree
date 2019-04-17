#!/usr/local/bin/Rscript

task <- dyncli::main()

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(monocle, warn.conflicts = FALSE)

#   ____________________________________________________________________________
#   Load data                                                               ####

parameters <- task$parameters
counts <- as.matrix(task$counts)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# just in case
if (is.factor(parameters$norm_method)) {
  parameters$norm_method <- as.character(parameters$norm_method)
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
if (parameters$filter_features) {
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= parameters$filter_features_mean_expression)
  cds <- setOrderingFilter(cds, ordering_genes)
  
  print(nrow(ordering_genes))
}

# if low # cells or features, do not select the number of cluster centers automatically -> https://github.com/cole-trapnell-lab/monocle-release/issues/26
# this avoids the error "initial centers are not distinct."
if (ncol(counts) < 500 || nrow(counts) < 500) {
  parameters$auto_param_selection <- FALSE
}

# reduce the dimensionality
cds <- monocle::reduceDimension(
  cds,
  max_components = parameters$max_components,
  reduction_method = parameters$reduction_method,
  norm_method = parameters$norm_method,
  auto_param_selection = parameters$auto_param_selection
)

# order the cells
cds <- monocle::orderCells(cds)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# extract the igraph and which cells are on the trajectory
gr <- monocle::minSpanningTree(cds)
to_keep <- setNames(rep(TRUE, nrow(counts)), rownames(counts))

# convert to milestone representation
cell_graph <- igraph::as_data_frame(gr, "edges") %>% mutate(directed = FALSE)

if ("weight" %in% colnames(cell_graph)) {
  cell_graph <- cell_graph %>% rename(length = weight)
} else {
  cell_graph <- cell_graph %>% mutate(length = 1)
}

cell_graph <- cell_graph %>% select(from, to, length, directed)

#   ____________________________________________________________________________
#   Save output                                                             ####

# create trajectory
output <- 
  dynwrap::wrap_data(cell_ids = colnames(cds)) %>%
  dynwrap::add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep
  ) %>%
  dynwrap::add_root(root_cell_id = cds@auxOrderingData[[parameters$reduction_method]]$root_cell) %>%
  dynwrap::add_timings(checkpoints)

# construct dimreds
dimred <- t(cds@reducedDimS)
colnames(dimred) <- paste0("Comp", seq_len(ncol(dimred)))


dimred_segments <-
  output$milestone_network %>%
  mutate(from_cell = gsub("milestone_", "", from), to_cell = gsub("milestone_", "", to)) %>% 
  as_tibble() %>% 
  rowwise() %>%
  mutate(
    path = igraph::shortest_paths(gr, from_cell, to_cell, mode = "out")$vpath %>% map(names),
    percentage = list(seq(0, 1, length.out = length(path)))
  ) %>%
  unnest(path, percentage) %>% 
  ungroup() 

dimred_segment_progressions <- 
  dimred_segments %>% 
  select(from, to, percentage)

dimred_segment_points <- 
  t(monocle::reducedDimK(cds))[dimred_segments$path, , drop = FALSE]
rownames(dimred_segment_points) <- NULL
colnames(dimred_segment_points) <- colnames(dimred)

dimred_milestones <- 
  t(monocle::reducedDimK(cds))[gsub("milestone_", "", output$milestone_ids), , drop = FALSE]
rownames(dimred_milestones) <- output$milestone_ids
colnames(dimred_milestones) <- colnames(dimred)

# add dimred
output <- 
  output %>% 
  dynwrap::add_dimred(
    dimred = dimred,
    dimred_milestones = dimred_milestones,
    dimred_segment_progressions = dimred_segment_progressions,
    dimred_segment_points = dimred_segment_points
  )

# write to file
dyncli::write_output(output, task$output)
