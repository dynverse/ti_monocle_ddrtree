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
gr <- cds@auxOrderingData[[parameters$reduction_method]]$cell_ordering_tree
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


reduced_dim_coords <- monocle::reducedDimK(cds)
space_df <- reduced_dim_coords %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("prin_graph_dim_1", "prin_graph_dim_2")) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))

edge_df <- cds %>%
  monocle::minSpanningTree() %>%
  igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
  left_join(space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")


#   ____________________________________________________________________________
#   Save output                                                             ####

output <- 
  dynwrap::wrap_data(cell_ids = rownames(dimred)) %>%
  dynwrap::add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep
  ) %>%
  dynwrap::add_dimred(dimred) %>%
  dynwrap::add_root(root_cell_id = cds@auxOrderingData[[parameters$reduction_method]]$root_cell) %>%
  dynwrap::add_timings(checkpoints)

dyncli::write_output(output, task$output)
