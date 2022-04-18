
library(PRROC)
library(gatom)
library(devtools)
library(igraph)
library(data.table)

simulate_data <- function(graph, alpha) {
  signals <- unique(E(graph)$signal)
  graph$signals <- setNames(rep(1, length(signals)), signals)
  module <- sample_subgraph(graph, subgraph_order=100, niter=1)
  modules_signals <- unique(E(graph)$signal[module])
  module_location <- signals %in% modules_signals
  graph$signals[module_location] <- rbeta(length(modules_signals), alpha, 1)
  graph$signals[!module_location] <- runif(length(signals) - length(modules_signals), 0, 1)
  for (i in seq_along(graph$signals)) {
    signal_name <- names(graph$signals)[i]
    location <- E(graph)$signal %in% signal_name
    E(graph)$pval[location] <- graph$signals[i]
  }
  graph <- delete_graph_attr(graph, "signals")
  return(list("module" = module, "graph" = graph, "signals" = modules_signals))
}

auc_roc_by_edges <- function(module, mcmc_score) {
  not_module <- !(seq_along(lengths(mcmc_score)) %in% module)
  roc <- roc.curve(mcmc_score[module], mcmc_score[not_module], curve = TRUE)
  return(roc)
}

auc_roc_by_signal <- function (graph, mcmc, module_signals) {
  groupsN <- as.integer(factor(E(graph)$signal))

  x <- Matrix::sparseMatrix(j=seq_along(groupsN),
                          i=groupsN,
                          x=rep(1, length(groupsN)))

  z <- x %*% t(mcmc$mat)
  z1 <- pmin(z, 1)
  rownames(z1) <- levels(factor(E(graph)$signal))

  prob <- rowMeans(z1)

  not_module <- names(graph$signals)[!(names(graph$signals) %in% module_signals)]
  roc <- roc.curve(prob[module_signals], prob[not_module], curve = TRUE)
  return(roc)
}

run_auc_roc <- function() {
  load("data/gatom_graph.rda")
  g <- simplify(gatom_graph, remove.multiple = TRUE, edge.attr.comb="min")
  simulated <- simulate_data(g, 0.2)
  gs <- score_graph(simulated$graph)
  V(gs)$likelihood <- 1


  mcs <- mcmc_sample(gs,  times=100, niter=30000, exp_lh = get_exp_lh(gs), edge_penalty = 0.025)
  mcmc_score <- colSums(mcs$mat)
  by_signal <- auc_roc_by_signal(gs, mcs, simulated$signals)
  by_edge <- auc_roc_by_edges(simulated$module, mcmc_score)
  return(list(by_edge=by_edge, by_signal=by_signal, simulated=simulated, mcmc_res=mcs))
}