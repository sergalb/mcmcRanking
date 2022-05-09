
library(PRROC)
library(gatom)
library(devtools)
library(igraph)
library(data.table)
library(mwcsr)
load_all()

simulate_data <- function(graph, alpha, subgraph_order) {
  if (is.null(E(graph)$signal)) {
    E(graph)$signal <- E(graph)$label
  }
  signals <- unique(unlist(E(graph)$signal))
  graph$signals <- setNames(rep(1, length(signals)), signals)
  module <- sample_subgraph(graph, subgraph_order=subgraph_order, niter=1)
  modules_signals <- lapply(unique(E(graph)$signal[module]), \(ms) sample(ms, 1))
  modules_signals <- unique(modules_signals)
  module_location <- signals %in% modules_signals
  graph$signals[module_location] <- rbeta(length(modules_signals), alpha, 1)
  graph$signals[!module_location] <- runif(length(signals) - length(modules_signals), 0, 1)
  E(graph)$pval <- lapply(E(graph)$signal,
                          function(s_name) as.numeric(graph$signals[s_name]))

  graph <- delete_graph_attr(graph, "signals")
  return(list("module" = module, "graph" = graph, "signals" = modules_signals))
}

auc_roc_by_edges <- function(module, mcmc_score) {
  not_module <- !(seq_along(lengths(mcmc_score)) %in% module)
  roc <- roc.curve(mcmc_score[module], mcmc_score[not_module], curve = TRUE)
  return(roc)
}


split_signals <- function(matrix, signals) {
  mat <- rbind(matrix, sapply(signals, length))
  split_row <- function(x, len) c(rep(0, max(0, x-1)), min(x, 1) , rep(0,len - max(x, 1)))
  splited <- apply(mat, 2, function(col)
    matrix(t(
      sapply(col,
             \(x) split_row(x,col[length(col)])
      ))
      , nrow=length(col)))
  return(do.call("cbind", splited)[seq_len(nrow(matrix)),])
}

auc_roc_by_signal <- function (graph, mcmc, module_signals) {
  groupsN <- as.integer(factor(unlist(E(graph)$signal)))

  x <- Matrix::sparseMatrix(j=seq_along(groupsN),
                          i=groupsN,
                          x=rep(1, length(groupsN)))

  z <- x %*% t(split_signals(mcmc$mat, E(graph)$signal))
  z1 <- pmin(z, 1)
  rownames(z1) <- levels(factor(unlist(E(graph)$signal)))

  prob <- rowMeans(z1)

  not_module <- names(graph$signals)[!(names(graph$signals) %in% module_signals)]
  roc <- roc.curve(prob[unlist(module_signals)], prob[not_module], curve = TRUE)
  return(roc)
}

run_auc_roc <- function(alpha, subgraph_order,  times, niter, edge_penalty, graph) {
  gs <- simplify(graph, remove.multiple = TRUE, edge.attr.comb="concat")
  simulated <- simulate_data(gs, alpha, subgraph_order)
  gs <- score_graph(simulated$graph)
  V(gs)$likelihood <- 1


  mcs <- mcmc_sample(gs, times=times, niter=niter, exp_lh = get_exp_lh(gs), edge_penalty = edge_penalty)
  mcmc_score <- colSums(pmin(mcs$mat, 1))
  by_signal <- auc_roc_by_signal(gs, mcs, simulated$signals)
  by_edge <- auc_roc_by_edges(simulated$module, mcmc_score)
  return(list(by_edge=by_edge, by_signal=by_signal, simulated=simulated, mcmc=mcs, graph=gs))
}

auc_roc_by_graph <- function (gs, res, signals) {
  in_res <- unique(E(res)$signal)
  not_module <- setdiff(unique(E(gs)), signals)
  roc <- roc.curve(as.numeric(signals %in% in_res), as.numeric(not_module %in% in_res) , curve = TRUE)
  return(roc)
}

rnc_auc_roc <- function() {
  load("data/gatom_graph.rda")
  gs <- simplify(gatom_simpl, remove.multiple = TRUE, edge.attr.comb="min")
  simulated <- simulate_data(gs, 0.1)
  gs <- scoreGraph(simulated$graph, k.gene=35, k.met=NULL)


  solver <- rnc_solver()
  solver_res <- solve_mwcsp(solver, gs)
  m <- solver_res$graph

  simulated$signals <- unique(E(gs)$signal[simulated$module])
  auc <- auc_roc_by_graph(gs, m, simulated$signals)

  return(list(auc=auc, simulated=simulated, solver_res=solver_res, graph=gs))
}

bench <- function(alphas, orders_penalties, times, niters, graphs) {
  res <- list()
  for (alpha in alphas) {
    for (order_penalty in orders_penalties) {
      for (time in times) {
        for (niter in niters) {
          for (graph_name in names(graphs)) {
            auc <- run_auc_roc(alpha, order_penalty[[1]], time, niter, order_penalty[[2]], graphs[[graph_name]])
            res <- c(res, setNames(list(auc), paste(graph_name, alpha, order_penalty[[1]], time, niter, order_penalty[[2]], sep="_")))
          }
        }
      }
    }
  }
  return(res)
}

experiments_table <- function(experiments) {
  res <- do.call(rbind, lapply(experiments, function(organism) {
    t(data.frame(
      lapply(names(organism), function(name) {
        experiment <- organism[[name]]
       row <- append(stri_split(name, regex="_")[[1]], c(experiment$by_edge$auc, experiment$by_signal$auc))
       return(row)
      }
    )))
  }))
  return(res)
}

# run_auc_roc(0.1, 100, graphs$kegg_mouse_metabolites)
# rnc_auc_roc()
# bench(c(0.1), c(50), c(1), c(20), list(kegg_hs_metabolites=graphs$kegg_hs_metabolites))
