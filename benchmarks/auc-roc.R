
library(PRROC)
library(gatom)
library(devtools)
library(igraph)
library(data.table)
library(mwcsr)
library(stringi)
load_all()

simulate_data_for_simpl_graphs <- function(graph, alpha, subgraph_order) {
  gs <- simplify(graph, remove.multiple = TRUE, edge.attr.comb="concat")
  gs_simple <- simplify(graph, remove.multiple = TRUE, edge.attr.comb="min")
  if (is.null(E(gs)$signal)) {
    E(gs)$signal <- E(gs)$label
    E(gs_simple)$signal <- E(gs_simple)$label

  }
  signals <- unique(unlist(E(gs)$signal))
  gs$signals <- setNames(rep(1, length(signals)), signals)
  module <- sample_subgraph(gs, subgraph_order=subgraph_order, niter=1)
  modules_signals <- lapply(unique(E(gs)$signal[module]), \(ms) sample(ms, 1))
  modules_signals <- unique(modules_signals)
  module_location <- signals %in% modules_signals
  gs$signals[module_location] <- rbeta(length(modules_signals), alpha, 1)
  gs$signals[!module_location] <- runif(length(signals) - length(modules_signals), 0, 1)

  E(gs_simple)$pval <- sapply(E(gs_simple)$signal,
                          function(s_name) as.numeric(gs$signals[s_name]))

  return(list("module" = module, "graph" = gs_simple, "signals" = modules_signals))
}

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

prob_by_mcmc <- function (graph, mcmc, signal_selector) {
  groupsN <- as.integer(factor(unlist(signal_selector(E(graph)))))

  x <- Matrix::sparseMatrix(j=seq_along(groupsN),
                          i=groupsN,
                          x=rep(1, length(groupsN)))

  z <- x %*% t(split_signals(mcmc$mat, signal_selector(E(graph))))
  z1 <- pmin(z, 1)
  rownames(z1) <- levels(factor(unlist(signal_selector(E(graph)))))

  return(rowMeans(z1))
}

auc_roc_by_signal <- function (graph, mcmc, module_signals) {
  prob <- prob_by_mcmc(graph, mcmc, \(edges) edges$signal)

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

roc_point <- function (gs, res, signals) {
  in_res <- unique(E(res)$signal)
  tp <- sum(in_res %in% signals)
  fp <- length(in_res) - tp
  p <- length(signals)
  n <- length(unique(E(gs)$signal)) - length(signals)
  return(list(x=fp/n, y=tp/p))
}

rnc_roc_point <- function(alpha, subgraph_order, k.gene, edge.threshold.min=1, graph) {
  gs <- simplify(graph, remove.multiple = TRUE, edge.attr.comb="min")
  simulated <- simulate_data(gs, alpha, subgraph_order)
  gs <- scoreGraph(simulated$graph, k.gene=k.gene, edge.threshold.min=edge.threshold.min, k.met=NULL)


  solver <- rnc_solver()
  solver_res <- solve_mwcsp(solver, gs)
  m <- solver_res$graph

  simulated$signals <- unique(E(gs)$signal[simulated$module])
  point <- roc_point(gs, m, simulated$signals)

  return(list(point=point, simulated=simulated, solver_res=solver_res, graph=gs))
}

virgo_cplex_roc_point <- function(alpha, subgraph_order, k.gene, threads = 4, edge.threshold.min=1, graph) {
  simulated <- simulate_data_for_simpl_graphs(graph, alpha, subgraph_order)
  gs <- scoreGraph(simulated$graph, k.gene=k.gene, edge.threshold.min=edge.threshold.min, k.met=NULL)


  vsolver <- virgo_solver(cplex_jar="/Applications/CPLEX_Studio221/cplex/lib/cplex.jar", cplex_bin="/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/libcplex2210.dylib",  threads=threads, penalty=0.001)
  sol <- solve_mwcsp(vsolver, gs)
  m <- sol$graph

  simulated$signals <- unique(E(gs)$signal[simulated$module])
  point <- roc_point(gs, m, simulated$signals)

  return(list(point=point, simulated=simulated, solver_res=sol, graph=gs))
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

table_from_experiment <- function(experiment) {
  t(data.frame(
    lapply(names(experiment), function(name) {
      experiment_iter <- experiment[[name]]
      append(stri_split(name, regex="_")[[1]], c(experiment_iter$by_edge$auc, experiment_iter$by_signal$auc, mean(rowSums(experiment_iter$mcmc$mat))))
    }
  )))
}

experiments_table <- function(experiments, postprocess=identity) {
  res <- do.call(rbind, lapply(experiments, table_from_experiment))
  return(postprocess(res))
}

remove_kegg_add_const_colnames <- function (table) {
  res <- table[,-1]
  colnames(res) <- c("organism", "topology", "alpha", "module_order", "mcmc_times", "mcmc_niter", "auc_by_edge", "auc_by_signal")
  rownames(res) <- seq_len(nrow(res))
  res
}

remove_kegg_add_pen_colnames <- function (table) {
  res <- table[,-1]
  colnames(res) <- c("organism", "topology", "alpha", "module_order", "mcmc_times", "mcmc_niter", "edge_penalty", "auc_by_edge", "auc_by_signal", "mean_result_module_order")
  rownames(res) <- seq_len(nrow(res))
  res
}


check_pathways <- function (graph, mcmc, size) {
  org.Mm.eg.gatom.anno <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
  genes_prob <- prob_by_mcmc(graph, mcmc, \(edges) edges$gene)
  module_genes <- names(genes_prob)[tail(order(genes_prob),size)]
  foraRes <- fgsea::fora(pathways=org.Mm.eg.gatom.anno$pathways,
                       genes=module_genes,
                       universe=unique(unlist(E(graph)$gene)),
                       minSize=5)
  return(foraRes)
}

plot_auc_roc_with_points <- function(auc_roc, points) {
  par(mar=c(6,7,6,3))
  plot(x=1,
       xlab="FPR",
       ylab="TPR",
       xlim = c(0,1),
       ylim=c(0,1),
       type="l",
       main=paste("ROC curve\nAUC = ", format(round(auc_roc$auc, 5), nsmall = 5), sep=""),
       cex.lab = 3,
     cex.axis = 3,
     cex.main = 3,
     cex.sub = 3)
  auc_points <- auc_roc$curve[,-3]
  lines(x=auc_points[,1], y=auc_points[,2], lwd=4, col=2, t="l", )
  points(lapply(points, \(x) x[[1]]), lapply(points, \(x) x[[2]]), pch =16, cex = 2, col = 3)
  legend("bottomright", legend = c("MCMC", "GATOM"),
       col = c(2, 3), lty = c(1,NA), lwd=4, pch=c(NA,16), cex=3, pt.cex = 3)

}
plot_auc_roc<- function(auc_roc) {
  par(mar=c(6,7,6,3))
  plot(x=1,
       xlab="FPR",
       ylab="TPR",
       xlim = c(0,1),
       ylim=c(0,1),
       type="l",
       main=paste("ROC curve\nAUC = ", format(round(auc_roc$auc, 5), nsmall = 5), sep=""),
       cex.lab = 3,
     cex.axis = 3,
     cex.main = 3,
     cex.sub = 3)
  auc_points <- auc_roc$curve[,-3]
  lines(x=auc_points[,1], y=auc_points[,2], lwd=4, col=2, t="l", )
  legend("bottomright", legend = c("MCMC"),
       col = c(2), lty = c(1), lwd=4, cex=3)

}
# run_auc_roc(0.05, 100, 30, 25000, 1, graphs$kegg_mouse_metabolites)
# run_auc_roc(0.1, 100, graphs$kegg_mouse_metabolites)
# rnc_auc_roc()
