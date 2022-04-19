library(data.table)


#' Number of runs of iterative likelihood increase of vertices.
#'
#' Running MCMC with likelihood weights on vertices can cause to get stuck in
#' locally good vertices, thereby finding locally significant solutions. So, we
#' will iteratively increase vertex weights and this function calculates number
#' of necessary iterations.
#'
#' @param x According to this value repetition depth calculates.
#' @return A number of runs.
#' @export
#' @examples
#' repetition_depth(16)
#' repetition_depth(81)
repetition_depth <- function(x) {
  d <- 0
  while (x > 4) {
    x <- sqrt(x)
    d <- d + 1
  }
  return(d - 2)
}

get_exp_lh <- function(graph) {
  depth <- repetition_depth(exp(max(graph$signals)) / exp(min(graph$signals)))
  return(1 / 2^depth)
}


#' Frequency of vertices.
#'
#' Calculates the frequency of occurrences of vertices in matrix object.
#'
#' @param mcmcObj Object of class MCMC.
#' @param inds Index numbers of rows involved in the calculation.
#' @param prob Logical scalar. Whether to return frequentist probability.
#' @return A named vector of frequency.
#' @seealso \code{\link{mcmc_sample}, \link{mcmc_onelong}}
#' @export
#' @examples
#' data(exampleGraph)
#' x <- mcmc_sample(exampleGraph, subgraph_order = 0, times = 1e3, niter = 100)
#' freq <- get_frequency(x)
#' tail(sort(freq))
#' p <- get_frequency(x, prob = TRUE)
#' tail(sort(p))
get_frequency <-
  function(mcmcObj,
           inds = seq_len(nrow(mcmcObj$mat)),
           prob = FALSE) {
    freq <- colSums(mcmcObj$mat[inds, , drop = FALSE])
    if (prob)
      return(freq / nrow(mcmcObj$mat))
    return(freq)
  }


#' Set likelihood.
#'
#' Set \emph{likelihood} vertex attribute to an \code{igraph} graph using pval
#' vertex attribute.
#'
#' @param graph An \code{igraph} graph.
#' @param fdr Numeric constant, from the false discovery rate a p-value
#'   threshold is calculated.
#' @return An \code{igraph} graph with \emph{likelihood} vertex attribute.
#' @importFrom igraph V V<-
#' @importFrom BioNet fitBumModel scoreFunction
#' @export
#' @examples
#' data(exampleGraph)
#' set_likelihood(exampleGraph, 1e-7)
set_likelihood <- function(graph, fdr) {
  pvals <- V(graph)$pval
  names(pvals) <- V(graph)$name
  fb <- fitBumModel(pvals, plot = FALSE)
  V(graph)$likelihood <- exp(scoreFunction(fb = fb, fdr = fdr))
  graph
}


#' score_graph.
#'
#' score graph.
#'
#' @param graph An \code{igraph} graph.
#' @return An \code{igraph} graph with \emph{likelihood} vertex attribute.
#' @importFrom data.table data.table
#' @importFrom BioNet fitBumModel scoreFunction
#' @export
#' @examples
#' data(exampleGraph)
#' score_graph(exampleGraph)
score_graph <- function(graph) {
  edge.table <- data.table(as_data_frame(graph, what = "edges"))

  pvalsToFit <- edge.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]

  edge.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F)
  if (edge.bum$a > 0.5) {
    E(graph)$score <- 0
    warning("Edge scores have been assigned to 0 due to an inappropriate p-value distribution")
  } else {
    edge.threshold <- BioNet::fdrThreshold(0.1, edge.bum)
    E(graph)$score <- with(edge.table,
                       (edge.bum$a - 1) *
                         (log(.replaceNA(pval, 1)) - log(edge.threshold)))
  }


  graph$signals <- setNames(E(graph)$score, E(graph)$signal)
  return(graph)

}


get_gatom_graph <- function() {
  network <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
  org.Mm.eg.gatom.anno <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))
  # org.Hs.eg.gatom.anno <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rds"))
  met.db <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))
  library(R.utils)
  library(data.table)
  met.de.raw <- fread("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.met.de.tsv.gz")
  gene.de.raw <- fread("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.gene.de.tsv.gz")
  g <- makeMetabolicGraph(network = network,
                          topology = "metabolites",
                          org.gatom.anno = org.Mm.eg.gatom.anno,
                          gene.de = gene.de.raw,
                          met.db = met.db,
                          met.de = met.de.raw)
  return(g)
}


.replaceNA <- function(x, y) { ifelse(is.na(x), y, x) }


