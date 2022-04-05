#ifndef UTILS_MCMC_H
#define UTILS_MCMC_H

#include <Rcpp.h>
#include <vector>
#include "edge.h"

using namespace Rcpp;
using namespace std;

/**
 * @brief convert list of all edges to vector of edges int incident to vertex of particular index.
 * result[0] contains all edges of vertex[0].
 * @param edgelist elements: edge[0] - int number of 'from' vertex; edge[1] - int number of 'to' vertex; edge[2] - StringVector of signals for edge.
 * @param signal map from signal to it weight signals[0] - name of signal; signal[1] - weight of signal
 * @param gorder count of vertex in graph
 * @return vector<vector<mcmc::Edge>>
 */
vector<mcmc::Edge> adj_list(List edgelist, size_t gorder);

vector<vector<unsigned>> adj_list(IntegerMatrix edgelist, size_t gorder);


#endif
