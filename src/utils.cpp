#include <Rcpp.h>
#include <vector>
#include <queue>
#include "utils.h"
#include "edge.h"

using namespace Rcpp;
using namespace std;

vector<mcmc::Edge> adj_list(List edgelist) {
    vector<mcmc::Edge> edges;
    for (size_t i = 0; i < edgelist.size(); ++i) {
        List edge = edgelist[i];
         edges.emplace_back(edge[0], edge[1], edge[2], i);
    }
    return edges;
}


vector<vector<unsigned>> adj_list(IntegerMatrix edgelist, size_t gorder) {
    vector<vector<unsigned>> edges(gorder);
    for (int i = 0; i < edgelist.nrow(); ++i) {
        edges[edgelist(i, 0)].push_back(edgelist(i, 1));
        edges[edgelist(i, 1)].push_back(edgelist(i, 0));
    }
    return edges;
}
