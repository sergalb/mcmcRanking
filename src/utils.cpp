#include <Rcpp.h>
#include <vector>
#include <queue>
#include "utils.h"
#include "edge.h"

using namespace Rcpp;
using namespace std;

vector<mcmc::Edge> adj_list(List edgelist, size_t gorder) {
    vector<mcmc::Edge> edges;
    for (size_t i = 0; i < edgelist.size(); ++i) {
      //  cout << "iter: " << i << endl;
        List edge = edgelist[i];
        string froms = edge[0];
        string tos = edge[1];
        int from = stoi(froms);
        int to = stoi(tos);
        string signalName = edge[2];
        edges.emplace_back(from, to, signalName, i);
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
