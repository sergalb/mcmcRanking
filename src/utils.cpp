#include <Rcpp.h>
#include <vector>
#include <queue>
#include "utils.h"
#include "edge.h"

using namespace Rcpp;
using namespace std;

vector<vector<mcmc::Edge>> adj_list(List edgelist, size_t gorder) {
    vector<vector<mcmc::Edge>> edges(gorder);
    cout << "start adj list" << endl;
    for (size_t i = 0; i < edgelist.size(); ++i) {
        cout << "iter: " << i << endl;
        List edge = edgelist[i];
        string froms = edge[0];
        string tos = edge[1];
        int from = stoi(froms);
        int to = stoi(tos);
        cout << "take edge, from: " << from << " to: " << to << endl;
        string signalName = edge[2];
        cout << "take signal: " << signalName << endl;
        edges[from].push_back(mcmc::Edge(from, to, 1, "static_cast<string>(signalName[0])"));
        edges[to].push_back(mcmc::Edge(to, from, 1, "static_cast<string>(signalName[0])"));
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
