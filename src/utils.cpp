#include <Rcpp.h>
#include <vector>
#include <queue>
#include <type_traits>
#include "utils.h"
#include "edge.h"

using namespace Rcpp;
using namespace std;

vector<mcmc::Edge> adj_list(List edgelist) {
    vector<mcmc::Edge> edges;
    // cout << "adj list" << endl;
    for (size_t i = 0; i < edgelist.size(); ++i) {
        List const& edge = edgelist[i];
        int from, to;
        if (TYPEOF(edge[0]) != TYPEOF(edge[1])) throw domain_error("from and to had different types in edgelist");
        if (TYPEOF(edge[0]) == STRSXP) {
            string froms = edge[0];
            string tos = edge[1];
            from = stoi(froms);
            to = stoi(tos);
        } else if (TYPEOF(edge[0]) == REALSXP) {
            from =edge[0];
            to =edge[1];
        }
        edges.emplace_back(from, to, edge[2], i);
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

unordered_map<string, pair<double, int>> convert_signals_to_map(DataFrame const&  signals) {
    unordered_map<string, pair<double, int>> result;
    StringVector const& signals_names = signals[0];
    NumericVector const& signals_weights = signals[1];
    for (size_t i = 0; i < signals_names.size(); ++i) {
        result.emplace(signals_names[i], pair<double, int>(signals_weights[i], 0));
    }
    return result;
}

vector<pair<string, double>> get_signals_value(DataFrame const& signals, size_t ind) {
    vector<pair<string, double>> result;
    StringVector const& signals_names = signals[0];
    NumericVector const& signals_weights = signals[ind];
    for (size_t i = 0; i < signals_names.size(); ++i) {
        result.emplace_back(signals_names[i], signals_weights[i]);
    }
    return result;
}