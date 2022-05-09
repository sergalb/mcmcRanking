#include <Rcpp.h>
#include <vector>
#include <RProgress.h>
#include <string>
#include "mcmc.h"
#include "utils.h"
#include <typeinfo>       // operator typeid
#include <functional>


using namespace Rcpp;
using namespace std;
using mcmc::Graph;

// [[Rcpp::export]]
LogicalVector sample_subgraph_internal(List edgelist, DataFrame signals, int gorder, int module_size, size_t niter, double edge_penalty) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    // cout << "start sample" << endl;
    vector<double> nodes(gorder, 1);
    Graph g = Graph(nodes, adj_list(edgelist), convert_signals_to_map(signals), true, edge_penalty);
    g.initialize_module(g.random_subgraph(module_size));
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    LogicalVector ret(edgelist.size(), false);
    for (size_t x : g.get_inner_edges()) {
       ret[x] = true;
    }
    return ret;
}

// [[Rcpp::export]]
NumericVector sample_llh_internal(List edgelist, DataFrame signals, NumericVector likelihood, size_t niter, bool fixed_size,
                                  LogicalMatrix start_module, double edge_penalty) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    Graph g = Graph(likelihood, adj_list(edgelist), convert_signals_to_map(signals), fixed_size, edge_penalty);
    vector<size_t> module;
    for (size_t j = 0; j < start_module.ncol(); ++j) {
        if (start_module(0, j)) {
            module.push_back(j);
        }
    }
    g.initialize_module(module);
    NumericVector llhs(niter, 0);
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
        for (auto const& x : g.get_signals()) {
            if (x.second.second != 0) {
                llhs[i] += log(x.second.first);
            }
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    return llhs;
}

// template<typename T>
void print_vec(vector<double> vec, string spliterator = " ") {
    for (auto e : vec) {
        cout << e << spliterator;
    }
    cout << endl;
}


// [[Rcpp::export]]
NumericMatrix mcmc_sample_internal(List edgelist, DataFrame signals, NumericMatrix likelihood, bool fixed_size, size_t niter,
                                   LogicalMatrix start_module, double edge_penalty) {
                                    //    cout << "start mcmc sample" << endl;
    Graph g = Graph((NumericVector) likelihood(_, 0), adj_list(edgelist), convert_signals_to_map(signals), fixed_size, edge_penalty);
    size_t order = edgelist.size();
    unsigned times = start_module.nrow();
    //LogicalVector ret(order * times, false);
    NumericMatrix ret(times, order);
    // colnames(ret) = rownames(edgelist);

    RProgress::RProgress pb;
    if (times * likelihood.ncol() * niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", times * likelihood.ncol() * niter);
        pb.tick(0);
    }
    for (int i = 0; i < times; ++i) {
        vector<size_t> module;
        for (size_t j = 0; j < order; ++j) {
            if (start_module(i, j)) {
                module.push_back(j);
            }
        }
        g.initialize_module(module);
        for (int k = 1; k < signals.ncol(); ++k) {
            // cout << "set signals, k: " << k << " i: " << i << endl;
            g.set_signals(get_signals_value(signals, k));
            for (size_t j = 0; j < niter; ++j) {
                g.next_iteration();
                if (j % 10000 == 9999) {
                    Rcpp::checkUserInterrupt();
                    pb.tick(10000);
                }
            }
            if (niter > 0)
                pb.tick(niter % 10000);
        }
        for (size_t x : g.get_inner_edges()) {
            //ret[x + i * order] = true;
            ret(i, x) = g.get_active_signal_by_edge(x) + 1;
        }
    }
    return ret;
}

// [[Rcpp::export]]
LogicalMatrix
mcmc_onelong_internal(List edgelist, DataFrame signals, NumericVector likelihood, bool fixed_size, int module_size, size_t start,
                      size_t niter, double edge_penalty) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    Graph g = Graph(likelihood, adj_list(edgelist), convert_signals_to_map(signals), fixed_size, edge_penalty);
    size_t order = likelihood.size();
    g.initialize_module(g.random_subgraph(module_size));
    //LogicalVector ret(order * (niter - start), false);
    LogicalMatrix ret(niter - start, order);
    colnames(ret) = CharacterVector(likelihood.names());
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
        if (i < start) {
            continue;
        }
        for (size_t x : g.get_inner_edges()) {
            //ret[x + (i - start) * order] = true;
            ret(i-start, x) = true;
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    return ret;
}

// [[Rcpp::export]]
IntegerVector
mcmc_onelong_frequency_internal(List edgelist, DataFrame signals, NumericVector likelihood, bool fixed_size, int module_size,
                                size_t start, size_t niter, double edge_penalty) {
    RProgress::RProgress pb;
    if (niter > 0) {
        pb = RProgress::RProgress("[:bar] ETA: :eta", niter);
        pb.tick(0);
    }
    Graph g = Graph(likelihood, adj_list(edgelist), convert_signals_to_map(signals), fixed_size, edge_penalty);
    size_t order = likelihood.size();
    g.initialize_module(g.random_subgraph(module_size));
    IntegerVector ret(order, 0);
    ret.names() = likelihood.names();
    for (size_t i = 0; i < niter; ++i) {
        g.next_iteration();
        if (i % 10000 == 9999) {
            Rcpp::checkUserInterrupt();
            pb.tick(10000);
        }
        if (i < start) {
            continue;
        }
        for (size_t x : g.get_inner_edges()) {
            ret[x]++;
        }
    }
    if (niter > 0)
        pb.tick(niter % 10000);
    return ret;
}
