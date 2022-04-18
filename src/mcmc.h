#ifndef MCMC_RANKING_MCMC_H
#define MCMC_RANKING_MCMC_H

#include <random>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string> 
#include <Rcpp.h>
#include "edge.h"
#include "hsa.h"

namespace mcmc {
    using namespace std;

    class Graph {
        const double edge_penalty;
        bool fixed_size;
        size_t order;
        vector<double> nodes;
        vector<Edge> list_edges;
        vector <vector<pair<size_t, size_t>>> edges;

        HSA inner;
        HSA outer;

        vector <size_t> in_nei_c;

        // neis[v][i].first - neighbour
        // neis[v][i].second - index of v in neis[neighbour]
        vector <vector<pair<size_t, size_t> >> neis;

        mt19937 gen;
        uniform_real_distribution<> unirealdis;

        // These vectors are not local in is_connected in order to avoid memory allocations
        vector<unsigned> bfsUsed;
        vector<unsigned> dsu, dsuCnt;
        unsigned bfsUsedIteration = 0;
        vector<size_t> bfsQueue;

        unordered_map<string, pair<double, int>> signals;

        void update_outer_edges(unsigned cand_in, unsigned cand_out);

        void update_neighbours(unsigned v, bool is_erased);

        void inner_update(unsigned v, bool is_erased);

        bool is_connected(Edge const& erased);

        void process_edge_for_random_subgraph(size_t ind, unordered_set<size_t> const & sg, HSA & candidates);

        size_t edge_neighbours_size(size_t edge);

        vector<pair<size_t, size_t>> get_neighbours_edge(Edge const& e);

        vector<pair<size_t, size_t>> get_neighbours_edge(size_t e);

        void remove_vertex_from_neis(size_t vertex, size_t position_to_erase);

        void update_signals(Edge const& edge, bool removed);

        double probability_on_change_vertex(size_t cand_in, size_t cand_out, size_t cur_size_outer, size_t new_size_outer);

        double probability_on_change_vertex(size_t cand, size_t cur_size, size_t new_size, bool erase);

    public:
        Graph(vector<double> nodes, vector<Edge> list_edges, unordered_map<string, pair<double, int>> signals, bool fixed_size, double edge_penalty);

        Graph(Rcpp::NumericVector nodes, vector<Edge> list_edges, unordered_map<string, pair<double, int>> signals, bool fixed_size, double edge_penalty);

        void set_nodes(Rcpp::NumericVector nodes);

        bool next_iteration();

        void initialize_module(vector<size_t> const & edges);

        vector<size_t> random_subgraph(size_t size);

        vector <size_t> get_inner_edges();

        vector <size_t> get_outer_edges();
        
        unordered_map<string, pair<double, int>> const& get_signals() const;

    };
    


}

#endif //MCMC_RANKING_MCMC_H
