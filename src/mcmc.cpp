#include <queue>
#include <cmath>
#include "mcmc.h"

namespace mcmc {
    Graph::Graph(vector<double> nodes, vector<Edge> list_edges, unordered_map<string, double> signals, bool fixed_size)
            : fixed_size(fixed_size), order(list_edges.size()), nodes(nodes), list_edges(list_edges), edges(nodes.size()), inner(order), outer(order),
              in_nei_c(order, 0), neis(nodes.size()),
              bfsUsed(nodes.size(), 0), signals(signals) {
        random_device rd;
        //todo back random
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
        bfsQueue.reserve(nodes.size());
        for (size_t i = 0; i < list_edges.size(); ++i) {
            auto edge = list_edges[i];
            edges[edge.first].emplace_back(edge.second, edge.id);
            edges[edge.second].emplace_back(edge.first, edge.id);
        }
    }

    string pair_to_str(pair<size_t, size_t> p) {
        return to_string(p.first) + " " + to_string(p.second);
    }

    template<typename T>
    void print_vec(vector<T> vec, function<string(T)> to_string, string spliterator = " ") {
        for (auto e : vec) {
            cout << to_string(e) << spliterator;
        }
        cout << endl;
    }

    Graph::Graph(Rcpp::NumericVector nodes, vector<Edge> list_edges, unordered_map<string, double> signals, bool fixed_size)
            : Graph(vector<double>(nodes.begin(), nodes.end()), list_edges, signals, fixed_size) {}

    void Graph::set_nodes(Rcpp::NumericVector nodes) {
        this->nodes = vector<double>(nodes.begin(), nodes.end());
    }
    
    void Graph::process_edge_for_random_subgraph(size_t ind, unordered_set<size_t> const & sg, HSA & candidates) {
        auto neighbours = get_neighbours_edge(ind);
        for (auto const& neighbour : neighbours) {
            if (!candidates.contains(neighbour.second)) {
                if (sg.count(neighbour.second) < 1){
                    candidates.insert(neighbour.second);
                }
            }
        }
    }

    vector<size_t> Graph::random_subgraph(size_t size) {
        cout << endl << "start random subgraph" << endl;
        if (size == 0) {
            return vector<size_t>();
        }
        unordered_set<size_t> sg;
        HSA candidates(order);
        size_t ind = uniform_int_distribution<>(0, order - 1)(gen);
        sg.insert(ind);
        process_edge_for_random_subgraph(ind, sg, candidates);
        while (sg.size() != size) {
            ind = uniform_int_distribution<>(0, candidates.size() - 1)(gen);
            unsigned new_edge = candidates.get(ind);
            process_edge_for_random_subgraph(new_edge, sg, candidates);
            sg.insert(new_edge);
            candidates.erase(new_edge);
        }
        return vector<size_t>(sg.begin(), sg.end());
    }

    vector<pair<size_t, size_t>> Graph::get_neighbours_edge(Edge e) {
        vector<pair<size_t, size_t>> neighbours(edges[e.first].size() - 1 + edges[e.second].size() - 1);
        auto it = copy_if(edges[e.first].begin(), edges[e.first].end(), neighbours.begin(), [&e](pair<size_t, size_t> const& edge) { return e.id != edge.second; } );
        copy_if(edges[e.second].begin(), edges[e.second].end(), it, [&e](pair<size_t, size_t> const& edge) { return e.id != edge.second; } );
        //todo check return optimization
        return neighbours;
    }


    vector<pair<size_t, size_t>> Graph::get_neighbours_edge(size_t e) {
        return get_neighbours_edge(list_edges[e]);
    }

    void Graph::initialize_module(vector<size_t> const & initial_edges) {
        cout << endl << "start initialize subgraph" << endl;
        inner.clear();
        outer.clear();
        std::fill(in_nei_c.begin(), in_nei_c.end(), 0);
        for (auto &x : neis) {
            x.clear();
        }
        for (size_t edge : initial_edges) {
            inner.insert(edge);
        }
        for (size_t i = 0; i < inner.size(); ++i) {
            Edge e = list_edges[inner.get(i)];
            vector<pair<size_t, size_t>> neighbours = get_neighbours_edge(e);
            for (auto neighbour: neighbours) {
                size_t neighbour_edge = neighbour.second;
                if (neighbour_edge == e.id) continue;
                in_nei_c[neighbour_edge]++;                    
                if (!inner.contains(neighbour_edge) && !outer.contains(neighbour_edge)) {
                    outer.insert(neighbour_edge);
                }
            }
            //todo understand how neis used (only for conection checking?)
            neis[e.first].emplace_back(e.second, neis[e.second].size());
            neis[e.second].emplace_back(e.first, neis[e.first].size() - 1);
        }
    }

    static unsigned dsuGetRoot(vector<unsigned> & dsu, unsigned v) {
        auto &p = dsu[v];
        return p == v ? v : p = dsuGetRoot(dsu, p);
    }

    bool Graph::is_connected(Edge erased) {
        cout << "is connected " << erased.to_string() << endl;
        inner_update(erased.id, true);

        bfsQueue.clear();
        
        fill(bfsUsed.begin(), bfsUsed.end(), 0);
        bfsUsed[erased.first] = 1;
        bfsQueue.push_back(erased.first);
    
        for (size_t i = 0; i < bfsQueue.size(); ++i) {
            auto v = bfsQueue[i];
            // cout << "v: " << v << endl;
            for (auto const& tov : neis[v]) {
                auto to = tov.first;
                // cout << "to: " << to << endl;
                if (bfsUsed[to]) continue;
                bfsQueue.push_back(to);
                bfsUsed[to] = 1;
                if (to == erased.second) break;
            }
            if (bfsUsed[erased.second]) break;
        }
        inner_update(erased.id, false);
        return bfsUsed[erased.second];
    }

    void Graph::update_outer_edges(unsigned cand_in, unsigned cand_out) {
        outer.swap(cand_out, cand_in);
        auto out_neighbour_edges = get_neighbours_edge(cand_out);
        for (auto const& edge : out_neighbour_edges) {
            size_t neighbour = edge.second;
            if (in_nei_c[neighbour]++ == 0 && neighbour != cand_in) {
                if (!inner.contains(neighbour)) {
                    outer.insert(neighbour);
                } else {
                    throw domain_error("inner edge " + to_string(neighbour) + " had no inner neighbours");
                }
            }
        }
        auto in_neighbour_edges = get_neighbours_edge(cand_in);
        for (auto const& edge : in_neighbour_edges) {
            size_t neighbour = edge.second;
            if (--in_nei_c[neighbour] == 0 && neighbour != cand_out) {
                if (!inner.contains(neighbour)) {
                    outer.erase(neighbour);
                } else {
                    auto n_of_n = get_neighbours_edge(neighbour);
                    for (auto const& n_e : n_of_n) {
                        cout << n_e.second << " inner: " << inner.contains(n_e.second) << "; ";
                    }
                    cout << endl;
                    throw domain_error("inner edge " + to_string(neighbour) + " have no inner neighbours");
                }
            }
        }
    }

    void Graph::remove_vertex_from_neis(size_t vertex, size_t position_to_erase) {
        auto &vec = neis[vertex];
        if (position_to_erase != vec.size() - 1) {
            auto newP = vec[position_to_erase] = neis[vertex].back();
            neis[newP.first][newP.second].second = position_to_erase;
        }
        vec.pop_back();
        // if (vertex == 169 || newP.first == 169 || vertex == 157 || newP.first == 157) {
        //     cout << "neis of 169 and 157" << endl;
        //     print_vec<pair<size_t, size_t>>(neis[169], &pair_to_str);
        //     print_vec<pair<size_t, size_t>>(neis[157], &pair_to_str);
        // }
        
    }

    void Graph::inner_update(unsigned e, bool is_erased) {
        // cout << "inner update, e: " << e << "; erased: " << is_erased << endl;
        if (is_erased) {
            inner.erase(e);
            Edge edge = list_edges[e];
            auto v1 = find_if(neis[edge.first].begin(), neis[edge.first].end(), [&edge](pair<size_t, size_t> const& nei) { return edge.second == nei.first; } );
            if (v1 == neis[edge.first].end()) {
                for (auto const& nei : neis[edge.first]) {
                    cout << nei.first << " " << nei.second << "; ";
                }
                cout << endl; 
                cout << edge.first << " " << edge.second << " " << edge.id << endl;
                throw domain_error("couldn't find eachother ends of edge in neis vector");
            }
            auto v2 = neis[edge.second][v1->second];
            if (e == 154) {
                cout << "v1: " << v1->first << " " << v1->second << endl;
                cout << "v2: " << v2.first << " " << v2.second << endl;
                print_vec<pair<size_t, size_t>>(neis[edge.first], &pair_to_str);
                print_vec<pair<size_t, size_t>>(neis[edge.second], &pair_to_str);
                print_vec<pair<size_t, size_t>>(neis[169], &pair_to_str);
                print_vec<pair<size_t, size_t>>(neis[157], &pair_to_str);
            }
            remove_vertex_from_neis(v1->first, v1->second);
            remove_vertex_from_neis(v2.first, v2.second);
        } else {
            if (inner.contains(e)) {
                throw domain_error("try to add in inner existing edge");
            } 
            inner.insert(e);
            Edge inserted_edge = list_edges[e];
            // size_t new_vertex;
            // if (neis[inserted_edge.first].empty()) {
            //     new_vertex = inserted_edge.first;
            // } else {
            //     new_vertex = inserted_edge.second;
            // } 
             // else {
                 // throw invalid_argument("both vertexes: " + to_string(inserted_edge.first) + " and " + to_string(inserted_edge.second) + " of add edge " + to_string(e) + " included to graph");
            // }
            // for (auto const& edge : edges[new_vertex]) {
            //     size_t v2 = edge.second;
            // if (!neis[v2].empty()) {
            neis[inserted_edge.first].emplace_back(inserted_edge.second, neis[inserted_edge.second].size());
            neis[inserted_edge.second].emplace_back(inserted_edge.first, neis[inserted_edge.first].size() - 1);
            // }
            // }
        }
    }

    void Graph::update_neighbours(unsigned v, bool is_erased) {
        if (is_erased) {
            if (inner.size() != 0)
                outer.insert(v);
            for (auto const& edge : edges[v]) {
                size_t neighbour = edge.second;
                if (--in_nei_c[neighbour] == 0) {
                    if (!inner.contains(neighbour)) {
                        outer.erase(neighbour);
                    }
                }
            }
        } else {
            if (inner.size() != 1)
                outer.erase(v);
            for (auto const& edge : edges[v]) {
                size_t neighbour = edge.second;
                if (in_nei_c[neighbour]++ == 0) {
                    if (!inner.contains(neighbour)) {
                        outer.insert(neighbour);
                    }
                }
            }
        }
    }

    double Graph::probability_on_change_vertex(size_t cand_in, size_t cand_out, size_t cur_size_outer, size_t new_size_outer) {
        double in_likelihood = signals[list_edges[cand_in].signal];
        double out_likelihood = signals[list_edges[cand_out].signal];
        return (out_likelihood * cur_size_outer) / (in_likelihood * new_size_outer);
    }

    bool Graph::next_iteration() {
        cout << endl << "start next iteration, fixed size: " << fixed_size << "; inner size: " << inner.size() << endl;
        if (fixed_size) {
            if (inner.size() == 0) {
                return true;
            }
            if (inner.size() == order) {
                return false;
            }
            unsigned cand_in = inner.get(uniform_int_distribution<>(0, inner.size() - 1)(gen));
            unsigned cand_out = outer.get(uniform_int_distribution<>(0, outer.size() - 1)(gen));
            double gen_p = unirealdis(gen);
            unsigned cur_size_outer = outer.size();
            unsigned new_size_outer;

            if (inner.size() == 1) {
                Edge const& in_e = list_edges[cand_in];
                new_size_outer = edges[in_e.first].size() - 1 + edges[in_e.second].size() - 1;
            } else {
                new_size_outer = outer.size();
                auto in_neighbour_edges = get_neighbours_edge(cand_in);
                auto out_neighbour_edges = get_neighbours_edge(cand_out);
                for (auto const& x : in_neighbour_edges) {
                    if (!--in_nei_c[x.second]) --new_size_outer;
                }
                for (auto const& x : out_neighbour_edges) {
                    if (!in_nei_c[x.second]) ++new_size_outer;
                }
                for (auto const& x : in_neighbour_edges) {
                    ++in_nei_c[x.second];
                }
            }
            double p = probability_on_change_vertex(cand_in, cand_out, cur_size_outer, new_size_outer);
            if (gen_p >= p)
                return false;
            inner_update(cand_out, false);
            Edge in_edge = list_edges[cand_in];
            if (!is_connected(in_edge)) {
                inner_update(cand_out, true);
                return false;
            }
            inner_update(cand_in, true);
            update_outer_edges(cand_in, cand_out);
            return true;
        } else {
            unsigned cur_in_out_size = inner.size() == 0 ? order : outer.size() + inner.size();
            bool erase = unirealdis(gen) < (1.0 * inner.size()) / cur_in_out_size;
            unsigned cand;
            if (erase) {
                cand = inner.get(uniform_int_distribution<>(0, inner.size() - 1)(gen));
            } else {
                cand = outer.size() == 0 ?
                       uniform_int_distribution<>(0, order - 1)(gen) :
                       outer.get(uniform_int_distribution<>(0, outer.size() - 1)(gen));
            }
            double gen_p = unirealdis(gen);
            if (erase) {
                int new_in_out_size = inner.size() + outer.size() - edges[cand].size() + in_nei_c[cand];
                if (inner.size() == 1) {
                    new_in_out_size = order;
                } else if (cur_in_out_size <= edges[cand].size() - in_nei_c[cand]) {
                    new_in_out_size = 1;
                }
                double p = cur_in_out_size / (nodes[cand] * new_in_out_size);
                if (gen_p >= p) {
                    return false;
                }
                if (!is_connected(list_edges[cand])) {
                    return false;
                }
                inner_update(cand, erase);
            } else {
                double p = nodes[cand] * (inner.size() == 0 ? 1.0 * order / (1 + edges[cand].size()) : 1);
                if (gen_p >= p) {
                    return false;
                }
                inner_update(cand, erase);
            }
            update_neighbours(cand, erase);
            unsigned new_in_out_size = inner.size() == 0 ? order : outer.size() + inner.size();
            double p = (erase ? 1 / nodes[cand] : nodes[cand]) * (1.0 * cur_in_out_size / new_in_out_size);
            if (gen_p < p) {
                return true;
            }
            inner_update(cand, !erase);
            update_neighbours(cand, !erase);
            return false;
        }
    }

    vector <size_t> Graph::get_inner_nodes() {
        
        return inner.get_all();
    }

    vector <size_t> Graph::get_outer_nodes() {
        return outer.get_all();
    }
};
