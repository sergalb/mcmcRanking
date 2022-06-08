#include <queue>
#include <cmath>
#include "mcmc.h"

namespace mcmc {
    Graph::Graph(vector<double> nodes, vector<Edge> list_edges, unordered_map<string, pair<double, int>> signals, bool fixed_size, double edge_penalty)
            : fixed_size(fixed_size), order(list_edges.size()), nodes(nodes), list_edges(list_edges), edges(nodes.size()), inner(order), outer(order),
              in_nei_c(order, 0), neis(nodes.size()),
              bfsUsed(nodes.size(), {false, false}), signals(signals), edge_penalty(edge_penalty){
        random_device rd;
        gen = mt19937(rd());
        unirealdis = uniform_real_distribution<>(0, 1);
        bfsQueue.reserve(nodes.size());
        for (size_t i = 0; i < list_edges.size(); ++i) {
            auto edge = list_edges[i];
            edges[edge.first].emplace_back(edge.second, edge.id);
            edges[edge.second].emplace_back(edge.first, edge.id);
        }
    }

    Graph::Graph(Rcpp::NumericVector nodes, vector<Edge> list_edges, unordered_map<string, pair<double, int>> signals, bool fixed_size, double edge_penalty)
            : Graph(vector<double>(nodes.begin(), nodes.end()), list_edges, signals, fixed_size, edge_penalty) {}

    string pair_to_str(pair<size_t, size_t> p) {
        return to_string(p.first) + " " + to_string(p.second);
    }

    template<typename T>
    void print_vec(vector<T> vec, function<string(T)> to_string, string spliterator = "; ") {
        for (auto e : vec) {
            cout << to_string(e) << spliterator;
        }
        cout << endl;
    }

    void Graph::set_signals(vector<pair<string, double>> signals_value) {
        for (auto const& sig_llh: signals_value) {
            auto it = signals.find(sig_llh.first);
            it->second = {sig_llh.second, it->second.second};
        }
    }

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
        //  cout << endl << "start random subgraph" << endl;
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

    size_t Graph::edge_neighbours_size(size_t edge) {
        Edge const& e = list_edges[edge];
        return edges[e.first].size() - 1 + edges[e.second].size() - 1;
    }

    vector<pair<size_t, size_t>> Graph::get_neighbours_edge(Edge const& e) {
        vector<pair<size_t, size_t>> neighbours(edge_neighbours_size(e.id));
        auto it = copy_if(edges[e.first].begin(), edges[e.first].end(), neighbours.begin(), [&e](pair<size_t, size_t> const& edge) { return e.id != edge.second; } );
        copy_if(edges[e.second].begin(), edges[e.second].end(), it, [&e](pair<size_t, size_t> const& edge) { return e.id != edge.second; } );
        //todo check return optimization or even make iterator instead of copying all edges
        return neighbours;
    }


    vector<pair<size_t, size_t>> Graph::get_neighbours_edge(size_t e) {
        return get_neighbours_edge(list_edges[e]);
    }

    void Graph::initialize_module(vector<size_t> const & initial_edges) {
        // cout << endl << "start initialize subgraph" << endl;
        inner.clear();
        outer.clear();
        std::fill(in_nei_c.begin(), in_nei_c.end(), 0);
        for (auto &x : neis) {
            x.clear();
        }
        for(auto & signal: signals) {
            signal.second.second = 0;
        }
        for (auto & edge: list_edges) {
            edge.active_signal = -1;
        }
        
        for (size_t edge : initial_edges) {
            inner.insert(edge);
        }
        for (size_t i = 0; i < inner.size(); ++i) {
            Edge &e = list_edges[inner.get(i)];
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

            int active_signal = uniform_int_distribution<>(0, e.signals.size() - 1)(gen);
            update_signals_on_add(e, active_signal);
        }
    }

    static unsigned dsuGetRoot(vector<unsigned> & dsu, unsigned v) {
        auto &p = dsu[v];
        return p == v ? v : p = dsuGetRoot(dsu, p);
    }

    bool Graph::is_connected(Edge const& erased) {
        // cout << "is connected " << erased.to_string() << endl;
        inner_update(erased.id, true);

        bfsQueue.clear();
        
        fill(bfsUsed.begin(), bfsUsed.end(), make_pair(false, false));
        bfsUsed[erased.first] = {true, 0};
        bfsUsed[erased.second] = {true, 1};
        bfsQueue.emplace_back(erased.first, 0);
        bfsQueue.emplace_back(erased.second, 1);
        int count[2] = {1,1};
        bool connected = false;
        for (size_t i = 0; i < bfsQueue.size(); ++i) {
            auto v = bfsQueue[i];
            count[v.second]--;
            // cout << "v: " << v << endl;
            for (auto const& tov : neis[v.first]) {
                auto to = tov.first;
                // cout << "to: " << to << endl;
                if (bfsUsed[to].first) {
                    if (bfsUsed[to].second != v.second) {
                        connected = true;
                        break;
                    } else continue;
                } 
                bfsQueue.emplace_back(to, v.second);
                count[v.second]++;
                bfsUsed[to] = {1, v.second};
                
            }
            if (connected || count[0] == 0 || count[1] == 0) break;
        }
        inner_update(erased.id, false);
        return connected;
    }

    void Graph::update_outer_edges(unsigned cand_in, unsigned cand_out) {
        // cout << "update outer edges" << endl;
        outer.swap(cand_out, cand_in);
        auto out_neighbour_edges = get_neighbours_edge(cand_out);
        for (auto const& edge : out_neighbour_edges) {
            size_t neighbour = edge.second;
            if (in_nei_c[neighbour]++ == 0 && neighbour != cand_in) {
                if (!inner.contains(neighbour)) {
                    outer.insert(neighbour);
                }
            }
        }
        auto in_neighbour_edges = get_neighbours_edge(cand_in);
        for (auto const& edge : in_neighbour_edges) {
            size_t neighbour = edge.second;
            if (--in_nei_c[neighbour] == 0 && neighbour != cand_out) {
                if (!inner.contains(neighbour)) {
                    outer.erase(neighbour);
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
    }

    void Graph::update_signals_on_erase(Edge & edge) {
        // cout << "update signals on erase" << endl;
        update_signals(edge.signals[edge.active_signal], -1);
        edge.active_signal = -1;
    }
    
    void Graph::update_signals_on_add(Edge & edge, int signal) {
        // cout << "update signals on add" << endl;
        edge.active_signal = signal;
        update_signals(edge.signals[signal], 1);
    }

    void Graph::update_signals(string const& signal, int diff) {
        auto it = signals.find(signal);
        it->second = {it->second.first, it->second.second + diff};
        if (it->second.second < 0) {
            throw domain_error("using signal less than 0 times");
        }
    }

    void Graph::inner_update(unsigned e, bool is_erased) {
        // cout << "inner update, e: " << e << "; erased: " << is_erased << endl;
        
        if (is_erased) {
            inner.erase(e);
            Edge const& edge = list_edges[e];
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
            remove_vertex_from_neis(v1->first, v1->second);
            remove_vertex_from_neis(v2.first, v2.second);
        
        } else {
            if (inner.contains(e)) {
                throw domain_error("try to add in inner existing edge");
            } 
            inner.insert(e);
            Edge inserted_edge = list_edges[e];
            neis[inserted_edge.first].emplace_back(inserted_edge.second, neis[inserted_edge.second].size());
            neis[inserted_edge.second].emplace_back(inserted_edge.first, neis[inserted_edge.first].size() - 1);
        }
    }

    void Graph::update_neighbours(unsigned e, bool is_erased) {
        auto neighbours = get_neighbours_edge(e);
        if (is_erased) {
            if (inner.size() != 0)
                outer.insert(e);
            for (auto const& neighbour : neighbours) {
                if (--in_nei_c[neighbour.second] == 0) {
                    if (!inner.contains(neighbour.second)) {
                        outer.erase(neighbour.second);
                    }
                }
            }
        } else {
            if (inner.size() != 1)
                outer.erase(e);
            for (auto const& neighbour : neighbours) {
                if (in_nei_c[neighbour.second]++ == 0) {
                    if (!inner.contains(neighbour.second)) {
                        outer.insert(neighbour.second);
                    }
                }
            }
        }
    }

    size_t Graph::peek_signal(vector<string> const& signals_cand) {
        if (signals_cand.size() > 1) {
            double sum_llh_out = 0;
            for (string const& signal: signals_cand) {
                sum_llh_out += signals[signal].first;
            }
            double peek_signal = unirealdis(gen) * sum_llh_out;
            double cur_peek = 0;
            for (size_t i = 0; i < signals_cand.size(); ++i) {
                cur_peek += signals[signals_cand[i]].first;
                if (peek_signal <= cur_peek) {
                    return i;
                }
            }
        }
        return 0;
    }

    pair<double, int> Graph::probability_on_change_vertex(size_t cand_in, size_t cand_out, size_t cur_size_outer, size_t new_size_outer) {
        vector<string> const& signals_out= list_edges[cand_out].signals;
        size_t signal_out_ind = uniform_int_distribution<>(0, signals_out.size() - 1)(gen);//peek_signal(signals_out);
        Edge const& in_edge = list_edges[cand_in];
        auto in_signal = signals[in_edge.signals[in_edge.active_signal]];
        auto out_signal = signals[signals_out[signal_out_ind]];

        double in_likelihood = (in_signal.second == 1) ? in_signal.first : 1;
        double out_likelihood = (out_signal.second == 0) ? out_signal.first : 1;

        double multisignals_penalty = log(exp(1) - 1 + signals_out.size());

        auto res = (out_likelihood * cur_size_outer) / (in_likelihood * new_size_outer * multisignals_penalty); 
        // cout << "in llh " << in_likelihood << " out llh " << out_likelihood  << " mult penalty " << multisignals_penalty << " cur size " << cur_size_outer << " new size " << new_size_outer << " prob " << res << endl;
        return {(out_likelihood * cur_size_outer) / (in_likelihood * new_size_outer * multisignals_penalty), signal_out_ind};
    }

    pair<double, int> Graph::probability_on_change_vertex(size_t cand, size_t cur_size, size_t new_size, bool erase) {
        // cout << "probability floating cand " << cand << " , cur_size: " << cur_size << " , new size: " << new_size << " , erase: " << erase << endl;
        size_t signal_ind;
        Edge const& cand_edge = list_edges[cand];
        if (erase) {
            signal_ind = cand_edge.active_signal;
            if (cand_edge.active_signal == -1) {
                throw domain_error("active signal is -1 while erasing edge: " + cand_edge.to_string());
            }
        } else {
            signal_ind = uniform_int_distribution<>(0, cand_edge.signals.size() - 1)(gen);
        }

        auto signal = signals[cand_edge.signals[signal_ind]];

        
        double multisignals_penalty = 1;
        double size_penalty = 1;
        double likelihood = 1;
        double sizes_relation = (double)cur_size/ (double)new_size;
        if (erase) {
            likelihood = (signal.second == 1) ? 1.0/signal.first : 1;
            size_t new_signal_ind = uniform_int_distribution<>(0, cand_edge.signals.size() - 1)(gen);
            if (signal_ind != new_signal_ind) {
                auto new_signal = signals[cand_edge.signals[new_signal_ind]];
                likelihood *= (new_signal.second == 0) ? new_signal.first : 1;
                sizes_relation = 1;
                signal_ind = new_signal_ind;
            } else {
                size_penalty = 1.0 / edge_penalty;
                signal_ind = -1;
            }
        } else {

            likelihood = (signal.second == 0) ? signal.first : 1;
            size_penalty = edge_penalty;
//            multisignals_penalty = 1.0 / cand_edge.signals.size();
        }
        double res = likelihood * sizes_relation * size_penalty * multisignals_penalty;
        if (cand_edge.signals.size() != 1) {
                // cout << "signal: " << signal_ind << " erase: " << erase << ", signals size: " << cand_edge.signals.size() << ", multisignals penalty: " << multisignals_penalty << ", res: " << res << endl;
            }
        
        return {likelihood * sizes_relation * size_penalty, signal_ind};
    }

    bool Graph::next_iteration() {
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
            size_t cur_size_outer = outer.size();
            size_t new_size_outer;
            if (inner.size() == 1) {            
                new_size_outer = edge_neighbours_size(cand_in);
            } else {
                new_size_outer = outer.size();
                auto in_neighbour_edges = get_neighbours_edge(cand_in);
                auto out_neighbour_edges = get_neighbours_edge(cand_out);
                for (auto const& x : in_neighbour_edges) {
                    if (in_nei_c[x.second] == 1) --new_size_outer;
                }
                for (auto const& x : out_neighbour_edges) {
                    if (!in_nei_c[x.second]) ++new_size_outer;
                }
            }
            pair<double, int> p_change_vertex = probability_on_change_vertex(cand_in, cand_out, cur_size_outer, new_size_outer);
            if (gen_p >= p_change_vertex.first)
                return false;
            Edge & in_edge = list_edges[cand_in];
            if (!is_connected(in_edge)) {
                return false;
            }
            inner_update(cand_in, true);
            update_outer_edges(cand_in, cand_out);
            update_signals_on_add(list_edges[cand_out], p_change_vertex.second);
            update_signals_on_erase(in_edge);
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
                int new_in_out_size = inner.size() + outer.size() - edge_neighbours_size(cand) + in_nei_c[cand];
                if (inner.size() == 1) {
                    new_in_out_size = order;
                } else if (cur_in_out_size <= edge_neighbours_size(cand) - in_nei_c[cand]) {
                    new_in_out_size = 1;
                    cout << "cand: " << cand << " cur in out: " << cur_in_out_size << endl;
                    print_vec<pair<size_t,size_t>>(get_neighbours_edge(cand), &pair_to_str);
                    throw domain_error("count of innner and outer edges less then count outer neighbours of some edge");
                }
                pair<double, int> p_change_vertex = probability_on_change_vertex(cand, cur_in_out_size, new_in_out_size, true);
                int peeked_signal = p_change_vertex.second;
                if (gen_p >= p_change_vertex.first) {
                    return false;
                }
                if (peeked_signal != -1) {
                    update_signals_on_erase(list_edges[cand]);
                    update_signals_on_add(list_edges[cand], peeked_signal);
                    return true;
                }
                if (!is_connected(list_edges[cand])) {
                    return false;
                }
            } else {
                double estimation = (inner.size() == 0 ? 1.0 * order / (1 + edge_neighbours_size(cand)) : 1);
                pair<double, int> p_change_vertex = probability_on_change_vertex(cand, estimation, 1, false);
                if (gen_p >= p_change_vertex.first) {
                    return false;
                }
            }
            inner_update(cand, erase);
            update_neighbours(cand, erase);
            unsigned new_in_out_size = inner.size() == 0 ? order : outer.size() + inner.size();
            pair<double, int> p_change_vertex = probability_on_change_vertex(cand, cur_in_out_size, new_in_out_size, erase);
            if (gen_p < p_change_vertex.first) {
                if (erase) {
                    update_signals_on_erase(list_edges[cand]);
                } else {
                    update_signals_on_add(list_edges[cand], p_change_vertex.second);
                }
                return true;
            }
            inner_update(cand, !erase);
            update_neighbours(cand, !erase);
            return false;
        }
    }

    vector <size_t> Graph::get_inner_edges() {
        return inner.get_all();
    }

    vector <size_t> Graph::get_outer_edges() {
        return outer.get_all();
    }

    unordered_map<string, pair<double, int>>  const& Graph::get_signals() const {
        return signals;
    }

    int Graph::get_active_signal_by_edge(size_t e) {
        return list_edges[e].active_signal;
    }
};
