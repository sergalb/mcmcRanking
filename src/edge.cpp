#include "mcmc.h"

namespace mcmc
{

    Edge::Edge(size_t first, size_t second, vector<string> signals, size_t id): first(first), second(second), signals(signals), id(id), active_signal(-1) {        
    }

    string Edge::to_string() const {
        return "first: " + std::to_string(first) + ", second: " + std::to_string(second) + " active signal: " + std::to_string(active_signal) + " , signal: " + signals[max(active_signal, 0)] + " signals size: "+ std::to_string(signals.size()) + ", id: " + std::to_string(id);
    }
};
