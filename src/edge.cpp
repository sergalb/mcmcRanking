#include "mcmc.h"

namespace mcmc
{

    Edge::Edge(size_t first, size_t second, vector<string> signals, size_t id): first(first), second(second), signals(signals), id(id) {

    }

    string Edge::to_string() {
        return "first: " + std::to_string(first) + ", second: " + std::to_string(second) + ", signals[0]: " + signals[0] + "signals size: "+ std::to_string(signals.size()) + ", id: " + std::to_string(id);
    }
};
