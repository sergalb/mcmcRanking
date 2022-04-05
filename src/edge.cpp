#include "mcmc.h"

namespace mcmc
{

    Edge::Edge(size_t first, size_t second, string signal, size_t id): first(first), second(second), signal(signal), id(id) {

    }

    string Edge::to_string() {
        return "first: " + std::to_string(first) + ", second: " + std::to_string(second) + ", signal: " + signal + ", id: " + std::to_string(id);
    }
};
