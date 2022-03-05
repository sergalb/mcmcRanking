#include "mcmc.h"

namespace mcmc
{

    Edge::Edge(size_t from, size_t to, double weight, string signal) : from(from), to(to), weight(weight), signal(signal) {

    }
};
