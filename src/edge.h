#ifndef MCMC_RANKING_EDGE
#define MCMC_RANKING_EDGE

namespace mcmc
{
    using namespace std;

    class Edge
    {
    public:
        size_t from;
        size_t to;
        double weight;
        string signal;
//    public:
        Edge(size_t from, size_t to, double weight, string signal);
    };
}

#endif //MCMC_RANKING_EDGE
