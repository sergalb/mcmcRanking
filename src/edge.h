#ifndef MCMC_RANKING_EDGE
#define MCMC_RANKING_EDGE

namespace mcmc
{
    using namespace std;

    class Edge
    {
    public:
        size_t first;
        size_t second;
        string signal;
        size_t id;   
        Edge(size_t first, size_t second, string signal, size_t id);
        string to_string();
    };
}

#endif //MCMC_RANKING_EDGE
