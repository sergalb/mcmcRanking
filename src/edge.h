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
        vector<string> signals;
        size_t id;
        int active_signal = -1;   
        Edge(size_t first, size_t second, vector<string> signals, size_t id);
        string to_string();
    };
}

#endif //MCMC_RANKING_EDGE