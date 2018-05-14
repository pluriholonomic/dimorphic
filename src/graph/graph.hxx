#ifndef __GRAPH_HXX__
#define __GRAPH_HXX__

#include <glog/logging.h>
#include <cmath>  // for std::round
#include <sstream>
#include <vector>
#include "../util/prng.hxx"

namespace Dimers {

// DimerGraph is a simple graph structure that is made to
// hold a dimer configuration and can perform swaps
//
// Since these graphs are usually regular graphs, we're having
// each node store a list of its edges

struct DimerGraph {
    struct DimerNode {
        bool color;  // black == true, white == false
        int degree;
        int uid;  // Unique ID for keeping track of diffusion of dimers/monomers
        int dimer_nbor;  // Index into edges; is a monomer if == -1, is
                         // uninitialized if == -2
        // FIXME: We could go back to int arrays, I was just lazy
        std::vector<int> edges;
    };

    DimerGraph(int V, int MAX_DEGREE, std::default_random_engine *gen,
               bool verbose);
    DimerGraph(const DimerGraph &) = delete;

    void reset_graph();
    void add_edge(int v, int w);
    void perform_swap(int v, int w);
    void generate_random_dimer_configuration(double monomer_density);

    int get_dimer_nbor(int v) {
        return (is_monomer(v)) ? -1 : _nodes[v].edges[_nodes[v].dimer_nbor];
    }
    int get_edge_idx(int v, int w);

    bool has_edge(int v, int w);
    bool is_monomer(int v) { return _nodes[v].dimer_nbor == -1; }
    bool is_uninitialized(int v) { return _nodes[v].dimer_nbor == -2; }

    void draw_random_neighbor(int v, int &nb, int &nb_idx);

    void get_monomers(std::vector<int> &mv);

    // Former template parameters
    int V;
    int MAX_DEGREE;

    std::vector<DimerNode> _nodes;

    random_int _nbor_prng;
    random_int _vertex_prng;
    random_dbl _real_prng;

    int _nmonomers;

    bool _verbose;
};

}  // namespace Dimers

#endif  // __GRAPH_HXX__
