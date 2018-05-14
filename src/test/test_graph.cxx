#include <random>
#include <set>
#include <sstream>
#include <unordered_map>
#include "../graph/graph.hxx"

using namespace Dimers;

void test_graph(DimerGraph &g, bool verbose) {
    auto V = g.V;
    auto MD = g.MAX_DEGREE;

    if (verbose)
        LOG(INFO)
            << "About to generate random dimer configuration with density 0.25";
    g.generate_random_dimer_configuration(0.25);
    if (verbose)
        LOG(INFO) << "Generated configuration with " << g._nmonomers
                  << " monomers";

    std::vector<int> monomers;
    g.get_monomers(monomers);

    if (monomers.size() != (unsigned)g._nmonomers) {
        LOG(FATAL) << "Mismatch between nmonomers " << g._nmonomers
                   << " and result of get_monomers().size(): "
                   << monomers.size();
    }

    if (verbose) {
        std::stringstream ss;

        for (auto m : monomers) {
            ss << m << " ";
        }

        LOG(INFO) << "Monomers: " << ss.str();
    }

    if (verbose) LOG(INFO) << "Testing non-adjacent swap";
    int v = -1, w = -1;
    int i = 0, j = 0, ctr = 0;
    while ((v < 0) && (w < 0)) {
        if (!g.has_edge(i, j)) {
            v = i;
            w = j;
            break;
        }

        ctr++;
        j++;
        j = j % MD;
        if (j == 0) {
            i++;
        }
        if (ctr > V * V) {
            LOG(FATAL) << "Graph is either complete or malformed";
        }
    }
    g.perform_swap(v, w);

    int non_monomer = 0;
    while (non_monomer == monomers[non_monomer]) {
        non_monomer++;
    }

    if (verbose)
        LOG(INFO) << "Testing non-monomer swap on index " << non_monomer;
    g.perform_swap(non_monomer, non_monomer + MD);

    if (verbose) LOG(INFO) << "Testing all admissible monomer-dimer swaps";
    int sites = 0;
    for (auto monomer : monomers) {
        for (int i = 0; i < MD; i++) {
            int nb = g._nodes[monomer].edges[i];
            if (nb < 0) {
                break;
            }

            if (g.is_monomer(nb)) {
                continue;
            }

            int nb_nb = g.get_dimer_nbor(nb);

#ifdef _DEBUG
            if (verbose) {
                LOG(INFO) << "Monomer: " << monomer << ", "
                          << "Neighbor: " << nb << ", "
                          << "Neighbor of Neighbor: " << nb_nb;
            }
#endif

            // If the move is a rotation (e.g. has_edge(monomer,nb_nb) is true),
            // Then monomer is swapped with nb
            // Otherwise monomer is swapped with nb_nb
            int swapped = (g.has_edge(monomer, nb_nb)) ? nb : nb_nb;
            int unswapped = (g.has_edge(monomer, nb_nb)) ? nb_nb : nb;
            g.perform_swap(monomer, nb);

            if (g.is_monomer(monomer) || !g.is_monomer(swapped)) {
                LOG(FATAL) << "Attempted to swap monomer at " << monomer
                           << " with dimer (" << nb << "," << nb_nb
                           << ") but after swap, neighbor is not a monomer; "
                           << "is_monomer(monomer)=" << g.is_monomer(monomer)
                           << "; is_monomer(swapped)=" << g.is_monomer(swapped);
            }

            if (g.get_dimer_nbor(monomer) != unswapped) {
                LOG(FATAL) << "After performing swap of monomer " << monomer
                           << " with dimer (" << nb << "," << nb_nb
                           << "), neighbor of monomer was not " << unswapped
                           << "; has_edge(monomer, nb_nb)="
                           << g.has_edge(monomer, nb_nb);
            }

            if (g.get_dimer_nbor(unswapped) != monomer) {
                LOG(FATAL) << "After performing swap of monomer " << monomer
                           << " with dimer (" << nb << "," << nb_nb
                           << "), neighbor of unswapped, " << unswapped
                           << "was not monomer"
                           << "; has_edge(monomer, nb_nb)="
                           << g.has_edge(monomer, nb_nb);
            }

            // Finally, reverse the swap and make sure that nothing changed
            int to_swap = (g.has_edge(monomer, nb_nb)) ? monomer : nb;
            g.perform_swap(swapped, to_swap);

            if (!g.is_monomer(monomer) || g.is_monomer(nb)) {
                LOG(FATAL) << "Attempted to unswap monomer at " << monomer
                           << " with neighbor " << nb
                           << " but after swap, neighbor is a monomer or "
                              "monomer is not a monomer; "
                           << "is_monomer(monomer)=" << g.is_monomer(monomer)
                           << "; is_monomer(nb)=" << g.is_monomer(nb);
            }

            if (g.get_dimer_nbor(nb) != nb_nb) {
                LOG(FATAL) << "After performing unswap of monomer " << monomer
                           << " with dimer (" << nb << "," << nb_nb
                           << "), neighbor of nb was not " << nb_nb;
            }

            if (g.get_dimer_nbor(nb_nb) != nb) {
                LOG(FATAL) << "After performing swap of monomer " << monomer
                           << " with dimer (" << nb << "," << nb_nb
                           << "), neighbor of nb_nb was not " << nb;
            }
        }
        LOG(INFO) << "Finished with " << ++sites << " sites";
    }

    if (verbose)
        LOG(INFO) << "Testing if the dimer covering is actually a matching";

    // Map stores (node, node_who_added_me)
    std::unordered_map<int, int> dimer_map;
    std::set<int> mmers;
    for (int i = 0; i < V; i++) dimer_map[i] = -1;

    for (int i = 0; i < V; i++) {
        auto n = g._nodes[i];
        switch (n.dimer_nbor) {
            case -2:
                LOG(FATAL) << "Found an unassigned node at " << i;
                break;
            case -1:
                mmers.insert(i);
                break;
            default:
                int nb = g.get_dimer_nbor(i);

                if (mmers.find(nb) != mmers.end()) {
                    LOG(FATAL) << "Dimer neighbor " << nb << " of " << i
                               << " is also a monomer?!?";
                }

                if ((dimer_map[i] == nb) && (dimer_map[nb] == nb)) {
                    continue;
                }

                if ((dimer_map[i] != nb) && (dimer_map[i] != -1)) {
                    LOG(FATAL)
                        << "Vertex " << i
                        << " has two edges incident on it with "
                        << "other vertices " << nb << " and " << dimer_map[i];
                }

                if ((dimer_map[nb] != nb) && (dimer_map[nb] != -1)) {
                    LOG(FATAL) << "Neighbor " << nb << " of vertex " << i
                               << " has another neighbor " << dimer_map[nb];
                }

                dimer_map[i] = dimer_map[nb] = i;
                break;
        }
    }
}

inline int shitty_mod(int q, int modulus) {
    return (q < 0) ? (q + modulus) : (q % modulus);
}

int main(void) {
    // Get std c++ gen
    std::default_random_engine generator;

    // Test 1: 2D 5x4 rectangle
    DimerGraph rectangle(20, 4, &generator, false);

    for (int i = 0; i < 20; i++) {
        // Right Edge
        if ((i + 1) % 4 != 0) {
            rectangle.add_edge(i, i + 1);
        }
        // Up Edge
        if (i < 16) {
            rectangle.add_edge(i, i + 4);
        }
            // To test dupe edges:
#ifdef __DEBUG
            // Left Edge
            // if ( (i>0) && ((i-1) % 4 != 3)) {
            //     g.add_edge(i,i-1);
            // }
            // // Down Edge
            // if (i>4) {
            //     g.add_edge(i, i-4);
            // }
#endif
    }
    test_graph(rectangle, true);

    // Test 2: Periodic Triangular Graph
    DimerGraph periodic_triangular(25, 6, &generator, true);
    int length = 5;
    for (int r = 0; r < length; r++) {
        for (int c = 0; c < length; c++) {
            int ru = shitty_mod(r + 1, length);  // row above
            // int rd = shitty_mod(r-1, length); // row below
            int cr = shitty_mod(c + 1, length);  // col right
            int cl = shitty_mod(c - 1, length);  // col left
            int cul = (r % 2 == 0) ? cl : c;
            int cur = (r % 2 == 0) ? c : cr;

            int v = length * r + c;

            int right = length * r + cr;
            periodic_triangular.add_edge(v, right);
            // LOG(INFO) << "Adding right edge from (" << r << "," << c << ") to
            // (" << r << "," << cr << ")";

            // int l = length*r+cl;
            // periodic_triangular.add_edge(v,l);
            // LOG(INFO) << "Adding left edge from (" << r << "," << c << ") to
            // (" << r << "," << cl << ")";

            int ul = length * ru + cul;
            periodic_triangular.add_edge(v, ul);
            // LOG(INFO) << "Adding upper left edge from (" << r << "," << c <<
            // ") to (" << ru << "," << cul << ")";

            int ur = length * ru + cur;
            periodic_triangular.add_edge(v, ur);
            // LOG(INFO) << "Adding upper right edge from (" << r << "," << c <<
            // ") to (" << ru << "," << cur << ")";

            // int ll = length*rd+cul;
            // periodic_triangular.add_edge(v,ll);
            // LOG(INFO) << "Adding lower left edge from (" << r << "," << c <<
            // ") to (" << rd << "," << cul << ")";

            // int lr = length*rd+cur;
            // periodic_triangular.add_edge(v,lr);
            // LOG(INFO) << "Adding lower right edge from (" << r << "," << c <<
            // ") to (" << rd << "," << cur << ")";
        }
    }
    test_graph(periodic_triangular, true);

    return 0;
}
