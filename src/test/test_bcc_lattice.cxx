#include <random>
#include <set>
#include <sstream>
#include "../graph/bcc_lattice.hxx"

using namespace Dimers;

#define BOX_SIZE 5

int main(void) {
    // Generate std gen
    std::default_random_engine gen;

    // 0. Check if we are adding the correct edges
    BccLattice bccL(BOX_SIZE, BOX_SIZE, BOX_SIZE, &gen, 0.0);
    for (int x = 0; x < BOX_SIZE; x++) {
        for (int y = 0; y < BOX_SIZE; y++) {
            for (int z = 0; z < BOX_SIZE; z++) {
                int v = bccL.serialize_cubic_index(x, y, z);

                // FIXME: This is obviously ugly and needs to be fixed ASAP
                int wn, wp;

                wn = bccL.serialize_cubic_index(mod_close(x - 1, BOX_SIZE), y,
                                                z);
                wp = bccL.serialize_cubic_index(mod_close(x + 1, BOX_SIZE), y,
                                                z);
                if (!bccL.has_edge(v, wn)) {
                    LOG(FATAL) << "Missing cubic-cubic edge [x-neg] (" << v
                               << "," << wn << ")";
                }

                if (!bccL.has_edge(v, wp)) {
                    LOG(FATAL) << "Missing cubic-cubic edge [x-pos] (" << v
                               << "," << wp << ")";
                }

                wn = bccL.serialize_cubic_index(x, mod_close(y - 1, BOX_SIZE),
                                                z);
                wp = bccL.serialize_cubic_index(x, mod_close(y + 1, BOX_SIZE),
                                                z);
                if (!bccL.has_edge(v, wn)) {
                    LOG(FATAL) << "Missing cubic-cubic edge [y-neg] (" << v
                               << "," << wn << ")";
                }

                if (!bccL.has_edge(v, wp)) {
                    LOG(FATAL) << "Missing cubic-cubic edge [y-pos] (" << v
                               << "," << wp << ")";
                }

                wn = bccL.serialize_cubic_index(x, y,
                                                mod_close(z - 1, BOX_SIZE));
                wp = bccL.serialize_cubic_index(x, y,
                                                mod_close(z + 1, BOX_SIZE));
                if (!bccL.has_edge(v, wn)) {
                    LOG(FATAL) << "Missing cubic-cubic edge [z-neg] (" << v
                               << "," << wn << ")";
                }
                if (!bccL.has_edge(v, wp)) {
                    LOG(FATAL) << "Missing cubic-cubic edge [z-pos] (" << v
                               << "," << wp << ")";
                }

                if ((x < BOX_SIZE - 1) && (y < BOX_SIZE - 1) &&
                    (z < BOX_SIZE - 1)) {
                    int v_bcc = bccL.serialize_bcc_index(x, y, z);
                    for (int dx = 0; dx < 2; dx++) {
                        for (int dy = 0; dy < 2; dy++) {
                            for (int dz = 0; dz < 2; dz++) {
                                // We don't need to do mod close since bcc
                                // vertices don't have edges across a gluing
                                // face
                                int xn = mod_close(x + dx, BOX_SIZE);
                                int yn = mod_close(y + dy, BOX_SIZE);
                                int zn = mod_close(z + dz, BOX_SIZE);

                                int v_cubic =
                                    bccL.serialize_cubic_index(xn, yn, zn);
                                if (!bccL.has_edge(v_bcc, v_cubic)) {
                                    LOG(WARNING) << "X: " << x << " Y: " << y
                                                 << " Z: " << z;
                                    LOG(FATAL)
                                        << "Missing bcc-cubic edge (" << v_bcc
                                        << "," << v_cubic << ")";
                                }
                            }
                        }
                    }
                    LOG(INFO) << "Body-Centered Vertex " << v_bcc << " passed";
                }
            }
        }
    }

    // Run the simulation and make sure that there is no overlap
    LOG(INFO) << "About to check other densities";
    for (int nmonomers = 1; nmonomers < BOX_SIZE * BOX_SIZE; nmonomers *= 2) {
        double density =
            static_cast<double>(nmonomers) / (BOX_SIZE * BOX_SIZE * BOX_SIZE);
        LOG(INFO) << "About to check density " << density;
        BccLattice bccLm(BOX_SIZE, BOX_SIZE, BOX_SIZE, &gen, density);

        DLOG(INFO) << "Initial monomers";
#ifndef NDEBUG
        bccLm.print_monomers();
#endif
        for (int ts = 0; ts < 1000; ts++) {
            if (ts % 50 == 0) {
                DLOG(INFO) << "Finished with time step " << ts
                           << " for density " << density;
            }
            bccLm.check_overlaps();
            bccLm.run_one_timestep();
#ifndef NDEBUG
            bccLm.print_monomers();
#endif
        }
    }
}
