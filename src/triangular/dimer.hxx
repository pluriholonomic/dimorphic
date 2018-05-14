#pragma once

#include <glog/logging.h>
#include <cassert>
#include <random>
#include <vector>

#include "../util/prng.hxx"
#include "../util/writer.hxx"

namespace Dimers {

//
// Static Defs
//
static const bool Black = true;
static const bool White = false;

// N.B. We store the directions around in the hexagonal neighborhood
//      in clockwise order, as this makes it easy to do some computations
//      NW      NE
//  W       x        E
//      SW      SE
//
//  A clockwise rotation is +1, counterclockwise is -1
//
enum class Neighbor { NW = 0, NE, E, SE, SW, W, None, Unoccupied };
Neighbor find_opposite_direction(const Neighbor nb);

// FIXME: Check if the compiler inlines these
static inline Neighbor rotate_neighbor(const Neighbor nb, int offset) {
    return static_cast<Neighbor>(((static_cast<int>(nb) + 6 + offset) % 6));
}

static inline Neighbor rotate_clockwise(const Neighbor nb) {
    return rotate_neighbor(nb, 1);
}

static inline Neighbor rotate_counterclockwise(const Neighbor nb) {
    return rotate_neighbor(nb, -1);
}

static inline int shitty_mod(int i, int offset, int modulus) {
    i += offset;
    if (i > modulus) {
        return i - modulus;
    }
    if (i < 0) {
        return i + modulus;
    }
    return i;
}

// RhombicLattice Operations
class RhombicLattice {
    struct RhombicLatticeSite {
        bool color; /* Black == True, White == False */
        Neighbor neighbor;
    };

    class TriangularWriter : public Writer<RhombicLatticeSite> {
       public:
        TriangularWriter(std::string outfile, int* dims, int ndims,
                         bool use_blosc = false)
            : Writer<RhombicLatticeSite>(outfile, dims, ndims, use_blosc, true),
              _nrows(dims[0]),
              _ncols(dims[1]) {}

        void write(RhombicLatticeSite* lattice, int* positions);

       private:
        int _nrows;
        int _ncols;
    };

   public:
    RhombicLattice(int nrows, int ncols, int ndefects, float alpha, float beta,
                   std::default_random_engine* _generator, std::string outfile,
                   bool use_writer = false);
    ~RhombicLattice() {
        if (_w) {
            delete _w;
        }
    }

    RhombicLattice(RhombicLattice&) = delete;

    const std::vector<RhombicLatticeSite>& get_lattice() const {
        return _lattice;
    }

    // Run various simulations
    void run_to_acceptance(long nsteps);
    void run_sequential(long nsteps);

    void write_lattice_to_text(FILE* lattice_file, FILE* defect_file);
    void write_lattice_to_hdf5();

   private:
    // MC moves --- assumes that (a,b) is a defect
    bool perform_single_metropolis_move(int a, int b);

    // Helper functions
    inline RhombicLatticeSite& get_site(int i, int j) {
        return _lattice[i * _nc + j];
    }
    inline const RhombicLatticeSite& get_site_const(int i, int j) const {
        return _lattice[i * _nc + j];
    }
    bool is_defect(int i, int j) const {
        return get_site_const(i, j).neighbor == Neighbor::None;
    }

    // Testing
    void check_neighbors();
    void check_defects() const;
    void test_swap(int i, int j);

    // DEPRECATED:
    // inline int&                get_defect(int i, int j) { return
    // _defects[i*_nc+j]; }
    inline Neighbor& get_neighbor(int i, int j) {
        return _lattice[i * _nc + j].neighbor;
    }
    inline Neighbor get_neighbor_const(int i, int j) const {
        return _lattice[i * _nc + j].neighbor;
    }
    inline void find_neighbor(int i, int j, int ni, int nj, Neighbor& nb) const;

    bool is_site_occupied(int i, int j) const {
        return (get_site_const(i, j).neighbor != Neighbor::Unoccupied);
    }
    void find_neighbor_index(int i, int j, int& nb_i, int& nb_j,
                             Neighbor nb) const;
    void find_neighbor_index(int i, int j, int& nb_i, int& nb_j) const;
    bool all_neighbors_defects(int i, int j) const;
    bool is_in_neighborhood(int i, int j, int pi, int pj) const;
    int count_defects_in_neighborhood(int i, int j) const;

    // Swap functions
    int swap(int defect_i, int defect_j, Neighbor nb, int& dest_a,
             int& dest_b);  // FIXME: Implement

    // perform the swap on lattice
    // FIXME: Overoptimization --- convert swap into xor-only
    void swap_site(int src_i, int src_j, int dst_i, int dst_j);

    // Data; Since a 1K x 1K lattice takes up <8MB, we might as well store
    // everything in memory
    std::vector<RhombicLatticeSite> _lattice;
    std::vector<int> _defects;
    std::vector<int> _positions;

    // PRNG
    random_int _row_dist;
    random_int _col_dist;
    random_int _dir_dist;
    random_dbl _dbl_dist;

    // Static;
    const int _nr;
    const int _nc;
    int _nd;

    const float _alpha;
    const float _beta;

    int ndefects;
    int nexcess_defects;

    bool _use_writer;
    TriangularWriter* _w;
};

}  // namespace Dimers
