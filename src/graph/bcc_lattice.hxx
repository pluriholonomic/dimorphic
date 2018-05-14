#ifndef __BCC_LATTICE_HXX__
#define __BCC_LATTICE_HXX__

// Utils
#include "../util/blosc/blosc_filter.h"
#include "../util/mod.hxx"
#include "../util/writer.hxx"

// Graph-based includes
#include "graph.hxx"
#include "ilattice.hxx"

#define BCC_MAX_DEGREE 14

// Helper functions
namespace {}

namespace Dimers {  //
// BccLattice:
// Wrapper for the periodic body-centered cubic (bcc) lattice
// for monomer-dimer simulations.
//
// A bcc lattice of size Nx * Ny * Nz contains:
// a. [Main] Cubic lattice with Nmain = Nx * Ny * Nz sites
// b. [Body-Centered] Cubic lattice with Nbcc = (Nx-1) * (Ny-1) * (Nz-1)
// That is, there are technically N = Nmain+Nbcc nodes, but we will use the
// number of standard cubic vertices to provide us with a linear size.
//
// We will use the convention that the first Nmain nodes are those of
// the main cubic part, while the remaining Nbcc nodes will be the body-centered
// portion of the lattice.
//
class BccLattice : public ILattice {
    class BccWriter : public Writer<DimerGraph::DimerNode> {
       public:
        BccWriter(std::string outfile, int *dims, int ndims, int nuids,
                  bool use_blosc = false)
            : Writer<DimerGraph::DimerNode>(outfile, dims, ndims, use_blosc) {}
        void write(DimerGraph::DimerNode *, int *) {
            LOG(WARNING)
                << "Attempting to use deprecated write mode in class BccWriter";
        }
        void write_node(BccLattice *lptr, DimerGraph &g);
    };

   public:
    BccLattice(int Nx, int Ny, int Nz, std::default_random_engine *gen,
               double monomer_density, std::string outfile = "");

    BccLattice(int Nx, int Ny, int Nz, std::default_random_engine *gen,
               int nmonomers, std::string outfile = "")
        : BccLattice(Nx, Ny, Nz, gen,
                     static_cast<double>(nmonomers) / (Nx * Ny * Nz),
                     outfile){};

    ~BccLattice() {
        delete _monomer_prng_ptr;
        if (_w) delete _w;
    }

    // Screw the rule of three for this --- this should never be copied!
    BccLattice(BccLattice &) = delete;

    // From interface, ILattice
    void run_one_timestep();
    void run_sequential(long time_to_run) {
        for (long i = 0; i < time_to_run; i++) {
            run_one_timestep();
        }
    }
    void write_lattice_to_hdf5();

    // DEPRECATED:
    // void write_lattice_to_hdf5() {
    //     DLOG(WARNING) << "Writing to HDF5 not implemented yet for
    //     BccLattice";
    // }

    void write_lattice_to_text(FILE *, FILE *) {
        DLOG(WARNING) << "Writing to text file not supported for BccLattice";
    }

    // deserialize: i -> (x,y,z)
    // serialize: (x,y,z) -> i
    void deserialize_cubic_index(int i, int &x, int &y, int &z) {
        x = i % Nx;
        y = ((i - x) / Nx) % Ny;
        z = i / (Nx * Ny);
    }

    int serialize_cubic_index(int x, int y, int z) {
        return Ny * Nx * z + Nx * y + x;
    }

    void deserialize_bcc_index(int i, int &x, int &y, int &z) {
        int j = i - Nmain;
        x = j % (Nx - 1);
        y = ((j - x) / (Nx - 1)) % (Ny - 1);
        z = j / ((Nx - 1) * (Ny - 1));
    }

    int serialize_bcc_index(int x, int y, int z) {
        return (Ny - 1) * (Nx - 1) * z + (Nx - 1) * y + x + Nmain;
    }

    void get_direction(int x, int y, int z, int nx, int ny, int nz, double *ret,
                       bool is_cubic, bool nbor_is_cubic);

    bool has_edge(int v, int w) { return _graph.has_edge(v, w); }

    void check_overlaps();

    // size accessors
    int get_nuids() const { return Nuids; }
    int get_nmain() const { return Nmain; }
    int get_nx() const { return Nx; }
    int get_ny() const { return Ny; }
    int get_nz() const { return Nz; }

#ifndef NDEBUG
    void print_monomers() {
        std::stringstream ss;
        ss << "We have the monomers: ";
        for (auto m : _monomers) {
            ss << m << "; ";
        }
        DLOG(INFO) << ss.str();
    }
#endif

   private:
    // Former template parameters
    int Nx, Ny, Nz;
    int Nmain, Nbcc, N, Nuids;

    void initialize_edges();

    // N.B. Each of the normal cubic edges have degree 14, whereas
    //      the body centered edges have degree 8.
    DimerGraph _graph;
    std::vector<int> _monomers;
    random_int *_monomer_prng_ptr;

    int _nmonomers, _ndimers;
    int _timestep;  // N.B. We're never going to run more than a billion steps,
                    // right? :p

    BccWriter *_w;
};

}  // namespace Dimers

#endif  // __BCC_LATTICE_HXX__
