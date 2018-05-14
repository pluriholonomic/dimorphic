#ifndef __ILATTICE_HXX__
#define __ILATTICE_HXX__

namespace Dimers {
struct ILattice {
    virtual ~ILattice(){};
    virtual void run_one_timestep() = 0;
    virtual void run_sequential(long) = 0;
    virtual void write_lattice_to_hdf5() = 0;
    virtual void write_lattice_to_text(FILE*, FILE*) = 0;
#ifndef NDEBUG
    virtual void print_monomers() = 0;
#endif
};
}  // namespace Dimers

#endif  // __ILATTICE_HXX__
