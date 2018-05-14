#include "bcc_lattice.hxx"

#include <algorithm>  // FIXME: Remove once we're done debugging
#include <set>
#include <unordered_map>

namespace Dimers {
BccLattice::BccLattice(int Nx, int Ny, int Nz, std::default_random_engine *gen,
                       double monomer_density, std::string outfile)
    : Nx(Nx),
      Ny(Ny),
      Nz(Nz),
      Nmain(Nx * Ny * Nz),
      Nbcc((Nx - 1) * (Ny - 1) * (Nz - 1)),
      N(Nmain + Nbcc),
      Nuids(-1),
      _graph(N, BCC_MAX_DEGREE, gen, false),
      _nmonomers(-1),
      _ndimers(-1),
      _timestep(0) {
    // Set up graph, tiling
    initialize_edges();
    _graph.generate_random_dimer_configuration(monomer_density);
    _nmonomers = _graph._nmonomers;
    if ((N - _nmonomers) % 2 != 0) {
        LOG(FATAL) << "BccLattice attempted to fill lattice with monomers, but "
                      "failed as N-_nmonomers="
                   << N - _nmonomers;
    }
    _ndimers = (N - _nmonomers) / 2;
    Nuids = _ndimers + _nmonomers;
    _monomer_prng_ptr = new random_int(_nmonomers, gen);

    // Set up writer
    if (outfile.length() > 0) {
        int dims[2] = {Nuids, 3};
        _w = new BccWriter(outfile, dims, 2, Nuids);
    } else {
        LOG(WARNING) << "Generating BccLattice object without an output file";
        _w = NULL;
    }

    LOG(INFO) << "Set up bcc lattice with " << _nmonomers << " monomers, "
              << _ndimers << " dimers, " << Nuids << " uids, " << Nmain
              << " cubic sites, " << Nbcc << " bcc sites, and " << N
              << " total sites";
}

void BccLattice::initialize_edges() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                int v = serialize_cubic_index(x, y, z);

                // Add normal cubic edges
                // FIXME:
                // I think that we only *technically* need to add the
                // positive octant of edges for each site --- make sure that
                // this is true
                int wx = serialize_cubic_index(mod_close(x + 1, Nx), y, z);
                int wy = serialize_cubic_index(x, mod_close(y + 1, Ny), z);
                int wz = serialize_cubic_index(x, y, mod_close(z + 1, Nz));

                // DLOG(INFO) << "About to add bcc edge (" << v << "," << wx <<
                // ")";
                _graph.add_edge(v, wx);
                // DLOG(INFO) << "About to add bcc edge (" << v << "," << wy <<
                // ")";
                _graph.add_edge(v, wy);
                // DLOG(INFO) << "About to add bcc edge (" << v << "," << wz <<
                // ")";
                _graph.add_edge(v, wz);

                if ((x < Nx - 1) && (y < Ny - 1) && (z < Nz - 1)) {
                    int v_bcc = serialize_bcc_index(x, y, z);
                    // Is a valid bcc vertex index
                    for (int dx = 0; dx < 2; dx++) {
                        for (int dy = 0; dy < 2; dy++) {
                            for (int dz = 0; dz < 2; dz++) {
                                // We don't need to do mod close since bcc
                                // vertices don't have edges across a gluing
                                // face
                                int xn = mod_close(x + dx, Nx);
                                int yn = mod_close(y + dy, Ny);
                                int zn = mod_close(z + dz, Nz);

                                int v_cubic = serialize_cubic_index(xn, yn, zn);
                                // DLOG(INFO) << "About to add edge (" << v_bcc
                                // << "," << v_cubic << ")";
                                _graph.add_edge(v_bcc, v_cubic);
                            }
                        }
                    }
                }
            }
        }
    }
}

void BccLattice::run_one_timestep() {
    // 0. Get monomer to swap
    // 1. Draw random neighbor
    // 2. Perform swap of monomer w/ random neighbor

    // FIXME: Inefficient a.f.
    _monomers.clear();
    _graph.get_monomers(_monomers);

    int monomer_idx = _monomer_prng_ptr->generate();
    int monomer_v = _monomers[monomer_idx];
    int rnd_nb, rnd_nb_idx;
    _graph.draw_random_neighbor(monomer_v, rnd_nb, rnd_nb_idx);

    DLOG(INFO) << "At timestep " << _timestep
               << " BEFORE PERFORMING swap between " << monomer_v
               << "(uid=" << _graph._nodes[monomer_v].uid
               << ", dimer_nbor=" << _graph._nodes[monomer_v].dimer_nbor
               << ") and " << rnd_nb << "(uid=" << _graph._nodes[rnd_nb].uid
               << ", dimer_nbor=" << _graph._nodes[rnd_nb].dimer_nbor
               << ") and "
               << ")";

    if (rnd_nb < 0 || rnd_nb_idx < 0) {
        std::stringstream ss;
        ss << "Unable to find admissible swap for proposal " << monomer_v << ":"
           << _graph.get_dimer_nbor(monomer_v) << "; neighbors: ";
        for (auto e : _graph._nodes[monomer_v].edges) {
            if (e < 0) break;
            ss << e << ":" << _graph.get_dimer_nbor(e) << " ";
        }
        LOG(FATAL) << ss.str();
    } else {
        _graph.perform_swap(monomer_v, rnd_nb);
        _monomers.clear();
        _graph.get_monomers(_monomers);  // FIXME: inefficient!
        DLOG(INFO) << "At timestep " << _timestep << " performed swap between "
                   << monomer_v << "(uid=" << _graph._nodes[monomer_v].uid
                   << ", dimer_nbor=" << _graph._nodes[monomer_v].dimer_nbor
                   << ") and " << rnd_nb << "(uid=" << _graph._nodes[rnd_nb].uid
                   << ", dimer_nbor=" << _graph._nodes[rnd_nb].dimer_nbor
                   << ") and "
                   << ")";
    }

    _timestep++;
}

void BccLattice::check_overlaps() {
    std::unordered_map<int, int> dimer_map;
    std::set<int> monomers;
    for (int i = 0; i < N; i++) dimer_map[i] = -1;

    for (int i = 0; i < N; i++) {
        auto n = _graph._nodes[i];
        switch (n.dimer_nbor) {
            case -2:
                LOG(FATAL) << "Found an unassigned node at " << i;
                break;
            case -1:
                monomers.insert(i);
                break;
            default:
                int nb = _graph.get_dimer_nbor(i);

                if (monomers.find(nb) != monomers.end()) {
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

void BccLattice::write_lattice_to_hdf5() {
    if (_w) {
        _w->write_node(this, _graph);
    } else {
        LOG(WARNING) << "Running BccLattice without a writer!";
    }
}

void BccLattice::get_direction(int x, int y, int z, int nx, int ny, int nz,
                               double *ret, bool is_cubic, bool nbor_is_cubic) {
    // Notes:
    // There are a total of 22 directions, 14 from the cubic vertices and 8 from
    // the bcc vertices.
    //
    //
    // Direction Vector == (nx,ny,nz)-(x,y,z)
    //
    // Each of the 22 direction will be stored as a double; each entry has
    // modulus either 0 or 1 [if cubic-cubic] or 1/sqrt(3) [if bcc-cubic,
    // cubic-bcc].
    //
    static const double noncubic_modulus = 1.0 / sqrt(3.0);
    int dx = nx - x;
    int dy = ny - y;
    int dz = nz - z;

    if ((dx != 0) && (abs(dx) != 1) && (abs(dx) != Nx - 1)) {
        LOG(FATAL) << "Incorrect x-direction of neighbor, dx=" << dx;
    }

    if ((dy != 0) && (abs(dy) != 1) && (abs(dy) != Ny - 1)) {
        LOG(FATAL) << "Incorrect x-direction of neighbor, dy=" << dy;
    }

    if ((dz != 0) && (abs(dz) != 1) && (abs(dz) != Nz - 1)) {
        LOG(FATAL) << "Incorrect x-direction of neighbor, dz=" << dz;
    }

    // If we are wrapping, force the deltas into {-1,1}
    dx = (abs(dx) == Nx - 1) ? (dx > 0) - (dx < 0) : dx;
    dy = (abs(dy) == Ny - 1) ? (dy > 0) - (dy < 0) : dy;
    dz = (abs(dz) == Nz - 1) ? (dz > 0) - (dz < 0) : dz;

    if (is_cubic && nbor_is_cubic) {
        // cubic-cubic edge
        if (!((dx && !dy && !dz) || (!dx && dy && !dz) || (!dx && !dy && dz))) {
            // FIXME: UGH. Things like this make me *almost* want to use
            // boost::format...
            LOG(FATAL) << "Attempting to compute a cubic-cubic direction "
                          "between points ("
                       << x << "," << y << "," << z << ") and nb (" << nx << ","
                       << ny << "," << nz << ") with dx=" << dx << " dy=" << dy
                       << " dz=" << dz;
        }

        ret[0] = dx;
        ret[1] = dy;
        ret[2] = dz;
    } else if (!is_cubic && nbor_is_cubic) {
        // bcc-cubic edge
        if (abs(dx) > 1 || abs(dy) > 1 || abs(dz) > 1) {
            LOG(FATAL) << "Attempting to compute a bcc-cubic direction between "
                          "points ("
                       << x << "," << y << "," << z << ") and nb (" << nx << ","
                       << ny << "," << nz << ") with dx=" << dx << " dy=" << dy
                       << " dz=" << dz;
        }

        ret[0] = (2 * dx - 1) * noncubic_modulus;
        ret[1] = (2 * dy - 1) * noncubic_modulus;
        ret[2] = (2 * dz - 1) * noncubic_modulus;
    } else if (is_cubic && !nbor_is_cubic) {
        // cubic-bcc edge
        if (abs(dx) > 1 || abs(dy) > 1 || abs(dz) > 1) {
            LOG(FATAL) << "Attempting to compute a cubic-bcc direction between "
                          "points ("
                       << x << "," << y << "," << z << ") and nb (" << nx << ","
                       << ny << "," << nz << ") with dx=" << dx << " dy=" << dy
                       << " dz=" << dz;
        }

        ret[0] = (2 * dx + 1) * noncubic_modulus;
        ret[1] = (2 * dy + 1) * noncubic_modulus;
        ret[2] = (2 * dz + 1) * noncubic_modulus;
    } else {
        // bcc-bcc edge (doesn't exist!)
        LOG(FATAL) << "Attempting to generate a direction for a bcc-bcc edge, "
                      "which cannot exist!";
    }
}

// FIXME FIXME FIXME: Needs a big refactoring, can be simplified a LOT
void BccLattice::BccWriter::write_node(BccLattice *lptr, DimerGraph &g) {
    int N = lptr->get_nuids();
    int Nmain = lptr->get_nmain();
    int Nx = lptr->get_nx();
    int Ny = lptr->get_ny();
    int Nz = lptr->get_nz();
    DLOG(INFO) << "BccWriter writing with N=" << N << ";Nx=" << Nx
               << ";Ny=" << Ny << ";Nz=" << Nz;

    // 0. Store the direction, position data in the format to be written
    //    Note that our output is of the format,:
    //    Position:  (time, uid, 3)
    //    Direction: (time, uid, 3)
    std::vector<bool> uids_written(N, false);
    std::vector<double> direction_data(3 * N, -2.0);
    std::vector<int> position_data(3 * N, -2);

    // TODO: Should this be a lambda? It is just kind of convenient...
    auto process_single_atom = [N, Nmain, Nx, Ny, Nz, lptr, &g, &direction_data,
                                &position_data, &uids_written](
                                   int x, int y, int z, bool is_cubic) {
        int serialized_id = (is_cubic) ? lptr->serialize_cubic_index(x, y, z)
                                       : lptr->serialize_bcc_index(x, y, z);
        auto node = g._nodes[serialized_id];
        int uid = node.uid;
        if (uids_written[uid]) {
            DLOG(INFO) << "[write_node:process_single_atom] Already wrote uid "
                       << uid << " at position (x=" << x << ",y=" << y
                       << ",z=" << z << ",is_cubic=" << is_cubic << ")";
            return;
        }

        bool dir_unwritten = true, pos_unwritten = true;
        for (int i = 0; i < 3; i++) {
            dir_unwritten &= direction_data[3 * uid + i] < -1.0;
            pos_unwritten &= position_data[3 * uid + i] < -1;
        }

        if (dir_unwritten && pos_unwritten) {
            // N.B.
            // 1. Position == position of black
            // 2. bcc point = cubic point + (0.5,0.5,0.5)
            // 3. We multiply all positions by 2, so that all points are
            // integral
            int nb = g.get_dimer_nbor(serialized_id);
            if ((nb >= Nmain) && (!is_cubic)) {
                LOG(FATAL)
                    << "Attempting to write a bcc-bcc edge, which is forbidden";
            }

            int nx, ny, nz;
            bool nb_is_cubic = nb < Nmain;

            if (nb_is_cubic) {
                lptr->deserialize_cubic_index(nb, nx, ny, nz);
            } else {
                lptr->deserialize_bcc_index(nb, nx, ny, nz);
            }

            if (node.color) {
                position_data[3 * uid + 0] = 2 * x + (!is_cubic);
                position_data[3 * uid + 1] = 2 * y + (!is_cubic);
                position_data[3 * uid + 2] = 2 * z + (!is_cubic);
            } else {
                position_data[3 * uid + 0] =
                    2 * nx + (is_cubic) * (!nb_is_cubic);
                position_data[3 * uid + 1] =
                    2 * ny + (is_cubic) * (!nb_is_cubic);
                position_data[3 * uid + 2] =
                    2 * nz + (is_cubic) * (!nb_is_cubic);
            }

            if (node.dimer_nbor >= 0) {
                int seridx = lptr->serialize_cubic_index(x, y, z);

                DLOG(INFO) << "is_cubic (" << is_cubic
                           << ") particle (uid=" << uid << ", seridx=" << seridx
                           << ") at (" << x << "," << y << "," << z
                           << ") has serialized neighbor " << nb
                           << " and deserialized neighbor (" << nx << "," << ny
                           << "," << nz << "); nb_is_cubic: " << nb_is_cubic;

                // N.B. (x,y,z) passed to get_direction is always black
                if (node.color) {
                    lptr->get_direction(x, y, z, nx, ny, nz,
                                        &direction_data[3 * uid], is_cubic,
                                        nb_is_cubic);
                } else {
                    lptr->get_direction(nx, ny, nz, x, y, z,
                                        &direction_data[3 * uid], nb_is_cubic,
                                        is_cubic);
                }
            } else {
                direction_data[3 * uid + 0] = -1.0;
                direction_data[3 * uid + 1] = -1.0;
                direction_data[3 * uid + 2] = -1.0;
            }

            uids_written[uid] = true;
        } else {
            LOG(FATAL) << "[write_node:process_single_atom] Already written "
                          "direction (dir_unwritten="
                       << dir_unwritten
                       << ") or position (pos_unwritten=" << pos_unwritten
                       << ") for uid " << uid;
        }
    };

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                process_single_atom(x, y, z, true);
                if ((x < Nx - 1) && (y < Ny - 1) && (z < Nz - 1))
                    process_single_atom(x, y, z, false);
            }
        }
    }

    // FIXME: Chomp off a lot of this, maybe make this a debug feature
    int all_positions_written = true;
    int npositions_written = 0;
    for (int i = 0; i < N; i++) {
        all_positions_written &= position_data[3 * i + 0] >= 0;
        all_positions_written &= position_data[3 * i + 1] >= 0;
        all_positions_written &= position_data[3 * i + 2] >= 0;
        auto all_written = (position_data[3 * i + 0] >= 0) *
                           (position_data[3 * i + 1] >= 0) *
                           (position_data[3 * i + 2] >= 0);
        npositions_written += all_written;
    }
    if (npositions_written != N) {
        LOG(FATAL) << "Only found " << npositions_written
                   << " positions that were written out of a total of " << N
                   << " required positions";
    }

    for (int i = 0; i < N; i++) {
        if (i < Nmain) {
            DLOG(INFO) << "cubic uid " << i << " is at position ("
                       << position_data[3 * i + 0] << ","
                       << position_data[3 * i + 1] << ","
                       << position_data[3 * i + 2] << ")";
        } else {
            DLOG(INFO) << "bcc uid " << i << " is at position ("
                       << position_data[3 * i + 0] << ","
                       << position_data[3 * i + 1] << ","
                       << position_data[3 * i + 2] << ")";
        }
    }

    // N.B. The bottom is mostly the same as what the triangular lattice uses
    // FIXME: If so, add this into the interfac,!
    // 1. Update offset, check extent
    offset[0] = _nwritten++;

    DLOG(INFO) << "Before write";
    if (_nwritten > current_dims[0]) {
        current_dims[0] += 32;
        status = H5Dset_extent(pos_dataset, current_dims);
        DLOG(INFO) << "[BccWriter] Status " << status
                   << " after set pos extent";
        status = H5Dset_extent(dir_dataset, current_dims);
        DLOG(INFO) << "[BccWriter] Status " << status
                   << " after set dir extent";
    }
#ifndef NDEBUG
    else {
        hsize_t extent[3] = {0, 0, 0};
        status = H5Sget_simple_extent_dims(pos_dataspace, extent, max_dims);
        printf(
            "[BccWriter] Current pos extent: (%llu,%llu,%llu), "
            "current_dims=(%llu,%llu,%llu), offset=(%llu,%llu,%llu)\n",
            extent[0], extent[1], extent[2], current_dims[0], current_dims[1],
            current_dims[2], offset[0], offset[1], offset[2]);
        status = H5Sget_simple_extent_dims(dir_dataspace, extent, max_dims);
        printf("[BccWriter] Current dir extent: (%llu,%llu,%llu)\n", extent[0],
               extent[1], extent[2]);
    }
#endif

    // 2. Write pos data
    DLOG(INFO) << "Before pos write; _ndims: " << _ndims;
    hid_t pos_filespace = H5Dget_space(pos_dataset);
    status = H5Sselect_hyperslab(pos_filespace, H5S_SELECT_SET, offset, NULL,
                                 frame_size, NULL);
    hid_t pos_memspace = H5Screate_simple(_ndims + 1, frame_size, max_dims);
    status = H5Dwrite(pos_dataset, H5T_NATIVE_INT, pos_memspace, pos_filespace,
                      H5P_DEFAULT, position_data.data());
    DLOG(INFO) << "[BccWriter] pos write status: " << status;
    H5Sclose(pos_filespace);
    H5Sclose(pos_memspace);

    // 3. Write dir data
    DLOG(INFO) << "Before dir write";
    hid_t dir_filespace = H5Dget_space(dir_dataset);
    status = H5Sselect_hyperslab(dir_filespace, H5S_SELECT_SET, offset, NULL,
                                 frame_size, NULL);
    hid_t dir_memspace = H5Screate_simple(_ndims + 1, frame_size, max_dims);
    status = H5Dwrite(dir_dataset, H5T_NATIVE_DOUBLE, dir_memspace,
                      dir_filespace, H5P_DEFAULT, direction_data.data());
    DLOG(INFO) << "[BccWriter] dir write status: " << status;
    H5Sclose(dir_filespace);
    H5Sclose(dir_memspace);
}
}  // namespace Dimers
