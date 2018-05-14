#include "dimer.hxx"

namespace Dimers {

// Helper
Neighbor find_opposite_direction(Neighbor nb) {
    return static_cast<Neighbor>((static_cast<int>(nb) + 3) % 6);
}

// RhombicLattice
RhombicLattice::RhombicLattice(int nrows, int ncols, int ndefects, float alpha,
                               float beta,
                               std::default_random_engine* generator,
                               std::string outfile, bool use_writer)
    : _row_dist(nrows, generator),
      _col_dist(ncols, generator),
      _dir_dist(6, generator),
      _dbl_dist(generator),
      _nr(nrows),
      _nc(ncols),
      _nd(ndefects),
      _alpha(alpha),
      _beta(beta),
      _use_writer(use_writer) {
    // Setup writer
    int dims[2] = {nrows, ncols};
    _w = new TriangularWriter(outfile, dims, 2);

    // Allocate memory, initialize
    // N.B.: Row Major
    _lattice.resize(nrows * ncols);
    _positions.resize(nrows * ncols);
    for (int i = 0; i < nrows * ncols; i++) {
        _lattice[i].neighbor = Neighbor::Unoccupied;
        _positions[i] = i;
    }

    // Setup defects array
    _defects.resize(nrows * ncols);
    memset(&_defects[0], 0, nrows * ncols);

    // Generate random defects
    int ndefects_generated = 0;
    while (ndefects_generated < ndefects) {
        printf("Generated %i defects out of %i defects\n", ndefects_generated,
               ndefects);
        int a = _row_dist.generate();
        int b = _col_dist.generate();

        if (!is_defect(a, b)) {
            printf("Place defect at site (%i, %i), is_defect(a,b)=%i\n", a, b,
                   is_defect(a, b));
            get_site(a, b).color = Black;
            get_site(a, b).neighbor = Neighbor::None;
            ndefects_generated++;
        }
    }

    // Generate Random Tiling among remaining tiles
    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            // If defect / occupied, continue
            if (is_defect(r, c) || is_site_occupied(r, c) ||
                all_neighbors_defects(r, c)) {
                if (is_defect(r, c)) {
                    printf(
                        "Found defect at (%i,%i) during initialization; nb: "
                        "%i\n",
                        r, c, get_neighbor(r, c));
                }
                continue;
            }

            int nb_i, nb_j;

            // FIXME: Handle the case in which all of a site's neighbors are
            // defects
            int ctr = 0;
            while (ctr < 36) {
                //        printf("ctr: %i\n", ctr);
                Neighbor nb = (Neighbor)_dir_dist.generate();
                find_neighbor_index(r, c, nb_i, nb_j, nb);

                // If we chose a defect or if the site is unoccupied try again
                if (nb_i < 0 || nb_j < 0 || is_site_occupied(nb_i, nb_j)) {
                    ctr++;
                    continue;
                }

                if (!is_defect(nb_i, nb_j)) {
                    // Draw a random bit to decide if we are black or write
                    // (FIXME: This is an expensive bit!)
                    bool my_color = _dir_dist.generate() % 2;
                    get_site(r, c).color = my_color;
                    get_site(r, c).neighbor = nb;
                    if (is_defect(r, c) || is_defect(nb_i, nb_j)) {
                        printf("WTF!!!\n");
                    }

                    get_site(nb_i, nb_j).color = !my_color;
                    get_neighbor(nb_i, nb_j) = find_opposite_direction(nb);

#ifdef __DEBUG
                    printf(
                        "Connecting (%i, %i, nb:%i, gd:%i) [%i] with dipole to "
                        "(%i, %i, nb:%i, gd:%i) [%i]\n",
                        r, c, get_neighbor(r, c), is_defect(r, c), my_color,
                        nb_i, nb_j, get_neighbor(nb_i, nb_j),
                        is_defect(nb_i, nb_j), !my_color);
#endif

                    break;
                }
#ifdef __DEBUG
                else {
                    printf("Found defect at (%i,%i) whose neighbor is %i\n",
                           nb_i, nb_j,
                           static_cast<int>(get_neighbor(nb_i, nb_j)));
                }
#endif
            }
            ctr++;
        }
    }

    // FIXME: Hack; just make any other unassigned sites equal to a defect
    nexcess_defects = 0;
    for (int i = 0; i < nrows * ncols; i++) {
        if (_lattice[i].neighbor == Neighbor::Unoccupied) {
            _lattice[i].neighbor = Neighbor::None;
            _lattice[i].color = Black;
            _defects[i] = 1;
            nexcess_defects++;
            printf("Adding excess defect at (%i,%i)\n", i / _nc, i % _nc);
        }
    }
    _nd = ndefects_generated + nexcess_defects;

    printf("Finished generating lattice with %i defects (%i excess)\n", _nd,
           nexcess_defects);
}

// FIXME: Compare performance of sites permanently static in memory to this
// implementation (which keeps reads local, w/in a cache line)
int RhombicLattice::swap(int defect_i, int defect_j, Neighbor nb, int& dest_i,
                         int& dest_j) {
    /* Swapping is a bit complicated, sadly.

       Define:
       a. N(i,j) = { RhombicLatticeSite(i,j) for (i,j) in Interaction(i,j) } and
       b. neighbor(i,j) = { (i',j') :
                            if RhombicLatticeSite(i,j).neighbor = (i',j')
                            and RhombicLatticeSite(i',j').neighbor=(i,j) }
       c. (neighbor_i, neighbor_j) = neighbor(proposal_i, proposal_j)

       We *assume* that (proposal_i, proposal_j) in N(defect_i, defect_j).

       Given this assumption, here is the swapping logic
       1. If neighbor(proposal_i, proposal_j) is in N(defect_i, defect_j),
          Then simply swap (proposal_i, proposal_j) with neighbor(proposal_i,
       proposal_j)
       2. Otherwise, we have to perform a cyclic rotation in the direction of
       the dimer. That is, apply the permutation (on positions, not
       connectivity, color),

          original                 ||    new
          =====================================================
          (defect_i, defect_j)     || (neighbor_i, neighbor_j)
          (proposal_i, proposal_j) || (defect_i, defect_j)
          (neighbor_i, neighbor_j) || (proposal_i, proposal_j)

          Why this choice? It preserves the orientation of dimer (relative to
       the defect) and is relatively intuitive. Here's a picture:

          original:                   new
               w'                        d'
           w b'                      w w'
          b d' d                    b b' d
           b-w                       b-w

       One cool thing about this type of swap is that it distinctly separates
       rotational and translational moves. Move 1 is strictly rotational,
       whereas move 2 is strictly translational; as such, we probably can look
       at move 1, move 2 correlators to determine coupling

     */

    // -1. If proposal is also a dimer, give up
    int proposal_i, proposal_j;
    find_neighbor_index(defect_i, defect_j, proposal_i, proposal_j, nb);
    if ((proposal_i < 0) && (proposal_j < 0)) {
        return -1;
    }

    // 0. Get Neighbor of proposal
    int pn_i, pn_j;
    find_neighbor_index(proposal_i, proposal_j, pn_i, pn_j,
                        get_neighbor(proposal_i, proposal_j));
    dest_i = -1;
    dest_j = -1;

    // 1. Handle defect swapping
    if (pn_i < 0 || pn_j < 0) {
#ifdef __DEBUG
        printf(
            "[Rejection] Attempting to swap a defect at (%i,%i,nb:%i) with "
            "defect at (%i,%i,nb:%i)\n",
            defect_i, defect_j,
            static_cast<int>(get_neighbor(defect_i, defect_j)), proposal_i,
            proposal_j, static_cast<int>(get_neighbor(proposal_i, proposal_j)));
#endif
    } else if ((pn_i >= 0) && (pn_j >= 0)) {
        // 2a. If both the proposal and it's neighbor are in the defect's
        // neighborhood,
        //     Then swap the proposal and the defect and patch up the
        //     directions/neighbors
        if (is_in_neighborhood(defect_i, defect_j, pn_i, pn_j)) {
            // We have to divide this up into two cases:
            // a. The direction that proposal's dimer points in corresponds to a
            // clockwise
            //    rotation around the defect (we will let the direction
            //    correspond to a vector from black to white)
            // b. Counter-clockwise rotation
            //
            // I'll draw the counter-clockwise directions
            //
            // Type  |Before| After
            // ------|------|-------
            // W->NW | w-b  | w x
            //       |  x   |  b
            // ---------------------
            // NW->NE|  w   |  w
            //       | x b  | b x
            // ---------------------
            // NE->E | x w  | b-w
            //       |  b   |  x
            // ---------------------
            // E->SE |  x   |  b
            //       | b-w  | x w
            // ---------------------
            // SE->SW| b x  | x b
            //       |  w   |  w
            // ---------------------
            // SW->W |  b   |  x
            //       | w x  | w-b
            //
            //
            // Note that this corresponds to simply rotating by one unit in the
            // *clockwise* direction
            //   NW NE            W  NW
            // W       E    ->  SW     NE
            //   SW SE            SE E
            //
            // This is kind of weird.
            //
            // As such, we need to do two things:
            // 1. Identify an orientation of the dimer
            //    a. Compute a difference direction, dp = (pn_i-proposal_i,
            //    pn_j-proposal_j) b. Classify based on direction. Note that
            //    there are only 5 directions in each loop
            //       If n is even, we have the following memory layout around a
            //       defect at (1,1)
            //         0  1  2
            //       0 NW NE
            //       1 W  X  E
            //       2 SW SE
            //
            //       If clockwise, dp in S_clock = [ (0,1), (1,0), (1,-1),
            //       (-1,-1), (-1,0) ] If counterclockwise in S_cclock = [
            //       (0,-1), (1,0), (1,1), (-1,1), (-1,0) ]
            //
            //       S_clock.intersection(S_cclock) = [ (1,0), (-1,0) ]
            //
            //       If n is odd, we have:
            //         0  1  2
            //       0    NW NE
            //       1 W  X  E
            //       2    SW SE
            //
            //       S_clock  = [ (1,1), (1,0), (0,-1), (-1,0), (-1, 1) ]
            //       S_cclock = [ (1,-1), (1,0), (0,1), (-1,0), (-1,-1) ]
            //       S_clock.intersection(S_cclock) = [ (1,0), (-1,0) ]
            //
            //
            // 2. Update Neighbor(pn_i, pn_j) as well
            //
            //
            // BIG TODO: Figure out how the hell to get rid of branches!!
            //

            auto shitty_mod = [](int i, int modulus) {
                if (i < -1) {
                    return i + modulus;
                } else if (i > 1) {
                    return i - modulus;
                } else {
                    return i;
                }
            };
            // Generate wrapped different vectors
            int dproposal_i = shitty_mod(proposal_i - defect_i, _nr);
            int dproposal_j = shitty_mod(proposal_j - defect_j, _nc);

            int dpn_i = shitty_mod(pn_i - defect_i, _nr);
            int dpn_j = shitty_mod(pn_j - defect_j, _nc);

            int cp_z = dproposal_i * dpn_j - dproposal_j * dpn_i;
            bool is_clockwise = cp_z < 0;

#ifdef __DEBUG
            printf(
                "[Rotation] About to swap defect at (%i,%i,nb:%i,gd:%i) with "
                "proposal (%i,%i,nb:%i,gd:%i)-(%i,%i,nb:%i,gd:%i), "
                "dproposal=(%i,%i), dpn=(%i,%i), is_clockwise=%i\n",
                defect_i, defect_j,
                static_cast<int>(get_neighbor(defect_i, defect_j)),
                is_defect(defect_i, defect_j), proposal_i, proposal_j,
                static_cast<int>(get_neighbor(proposal_i, proposal_j)),
                is_defect(proposal_i, proposal_j), pn_i, pn_j,
                static_cast<int>(get_neighbor(pn_i, pn_j)),
                is_defect(pn_i, pn_j), dproposal_i, dproposal_j, dpn_i, dpn_j,
                is_clockwise);
#endif

            // Now update proposal's direction.
            //
            // Note that if the proposed move is clockwise (relative to the
            // defect), then we select the new neighbor in an *anti-clockwise*
            // manner (see the diagrams above)
            get_neighbor(proposal_i, proposal_j) =
                (is_clockwise)
                    ? rotate_counterclockwise(
                          get_neighbor(proposal_i, proposal_j))
                    : rotate_clockwise(get_neighbor(proposal_i, proposal_j));

            // Perform the swap
            swap_site(defect_i, defect_j, proposal_i, proposal_j);
            dest_i = proposal_i;
            dest_j = proposal_j;

            // Update the neighbor's orientation
            get_neighbor(pn_i, pn_j) =
                find_opposite_direction(get_neighbor(defect_i, defect_j));

#ifdef __DEBUG
            printf(
                "[Rotation] After swapping (former) defect at "
                "(%i,%i,nb:%i,gd:%i) with proposal (%i,%i,nb:%i,gd:%i), pn "
                "(%i,%i,nb:%i,gd:%i)\n",
                defect_i, defect_j,
                static_cast<int>(get_neighbor(defect_i, defect_j)),
                is_defect(defect_i, defect_j), proposal_i, proposal_j,
                static_cast<int>(get_neighbor(proposal_i, proposal_j)),
                is_defect(proposal_i, proposal_j), pn_i, pn_j,
                static_cast<int>(get_neighbor(pn_i, pn_j)),
                is_defect(pn_i, pn_j));
#endif
            return 0;
        } else {
        // 2b. Otherwise, simpler perform the permutation as indicated in the
        // earlier comment;
        //     note that this doesn't require us to update orientation
#ifdef __DEBUG
            printf(
                "[Translation] about to swap defect (%i,%i,nb:%i,gd:%i) with "
                "dimer (%i,%i,nb:%i,gd:%i)-(%i,%i,nb:%i,gd:%i)\n",
                defect_i, defect_j,
                static_cast<int>(get_neighbor(defect_i, defect_j)),
                is_defect(defect_i, defect_j), proposal_i, proposal_j,
                static_cast<int>(get_neighbor(proposal_i, proposal_j)),
                is_defect(proposal_i, proposal_j), pn_i, pn_j,
                static_cast<int>(get_neighbor(pn_i, pn_j)),
                is_defect(pn_i, pn_j));
#endif

            swap_site(defect_i, defect_j, pn_i, pn_j);
            dest_i = pn_i;
            dest_j = pn_j;

#ifdef __DEBUG
            printf(
                "[Translation] swapped defect (%i,%i,nb:%i,gd:%i) with "
                "neighbor (%i,%i,nb:%i,gd:%i)\n",
                defect_i, defect_j,
                static_cast<int>(get_neighbor(defect_i, defect_j)),
                is_defect(defect_i, defect_j), pn_i, pn_j,
                static_cast<int>(get_neighbor(pn_i, pn_j)),
                is_defect(pn_i, pn_j));

            swap_site(defect_i, defect_j, proposal_i, proposal_j);
            printf(
                "[Translation] swapped proposal (%i,%i,nb:%i,gd:%i) with "
                "neighbor (%i,%i,nb:%i,gd:%i)\n",
                defect_i, defect_j,
                static_cast<int>(get_neighbor(defect_i, defect_j)),
                is_defect(defect_i, defect_j), proposal_i, proposal_j,
                static_cast<int>(get_neighbor(proposal_i, proposal_j)),
                is_defect(proposal_i, proposal_j));
#endif

            get_neighbor(proposal_i, proposal_j) = find_opposite_direction(nb);
            get_neighbor(defect_i, defect_j) = nb;

#ifdef __DEBUG
            printf(
                "[Translation] defect (%i,%i,nb:%i,gd:%i), proposal "
                "(%i,%i,nb:%i,gd:%i), pn (%i,%i,nb:%i,gd:%i)\n",
                defect_i, defect_j,
                static_cast<int>(get_neighbor(defect_i, defect_j)),
                is_defect(defect_i, defect_j), proposal_i, proposal_j,
                static_cast<int>(get_neighbor(proposal_i, proposal_j)),
                is_defect(proposal_i, proposal_j), pn_i, pn_j,
                static_cast<int>(get_neighbor(pn_i, pn_j)),
                is_defect(pn_i, pn_j));
#endif
            return 1;
        }
    } else {
        LOG(FATAL) << "Bad neighbor of proposal (" << proposal_i << ","
                   << proposal_j << "), proposal_neighbor (" << pn_i << ","
                   << pn_j << ")";
    }
    return -1;
}

// N.B.: We assume that the *even* rows are the rows that are furthest to the
// left!
void RhombicLattice::find_neighbor_index(int i, int j, int& nb_i, int& nb_j,
                                         Neighbor nb) const {
    switch (nb) {
        case Neighbor::NW:
            nb_i = (i > 0) ? i - 1 : _nr - 1;
            nb_j = (i % 2) ? j : (j > 0) ? j - 1 : _nc - 1;
            break;
        case Neighbor::NE:
            nb_i = (i > 0) ? i - 1 : _nr - 1;
            nb_j = (i % 2) ? (j < _nc - 1) ? j + 1 : 0 : j;
            break;
        case Neighbor::E:
            nb_i = i;
            nb_j = (j < _nc - 1) ? j + 1 : 0;
            break;
        case Neighbor::SE:
            nb_i = (i < _nr - 1) ? i + 1 : 0;
            nb_j = (i % 2) ? (j < _nc - 1) ? j + 1 : 0 : j;
            break;
        case Neighbor::SW:
            nb_i = (i < _nr - 1) ? i + 1 : 0;
            nb_j = (i % 2) ? j : (j > 0) ? j - 1 : _nc - 1;
            break;
        case Neighbor::W:
            nb_i = i;
            nb_j = (j > 0) ? j - 1 : _nc - 1;
            break;
        default:
            nb_i = nb_j = -1;
            break;
    }
}

int RhombicLattice::count_defects_in_neighborhood(int i, int j) const {
    int nb_i, nb_j;
    int ret = 0;
    for (int nb = 0; nb < 6; nb++) {
        find_neighbor_index(i, j, nb_i, nb_j, static_cast<Neighbor>(nb));
        if (is_defect(nb_i, nb_j)) ret++;
    }
    return ret;
}

bool RhombicLattice::all_neighbors_defects(int i, int j) const {
    return (count_defects_in_neighborhood(i, j) == 6);
}

void RhombicLattice::swap_site(int src_i, int src_j, int dst_i, int dst_j) {
// FIXME: Maybe add a copy ctor to LatticeSite?
//        I'm just lazy to do the rule-of-three.
//        This is also not theoretically thread-safe
//        --- although we can do an odd/even partition
#ifdef __DEBUG
    if (is_defect(src_i, src_j)) {
        printf(
            "[src] Swapping defect at (%i,%i,nb:%i) with site (%i,%i,nb:%i)\n",
            dst_i, dst_j, get_neighbor(dst_i, dst_j), src_i, src_j,
            get_neighbor(src_i, src_j));
    }
    if (is_defect(dst_i, dst_j)) {
        printf(
            "[dst] Swapping defect at (%i,%i,nb:%i) with site (%i,%i,nb:%i)\n",
            dst_i, dst_j, get_neighbor(dst_i, dst_j), src_i, src_j,
            get_neighbor(src_i, src_j));
    }
#endif

    RhombicLatticeSite tmp;
    tmp.color = get_site(src_i, src_j).color;
    tmp.neighbor = get_neighbor(src_i, src_j);

    get_site(src_i, src_j).color = get_site(dst_i, dst_j).color;
    get_neighbor(src_i, src_j) = get_neighbor(dst_i, dst_j);

    get_site(dst_i, dst_j).color = tmp.color;
    get_neighbor(dst_i, dst_j) = tmp.neighbor;

    // Update positions for writing to disk
    std::swap(_positions[src_i * _nc + src_j], _positions[dst_i * _nc + dst_j]);
}

bool RhombicLattice::is_in_neighborhood(int i, int j, int pi, int pj) const {
    Neighbor nb;
    find_neighbor(i, j, pi, pj, nb);
    return nb != Neighbor::None;
}

void RhombicLattice::check_neighbors() {
    for (int i = 0; i < _nr; i++) {
        for (int j = 0; j < _nc; j++) {
            if (is_defect(i, j)) {
                test_swap(i, j);
            } else {
                int n_i, n_j;
                int nn_i, nn_j;
                find_neighbor_index(i, j, n_i, n_j, get_neighbor_const(i, j));
                find_neighbor_index(n_i, n_j, nn_i, nn_j,
                                    get_neighbor_const(n_i, n_j));
                if ((nn_i != i) || (nn_j != j)) {
                    LOG(FATAL)
                        << "Site (" << i << "," << j << ") has neighbor ("
                        << n_i << "," << n_j << ") which has neighbor (" << nn_i
                        << "," << nn_j << ")";
                }
            }
        }
    }
}

void RhombicLattice::check_defects() const {
    int ndefects_in_arr = 0;
    for (int i = 0; i < _nr * _nc; i++) {
        ndefects_in_arr += is_defect(i / _nc, i % _nc);
    }
    if (ndefects_in_arr != _nd) {
        LOG(FATAL) << "Found only " << ndefects_in_arr
                   << " but was supposed to have " << _nd << " defects";
    }
}

void RhombicLattice::test_swap(int i, int j) {
    for (int nb = 0; nb < 6; nb++) {
        int nb_i = 0, nb_j = 0;
        int nnb_i, nnb_j;
        int tmp_i, tmp_j;
        int tmp2_i, tmp2_j;
        find_neighbor_index(i, j, nb_i, nb_j, static_cast<Neighbor>(nb));
        Neighbor nnb = get_neighbor_const(nb_i, nb_j);
        if (get_neighbor(nb_i, nb_j) == Neighbor::None) {
            printf("[test_swap] %ith neighbor of (%i,%i) is also a defect\n",
                   nb, i, j);
            continue;
        }
        find_neighbor_index(nb_i, nb_j, nnb_i, nnb_j, nnb);
        printf(
            "[test_swap] Before swaps; initial (%i,%i,n:%i,rn:%i), "
            "nb:(%i,%i,n:%i), nnb(%i,%i,n:%i)\n",
            i, j, static_cast<int>(get_neighbor(i, j)), static_cast<int>(nb),
            nb_i, nb_j, static_cast<int>(nnb), nnb_i, nnb_j,
            static_cast<int>(get_neighbor(nnb_i, nnb_j)));

        int isTrans = swap(i, j, static_cast<Neighbor>(nb), tmp_i, tmp_j);
        Neighbor new_dir;
        if (isTrans) {
            find_neighbor(tmp_i, tmp_j, nb_i, nb_j,
                          new_dir);  // nb = Neighbor(nb-tmp)
        } else {
            new_dir = find_opposite_direction(static_cast<Neighbor>(nb));
        }
        printf(
            "[test_swap]: After first swap; tmp: (%i,%i,n:%i), nbor "
            "(%i,%i,n:%i), relative_nbor: %i, isTrans: %i, inn: %i, %i\n",
            tmp_i, tmp_j, static_cast<int>(get_neighbor(tmp_i, tmp_j)), nb_i,
            nb_j, static_cast<int>(get_neighbor(nb_i, nb_j)),
            static_cast<int>(new_dir), isTrans,
            is_in_neighborhood(tmp_i, tmp_j, nb_i, nb_j),
            is_in_neighborhood(tmp_i, tmp_j, nnb_i, nnb_j));

        assert(new_dir != Neighbor::None);

        swap(tmp_i, tmp_j, new_dir, tmp2_i, tmp2_j);

        int tmp3_i, tmp3_j;
        find_neighbor_index(tmp_i, tmp_j, tmp3_i, tmp3_j,
                            find_opposite_direction(static_cast<Neighbor>(nb)));

        printf(
            "[test_swap] After swaps; opp_nb: (%i,%i,n:%i), tmp: (%i,%i,n:%i), "
            "tmp2: (%i,%i,n:%i)\n",
            tmp3_i, tmp3_j,
            static_cast<int>(
                find_opposite_direction(static_cast<Neighbor>(nb))),
            tmp_i, tmp_j, static_cast<int>(get_neighbor(tmp_i, tmp_j)), tmp2_i,
            tmp2_j, static_cast<int>(get_neighbor(tmp2_i, tmp2_j)));

        if ((i != tmp2_i) || (j != tmp2_j)) {
            LOG(FATAL) << "Initial point (" << i << "," << j
                       << ") is different from double swapped point (" << tmp2_i
                       << "," << tmp2_j << ")";
        }

        if (get_neighbor(nb_i, nb_j) != nnb) {
            LOG(FATAL) << "Failing swap test for (" << static_cast<int>(nb_i)
                       << "," << static_cast<int>(nb_j)
                       << "); previous nb: " << static_cast<int>(nnb)
                       << "; current nb: "
                       << static_cast<int>(get_neighbor(nb_i, nb_j));
        }
        if (!is_defect(i, j)) {
            LOG(FATAL) << "Failing swap test as defect at (" << i << "," << j
                       << ") didn't swap correctly";
        }

        int nnb2_i, nnb2_j;
        find_neighbor_index(nb_i, nb_j, nnb2_i, nnb2_j, nnb);
        if ((nnb2_i != nnb_i) || (nnb2_j != nnb_j)) {
            LOG(FATAL) << "Failed swap tests as next original neighbor ("
                       << static_cast<int>(nnb_i) << ","
                       << static_cast<int>(nnb_j)
                       << ") did not match post-swap neighbor ("
                       << static_cast<int>(nnb2_i) << ","
                       << static_cast<int>(nnb2_j) << ")";
        }
    }
}

bool RhombicLattice::perform_single_metropolis_move(int a, int b) {
    // Only perform metropolis move on defects
    if (!is_defect(a, b)) return false;

    // Find current number of defected nbors
    int old_defects = count_defects_in_neighborhood(a, b);

    // Generate a random neighbor
    Neighbor nb = static_cast<Neighbor>(_dir_dist.generate());
    int nb_i, nb_j;
    find_neighbor_index(a, b, nb_i, nb_j, nb);

    // Swap with neighbor, count defects
    int a_new, b_new;
    int isTrans = swap(a, b, nb, a_new, b_new);
    if ((a_new < 0) && (b_new < 0)) {
        return false;  // Attempting to swap w/ defect
    }
    int new_defects = count_defects_in_neighborhood(a_new, b_new);

#ifdef __DEBUG
    if (old_defects == new_defects) {
        printf(
            "[perform_single_metropolis_move] defect (%i, %i), neighbor (%i, "
            "%i), neighbor's neighbor (%i, %i), old_defects (%i), new_defects "
            "(%i)\n",
            a, b, nb_i, nb_j, a_new, b_new, old_defects, new_defects);
    }
#endif

    // Compute bias probability
    double p = exp(_beta * (old_defects - _alpha * new_defects));
    if (_dbl_dist.generate() >= p) {
#ifdef __DEBUG
        printf("[perform_single_metropolis_move] REJECTED move at (%i, %i)\n",
               a, b);
#endif

        // See test_swap() for why this works
        Neighbor new_dir;
        int nb_i, nb_j;
        if (isTrans) {
            find_neighbor_index(a, b, nb_i, nb_j, nb);
            find_neighbor(a_new, b_new, nb_i, nb_j, new_dir);
        } else {
            new_dir = find_opposite_direction(nb);
        }

        int a_rev, b_rev;
        swap(a_new, b_new, new_dir, a_rev, b_rev);
        if ((a_rev != a) || (b_rev != b)) {
            LOG(FATAL) << "Attempted to reverse a swap from (" << a << "," << b
                       << ") to (" << a_rev << "," << b_rev
                       << "), but we failed";
        }
        return false;
    } else {
#ifdef __DEBUG
        printf(
            "[perform_single_metropolis_move] ACCEPTED move at (%i, %i) with "
            "Boltzmann weight %g @ (T:%g, alpha:%g)\n",
            a, b, p, 1.0 / _beta, _alpha);
#endif
        return true;
    }
}

// Each "step" is simply the first hitting time until an acceptance
void RhombicLattice::run_to_acceptance(long nsteps) {
    long current_step = 0;
    while (current_step < nsteps) {
#ifdef __DEBUG
        check_defects();
        check_neighbors();
#endif

        // a. Draw row, col
        int a = _row_dist.generate();
        int b = _col_dist.generate();

        // FIXME: We should also store a list of defects that we can iterate
        // over and/or just go to something like an
        //        unordered_map and incur the O(log(ndefects)) penalty.
        bool success = perform_single_metropolis_move(a, b);
        if (success) current_step++;
    }
}

// Each "step" corresponds to performing an accept+reject @ each site
void RhombicLattice::run_sequential(long nsteps) {
    for (int r = 0; r < _nc; r++) {
        for (int c = 0; c < _nc; c++) {
            perform_single_metropolis_move(r, c);
        }
    }
}

void RhombicLattice::find_neighbor(int i, int j, int ni, int nj,
                                   Neighbor& nb) const {
    int dx = ni - i;
    int dy = nj - j;

    // Ugh, so ugly.
    if (dx == _nr - 1) {
        dx = -1;
    }
    if (dx == -_nr + 1) {
        dx = 1;
    }
    if (dy == _nc - 1) {
        dy = -1;
    }
    if (dy == -_nc + 1) {
        dy = 1;
    }

#ifdef __DEBUG
    printf("(%i,%i) - (%i,%i), dx (%i), dy (%i)\n", i, j, ni, nj, dx, dy);
#endif

    bool isNW = (i % 2) ? (dx == -1) && (dy == 0) : (dx == -1) && (dy == -1);
    bool isNE = (i % 2) ? (dx == -1) && (dy == 1) : (dx == -1) && (dy == 0);
    bool isE = (dx == 0) && (dy == 1);
    bool isSE = (i % 2) ? (dx == 1) && (dy == 1) : (dx == 1) && (dy == 0);
    bool isSW = (i % 2) ? (dx == 1) && (dy == 0) : (dx == 1) && (dy == -1);
    bool isW = (dx == 0) && (dy == -1);

    if (isNW) {
        nb = Neighbor::NW;
    } else if (isNE) {
        nb = Neighbor::NE;
    } else if (isW) {
        nb = Neighbor::W;
    } else if (isSE) {
        nb = Neighbor::SE;
    } else if (isSW) {
        nb = Neighbor::SW;
    } else if (isE) {
        nb = Neighbor::E;
    } else {
        nb = Neighbor::None;
    }
}

void RhombicLattice::write_lattice_to_text(FILE* lattice_file,
                                           FILE* defect_file) {
    for (int i = 0; i < _nr * _nc; i++) {
        fprintf(lattice_file, "%i %i,", _positions[i],
                static_cast<int>(_lattice[i].neighbor));
#ifdef __DEBUG
        if (is_defect(i / _nc, i % _nc)) {
            printf("Defect @ %i w/ neighbor %i\n", _positions[i],
                   static_cast<int>(_lattice[i].neighbor));
        }
#endif
        fprintf(defect_file, "%i,", _defects[i]);
    }
    fprintf(lattice_file, "\n");
    fprintf(defect_file, "\n");
}

void RhombicLattice::write_lattice_to_hdf5() {
    if (_use_writer) {
        _w->write(_lattice.data(), _positions.data());
    }
}

void RhombicLattice::TriangularWriter::write(RhombicLatticeSite* lattice,
                                             int* positions) {
    // 0. Invert lattice-centered arrays to make them position-centered
    // FIXME: Ugly because positions is stored "backwards"
    int* direction_data = (int*)malloc(_nrows * _ncols * sizeof(int));
    int* position_data = (int*)malloc(_nrows * _ncols * sizeof(int));

    for (int i = 0; i < _nrows; i++) {
        for (int j = 0; j < _ncols; j++) {
            // pos[i*ncols+j] contains the particle id at lattice site (i,j)
            direction_data[positions[i * _ncols + j]] =
                static_cast<int>(lattice[i * _ncols + j].neighbor);
            position_data[positions[i * _ncols + j]] = i * _ncols + j;
        }
    }

    // 1. Update offset, check extent
    offset[0] = _nwritten++;

    DLOG(INFO) << "Before write";
    if (_nwritten > current_dims[0]) {
        current_dims[0] += 32;
        status = H5Dset_extent(pos_dataset, current_dims);
        // printf("Status %i after set pos extent\n", status);
        status = H5Dset_extent(dir_dataset, current_dims);
        // printf("Status %i after set dir extent\n", status);
    }
#ifdef __DEBUG
    else {
        hsize_t extent[3] = {0, 0, 0};
        status = H5Sget_simple_extent_dims(pos_dataspace, extent, max_dims);
        printf(
            "[TriangularWriter] Current pos extent: (%llu,%llu,%llu), "
            "current_dims=(%llu,%llu,%llu), offset=(%llu,%llu,%llu)\n",
            extent[0], extent[1], extent[2], current_dims[0], current_dims[1],
            current_dims[2], offset[0], offset[1], offset[2]);
        status = H5Sget_simple_extent_dims(dir_dataspace, extent, max_dims);
        printf("[TriangularWriter] Current dir extent: (%llu,%llu,%llu)\n",
               extent[0], extent[1], extent[2]);
    }
#endif

    // 2. Write pos data
    DLOG(INFO) << "Before pos write; _ndims: " << _ndims;
    hid_t pos_filespace = H5Dget_space(pos_dataset);
    status = H5Sselect_hyperslab(pos_filespace, H5S_SELECT_SET, offset, NULL,
                                 frame_size, NULL);
    hid_t pos_memspace = H5Screate_simple(3, frame_size, max_dims);
    status = H5Dwrite(pos_dataset, H5T_NATIVE_INT, pos_memspace, pos_filespace,
                      H5P_DEFAULT, position_data);
    H5Sclose(pos_filespace);
    H5Sclose(pos_memspace);

    // 3. Write dir data
    DLOG(INFO) << "Before dir write";
    hid_t dir_filespace = H5Dget_space(dir_dataset);
    status = H5Sselect_hyperslab(dir_filespace, H5S_SELECT_SET, offset, NULL,
                                 frame_size, NULL);
    hid_t dir_memspace = H5Screate_simple(3, frame_size, max_dims);
    status = H5Dwrite(dir_dataset, H5T_NATIVE_INT, dir_memspace, dir_filespace,
                      H5P_DEFAULT, direction_data);
    H5Sclose(dir_filespace);
    H5Sclose(dir_memspace);

    // 4. Cleanup
    free(direction_data);
    free(position_data);
}

}  // namespace Dimers
