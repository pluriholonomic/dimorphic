#include "graph.hxx"

namespace Dimers {

DimerGraph::DimerGraph(int V, int MAX_DEGREE, std::default_random_engine *gen,
                       bool verbose)
    : V(V),
      MAX_DEGREE(MAX_DEGREE),
      _nbor_prng(MAX_DEGREE, gen),
      _vertex_prng(V, gen),
      _real_prng(gen),
      _nmonomers(0),
      _verbose(verbose) {
    // Resize nodes, make sure edge array is the correct size
    _nodes.resize(V);
    for (auto &node : _nodes) {
        node.edges.resize(MAX_DEGREE);
    }

    reset_graph();
}

void DimerGraph::reset_graph() {
    for (int i = 0; i < V; i++) {
        _nodes[i].degree = 0;
        _nodes[i].uid = -1;

        _nodes[i].dimer_nbor = -2;
        _nodes[i].color = true;
        for (int d = 0; d < MAX_DEGREE; d++) {
            _nodes[i].edges[d] = -1;
        }
    }
}

void DimerGraph::add_edge(int v, int w) {
    if ((v > V - 1) || (w > V - 1)) {
        LOG(FATAL) << "Attempting to add edge (" << v << "," << w
                   << ") but we only have " << V << " nodes";
    }

    if (has_edge(v, w)) {
        LOG(WARNING) << "Attempting to re-add edge (" << v << "," << w << ")";
        return;
    }

    int dv = _nodes[v].degree;
    int dw = _nodes[w].degree;
    if (dv > MAX_DEGREE) {
        LOG(FATAL) << "Attempting to add an edge to node " << v
                   << " with degree " << dv
                   << " but the maximum degree per vertex is " << MAX_DEGREE;
    }

    if (dw > MAX_DEGREE) {
        LOG(FATAL) << "Attempting to add an edge to node " << w
                   << " with degree " << dw
                   << " but the maximum degree per vertex is " << MAX_DEGREE;
    }
    _nodes[v].edges[dv] = w;
    _nodes[w].edges[dw] = v;

    _nodes[v].degree++;
    _nodes[w].degree++;
    if (v == 63) {
        DLOG(INFO) << "Edge for v=63: " << w << "; deg(v)=" << dv;
        std::stringstream ss;
        ss << "Neighbors of v=63: ";
        for (int i = 0; i < MAX_DEGREE; i++) {
            ss << _nodes[v].edges[i] << " ";
        }
        DLOG(INFO) << ss.str();
    }
    if (w == 63) {
        DLOG(INFO) << "Edge for w=63: " << v << "; deg(w)=" << dw;
        std::stringstream ss;
        ss << "Neighbors of w=63: ";
        for (int i = 0; i < MAX_DEGREE; i++) {
            ss << _nodes[w].edges[i] << " ";
        }
        DLOG(INFO) << ss.str();
    }
}

void DimerGraph::generate_random_dimer_configuration(double monomer_density) {
    int current_uid = 0;
    if (monomer_density > 1.0) {
        LOG(FATAL) << "Attempting to generate monomer configuration with "
                      "density > 1.0 [density: "
                   << monomer_density << "]";
    }

    // Generate monomer samples via reservoir sampling
    int target_nmonomers = std::round(monomer_density * V);

    // FIXME: Was running into a weird bug here, hence the allocation of maximum
    // size..
    int sample[V];

    for (int i = 0; i < target_nmonomers; i++) {
        sample[i] = i;
    }

    for (int i = target_nmonomers; i < V; i++) {
        int j = std::round(i * _real_prng.generate());
        if (j <= target_nmonomers) {
            sample[j] = i;
        }
    }

    for (int j = 0; j < target_nmonomers; j++) {
        int s = sample[j];
        if (_nodes[s].dimer_nbor != -2) {
            LOG(FATAL) << "Attempting to add random dimer configuration to "
                          "non-reset DimerGraph";
        }
        _nodes[s].color = true;  // monomers are black
        _nodes[s].dimer_nbor = -1;
        _nodes[s].uid = current_uid++;
        _nmonomers++;
        DLOG(INFO) << "Added monomer at site " << s
                   << ", current_uid: " << current_uid;
    }

    if (_verbose)
        LOG(INFO) << "Added " << _nmonomers << " monomers before adding dimers";

    // For each vertex:
    // 1. Draw a random neighbor
    // 2. Check if it is a monomer and/or assigned
    // 3. If not, assign it to be a monomer neighbor
    //
    for (int v = 0; v < V; v++) {
        if (_nodes[v].degree <= 0) {
            if (_verbose)
                LOG(INFO) << "Skipping vertex " << v
                          << " because it has degree 0";
            continue;
        }

        if (_nodes[v].dimer_nbor != -2) {
            if (_verbose)
                LOG(INFO) << "Skipping vertex " << v
                          << " because it has already been assigned";
            continue;
        }

        // 0. Find number of non-occupied, non-monomer sites
        int num_nbors_accessible = 0;
        for (int i = 0; i < _nodes[v].degree; i++) {
            if (_nodes[_nodes[v].edges[i]].dimer_nbor == -2) {
                num_nbors_accessible++;
            }
        }

        // 1. If we don't have any accessible neighbors, we become a monomer
        if (num_nbors_accessible == 0) {
            if (_verbose)
                LOG(INFO) << "Added monomer at site " << v
                          << ", current_uid: " << current_uid
                          << " [No accessible neighbors]";
            _nodes[v].color = true;
            _nodes[v].dimer_nbor = -1;
            _nodes[v].uid = current_uid++;
            _nmonomers++;
            continue;
        }

        // 2. Draw random accessible neighbor via reservoir sampling
        int nb_idx = 0, naccessible_seen = 1;
        for (int i = 0; i < _nodes[v].degree; i++) {
            if (_nodes[_nodes[v].edges[i]].dimer_nbor == -2) {
                if (_real_prng.generate() <= 1.0 / naccessible_seen) {
                    nb_idx = i;
                }
                naccessible_seen++;
            }
        }

        if (naccessible_seen - 1 != num_nbors_accessible) {
            LOG(FATAL) << "The number of accessible neighbors for vertex " << v
                       << " is " << num_nbors_accessible << " but we performed "
                       << "reservoir sampling on only " << naccessible_seen - 1;
        }

        // 3. Assign uid, neighbor index to v
        _nodes[v].uid = current_uid;
        _nodes[v].dimer_nbor = nb_idx;
        _nodes[v].color = true;  // first node is black

        // 4. Do the same for the neighbor
        int nb = _nodes[v].edges[nb_idx];
        _nodes[nb].dimer_nbor = get_edge_idx(nb, v);
        _nodes[nb].uid = current_uid++;
        _nodes[nb].color = false;  // nbor is white

        if (_verbose) {
            LOG(INFO) << "Added dimer (" << v << "," << nb
                      << "), current_uid: " << current_uid;
        }
    }

    if (2 * current_uid - _nmonomers != V) {
        LOG(FATAL) << "Unable to generate configuration at density "
                   << monomer_density << "; Generated " << _nmonomers
                   << " monomers, " << current_uid - _nmonomers
                   << " dimers, but V=" << V << " and "
                   << 2 * (current_uid - _nmonomers) + _nmonomers
                   << " = 2*num_dimers+nmonomers != V";
    }

    DLOG(INFO) << "Generated monomer-dimer configuration with " << current_uid
               << " uids";
}

void DimerGraph::perform_swap(int monomer_v, int proposal_v) {
    DLOG(INFO) << "Perform swap requested between " << monomer_v << " and "
               << proposal_v;
    if (monomer_v == proposal_v) {
        LOG(WARNING) << "Attempting to swap element at index " << monomer_v
                     << " with itself";
        return;
    }

    if (!has_edge(monomer_v, proposal_v)) {
        LOG(WARNING)
            << "Attempting to perform a swap between non-adjacent vertices "
            << monomer_v << " (monomer) and " << proposal_v;
        return;
    }

    if (!is_monomer(monomer_v)) {
        LOG(WARNING) << "Attempting to swap a non-monomer at vertex "
                     << monomer_v;
        return;
    }

    // N.B. (maybe FIXME?) Should we do anything special for monomer-monomer
    // swaps? My intuition is that the client class that inherits from this
    // class should really do that type of work, so we should basically let them
    // handle it and if it comes up here, we simply ignore monomer-monomer swaps
    if (is_monomer(proposal_v)) {
        LOG(WARNING)
            << "Attempting to swap a monomer with a monomer; vertices: "
            << monomer_v << "," << proposal_v;
        return;
    }

    int proposal_nb = get_dimer_nbor(proposal_v);
    DLOG(INFO) << "Proposal " << proposal_v << " has dimer neighbor "
               << proposal_nb;

    if (has_edge(monomer_v, proposal_nb)) {
        // If the monomer and the proposal's neighbor are adjacent,
        // perform a simple swap of proposal, monomer

        // 0. swap uids
        int tmp_uid = _nodes[monomer_v].uid;
        _nodes[monomer_v].uid = _nodes[proposal_v].uid;
        _nodes[proposal_v].uid = tmp_uid;

        // 1. swap colors
        _nodes[monomer_v].color =
            _nodes[monomer_v].color ^ _nodes[proposal_v].color;
        _nodes[proposal_v].color =
            _nodes[monomer_v].color ^ _nodes[proposal_v].color;
        _nodes[monomer_v].color =
            _nodes[monomer_v].color ^ _nodes[proposal_v].color;

        // 2. update monomer_v's neighbor to proposal_nb (and vice versa)
        _nodes[monomer_v].dimer_nbor = get_edge_idx(monomer_v, proposal_nb);
        _nodes[proposal_nb].dimer_nbor = get_edge_idx(proposal_nb, monomer_v);

        // 3. set proposal to be a monomer
        _nodes[proposal_v].dimer_nbor = -1;
    } else {
        // N.B. Ansatz:
        // Any swap of a monomer with a dimer vertex whose neighbor is not in
        // the monomer's neighborhood is direction preserving. This means that
        // we *DON'T* just swap the monomer and the neighbor of the neighbor,
        // but instead, we cyclically swap vertices. For instance, if we have:
        //
        // x - b -w
        //
        // then we go to
        //
        // b - w - x
        //
        // It is of course easier to just swap x, w, if orientations didn't
        // matter, but we would like to separate rotation and translation as
        // cleanly as possible.

        // 0. swap uids of monomer, neighbor of proposal
        int tmp_uid = _nodes[monomer_v].uid;
        _nodes[monomer_v].uid = _nodes[proposal_nb].uid;
        _nodes[proposal_nb].uid = tmp_uid;

        // 1. swap colors of proposal, neighbor of proposal
        //    [ Note: This is the only thing that changes when we swap b, w!]
        _nodes[proposal_v].color =
            _nodes[monomer_v].color ^ _nodes[proposal_nb].color;
        _nodes[proposal_nb].color =
            _nodes[proposal_v].color ^ _nodes[proposal_nb].color;
        _nodes[proposal_v].color =
            _nodes[monomer_v].color ^ _nodes[proposal_nb].color;

        // 2. swap colors of monomer, neighbor of proposal
        //    [ Note: This is the only thing that changes when we swap b, w!]
        _nodes[monomer_v].color =
            _nodes[monomer_v].color ^ _nodes[proposal_nb].color;
        _nodes[proposal_nb].color =
            _nodes[monomer_v].color ^ _nodes[proposal_nb].color;
        _nodes[monomer_v].color =
            _nodes[monomer_v].color ^ _nodes[proposal_nb].color;

        // 3. update monomer_v's neighbor to proposal (and vice-versa)
        _nodes[monomer_v].dimer_nbor = get_edge_idx(monomer_v, proposal_v);
        _nodes[proposal_v].dimer_nbor = get_edge_idx(proposal_v, monomer_v);

        // 4. Set proposal_nb to be a monomer
        _nodes[proposal_nb].dimer_nbor = -1;
    }
}

bool DimerGraph::has_edge(int v, int w) {
    bool v_has_w = false, w_has_v = false;
    for (int d = 0; d < MAX_DEGREE; d++) {
        v_has_w |= (_nodes[v].edges[d] == w);
        w_has_v |= (_nodes[w].edges[d] == v);
    }
    return v_has_w && w_has_v;
}

int DimerGraph::get_edge_idx(int v, int w) {
    int ret = -1;
    for (int d = 0; d < MAX_DEGREE; d++) {
        if (_nodes[v].edges[d] == w) {
            ret = d;
            break;
        }
    }
    return ret;
}

void DimerGraph::get_monomers(std::vector<int> &mv) {
    for (int v = 0; v < V; v++) {
        if (_nodes[v].dimer_nbor == -1) mv.push_back(v);
    }
}

void DimerGraph::draw_random_neighbor(int v, int &rnd_nb, int &rnd_nb_idx) {
    // Generate a random neighbor of vertex v and store the neighbor and index
    // into edge array
    auto edges = _nodes[v].edges;
    int num_admissible_neighbors = 0;
    rnd_nb = rnd_nb_idx = -1;
    for (int i = 0; i < MAX_DEGREE; i++) {
        int nb = edges[i];
        if ((nb < 0) || (nb > V - 1)) {
            DLOG(INFO)
                << "Found inadmissible neighbor in draw_random_neighbor: " << nb
                << " which was the " << i << "th neighbor of " << v;
            break;
        }

        if (_nodes[nb].dimer_nbor >= 0) {
            DLOG(INFO) << "Found admissible neighbor " << nb << " of " << v;
            num_admissible_neighbors++;
            if (rnd_nb_idx < 0) {
                rnd_nb = nb;
                rnd_nb_idx = i;
            }
        }
    }

    if (num_admissible_neighbors == 0) {
        std::stringstream ss;
        ss << "Found no admissible neighbors for vertex " << v
           << "; Neighbors of v: ";
        for (int i = 0; i < MAX_DEGREE; i++) {
            ss << edges[i] << " ";
        }
        LOG(WARNING) << ss.str();
        return;
    }

    int nseen = 1;
    for (int i = rnd_nb_idx + 1; i < MAX_DEGREE; i++) {
        int nb = edges[i];
        if ((nb < 0) || (nb > V - 1)) {
            DLOG(INFO) << "About to exit loop; nseen: " << nseen
                       << "; num_admissible_neighbors: "
                       << num_admissible_neighbors;
            break;
        }

        if (_nodes[nb].dimer_nbor >= 0) {
            if (_real_prng.generate() < 1.0 / (++nseen)) {
                rnd_nb = nb;
                rnd_nb_idx = i;
            }
        }
    }
    DLOG(INFO) << "Drawing neighbor " << rnd_nb << " of vertex " << v
               << " with index " << rnd_nb_idx;
}

}  // namespace Dimers
