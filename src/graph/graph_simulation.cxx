#include <getopt.h>
#include <ctime>
#include <iostream>
#include "bcc_lattice.hxx"

int main(int argc, char **argv) {
    // Options to collect
    unsigned seed = 0;
    int nx, ny, nz, nmonomers;
    long max_steps, steps_per_write = 1;
    float alpha = 1.0, beta = 1.0;
    char lattice_type;

    // IO
    std::string hdf5_file;

    int c;

    while (true) {
        static struct option long_options[] = {
            {"lattice", required_argument, NULL, 'l'},
            {"nx", required_argument, NULL, 'x'},
            {"ny", required_argument, NULL, 'y'},
            {"nz", required_argument, NULL, 'z'},
            {"seed", optional_argument, NULL, 's'},
            {"nmonomers", required_argument, NULL, 'm'},
            {"max_steps", required_argument, NULL, 't'},
            {"out_freq", optional_argument, NULL, 'w'},
            {"alpha", required_argument, NULL, 'a'},
            {"beta", required_argument, NULL, 'b'},
            {"outfile", required_argument, NULL, 'o'},
            {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "l:x:y:z:s:m:t:w:a:b:o:", long_options,
                        &option_index);

        if (c == -1) break;

        switch (c) {
            case 'l':
                lattice_type = *optarg;
                break;
            case 'x':
                nx = atoi(optarg);
                break;
            case 'y':
                ny = atoi(optarg);
                break;
            case 'z':
                nz = atoi(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'm':
                nmonomers = atoi(optarg);
                break;
            case 't':
                max_steps = atol(optarg);
                break;
            case 'w':
                steps_per_write = atol(optarg);
                break;
            case 'a':
                alpha = atof(optarg);
                break;
            case 'b':
                beta = atof(optarg);
                break;
            case 'o':
                hdf5_file = optarg;
                hdf5_file += ".h5";
                break;
                break;
            case '?':
                break;
            default:
                printf("Invalid option found: %s\n", optarg);
                exit(1);
                break;
        }
    }

    // FIXME: Add in seed
    std::default_random_engine generator;  //(seed);
    Dimers::ILattice *lattice_ptr;

    if (lattice_type != 'b') {
        printf("Only the bcc lattice is currently implemented. Exiting. \n");
        exit(1);
    } else {
        printf(
            "About to generate bcc lattice of size (%i, %i, %i) with %i "
            "monomers, alpha (%f), beta (%f); interval (%li), max (%li), "
            "hdf5_file (%s), seed (%i)\n",
            nx, ny, nz, nmonomers, alpha, beta, steps_per_write, max_steps,
            hdf5_file.c_str(), seed);
        lattice_ptr = new Dimers::BccLattice(nx, ny, nz, &generator, nmonomers,
                                             hdf5_file);
    }

    long total_steps_run = 0;
    long steps_to_run = 0;
    long nintervals_run = 0;

    printf("Before main loop, max_steps: %li\n", max_steps);
    while (total_steps_run < max_steps) {
        clock_t t1, t2;
        steps_to_run = std::min(steps_per_write, max_steps - steps_to_run);

        // Run a interval and time
        t1 = clock();
        lattice_ptr->run_sequential(steps_to_run);
        t2 = clock();

        double dt = (double)(t2 - t1) / CLOCKS_PER_SEC * 1000;  // ms
        total_steps_run += steps_per_write;

        if (nintervals_run % 100 == 0) {
            printf(
                "At step %li out of %li steps, took %g for the last %li "
                "steps\n",
                total_steps_run, max_steps, dt, steps_per_write);
        }

        lattice_ptr->write_lattice_to_hdf5();

        nintervals_run++;
    }

    if (lattice_type == 'b') {
        delete lattice_ptr;
    }

    return 0;
}
