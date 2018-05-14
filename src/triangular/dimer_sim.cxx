#include <getopt.h>
#include <cstdio>
#include <ctime>
#include <iostream>

#include "dimer.hxx"

/* Dimer system

   === Memory Layout ===

   Our dimers will be made up of staggered linear arrays; by convention, the
   even-indexed rows will be staggered ahead of the odd-indexed rows, e.g.

   Row
   0   x x x x x x ... x
   1  x x x x x x ... x
   2   x x x x x x ... x
   3  x x x x x x ... x

   where x marks a spot that can be occupied by either a black or white dimer
   cover. Note that the ith element in row j will interact with elements in the
   following set:

   Interaction(i, j) = (i % 2)
                     ? { (i-1, j-1), (i-1, j), (i, j-1), (i, j+1), (i+1, j-1),
   (i+1, j) } : { (i-1, j), (i-1, j+1), (i, j-1), (i, j+1), (i+1, j), (i+1, j+1)
   }

   Each element is connected to either no other element or one of the six
   elements in Interactions(i, j). As such, each lattice site is represented by
   a color (stored as a bool) and an edge type, which is an integer in the range
   0 to 6 (inclusive; stored as an enum). We will number the sites as one
   numbers a benzene ring:

         1 2          NW  NE
        6 * 3  or   W        E
         5 4          SW  SE

   TODO: Optimize the lattice site data structure

   === Monte Carlo Dynamics ===

   Acceptable moves involve swapping any dimer tile (i,j) with any element in
   Interaction(i,j). We will have:

   Delta_E(i,j,a,b) = 1[(a,b) in Interactions(i,j)] 
                    * ( 1[ Configuration that swaps (i,j) with (a,b) separates defects ] 
                    - alpha * 1[ Configuration that swaps (i,j) with (a,b) merges defect] )

   Let's reason about at what happens for different values of alpha:

   Case 0: alpha == 1

   We will effectively have infinite-time acceptor states that involve clusters
   of defects. This is because Delta_E will be <= 0 when we merge zero or more
   defects. As such, if we find a local hexagonal patch with every element equal
   to a defect, then we will never leave that state. Presumably, the system will
   evolve until it finds a state in which it has collected lots of defects into
   one place.

   Case 1: alpha < 1:

   Delta_E <= 0 iff alpha * num_merged_defects > num_separated_defects, so we
   need to merge 1/alpha * num_separated defects in order to move defects closer
   together. This is more akin to the behavior that we hope to see, since


   The second pair of indicator function condition needs to be drawn out; let b
   be a black vertex and w a white vertex

           b-w w-b
          w b w-b
           b b-w

   Note that diagonal lines are implied if none of the other admissible (b-w,
   w-b) lines are drawn. If b swaps with NW (b) or NE (w), we only need to swap
   b with one of the members and update their connectivity. If, however, b swaps
   with E or SE, we need to swap the defect b with the dimer b-w. This is how
   such a swap would look:

            b-w w-b
           w b w-b
            b w b
   (Again, diagonal bonds are implied). We will have to handle all of the
   various ways that this can happen.

   Once we perform the swap, we check whether the number of defects that b was
   touching perviously is less than or greater than the number after. This will
   give us Delta_E.

   === Constraints ===

   Our system will specify some constraints, which will take the form of a fixed
   number of defects. Clement suggests that defects all need to be of the same
   species due to Dubedat's bosonizations results

   -- Tarun (05/04/2016)
*/

int main(int argc, char **argv) {
    // Options to collect
    int nrows, ncols, ndefects;
    long max_steps, steps_per_write = 1;
    float alpha = 1.0, beta = 1.0;

    // IO
    FILE *lattice_outfile = NULL;
    FILE *defects_outfile = NULL;
    char *lattice_file;
    char *defects_file;
    std::string hdf5_file;

    int c;

    while (true) {
        static struct option long_options[] = {
            {"nrows", required_argument, NULL, 'r'},
            {"ncols", required_argument, NULL, 'c'},
            {"ndefects", required_argument, NULL, 'd'},
            {"max_steps", required_argument, NULL, 'm'},
            {"out_freq", optional_argument, NULL, 'w'},
            {"alpha", required_argument, NULL, 'a'},
            {"beta", required_argument, NULL, 'b'},
            {"outfile", required_argument, NULL, 'o'},
            {"textfile", optional_argument, NULL, 't'},
            {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "r:c:d:m:w:a:b:o:t:", long_options,
                        &option_index);

        if (c == -1) break;

        printf("c=%c\n", (char)c);

        switch (c) {
            case 'r':
                nrows = atoi(optarg);
                break;
            case 'c':
                ncols = atoi(optarg);
                break;
            case 'd':
                ndefects = atoi(optarg);
                break;
            case 'm':
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
                printf("in outfile\n");
                hdf5_file = optarg;
                hdf5_file += ".h5";
                break;
            case 't':
                printf("in textfile\n");
                lattice_file =
                    (char *)malloc(sizeof(char) * (strlen(optarg) + 5));
                defects_file =
                    (char *)malloc(sizeof(char) * (strlen(optarg) + 5));
                sprintf(lattice_file, "%s.lat", optarg);
                sprintf(defects_file, "%s.def", optarg);

                lattice_outfile = fopen(lattice_file, "w");
                defects_outfile = fopen(defects_file, "w");

                free(lattice_file);
                free(defects_file);

                break;
            case '?':
                break;
            default:
                printf("Invalid option found: %s\n", optarg);
                exit(1);
                break;
        }
    }

    printf(
        "About to generate lattice of size (%i, %i) with %i defects at alpha "
        "(%f), beta (%f); interval (%li), max (%li), hdf5_file (%s)\n",
        nrows, ncols, ndefects, alpha, beta, steps_per_write, max_steps,
        hdf5_file.c_str());
    std::default_random_engine generator;
    Dimers::RhombicLattice L(nrows, ncols, ndefects, alpha, beta, &generator,
                             hdf5_file, true);

    long total_steps_run = 0;
    long steps_to_run = 0;
    long nintervals_run = 0;

    printf("Before main loop\n");
    while (total_steps_run < max_steps) {
        clock_t t1, t2;
        steps_to_run = std::min(steps_per_write, max_steps - steps_to_run);

        // Run a interval and time
        t1 = clock();
        L.run_sequential(steps_to_run);
        t2 = clock();

        double dt = (double)(t2 - t1) / CLOCKS_PER_SEC * 1000;  // ms
        total_steps_run += steps_per_write;

        if (nintervals_run % 100 == 0) {
            printf(
                "At step %li out of %li steps, took %g for the last %li "
                "steps\n",
                total_steps_run, max_steps, dt, steps_per_write);
        }

        L.write_lattice_to_hdf5();

        if ((lattice_outfile != NULL) && (defects_outfile != NULL)) {
            L.write_lattice_to_text(lattice_outfile, defects_outfile);
        }

        nintervals_run++;
    }

    return 0;
}
