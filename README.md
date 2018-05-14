# Dimorphic: Simulations of Monomer-Dimer Systems


# Background 

This is work that I did in 2016 while trying to (numerically) show that monomer-dimers systems
exhibit similar statistical properties to fragile glass formers. These simulations move the monomers
throughout the graph and use an energy function that promotes clustering among monomers. The
hypothesis is that as one gets closer to the glass transition, one will find that the dynamics of
the system arrest in such a way that dimers (which represent physical molecules) will have to relax
by coordinating in a way that clusters monomers. These simulations (and associated media, such as
videos) aim to show this. 


This repo contains all of the code used to generate simulations of monomer-dimer systems. At the
moment, it only contains data from 2D simulations of monomer-dimer dynamics on the triangular
lattice; however, we hope to add 3D simulations in the near future. There is some code for handling
BCC lattices, although simulations have not been done. 

The results of this simulator have been posted on arXiv and I will update the identifier once I have
received it (ETA: 05/15/2018).

I will also provide some scripts to generate the figures in the paper.

# Installation Requirements

0. C++11 (or greater)
1. SCons
2. Google Log 
3. Google Flags
4. HDF5

You will need to modify the `SConscript` file to point to the paths of your
libraries.

# Triangular lattice

The binary constructed from `dimer_sim.cxx` is used to perform the simulation described in the
paper. The simulations output positions (of monomers and dimers) and orientations (of dimers) to an
HDF5 file. At the moment, blosc compression is turned off by default because it results in a big i/o
slowdown; however, we probably will have to use it in the future.   

# Graph 

In order to handle 3D graphs and to test the monomer-dimer hypothesis on other graphs, I abstracted
some of the lessons that I learned from the first implementation and made a new simulation system.
It has the following components:

1. `DimerGraph`: Abstract Graph that does the following a. Stores the lattice b. Keeps track of the
   monomer/dimer identities c. Generates a random initial condition d. Performs the swap

2. `ILattice`: A very simple interface/abstract class that handles what a lattice needs in order to
   run a simulation

3. `BccLattice`: A wrapper around a DimerGraph object that handles all of the Bcc specific tasks that
   a simulation needs

In addition, I separated out the prng, writer classes and made them a bit more generic so that the
graph-based simulations could share some code w/ the original simulation. 


# TODO

0. Figure out how to store direction/orientation for elements in BccLattice (annoying because the
   graph is not d-regular anymore) 
   EDIT [8/26/2016]: Made a first attempt at this; need to fix some bugs in the tests
1. Add a writer for the BccLattice
2. Make an FccLattice
