# advection


## Context
The purpose of this program is to study the flow of advection by building a multi-threaded simulation given a set of parameters. This phenomenon is modelled using a hyperbolic partial-differential equation. OpenMP is used to achieve multi-core efficiency when applying boundary conditions. The timed tests for this program were run on an 8 core Apple MacBook Pro with the Apple M1 Pro Chip. This contains up to 32GB of unified memory and offers up to 200GB/s memory bandwidth. The hard scaling tests were run iteratively from 1 to 12 cores and the speedup graphs can be found in /out/speedup_.png.

This code is written in C++ and was compiled using g++11.

## Usage
g++-11 -fopenmp advection_boost.cpp<br/>
./a.out N> NT L T u v

The optional flag for number of threads has not been added. To change this, edit numThreads in advection_boost.cpp.

## Accuracy
The serial and parallel implementations were compared using bitwise reproducability. The source code for this program can be found in /src/bitcompare/bitcompare.cpp.
