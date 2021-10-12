# CVTree-Parallelised
This repository houses the parallelisation of the sequential program, "Bioinformatics - Genome Similarity Using Frequency Vectors (C++)". The program takes in a file through the command-line containing the number of strings in the first line and following that several bacteria names. These names are then concatenated with  .faa matching to files in the /data/ directory with gene names and their respective gene represented as DNA. The program collects the "kmer" or subset of an entire gene and then compares that against other bacteria "kmers" and finds the correlation between any given bacteria. The parallelization language and frameworks in consideration for the project are as follows. As the language of the program is C++, the implementation of frameworks such as OpenMP and/or C++17 Parallel Algorithms. In addition to this the implementation of NUMA is also under consideration. Regarding the hardware, it will be targeted at an Intel Core i7-6700 Skylake with four (4) physical cores and eight (8) virtual cores. Moreover, Intel VTune is a possible profiling tool that may also be used due to the use of the program on an Intel CPU.
