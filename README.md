# motif_prob
Developed and maintained by: Roberto di Castro (robdic@libero.it); Luciano Prosperi (prosperi47@gmail.com); Mattia Prosperi (m.prosperi@ufl.edu)
The software is released under the MIT license.

# summary
Identification of DNA/RNA motifs and quantification of their occurrences are important for the study of genetic diseases, gene evolution, transcription sites, and other biological mechanisms. Although several algorithms for motif characterization are available, most of them are quasi-exact, and correct p-value calculation remains challenging. Exact formulae for motif occurrence, under Bernoullian or Markovian models, have exponential complexity, thus can be cumbersome to be implemented efficiently, but approximations can be calculated with constant cost. Prosperi et al. (2012) provided an exact formula for counting the distribution of strings that do not overlap with themselves (i.e. non-clumpable), coupled with a mathematical demonstration of its validity, under both Bernoullian and Markovian models.

# implementation notes
This software implements the exact formula by Prosperi et al. (2012). Two different implementations have been produced: one in Perl ("strperl1e") and another in C++ ("formlp03"). Both programs take the same input and parameters, namely: (1) a query string or multiple strings to be analyzed; (2) the length of the reference string (i.e. the genome); and (3) the nucleotide frequencies of the genome.

# source code and binaries
Each program consists of a unique source code file ("strperl1e"  Perl 5.3 and "formlp03.cpp" ISO C++ v14.0, respectively) and no other dependencies/libraries are required. The programs were compiled under Strawberry Perl and Microsoft Visual C++ suites.
Executable files for MS Windows are available on this repository ("strperl1e.exe", "formlp3.exe").

# input specifications
A single input file

