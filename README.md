# motif_prob
Developed and maintained by: Roberto di Castro (robdic@libero.it); Luciano Prosperi (prosperi47@gmail.com); Mattia Prosperi (m.prosperi@ufl.edu)
The software is released under the MIT license.

# summary
Identification of DNA/RNA motifs and quantification of their occurrences are important for the study of genetic diseases, gene evolution, transcription sites, and other biological mechanisms. Although several algorithms for motif characterization are available, most of them are quasi-exact, and correct p-value calculation remains challenging. Exact formulae for motif occurrence, under Bernoullian or Markovian models, have exponential complexity, thus can be cumbersome to be implemented efficiently, but approximations can be calculated with constant cost. Prosperi et al. (2012) provided an exact formula for counting the distribution of strings that do not overlap with themselves (i.e. non-clumpable), coupled with a mathematical demonstration of its validity, under both Bernoullian and Markovian models.

# implementation notes
This software implements the exact formula by Prosperi et al. (2012). Two different implementations have been produced: one in Perl ("strperl1e") and another in C++ ("formlp03"). Both programs take the same input and parameters to the formula, namely: (1) a query string or multiple strings to be analyzed; (2) the length of the reference string (i.e. the genome); and (3) the nucleotide frequencies of the genome.

The formula is obviously consistent only for motifs whose length is smaller than the genome length. If a clumped motif is present in the input, the string is flagged.
Since the computational complexity of the formula is exponential, motif occurrences are calculated at increasing counts until the occurrence probability becomes lower than 10e-7, or the upper limit of 500 counts is reached. These values can be changed in the source code, but can be considered appropriate for most of realistic ratios between motif and genome lengths.

In order to avoid issues with floating point operations when frequency/length ratios diverge, and to provide comprehensive estimations for relatively ill-posed configurations, we have further implemented the calculation of the expected number of strings and the motif's (stationary) occurrence probability at any text position, according to Robin et al. (2005) and Marschall & Rahmann (2008).

# source code and binaries
Each program consists of a unique source code file ("strperl1e" written in Perl 5.3 and "formlp03.cpp" written in ISO C++ v14.0, respectively) and no other dependencies/libraries are required. The programs were compiled under Strawberry Perl and Microsoft Visual C++ suites.
Executable files for MS Windows are available on this repository ("strperl1e.exe", "formlp3.exe").

# input specifications
A single input file (named INPPERL.TXT and INPCPLUSPLUS.TXT, which can be changed in the source code) where the first four lines include the relative nucleotide frequencies in the reference genome sequence, the fifth line is the length of the genome, and then each of all the following lines contains a motif (only ACGT characters allowed).

Example:

0.15

0.25

0.25

0.35

4500000

CAGATA

GATTACA

GAGGCGGCGTGC

AGCTGTCGA

CGAGGCTGGCG

GCGCGC

GCACTGC

CGGTCAAA

# references
Marschall T, Rahmann S. Efficient exact motif discovery, Bioinformatics 2009;25:i356–i364.

Prosperi MCF, Prosperi L, Gray RR, Salemi M. On counting the frequency distribution of string motifs in molecular sequences. International Journal of Biomathematics 2012;5(6):1250055.

Robin S, Rodolphe F, Schbath S. DNA, Words and Model: Statistics of Exceptional Words. Cambridge University Press 2005; Cambridge, UK.
