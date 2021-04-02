# motif_prob
Developed and maintained by: Roberto di Castro (robdic@libero.it); Luciano Prosperi (prosperi47@gmail.com); Mattia Prosperi (m.prosperi@ufl.edu); Simone Marini (smarini@ufl.edu). The software is released under the MIT license.

# Summary
Identification of DNA/RNA motifs and quantification of their occurrences are important for the study of genetic diseases, gene evolution, transcription sites, and other biological mechanisms. Although several algorithms for motif characterization are available, most of them are quasi-exact, and correct p-value calculation remains challenging. Exact formulae for motif occurrence, under Bernoullian or Markovian models, have exponential complexity, thus can be cumbersome to be implemented efficiently, but approximations can be calculated with constant cost. Prosperi et al. (2012) provided an exact formula for counting the distribution of strings that do not overlap with themselves (i.e. non-clumpable), coupled with a mathematical demonstration of its validity, under both Bernoullian and Markovian models.

# Implementation notes
This software implements the exact formula by Prosperi et al. (2012) under the Markovian model. Two different implementations have been produced: one in Perl ("strperl1e") and another in C++ ("formlp03"). Both programs take the same input and parameters to the formula, namely: (1) a query string or multiple strings to be analyzed; (2) the length of the reference string (i.e. the genome); and (3) the nucleotide frequencies of the genome.

The formula is obviously consistent only for motifs whose length is smaller than the genome length. If a clumped motif is present in the input, the string is flagged.
Since the computational complexity of the formula is exponential, motif occurrences are calculated at increasing counts until the occurrence probability becomes lower than 10e-7, or the upper limit of 500 counts is reached. These values can be changed in the source code, but can be considered appropriate for most of realistic ratios between motif and genome lengths.

In order to avoid issues with floating point operations when frequency/length ratios diverge, and to provide comprehensive estimations for relatively ill-posed configurations, we have further implemented the calculation of the expected number of strings and the motif's (stationary) occurrence probability at any text position, according to Robin et al. (2005) and Marschall & Rahmann (2008).

# Installation Notes

## Executables
Executable files for MS Windows are available on this repository ("strperl1e.exe", "formlp3.exe"). The programs were compiled under Padre Perl and MinGW (32bit) C++ suites.

## Source code and binaries
Each program consists of a unique source code file ("strperl1e" written in Perl 5.3 and "formlp03.cpp" written in ISO C++ v14.0, respectively) and no other dependencies/libraries are required. The Perl script can be invoked directly without compilation.

## Makefile
Binaries for UNIX are also available on the repository, along with the Makefile for the C++ version (using g++ compiler). To compile, just run:
```
make
```

# Input specifications
## C++
```
./formlp03 -i INPUT_FILE -o OUTPUT_FILE {-a freq_a -c freq_c -g freq_g -t freq_t -n genome_length | -f GENOME_FILE.fasta}

-i input file
-o output file
-a A frequency
-c C frequency
-g G frequency
-t T frequency
-n total genome length
-f genome file. If specfied, nucleotide frequencies will be calculated on this file,
and other options will be ignored. Fasta format. Note that the total genome length will
be calculated on the fasta full sequence lenght, while nucleotide frequencies will only
consider standard (A, C, G, T) characters.
```

### Examples:
User adds genome file (to be used to calculated the nucleotide frequencies, and length):
```
formlp03 -i INPCPLUSPLUS_unix.TXT -o OUTCPLUSPLUS.TXT -f mygenome.fasta
```
User inputs nucleotide frequencies and genome length from command line:
```
formlp03 -i INPCPLUSPLUS_unix.TXT -o OUTCPLUSPLUS.TXT -a 0.1 -c 0.3 -g 0.1 -t 0.5 -n 100000
```

Example of C++ input file:
Example of Perl input file:
```
CAGATA
GATTACA
GAGGCGGCGTGC
AGCTGTCGA
CGAGGCTGGCG
GCGCGC
GCACTGC
CGGTCAAA
```
## Perl
A single input file (named INPPERL.TXT, which can be changed in the source code) where the first four lines include the relative nucleotide frequencies in the reference genome sequence, the fifth line is the length of the genome, and then each of all the following lines contains a motif (only ACGT characters allowed).

Example of Perl input file:
```
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
```
# Motif discovery and motif data bases
* [Jasper](http://jaspar.genereg.net/)
* [The MEME Suite](https://meme-suite.org/meme/)
* [PROSITE](https://prosite.expasy.org/)
* [MOTIF Search](https://www.genome.jp/tools/motif/)
* [MolBiol motifs](https://molbiol-tools.ca/Motifs.htm)

# References
Marschall T, Rahmann S. Efficient exact motif discovery, Bioinformatics 2009;25:i356–i364.

Prosperi MCF, Prosperi L, Gray RR, Salemi M. On counting the frequency distribution of string motifs in molecular sequences. International Journal of Biomathematics 2012;5(6):1250055.

Robin S, Rodolphe F, Schbath S. DNA, Words and Model: Statistics of Exceptional Words. Cambridge University Press 2005; Cambridge, UK.


# MIT license

Copyright (c) [2021] [Mattia Prosperi, Simone Marini, et al.]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
