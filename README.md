# long_range_LD
Sarch for Long Range LD using the TAPER algorithm.

This is an experimental program to find long range LD using the TAPER algorithm ([Xiong <i>et al</i>., 2006](http://dx.doi.org/10.1109/TKDE.2006.1599388)).
In brief, this program uses an upper bound on the correlation coefficient to ignore pairwise SNP comparisons that could not have a correlation coefficient above a given threshold. 
In practice, the number of pairwise comparisons still remains too large for this program to be used in earnest, unless the threshold is set quite high. In addition, the use of the pairwise correlation coefficient is perhaps inappropriate for looking for long range LD across chromosomes. 
