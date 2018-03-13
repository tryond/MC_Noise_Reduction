Console Compile Command: 

g++ -std=c++11 -O3 -fopenmp -o simplept RNG.cpp Vec.cpp SimpleFilter.cpp ObjectFilter.cpp BrdfFilter.cpp NormalFilter.cpp DofFilter.cpp simplept.cpp

Console Run Command:

simplept 256 1 2 3 4 5 6 7 8 9 10 11 12 13 14 

1: simple mean parameter

2: simple median parameter

3: simple gaussian parameter

4: object mean parameter

5: object median parameter

6: object gaussian parameter

7: brdf mean parameter

8: brdf gaussian parameter

9: normal mean parameter

10: normal gaussian parameter

11: dof mean level parameter

12: dof mean focus level parameter

13: dof gray level parameter

14: dof gray focus level parameter


Note: edit frame parameter (line 334) in order to produce multiple frames
- if multiple frames, uncomment line 388 and comment out line 389