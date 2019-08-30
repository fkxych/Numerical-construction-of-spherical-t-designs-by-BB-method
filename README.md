# Numerical-construction-of-spherical-t-designs-by-Barzilai-Borwein-method
-----Computing of spherical t-designs-----
#
Author: Congpei An and Yuchen Xiao
#
Thanks for Professor Rob Womersley for the codes in [6] and Professor Mark Schmidt for the codes in [18].
#
References: 

[6] I. H. Sloan, R. S. Womersley, A variational characterisation of spherical designs, Journal of Approximation Theory 159 (2) (2009) 308-318.

[18] M. Schmidt, minfunc: unconstrained differentiable multivariate optimization in matlab, http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html (2005).
#
Time: Nov 18, 2018
##
The datasets including all results of spherical $t$-design have been shown in the file named "Results of Paper", the main characters are illustrated as follows
#
t is the degree of spherical polynomials;
#
N is the number of points, N=(t+1)^2;
#
X0 is the initial points set, we use extremal points, you can find in https://web.maths.unsw.edu.au/~rsw/Sphere/;
#
XX is the terminal points set, which is the result of BB method or QN method in computing spherical $t$-design;
#
minY is the minimal singular value, which is used to verify XX is a spherical $t$-design.
##
File named "BBstd" includes all codes for numerical constructing spherical $t$-design in this paper.
#
The objective function is $A_{N,t}$.
#
The original code of BB method comes from Package of "minFunc", which is made by Professor Schmidt. We modified it for a numerical constructing spherical $t$-design.
