# Theoretical error estimates for computing the matrix logarithm by Padé-type approximants

This is the companion repository to the paper
- L. Aceto, F. Durastante. Theoretical error estimates for computing the matrix logarithm by Pade-type approximants.

**Abstract** In this article, we focus on the error that is committed when computing
the matrix logarithm using the Gauss–Legendre quadrature rules. These formulas can be
interpreted as Pade approximants of a suitable Gauss hypergeometric function.
Empirical observation tells us that the convergence of these quadratures becomes slow
when the matrix is not close to the matrix identity, thus suggesting the usage of
an inverse scaling and squaring approach for obtaining a matrix with this property.
The novelty of this work is the introduction of error estimates that can be used to
select *a priori* both the number of Legendre points needed to obtain a given accuracy
and the number of inverse scaling and squaring to be performed. We include
some numerical experiments to show the reliability of the estimates introduced.

## Obtain the code

The code can be obtained by doing
```
git clone https://github.com/Cirdans-Home/padelogarithm.git
```
It use code from the [Chebfun project](https://www.chebfun.org/), that is included here as a [Git submodule](https://www.git-scm.com/book/en/v2/Git-Tools-Submodules).
You can install Chebfun on your machine by doing
```
unzip('https://github.com/chebfun/chebfun/archive/master.zip')
movefile('chebfun-master', 'chebfun'), addpath(fullfile(cd,'chebfun')), savepath
```
The examples in the paper were made with the 5.7.0 version, to be sure of using it
you can check out the appropriate commit by doing
```
git submodule init
```
and then add to your Matlab path the version in the cloned repository folder.

## MATLAB's logm

The folder `scalingandsquaring` contains a copy **as is** of MATLAB's `logm`
function (Higham & Relton) with an added print statement that tells the
user what choices have been made for the number of scaling and squaring
step and the number of terms in the diagonal Padé expansion. This has been
used to make the comparison reported in the paper.

For information about that function see
- A. H. Al-Mohy and Nicholas J. Higham, Improved inverse scaling and
       squaring algorithms for the matrix logarithm, SIAM J. Sci. Comput.,
       34(4), (2012), pp. C153-C169.
-   A. H. Al-Mohy, Nicholas J. Higham and Samuel D. Relton, Computing the
       Frechet derivative of the matrix logarithm and estimating the
       condition number, SIAM J. Sci. Comput., 35(4), (2013), C394-C410.

## Advanpix

To focus on the error of the rational approximation more than to the numerical
conditioning of the computation of the inverses, one of the example using the
bound based on pseudospectra uses computation in quadruple precision. These
are made available by the [advanpix library](https://www.advanpix.com/).
