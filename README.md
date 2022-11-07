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
