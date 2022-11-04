# Theoretical error estimates for computing the matrix logarithm by Pade-type approximants

This is the companion repository to the paper
- L. Aceto, F. Durastante. Theoretical error estimates for computing the matrix logarithm by Pade-type approximants.

**Abstract** In this article we focus on the computation of the matrix logarithm by a Gauss-Legendre rule. This leads to a Padé-type approximant.
Results on Padé approximants to Gauss hypergeometric functions help in finding the corresponding results to the logarithm.

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
