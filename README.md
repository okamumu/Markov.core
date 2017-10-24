# Markov.src

Fortran programs for Markov chains.

## Modules

### matrix directory

#### Fundamentals of matrix operation

- sparse (sparse_mod.f90): utilities for sparse matrix. The first index (0 or 1) is defined in this module.
- spblas (spblas_mod.f90): BLAS for sparse matrix. This module depends on `sparse`.
- kron (kron_mod.f90): The fundamental functions for Kronecker product and sum for dense matrix.
- spkron (spkron_mod.f90): The fundamental functions for Kronecker product and sum for sparse matrix. This module depends on `spblas`.

#### Krylov subspace

- arnoldi (arnoldi_mod.f90): Arnoldi processes for dense and sparse matrices. The functions for sparse matrix depends on `spblas`.

#### Gauss-Seidel-type algorithm

- gsstep (gsstep_mod.f90): The fundamental functions of Gauss-Seidel algorithm for dense and sparse matrices. The functions for sparse matrix depends on `sparse`.
- gsstep_mm (gstep_mm_mod.f90): The fundamental functions of a variant of Gauss-Seidel algorithm (experimental).

### math directory

#### Math & Distribution

- gamma (gamma_mod.f90): Mathematical functions of the gamma family.
- gamma_dist (gamma_dist_mod.f90): pdf and cdf of the gamma distribution. This depends on `gamma`.
- erlang_dist (erlang_dist_mod.f90): pdf and cdf of the Erlang distribution. This depends on `gamma_dist`.
- poisson (poisson_mod.f90): pmf and rightbound of the Poisson distribution.
- accerelation (quad_extrea_mod.f90): Acceleration algorithms (experimental).

#### Integral

- gaussinte (gaussinte_mod.f90): Integral with Gauss-Legendre quadrature.

### mexp directory

#### Matrix exponential (general)

- mpow (mpow_mod.f90): Matrix power for dense and sparse matrices. The functions for sparse matrix depends on `sparse` and call the matrix power function for dense matrix after converting the sparse matrix to a dense matrix.
- mexp_pade (mexp_pade_mod.f90): Matrix exponential with Pade approximation for dense and sparse matrices. The functions for sparse matrix depends on `sparse` and call the Pade approximation for dense matrix after converting the sparse matrix to a dense matrix.

#### Matrix exponential for continuous-time Markov chain (CTMC)

- unif_matrix (unif_matrix_mod.f90): The uniformed CTMC kernel for dense and sparse CTMC kernel. The functions for sparse kernel depends on `sparse`.
- mexp_unif (mexp_unif_mod.f90): The uniformization of dense and sparse CTMC kernels. This supports both cases of vector-by-matrix and matrix-by-matrix in the computation of matrix exponential. The functions for sparse kernel depends on `spblas`.
- mexpint_unif (mexpint_unif_mod.f90): The uniformization for the integral of dense and sparse CTMC kernels. This supports both cases of vector-by-matrix and matrix-by-matrix in the computation of matrix exponential. The functions for sparse kernel depends on `spblas`.
- mexpconv_unif (mexpconv_unif_mod.f90): The uniformization for the convolution integral of dense and sparse CTMC kernels. This supports only the case of vector-by-matrix. The functions for sparse kernel depends on `spblas`.

#### Matrix exponential for continuous-time Markov chain (CTMC) with Markovian arrival process (MAP) structure

- map_unif_matrix (map_unif_matrix_mod.f90): The uniformed CTMC kernels for dense and sparse CTMC kernel with MAP structure. The functions for sparse kernel depends on `sparse`.
- map_mexp_unif (map_mexp_unif_mod.f90): The uniformization of dense and sparse CTMC kernels with MAP structure. This supports only the case of vector-by-matrix. The functions for sparse kernel depends on `spblas`.
- map_mexpconv_unif (map_mexpconv_unif_mod.f90): The uniformization for the convolution integral of dense and sparse CTMC kernels with MAP structure. This supports only the case of vector-by-matrix. The functions for sparse kernel depends on `spblas`.

### st directory

#### Stationary analysis for continuous-time Markov chain (CTMC)

- markovst_dgesv (markovst_dgesv_mod.f90): The stationary vector of dense and sparse CTMC kernels with the dgesv routine for a dense matrix in the common BLAS. This depends on `sparse` and is an experimental module.
- markovst_gth (markovst_gth_mod.f90): The stationary vector of dense and sparse CTMC kernels with the GTH algorithm for a dense matrix. This depends on `sparse` in the functions for sparse matrix.
- markovst_gs (markovst_gth_mod.f90): The stationary vector of dense and sparse CTMC kernels with Gauss-Seidel algorithm. This module depends on `spblas` and `gsstep`.
- markovst_sor (markovst_sor_mod.f90): The stationary vector of dense and sparse CTMC kernels with SOR (successive over-relaxation) algorithm. This module depends on `spblas` and `gsstep`. (experimental)

### inv directory

#### Inverse matrix of the kernels of continuous-time Markov chain (CTMC)

- markovinv_dgesv (markovinv_dgesv_mod.f90): The inverse of dense and sparse CTMC kernels with the dgesv routine for a dense matrix in the common BLAS. This depends on `sparse`.
- markovinv_gs (markovinv_gs_mod.f90): The inverse of dense and sparse CTMC kernels with Gauss-Seidel algorithm. This module depends on `gsstep`.

### stsen directory

#### Sensitivity function of stationary vector of continuous-time Markov chain (CTMC)

- markovstsen_dgesv (markovstsen_dgesv_mod.f90): The sensitivity function (the first derivative) of stationary vector for dense and sparse CTMC kernels with the dgesv routine for a dense matrix in the common BLAS. This depends on `sparse`.
- markovstsen_gs (markovstsen_gs_mod.f90): The sensitivity function (the first derivative) of stationary vector for dense and sparse CTMC kernels with Gauss-Seidel algorithm. This module depends on `spblas` and `gsstep`.

### qst directory

#### Quasi-stationary vector of continuous-time Markov chain (CTMC): Experimental

- markovqst_gs (markovqst_gs_mod.f90): The quasi-stationary distribution of dense and sparse CTMC kernels with Gauss-Seidel-type algorithm. This depends on `gsstep`. (experimental)
- markovqst2_gs (markovqst2_gs_mod.f90): The quasi-stationary distribution and right eigenvector of dense and sparse CTMC kernels with Gauss-Seidel-type algorithm. This depends on `gsstep`. (experimental)

### em directory

#### EM algorithms for phase-type (PH) distribution

- gph_estep_wtime (gph_estep_wtime_mod.f90): E-step for general PH distribution (dense and sparse kernels) with point samples and their weights. This module depends on `poisson`, `sparse`, `spblas`, `mexp_unif` and `mexpconv_unif`.
- gph_estep_group (gph_estep_group_mod.f90): E-step for general PH distribution with generalized grouped data which includes simple point samples, grouped data and grouped data with missing samples. This module depends on `poisson`, `sparse`, `spblas`, `gamma`, `mexp_unif` and `mexpconv_unif`.
- gph_estep_grouppoi (gph_estep_grouppoi_mod.f90): E-step for general PH distribution (dense and sparse kernels) with generalized grouped data which includes simple point samples, grouped data and grouped data with missing samples, provided that the total number of samples follows a Poisson random variable. This module depends on `poisson`, `sparse`, `spblas`, `gamma`, `mexp_unif` and `mexpconv_unif`.
- gph_mstep (gph_mstep_mod.f90): M-step for general PH distribution (dense and sparse kernels). This depends on `sparse`.
- cf1_mstep (cf1_mstep_mod.f90): M-step for canonical form 1 (CF1) with both dense and sparse kernels. This depends on `sparse` and `gph_mstep`.
- herlang_emstep (herlang_emstep_mod.f90): E-step and M-step for hyper-Erlang distribution with a fixed shape parameter vector. The weighted samples, generalized grouped data and generalized grouped data with a Poisson total number are available. This depends on `gamma` and `erlang_dist`.

#### EM algorithms for Markovian arrival process (MAP)

- map_estep_group (map_estep_group_mod.f90): E-step for general MAP with dense and sparse kernels when the generalized grouped data is given. This module depends on `sparse`, `spblas`, `poisson`, `map_mexp_unif`, and `map_mexpconv_unif`.
- map_mstep (map_mstep_mod.f90): M-step for general MAP with dense and sparse kernels. This depends on `sparse`.
- erhmm_estep_time (erhmm_estep_time_mod.f90): E-step for hidden Markov model with Erlang outputs (ERHMM). Only the time data is available. This module depends on `erlang_dist`.
- erhmm_mstep (erhmm_mstep_mod.f90): M-step for hidden Markov model with Erlang outputs (ERHMM).
