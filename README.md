# Approximate Joint Diagonalization (AJD)

This provides functions to approximately diagonalize a set of 2x2 complex, hermitian matrices.
The function `sweepAJD!(ops,U)` takes an array of matrices `ops` and a matrix `U` and updates both based on the method of Jacobi angles as described in [Gygi et al.][1].

Additionally, if the package `Optim` is used, `step_SteepestDescent!(ops,U)` provides the same function with an alternative algorithm as described in [Marzari et al.][2].

[1]: https://doi.org/10.1016/S0010-4655(03)00315-1
[2]: https://link.aps.org/doi/10.1103/RevModPhys.84.1419
