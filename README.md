# hmcs-qss
Hamiltonian Monte Carlo (HMC) sampling from quantum state space

------
This repository contains the MATLAB code that I wrote during my M.Sc. thesis project "Hamiltonian Monte Carlo sampling from quantum state space" from August 2014 to May 2015, under the supervision of Prof. Berthold-Georg Englert in Centre for Quantum Technologies and Department of Physics, National University of Singapore.

###### MATLAB code for two-qubit states:
* Cholesky Decomposition with Primitive Prior: `code/cholesky_2qb_flat.m`
* Cholesky Decomposition with Jeffreys Prior or Hedged Prior: `code/cholesky_2qb_non_flat.m`
* Spectral Decomposition with Primitive Prior: `code/spect_2qb_flat.m`
* Fidelity and Distance: `code/Fidelity_Distance.m`

###### Note:
* A main HMC program is not included in this repo.
* For two-qubit states only.

###### Usage:
`code/cholesky_2qb_flat.m`, `code/cholesky_2qb_non_flat.m` and `code/spect_2qb_flat.m` are functions to be called by a main HMC program, in order to compute numerically the density matrix, probabilities, Jacobian determinant and potential gradients based on the angle variables, density matrix dimension and POMs provided by the main HMC program. 

Code is released under the MIT license.
