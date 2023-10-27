# All electromagnetic scattering bodies are matrix-valued oscillators
Supporting data and scripts for simulation and plotting for the [paper](https://arxiv.org/abs/2211.04457).

## Prerequisites
1. Integration in frequency domain is done with a Legendre-Gaussian quadrature. The Matlab function _legpts.m_ provides Legendre points and Gauss-Legendre quadrature weights. See http://www.chebfun.org/ for Chebfun information.
   __For testing purposes, a total of 101 frequency sampling points is assumed. For the results as reported on the paper, 2001 frequencies are required for good convergence.__

2. Integration in 2D spatial domain is done with the trapezoidal rule. The Matlab function _trapez_mat.m_ generates a 2D map of trapezoidal weights.

## Figure 1
Figure 1 demonstates the properties of T matrix using the scattering of an elliptical cylinder as an example.

The simulation is done with a integral-equation based fast direct solver. This technique is based on the Discrete Dipole Approximation (DDA) method, with a distretization of Duan-Rokhlin method and a fast direct solver method enabled by Fast Multipople Method. 
For a more in-depth discussion of this techniqure, we recommend reading [Fullwave design of cm-scale cylindrical metasurfaces via fast direct solvers](https://arxiv.org/abs/2308.08569).

In this problem, we employed a simple implementation of this simulation technique, with the following 3 Matlab functions
1. _planck_win_get_D.m_ computes the D0, D1 factors used to obtain the Duan-Rokhlin correction terms.
2. _compute_T.m_ returns the T matrix with Duan-Rokhlin corrections for the given gemetry and the given frequency.
3. _BTTB_matvec.m_ performs fast BTTB matrix-vector mutiplication y = Tx where T is a nm by nm block toeplitz matrix with toeplitz blocks (BTTB) and x is an nm by 1 vector.

With these 3 function files, _fig1b.m_, _fig1c.m_, _fig1d.m_, _fig1e.m_ are sufficienct in producing all the data and figures for the corresponding subplots. For accurate results, remember to set number of frequencies _nl = 2001_ or larger.

## Figure 2
Figure 2 reports the key findings of the matrix oscillator bound for Near-Field Radiative Heat Transfer (NFRHT) and the comparisons to previous studies.

Fig2b can be produced by simply running _fig2b.m_, which includes the precomputed data of planar heat transfer coefficients (HTC) and bounds for different materials.

Fig2c can be produced with the plotting script _fig2c.m_ and the precomputed spectral HTC data in _spectral_htc.mat_.

Fig3d can be produced by simply running _fig2d.m_, which includes the precomputed data of optimal frequencies for different materials.

## Figure 3
Figure 3 presents the data supporting the domain monotonicity property.

The static scattered field from the double-tip structure is modeled in COMSOL Multiphyics, and the relevant data for plotting are summarized in _Escat_tip.m_. One can run _fig3.m_ to plot the data.

## Contact
If you have questions or comments regarding the paper or using the code, please feel free to reach out to Lang (lang.zhang@yale.edu) or Owen (owen.miller@yale.edu).
