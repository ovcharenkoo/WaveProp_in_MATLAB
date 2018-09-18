# **Simple FDTD wave propagation in 2D elastic VTI medium**

Single-file vectorized implementation of elastic wave propagation in vertical transversely-isotropic(VTI) medium in MATLAB. The program is solving second-order wave equation in displacement formulation. We don't account for derivatives of elastic parameters. Modelling area is surrounded by simple absorbing sponge boundaries with exponential decay (Cerjan, 1985). 

![Wavefield example](img/snap.jpg)

### **DISCRETIZATION DETAILS**:
* Finite-Differences in Time Domain (FDTD)
* Explicit time stepping
* O(2,2)
* Conventional stencils derived from Taylor series: 
    * in space [1: -2 :1]/dx2 and [1: -1: -1:1]/4dxdz
    * in time [1: -2 :1]/dt2

### **MODEL DETAILS**
* VTI anisotropy (c11, c13, c33, c44, rho)
* Sponge or reflecting boundaries

### **HOW TO USE**: 
Run `elastic_2D_FDTD_O22_VTI.m` in MATLAB

oleg.ovcharenko@kaust.edu.sa

vladimir.kazei@kaust.edu.sa