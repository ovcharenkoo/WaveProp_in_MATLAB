# **Simple FDTD wave propagation in 2D elastic isotropic medium**

Single-file vectorized implementation of elastic wave propagation in MATLAB. The program is solving second-order wave equation in displacement formulation. We don't account for derivatives of elastic parameters. Modelling area is surrounded by simple absorbing sponge boundaries with exponential decay (Cerjan, 1985). 

![Wavefield example](doc/snap.jpg)

### **DISCRETIZATION DETAILS**:
* Finite-Differences in Time Domain (FDTD)
* Explicit time stepping
* O(2,2)
* Conventional stencils derived from Taylor series: 
    * in space [1: -2 :1]/dx2 and [1: -1: -1:1]/4dxdz
    * in time [1: -2 :1]/dt2

### **MODEL DETAILS**
* Isotropic (vp, vs, rho)
* Sponge or reflecting boundaries

### **HOW TO USE**: 
Run `elastic_2D_FDTD_O22.m` in MATLAB

oleg.ovcahrenko@kaust.edu.sa

vladimir.kazei@kaust.edu.sa