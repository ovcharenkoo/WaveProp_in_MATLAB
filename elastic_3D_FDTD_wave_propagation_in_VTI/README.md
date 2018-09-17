# **Simple FDTD wave propagation in 3D elastic VTI medium**

[NOT DEBUGGED YET]

Single-file vectorized implementation of elastic wave propagation a medium with vertical transverse isotropy (VTI) in MATLAB. The program is solving second-order wave equation in displacement formulation. We don't account for derivatives of elastic parameters. Modelling area is surrounded by simple absorbing sponge boundaries with exponential decay (Cerjan, 1985). 

![Wavefield example](doc/snap.jpg)

### **DISCRETIZATION DETAILS**:
* Finite-Differences in Time Domain (FDTD)
* Explicit time stepping
* O(2,2)
* Conventional stencils derived from Taylor series: 
    * in space [1: -2 :1]/dx2 and [1: -1: -1:1]/4dxdz
    * in time [1: -2 :1]/dt2

### **MODEL DETAILS**
* VTI anisotropy (c11, c13, c33, c44, c66, rho)
* Spounge of reflecting boundaries

### **HOW TO USE**: 
Run `elastic_3D_FDTD_O22_VTI.m` in MATLAB

oleg.ovcahrenko@kaust.edu.sa

vladimir.kazei@kaust.edu.sa
