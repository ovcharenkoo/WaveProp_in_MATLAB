# **SIMPLE FDTD WAVE PROPAGATION IN MATLAB**

Single-file vectorized implementations of wave propagation in MATLAB. We solve second-order wave equation in displacement formulation in time domain(FDTD). We don't account for derivatives of elastic parameters to keep it simple. Medium is surrounded by simple absorbing sponge boundaries with exponential decay (Cerjan, 1985). 

![Wavefield examples](docs/assets/snap_all.jpg)

### **DISCRETIZATION DETAILS**:
* Finite-Differences in Time Domain (FDTD)
* Explicit time stepping
* O(2,2)
* Conventional stencils derived from Taylor series: 
    * in space [1: -2 :1]/dx2 and [1: -1: -1:1]/4dxdz
    * in time [1: -2 :1]/dt2

oleg.ovcharenko@kaust.edu.sa

vladimir.kazei@kaust.edu.sa
