# DFSR
DFSR is a divergence-free inflow/inlet turbulence generation method developed for large-eddy simulation of the atmospheric boundary layer(ABL) flows. The technique uses the spectral representation method to generate ABL-like turbulence with prescribed two-point flow statistics. DFSR  offers the following key features:
Capable of running on both structured and unstructured meshes
- Uses LAPACK subroutine for matrix factorization 
- Employes the FFT technique for efficient simulation of the velocity-time series 
- Parallel implementation using MPI    

Divergence-free condition is imposed on the generated turbulence by adding a gradient of a velocity potential. The final velocity field is written for each time step in the corresponding case directory.
