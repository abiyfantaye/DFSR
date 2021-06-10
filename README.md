# DFSR
The DFSR (Divergence-free Spectral Representation) method is an inflow/inlet turbulence generation method developed for large-eddy simulation of the atmospheric boundary layer(ABL) flows. The technique uses the spectral representation method to generate ABL-like turbulence with prescribed two-point flow statistics. The code is developed as a utility application using OpenFOAM as a framwork. DFSR  offers the following key features:

- Capable of running on **unstructured** meshes
- Uses **LAPACK** subroutine for matrix factorization 
- Employes the **FFT** technique for efficient time-series computation 
- Parallel implementation using **OpenMPI**    

The divergence-free condition is imposed on the generated turbulence by adding a gradient of a velocity potential instead of using the original procedure reported in the paper. The final velocity field is written for each time step in the corresponding case directory. Mathematical details of the implemented procedure can be found in [Melaku and Bitsuamlak](https://www.sciencedirect.com/science/article/pii/S0167610521000660)

## Installing DFSR
### Requirements
Before compiling the source code, you need to have the following packages already installed. 
- **OpenFOAM** preferably v8  
- Intel's oneAPI MKL package or **LAPACK** library    
- Fastest Fourier Transform in the West(**FFTW**) package

### Build from source  

Follow the following procedure to compile the source code. 

- Clone the code to your computer:

   `$ git clone https://github.com/abiyfantaye/DFSR.git DFSR` 

- Copy the `DFSR/src/DFSRTurb` directory to `$FOAM_UTILITIES/preProcessing/`:

    `$ cp -r DFSR/src/DFSRTurb $FOAM_UTILITIES/preProcessing/`

- Compile the code in `$FOAM_UTILITIES/preProcessing/DFSRTurb` directory:

    `$ cd $FOAM_UTILITIES/preProcessing/DFSRTurb`\
    `$ ./Allwclean`

The turbulence generator is now fully embedded in OpenFOAM and can be launched as `DFSRTurb` in serial or parallel.  

## Running the test case

### Input parameters 

The directory `example/openTerrain` contains an empty domain LES case configured for open terrain condition. The basic parameters for the DFSR method are defined in `constant/DFSRTurbDict`. The method reads the ABL profile from `constant/boundaryData/windProfile/profile` file. This file specifies the variation of the mean velocity, turbulence intesities and legth scales over the hieght in the following tabular format: 

| z[m]  | Uav[m/s] |  Iu   |   Iv   |   Iw   | Lu [m] | Lv[m]  | Lw[m]  |
| ----- | -------- | ----- | ------ | ------ | ------ | ------ | ------ |
| 0.02  | 14.0056  | 0.160 | 0.120  | 0.120  | 0.750  | 0.500  | 0.300  |
| ..... | .......  | ..... | .....  | .....  | .....  | .....  | .....  |
| ..... | .......  | ..... | .....  | .....  | .....  | .....  | .....  |

The parameters that need to specified in `constant/DFSRTurbDict` file area: 

| Parameter | Description | Optional | Default |
| --- | --- | --- | --- |
| `patchName` | Name of the patch to apply the inflow | no | -- |
| `startTime` | Starting time of the generated terbulence | yes | 0.0|
| `endTime` | Ending time of the generated terbulence  | no | --|
| `timeStep` | Time step of the inflow. Can be different from the time step of the solver. Should be 1/(2*fmax) to take advantage of the FFT algorithm. | no | --|
| `fMax` | The cutt of frequency | no | --|
| `nFreq` | The number of frequency steps/segments | no | 2048 |
| `nInterp` | Number of spectral interpolation frequencies | no | 50 |
| `cohUav` | Constant mean velocity used to calculate the coherency function | no | --|
| `seed` | Seed of the random generator used to sample random phase angles. Specify any random positive integer | no | --|
| `C` | Coherency decay coefficients.  | no | --|
| `scaleI` | Factors to scale turbulence intensities and length scale profiles in each direction helps to account for decay of turbulence | no | --|
| `scaleL` | Factors to scale turbulence intensities and length scale profiles in each direction helps to account for decay of turbulence  | no | --|
| `correctDivergence` | Apply div-free corrections or not  | no | --|
| `writeSamples` | Write the time-series for sampled points | no | --|
| `writeInflow` | Write the final inflow data | no | --|



