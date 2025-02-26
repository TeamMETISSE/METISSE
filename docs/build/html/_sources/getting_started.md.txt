
# Installation

## Prequisites

### Fortran compiler

METISSE requires a Fortran installation specifically gfortran that comes with gcc/6.4.0 and more recent versions.
Check out [this page](https://fortran-lang.org/learn/os_setup/install_gfortran/) if you need help installing gfortran. 

### Input tracks

The user also needs to provide stellar tracks for METISSE to work. 

A set of sample tracks for testing METISSE can be downloaded from [zenodo](https://zenodo.org/records/14918163). The tracks have been computed using the stellar evolution code [MESA](https://docs.mesastar.org/en/release-r24.03.1/) and are ready for use with METISSE. 
The folder `hydrogen` contains stellar tracks with initial masses between 0.1 and 300 M<sub>$_\odot$</sub> for each metallicity. Similarly the folder `helium` contains tracks of naked helium stars in the mass range 0.3 and 150 M<sub>$_\odot$</sub>. The tracks are for non-rotating stars at solar metallicity (Z=0.02).
A pre-packaged grid of MESA tracks for hydrogen and helium stars with metallicity values ranging from 10<sup>-5</sup> to 10<sup>-1</sup> will be available soon (a future Zenodo entry). 

One can also use their custom set of hydrogen and helium stellar tracks, computed using different input parameters, or even different stellar evolution codes with METISSE. For instructions on how to use a custom set of input stellar tracks with METISSE, refer to [](using_custom_input_tracks.md).

## Code
The code package for METISSE is available at this [GitHub Repository](https://github.com/TeamMETISSE/METISSE).

## Running METISSE

METISSE can be run in two different ways:

1. [Directly as a standalone code](usage_standalone.md) for evolving populations of single stars.
2. [In conjunction with binary evolution codes](usage_other.md) for evolving populations of single and binary stars. 



