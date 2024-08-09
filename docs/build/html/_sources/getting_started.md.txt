
# Getting started

## Prequisites

METISSE requires a fortran installation specifically gfortran that comes with gcc/6.4.0 and more recent versions.
Checkout [this page](https://fortran-lang.org/learn/os_setup/install_gfortran/) if you need help installing gfortran. 


## Code
The code package for METISSE is available at this [GitHub Repository](https://github.com/TeamMETISSE/METISSE).


## Input tracks

The user also needs to provide evolutionary stellar tracks for METISSE to work. 

A pre-packaged grid of evolutionary tracks for hydrogen and helium stars can be downloaded from (a future zenodo entry). These tracks have been computed using the stellar evolution code [MESA](https://docs.mesastar.org/en/release-r24.03.1/) and are ready for use with METISSE. 
The tracks are for non-rotating stars with metallicity values ranging from 10<sup>-5</sup> to 10<sup>-1</sup>. 
The folder `hydrogen` contains stellar tracks with initial masses between 0.1 and 300 M<sub>$_\odot$</sub> for each metallicity. Similarily the folder `helium` contains tracks of naked helium stars in the mass range 0.3 and 150 M<sub>$_\odot$</sub>.  


METISSE can also use stellar tracks computed using different input parameters, or even different stellar evolution codes. [](using_custom_input_tracks.md) describes how to use a custom set of input stellar tracks with METISSE.


## Running METISSE

METISSE can be run in two different ways:

1. [Directly as a standalone code](usage_standalone.md) for evolving populations of single stars.
2. [In conjunction with binary evolution codes](usage_other.md) for evolving populations of single and binary stars. 



