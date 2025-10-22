
# Installation

METISSE can be used in two different ways:
1. **Directly as a standalone code** for evolving populations of single stars.
2. **In conjunction with binary evolution codes** for evolving populations of single and binary stars. 

Depending on how you intend to use METISSE, the instructions for downloading and running it may vary.

## Get the code
The code package for METISSE is available at the [GitHub Repository](https://github.com/TeamMETISSE/METISSE).
If you plan to run METISSE by itself (standalone mode), obtain the code using one of the following methods:
### Clone the repository 

If you have Git installed, open a terminal and run:

```console

git clone https://github.com/TeamMETISSE/METISSE.git

```
This will create a local copy of the METISSE source code in a folder named METISSE. 


### Download as a ZIP archive

If you prefer not to use Git:

- Go to the METISSE GitHub repository
- Click the green “Code” button.
- Select “Download ZIP”.
- Extract the downloaded archive to your desired location.

To download METISSE with other codes, see code-specific instructions on [this page](usage_other.md)

## Get the compiler

Regardless of how you use it, METISSE must be compiled before use.
It requires a Fortran compiler — specifically gfortran, which comes with GCC 6.4.0 or newer version.

Check your compiler version with:

```console

gcc --version

```

Check out [this page](https://fortran-lang.org/learn/os_setup/install_gfortran/) if you need help installing gfortran. 

## Input tracks

METISSE also requires a set of stellar tracks to work. 

A set of sample tracks for testing METISSE can be downloaded from [zenodo](https://zenodo.org/records/14918163).
The tracks have been computed using the stellar evolution code [MESA](https://docs.mesastar.org/en/release-r24.03.1/) and are ready for use with METISSE. 
- The folder `hydrogen` contains stellar tracks with initial masses between 0.1 and 300 M<sub>$_\odot$</sub> for each metallicity.
- The folder `helium` contains tracks of naked helium stars in the mass range 0.3 and 150 M<sub>$_\odot$</sub>. 
- The tracks are for non-rotating stars at solar metallicity (Z=0.02).
<!-- A pre-packaged grid of MESA tracks for hydrogen and helium stars with metallicity values ranging from 10<sup>-5</sup> to 10<sup>-1</sup> will be available soon (a future Zenodo entry).  -->

One can also use their custom set of hydrogen and helium stellar tracks, computed using different input parameters, or even different stellar evolution codes with METISSE. For instructions on how to use a custom set of input stellar tracks with METISSE, refer to [](using_custom_input_tracks.md).


## Run METISSE

Refer to the following pages for build instructions and examples for METISSE:
 
1. [Using METISSE directly](usage_standalone.md).
2. [Using METISSE with binary evolution codes](usage_other.md).

