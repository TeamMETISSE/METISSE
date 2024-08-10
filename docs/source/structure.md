# Code structure


METISSE has been written completely using Modern Fortran. 
The source code contains two types of files located in the *src* directory: 


## SSE specific subroutines 


METISSE has been developed as an alternative to SSE [(Hurley et al. 2000)](https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract) code and, therefore, contains similar subroutines as SSE. The following subroutines replicate the behaviour of the SSE subroutines externally, with similar names and input/output variables.

| File Name           | Description                                                                                             |
|---------------------|---------------------------------------------------------------------------------------------------|
| METISSE_zcnsts.f90    | Controls metallicity (Z) specific functions.                                         |
| METISSE_star.f90      | Find relevant tracks from the input set and interpolate in mass to <br>get a track of a given mass. |
| METISSE_hrdiag.f90    | Interpolate within the new track to determine stellar parameters at a given age. <br>Also, compute the evolutionary phase of the star including the remnant phases and their properties. |
| METISSE_deltat.f90    | Calculate the time steps depending on the stage of evolution.                                   |
| METISSE_mlwind.f90    | Derive the mass loss through stellar winds.                                                    |
| METISSE_gntage.f90    | Determine the age of a giant star after merger or rejuvenation.                                 |



## Modules and other files

The following Fortran modules contain supporting data structures and utilities. 

| File Name           | Description                                                                                             |
|---------------------|---------------------------------------------------------------------------------------------------|
| track_support.f90   | Contains general data structures and functions needed throughout the program.                     |
| interp_support.f90  | Contains functions required for interpolation.                                                    |
| remnant_support.f90 | Contains functions needed to calculate properties of remnant phases.                               |
| z_support.f90       | Together with METISSE_zcnsts, it reads all input namelists and files, <br>including EEP files, and sets Z parameters and other metallicity-based functions. |
| sse_support.f90     | Contains subroutines to calculate SSE-specific quantities.                                         |



In addition to the above modules, there is an additional file called *METISSE_miscellaneous.f90*. This file contains miscellaneous subroutines that can be called by the binary codes. Ideally, these subroutines should also be organized into a module. However, this is not feasible because the overlying binary evolution codes often have segments written in Fortran 77, and Fortran 77 does not support the use of modules.

It contains:  
| File Name           | Description                                                                                             |
|---------------------|---------------------------------------------------------------------------------------------------------|
| alloc_track.f90     | Allocate the track object.                                                                              |
| dealloc_track.f90   | Deallocate the track object and arrays within.                                                          |
| initialize_front_end| Inform METISSE what code is using it.                                                                   |
| set_star_type      | Set star type to `rejuvenated` before calling star.   


One of these sets of files is used depending on how METISSE is being used.

**In the standalone mode:**
| File Name           | Description                                                                                             |
|---------------------|---------------------------------------------------------------------------------------------------------|
| main_metisse.f90    | Main program for running METISSE. Can only evolve single stars. <br> Reads the input files and sets up relevant parameters and data structures <br>before evolving stars of given masses. |
| evolv_metisse.f90   | Controls the evolution of each star and writes output to files.                                         |
| assign_commons_main.f90 | Assign values for variables used in METISSE from SSE_input_controls.                                    |

**As part of other codes:**
| File Name           | Description                                                                                             |
|---------------------|---------------------------------------------------------------------------------------------------------|
| assign_commons_xyz.f90  | To assign common variables when METISSE is used with other code say 'xyz'.                               |
| comenv_lambda.f90   | Get the appropriate ZAMS radius and calculate common envelope lambda.                                       |

## Workflow

Here is a flowchart describing the workflow of METISSE:

![METISSE's flowchart](METISSE_flowchart.png)
   


<!-- # module details

All kinds tracks in METISSE are stored in Fortran's derived data type (think class in python/c++) called track  which is defined in track_support. It mostly has two components : a header part (e.g. initial_mass, ntrack etc.)  and an array tr which stores file data as an array.  xa, sa, sa_he and tarr  all belong to this type.
xa is a temporary variable. We read all the input data through it and, filter whatever we need and store it either in sa (hydrogen tracks) or sa_he (helium tracks).
tarr is special; it is where we store interpolated tracks. It also makes use some other feature of track such as pars which stores stellar parameters at any given age.

First index of tr is denoters the column, and second denotes row.  Different EEPs actually are row numbers. so in tr, each row is like a time, which corresponds to an evolutionary stage. TPAGB_EEP and cCBurn_EEP let us know from which row/time it enters these stages. 


For reading input data in xa, we have two kind of functions. One is specific to filetype used by Choi et al. 2016 (MIST-EEP files) and other is more general. However they both do the same thing; read file and store header data and tr for each element of xa.  However, at first call they also identify various columns and their locations, controlled by the variable get_columns. These columns are identified by calling the subroutine get_named_columns. After xa is read, tracks are checked for errors through check_tracks . For complete tracks, we check here if any column needs to be converted to log within this subroutine.  We next find zparameters/mcrit/mcutoff using xa (See METISSE's 2020 paper for their definition). Finally we use copy_and_deallocatex to store only the relevant information in sa or sa_he.  We also do other some other checks and set what we call key_cols , where essential columns get reassigned to match to reduced array format of sa/sa_he. -->

<!-- The interpolation is either linear or monotonic piece-wise cubic, depending on the number of tracks available in the neighborhood. The resulting track is a collection of stellar parameters at each EEP. 
METISSE further interpolates within this mass-interpolated track to calculate stellar parameters at any instant. -->

