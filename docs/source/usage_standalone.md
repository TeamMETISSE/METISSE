# Using METISSE as a standalone code 

## Supplying input 

In the standalone mode values of mass, metallicity and other input parameters
are supplied using `SSE_input_controls` namelist contained in *evolve_metisse.in* file. 

```
&SSE_input_controls

  ! EVOLUTION CONTROLS
   
    ! For reading input masses for new tracks
    read_mass_from_file = .false.      

    input_mass_file = ''

    number_of_tracks = 0

    max_mass = -1.0

    min_mass = -1.0
    
    initial_Z = -1.0

    ! Maximum age in Myrs
    max_age = -1.0   

  ! REMNANT CONTROLS


    ! Options - "Mestel", "Modified_mestel"
    WD_mass_scheme = 'Modified_mestel'

    
    ! Only for White Dwarfs
    Use_initial_final_mass_relation = .false.       


    ! Options - "original_SSE", "Belczynski2002", "Belczynski2008", "Eldridge_Tout2004"
    BHNS_mass_scheme = 'Belczynski2008'

    
    ! Maximum neutron star mass 
    ! Suggested 1.8 for BHNS_mass_scheme="original_SSE", 3.0 otherwise

    Max_NS_mass = 3.d0
   
    ! Allow electron capture supernovae

    allow_electron_capture = .true.       
    
  ! TIMESCALE CONTROLS

    pts_1 = 0.05
    pts_2 = 0.01
    pts_3 = 0.02

  !OUTPUT CONTROLS

    ! 'write_track_to_file' generates a SSE-style output file 
    ! only at the END of the evolution

    write_track_to_file = .true.

/

/ ! end of SSE_controls inlist

```

`SSE_input_controls` is only used in standalone mode, it is ignored otherwise and input parameters provided by the overlying code are used.  

*evolve_metisse.in* also contains another Fortran namelist `METISSE_input_controls`. This namelist contains input parameters specific to METISSE, including `metallicity_file_list`.


```
&METISSE_input_controls

    ! A metallicity file contains details about 
    ! the set of input tracks for a given metallicity,
    ! such as the path to the folder, their metallicity value
    ! and other information/metadata (see metallicity_defaults.in)

    ! The option 'metallicity_file_list' is used for providing 
    ! path/location of that metallicity file.
    ! In the case of a grid of stellar tracks,
    ! with folders containing tracks of various metallicities,
    ! location of the metallicity file for each folder/metallicity
    ! can be provided as a list of comma-separated strings 
    ! for up to 20 metallicities.
    ! For example: metallicity_file_list = 'path1',
    !                                      'path2',
    !                                       ...
    !                                      'path20'

    
    metallicity_file_list = ''


    ! if (abs(Z_input-Z_required)/MIN(Z_input,Z_required)) > Z_accuracy_limit
    Z_accuracy_limit = 1d-2

    ! INTERPOLATION CONTROLS

    ! Skip interpolation in mass if there is already
    ! an input track with initial_mass within the 'mass_accuracy_limit'

    mass_accuracy_limit = 1d-4

    ! OTHER REMNANT CONTROLS
    ! If true, 'construct_wd_track' is used (for low-mass stars) to construct the track between 
    ! Thermally-Pulsating AGB phase or tip of the AGB to the white dwarf cooling track
    ! It is useful if input tracks do not contain this phase
    ! but can be used otherwise too.

    construct_wd_track = .true.

    
    ! OUTPUT CONTROLS

    ! if true, 'verbose' prints useful details when reading the files

    verbose = .false. 

    ! 'write_eep_file' generates MIST style output file 
    ! at EVERY step of mass interpolation
    ! useful for debugging and single-star evolution calculations with implicit mass loss

    write_eep_file = .false.    

    !OUTPUT_DIR = ''

/

```

Both are Fortran namelists, so comments (!) and blank lines can be used freely. Characters are **case-insensitive**. Although make sure to leave a blank line at the end of the file (after the `/` symbol)

Refer to `src/defaults/evolve_metisse_defaults.inc` for variable names and their default values. **Do not modify any file inside the defaults folder**.

## Running METISSE 

*makefile* in folder *make* contains all necessary instructions to compile METISSE in the standalone mode. 

Just do `./mk` to compile the package.

To run METISSE in the standalone mode, simply do:

`./metisse`

<!-- METISSE will evolve a 1 M$_\odot$ star, of input metallicity upto 10 Gyr for you. :blush: -->

Check *output* directory for output data files.

## Output 

METISSE can produce two types of output files in the standalone mode:

**1. files ending with .dat :**

SSE-like output files, controlled by `write_track_to_file` in SSE_input_controls.
These contain the following stellar parameters until max_age. Time and age at hydrogen ZAMS are assumed to be zero.

| Column Header | Description |
|-----------------|-----------------|
| time | Physical time [Myr] |
| age | Age of star [Myr] |
| mass | Current mass of the star [M$_\odot$] |
| core_mass | Mass of dominant core [M$_\odot$] |
| He_core | Mass of helium core [M$_\odot$] |
| CO_core | Mass of carbon-oxygen core [M$_\odot$] |
| log_L | Log of surface luminosity [L$_\odot$] |
| log_Teff | Log of effective temperature [K] |
| log_radius | Log of radius [R$_\odot$] |
| phase | SSE stellar type/phase |
| e | Extra information about WD track (to be removed in future) |



**2. files ending with .eep :**

For debugging purposes, METISSE can write mass-interpolated track to output file with the same columns as input files (plus phase column) and a MIST-style file structure. 
This output file only contains data from ZAMS to the end of nuclear burning phases, i.e., it does not contain information about remnant phase. 
It is controlled by `write_eep_file` in METISSE_input_controls.


# Examples


