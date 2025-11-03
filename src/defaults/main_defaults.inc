! Defaults for SSE_input_controls
! Do NOT make any changes to this file

! NOTE: SSE_input_controls are read only when METISSE is called directly though its main unit.
! Otherwise, values provided by the overlying code are used.
    
!SSE_input_controls

  ! EVOLUTION CONTROLS
   
    ! Initial metallicity
    initial_Z = -1.0

    ! Maximum age in Myrs
    max_age = -1.0
    
    ! There are two ways stellar mass values can be provided to METISSE:

    !1. `read_mass_from_file = .false.`
    !Use `min_mass` to specify the lower limit, `max_mass` for upper limit and `number_of_tracks` to be evolved distributed in mass between the two limits according to the `sampling_scheme`
    !To evolve just one star, specify `number_of_tracks = 1` and input mass value in `min_mass`. Ignore max_mass.
    !OR
    !2. `read_mass_from_file = .true.`

    !Provide input masses in a text file (one mass value in a line) and specify the location of that file in `input_mass_file`.
    !Use `number_of_tracks` to tell how many mass values are to be read from the file.

    read_mass_from_file = .false.      

    input_mass_file = ''

    number_of_tracks = 0

    max_mass = -1.0

    min_mass = -1.0
    
    !IMF sampling scheme
    ! (1) `Uniform` = uniform distribution
    ! (2) `Kroupa2001` = Kroupa, P., 2001 distribution
    sampling_scheme = 'Uniform'
    

  ! REMNANT CONTROLS

    ! For White dwarfs, two prescriptions are available for calculating luminosity:

    !(1) “Mestel” - Shapiro S. L., Teukolsky S. A., 1983
    !(2) “Modified_mestel” - Hurley J. R., Shara M. M., 2003
    WD_mass_scheme = 'Modified_mestel'

    ! Only for White Dwarfs
    !If .true. use the initial final mass relation for white dwarfs
    ! from Han, Z, Podsiadlowski, P., Eggleton, P. P. 1995.
    Use_initial_final_mass_relation = .false.
    
    ! The type and mass of the NS/BH can  be  calculated  from  one  of  the  following  prescriptions:

    !(1) “original_SSE” - Hurley et al. 2000
    !(2) “Belczynski2002” - Belczynski et al. 2002
    !(3) “Belczynski2008” - Belczynski et al. 2008
    !(4) “Eldridge_Tout2004” - Eldridge J. J., Tout C. A., 2004

    BHNS_mass_scheme = 'Belczynski2008'
    
    ! Maximum neutron star mass.
    ! Recommended 1.8 for BHNS_mass_scheme = “original_SSE”, 3.0 otherwise.

    Max_NS_mass = 3.d0
   
    ! Allow electron capture supernovae

    allow_electron_capture = .true.       
    
  ! TIMESTEP CONTROLS
    
    !Like SSE, METISSE uses `pts` variables to determine the timestep
    ! as the decimal fractions of the time spent in a phase.
    !pts_1 - 95% of MS, HeMS
    !pts_2 - last 5% of MS, cHeBurn, HeHG, HeGB
    !pts_3 - HG, RGB, EAGB, TPAGB
    
    pts_1 = 0.05
    pts_2 = 0.01
    pts_3 = 0.02

  !OUTPUT CONTROLS

    ! Generate a SSE-style output file (file ending in .dat) at the END of the evolution.

    write_output_to_file = .false.
    
  ! COLUMN CONTROLS
    
    ! By default METISSE only interpolates in the essential columns
    ! and the additional columns specified in the format_file,
    ! as interpolating in all columns of input files can slow down computations
    ! Using fewer columns means fewer calculations, therefore faster runs.
    ! If read_all_columns is true then all columns are used.
    ! The interpolated quantities are printed in MIST-style files
    ! if write_eep_file is true.
    
    read_all_columns = .false.

    ! Extra columns to be used for interpolation if read_all_columns is false
    ! Up to 100 column names can be listed as strings, separated by commas.
    ! (irrespective of the order or whitespace between the strings).

    extra_columns = ''

    ! Alternative way for listing extra column names.
    ! List one column name per line in the extra_columns_file
    ! and provide path of the file here.
    
    extra_columns_file = ''
