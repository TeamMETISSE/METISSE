! Defaults control for METISSE
! Do NOT make any changes to this file
! Use evolve_metisse.in for input
! Make sure to ADD a new line at the end (after / symbol)in the inlist


! NOTE: SSE_input_controls in evolve_metisse.in are read 
! only when METISSE is called directly through its main unit/ standalone mode.
! Otherwise, values provided by the overlying code are used.
    
!SSE_input_controls

  ! EVOLUTION CONTROLS
   
    ! There are two ways stellar mass values can be supplied to METISSE:

    !1. `read_mass_from_file = .false.`
    !Use `min_mass` to specify the lower limit, `max_mass` for upper limit and `number_of_tracks` to be evolved uniformly distributed in mass between the two limits.
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
    
    ! Initial metallicity
    initial_Z = -1.0

    ! Maximum age in Myrs
    max_age = -1.0   

  ! REMNANT CONTROLS


    ! For White dwarfs, two prescriptions are available for calculating luminosity:

    ! (1) p. 85 of [Shapiro S. L., Teukolsky S. A., 1983](https://ui.adsabs.harvard.edu/abs/1983bhwd.book.....S/abstract)

    ! (2) [Hurley J. R., Shara M. M., 2003, ApJ, 589, 179](https://iopscience.iop.org/article/10.1086/374637)

    ! Options - "Mestel", "Modified_mestel"
    WD_mass_scheme = 'Modified_mestel'

    
    ! Only for White Dwarfs
    Use_initial_final_mass_relation = .false.       

    
    ! The type and mass of the NS/BH can  be  calculated  from  one  of  the  following  prescriptions:

    !(1) Default SSE

    !(2) [Belczynski et al. 2002, ApJ, 572, 407](https://iopscience.iop.org/article/10.1086/340304)

    !(3) [Eldridge J. J., Tout C. A., 2004, MNRAS, 353, 87](https://ui.adsabs.harvard.edu/abs/2004MNRAS.353...87E/abstract)

    !(4) [Belczynski et al. 2008, ApJS, 174, 223](https://iopscience.iop.org/article/10.1086/521026)


    ! Options - "original_SSE", "Belczynski2002", "Belczynski2008", "Eldridge_Tout2004"
    BHNS_mass_scheme = 'Belczynski2008'

    
    ! Maximum neutron star mass 
    ! Suggested 1.8 for BHNS_mass_scheme="original_SSE", 3.0 otherwise

    Max_NS_mass = 3.d0
   
    ! Allow electron capture supernovae

    allow_electron_capture = .true.       
    
  ! TIMESCALE CONTROLS

    
    !Similar to SSE, METISSE determines timesteps as the fractions of the time spent in a phase
    !pts_1 - 95% of MS, HeMS
    !pts_2 - last 5% of MS, cHeBurn, HeHG, HeGB
    !pts_3 - HG, RGB, EAGB, TPAGB
    
    pts_1 = 0.05
    pts_2 = 0.01
    pts_3 = 0.02

  !OUTPUT CONTROLS

    ! 'write_track_to_file' generates a SSE style output file 
    ! only at the END of the evolution

    write_track_to_file = .true.
