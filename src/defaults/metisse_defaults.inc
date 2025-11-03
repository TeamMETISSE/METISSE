! Defaults for METISSE_input_controls
! Do NOT make any changes to this file

!METISSE_input_controls
    
  ! TRACK CONTROLS

    
    ! METISSE looks for files(s) ending with *metallicity.in*
    ! in the paths provided through 'METALLICITY_DIR' and 'METALLICITY_DIR_HE'
    ! where 'METALLICITY_DIR' is for providing the location of the metallicity file(s)
    ! for hydrogen star and  'METALLICITY_DIR_HE' is for providing the location
    ! of the metallicity file(s) for the naked helium/stripped stars.
    
    
    ! Absolute path to the folder containing the metallicity files
    ! for the hydrogen stars.
    ! METISSE will stop if METALLICITY_DIR is an empty string.
    
    METALLICITY_DIR = ''

    ! Absolute path to the folder containing the metallicity files
    ! for the naked helium stars.
    ! If METALLICITY_DIR_HE is an empty string,
    ! SSE fitting formulae for naked helium stars are used.
    
    METALLICITY_DIR_HE = ''
    
    
    ! METISSE checks for a match in metallicity based on the condition:
    !`(abs(Z_files-Z_input)/MIN(Z_files,Z_input)) > Z_accuracy_limit`,
    ! where Z_files is the metallicity value of the tracks (from the metallicity file) and Z_input is the value we want to use.
    
    Z_accuracy_limit = 1d-2

    ! METISSE skips mass-interpolation and uses a neighboring mass track
    ! if the absoulte difference between the initial mass of that track
    ! and the input mass is less than the 'mass_accuracy_limit'

    mass_accuracy_limit = 1d-4

    
  ! MISCELLANEOUS CONTROLS

    ! if true, print important details when reading the input tracks on the terminal

    verbose = .false. 

    ! Generate a MIST-style output file at every step of mass interpolation.
    ! Useful when computing evolution of isolated single stars
    ! with implicit mass loss and for debugging purposes.

    write_eep_file = .false.		

    ! If true error messages are written to a file (usually fort.99)
    ! otherwise they are printed on screen
    ! (does NOT apply to errors during reading input files)
    write_error_to_file = .true.

    ! If .true., artificially construct the HRD (Hertzsprung-Russell diagram)
    ! track between tip of the AGB to the WD cooling track (for low-mass stars).
    ! It is important if input tracks do not contain this phase
    ! but can be used regardless.

    construct_postagb_track = .false.

