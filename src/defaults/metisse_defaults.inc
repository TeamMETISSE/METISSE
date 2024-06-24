! Defaults control for METISSE
! Do NOT make any changes to this file
! Use evolve_metisse.in for input
! Make sure to ADD a new line at the end (after / symbol)in the inlist


! NOTE: SSE_input_controls in evolve_metisse.in are read 
! only when METISSE is called directly through its main unit/ standalone mode.
! Otherwise, values provided by the overlying code are used.

    !path is relative to the executable
    METISSE_DIR = '.'

!METISSE_input_controls
    
    ! A metallicity file contains details about
    ! the set of input tracks for a given metallicity,
    ! such as the path to the folder, their metallicity value
    ! and other information/metadata.

    ! METISSE looks for files(s) ending with *metallicity.in*
    ! in the paths provided through 'tracks_dir' and 'tracks_dir_he'
    ! where 'tracks_dir' is for providing the location of the metallicity file(s)
    ! for hydrogen star and  'tracks_dir_he' is for providing the location
    ! of the metallicity file(s) for the naked helium/stripped stars.

    ! It checks for the following condition to find a match in metallicity
    !`(abs(Z_input-Z_required)/MIN(Z_input,Z_required)) > Z_accuracy_limit`,
    ! where Z_input is the metallicity value of the tracks (contained in the metallicity file) and Z_required is the value we want to use.

    ! For a grid of stellar tracks, containing sets of stellar tracks with different metallicities,
    ! the above paths can contain a list of metallicity files,
    ! and METISSE will select the one that is closest to the input metallicity
    ! according to the aforementioned condition.
    
    tracks_dir = ''

    tracks_dir_he = ''
    
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
    ! useful for debugging and single star evolution calculations with implicit mass loss

    write_eep_file = .false.		

    ! If true error messages are written to a file (usually fort.99)
    ! otherwise they are printed on screen
    ! does NOT apply to errors during reading input files
    write_error_to_file = .true.




