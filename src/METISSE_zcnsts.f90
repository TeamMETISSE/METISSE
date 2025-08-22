subroutine METISSE_zcnsts(z,zpars,ierr)
    use track_support
    use z_support
    use c_m_interface

    real(dp), intent(in) :: z
    real(dp), intent(out) :: zpars(20)
    integer, intent(out) :: ierr
    
    character(LEN=strlen), allocatable :: track_list(:)
    character(LEN=strlen) :: USE_DIR, find_cmd, rnd, infile, temp_filename
    integer :: i,j,nloop, num_tracks
    logical :: load_tracks, debug
    
    debug = .false.
    ierr = 0
    ! At this point in the code front_end might not be assigned
    ! So we return ierr and let zcnsts.f of the overlying code
    ! decide how to deal with errors.

    code_error = .false.
    
    if (front_end <0) then
        print*, 'METISSE error: front_end is not initialized'
        ierr = 1; return
    endif
    
    if (debug) print*, 'in METISSE_zcsnts',z
    
    ! read one set of stellar tracks (of input Z)
    load_tracks = .false.
    use_sse_NHe = .false.
    
    if (allocated(sa) .eqv. .true.) then
        ! tracks have been loaded at least once, for initial_Z
        
        ! New tracks need to be loaded
        ! if input metallicity 'z' has changed significantly from the old 'initial_z'
        if (relative_diff(initial_Z,z) .ge. Z_accuracy_limit) load_tracks = .true.
        if (all(abs(zpars) == 0.d0)) load_tracks = .true.

        ! or maybe metallicity is the same, but paths may have changed
        ! (for example, for sets of tracks computed with different stellar parameters)
        ! Currently only for cosmic, as it can change path_to_tracks mid-computation
        ! through its python wrapper
        
        if (front_end == COSMIC) call check_path_change(load_tracks)
    
        if (load_tracks) then
            if (mode/=0) then
              print*, 'METISSE error: cannot change path or metallicity mid-run when using mpi'
              ierr = 1
              return
            endif
        else
            if (debug) print*, 'No change in metallicity or paths, exiting METISSE_zcnsts',initial_Z,z
            return
        endif
    else
        !first entry, read inputs and setup variables
        if (debug) print*, 'Initializing METISSE_zcnsts'

        load_tracks = .true.
        
        !path is relative to the executable
        METISSE_DIR = '.'

        !read default options first
        include 'defaults/metisse_defaults.inc'

        !read user inputs
        
        select case(front_end)
        case(test)
             call get_test_inputs()
        case(main)
            infile = trim(METISSE_DIR)// '/main.input'
            call read_main_input(infile,ierr)
            if (.not. defined(initial_Z ))then
                print*,"METISSE error: initial_Z is not defined in ",trim(infile)
                ierr = 1
            endif
            if (ierr/=0) call stop_code
            
            infile = trim(METISSE_DIR)// '/metisse.input'
            call read_metisse_input(infile,ierr)
            if (ierr/=0) call stop_code
        case(bse)
            infile = 'evolve_metisse.in'
            call read_metisse_input(infile,ierr)
            if (ierr/=0) call stop_code
        case(COSMIC)
             call get_COSMIC_input()
        case(AMUSE)
            ! If AMUSE has not set METALLICITY_DIR, it will use the defaults from METISSE
            if (len(trim(amuse_metallicity_dir)) > 0) METALLICITY_DIR = amuse_metallicity_dir
            if (len(trim(amuse_metallicity_dir_he)) > 0) METALLICITY_DIR_HE = amuse_metallicity_dir_he
            !write(*,*) "path_to_tracks", path_to_tracks
        case default
            print*, "METISSE error: reading inputs; unrecognized front_end_name"
            ierr = 1; return
        end select
        if (verbose) then
            write(*,*) "METALLICITY_DIR: ", METALLICITY_DIR
            write(*,*) "METALLICITY_DIR_HE: ", METALLICITY_DIR_HE
        endif
        
        !Some unit numbers are reserved: 5 is standard input, 6 is standard output.
        if (verbose) then
            ! write output to screen
            out_unit = 6
        else
            out_unit = alloc_iounit(ierr)
            open(out_unit,file='tracks_log.txt',action='write',status='unknown')
        endif
        
        if (write_error_to_file) then
            err_unit = 99   !will write to fort.99
        else
            err_unit = 6      !will write to screen
        endif
        
        ! use input file/path to locate list of *metallicity.in files
        ! these file contain information about eep tracks, their metallicity
        ! and the format file
        
        select case(front_end)
        case(COSMIC)
            ! hold for now since we read eeps directly
        case default    
            if (len(trim(METALLICITY_DIR))< 1) then
                write(*,*) "METISSE error: METALLICITY_DIR/path_to_tracks is an empty string"
                ierr = 1
                return
            else
                temp_filename = '.Zfilenames_H.txt'
                call get_metallicity_file_list(METALLICITY_DIR,metallicity_file_list,temp_filename)
                    
                if (.not. allocated(metallicity_file_list)) then
                    write(*,*) "METISSE error: metallicity file(s) not found in ", trim(METALLICITY_DIR)
                    write(*,*) "check if METALLICITY_DIR/path_to_tracks is correct"
                    ierr = 1
                    return
                else
                    if(debug) print*,'metallicity files: ',metallicity_file_list
                    call get_metallicity_list(metallicity_file_list,Z_H)
                endif
            endif
            
            if (len(trim(METALLICITY_DIR_HE))< 1) then
                write(out_unit,*) "Warning: METALLICITY_DIR_HE/path_to_he_tracks is an empty string"
                write(out_unit,*) "Switching to SSE formulae for helium stars "
                use_sse_NHe = .true.
            else
                temp_filename  = '.Zfilenames_He.txt'
                call get_metallicity_file_list(METALLICITY_DIR_HE,metallicity_file_list_he,temp_filename)
    
                if (.not. allocated(metallicity_file_list_he)) then
                    write(*,*) "METISSE error: metallicity file(s) not found in ", trim(METALLICITY_DIR_HE)
                    write(*,*) "check if METALLICITY_DIR_HE/path_to_he_tracks is correct"
                    ierr = 1
                    return
                else
                    if(debug) print*,'metallicity files he : ', metallicity_file_list_he
                    call get_metallicity_list(metallicity_file_list_he,Z_He)
                endif
            endif
        end select
    endif
    
    if (front_end > main) initial_Z = z
    write(out_unit,'(a,1p1e13.5)') ' Input Z is :', z

    if (use_sse_NHe)then
        ! only read hydrogen tracks
        nloop = 1
    else
        ! read both hydrogen and helium tracks
        nloop = 2
    endif

    ! need to intialize these seperately as they may be
    ! used uninitialized if he tracks are not present
    i_he_RCO = -1
    i_he_mcenv = -1
    i_he_Rcenv = -1
    i_he_age = -1
    do i = nloop,1, -1
        select case(front_end)
        case(COSMIC)
            ! we won't read in any metallicity files
            ! instead we will read in the format files and eeps 
            ! directly with COSMIC and pass them to METISSE

            if (i == 2) then
                ! Naked helium stars
                if (allocated(py_track_list_he)) then
                    ! First get info on the properties of the track_lis
                    if (allocated(track_list)) deallocate(track_list)
                    allocate(track_list(size(py_track_list_he)))
                    track_list = py_track_list_he
                    USE_DIR = METALLICITY_DIR_HE
                    num_tracks = size(track_list)
                    call apply_cosmic_format_controls('He')
                    call read_key_eeps_he()
                    if (debug) print*, "key eeps for he stars", key_eeps_he
                    call set_tracks_from_python_inputs(.true.)
                end if
            else
                ! Hydrogen rich stars
                if (allocated(py_track_list)) then
                    ! First get info on the properties of the track_list
                    if (allocated(track_list)) deallocate(track_list)
                    allocate(track_list(size(py_track_list)))
                    track_list = py_track_list
                    USE_DIR = METALLICITY_DIR
                    num_tracks = size(track_list)
                    call apply_cosmic_format_controls('H')
                    call read_key_eeps()
                    if (debug) print*, "key eeps", key_eeps   
                    call set_tracks_from_python_inputs(.false.)
                end if
            endif
            
            if (debug) print*, "num_tracks", num_tracks
            if (debug) print*, "tracks set by COSMIC"

        case default
            !read the files explicitly
            !read metallicity related variables
            if (i == 2) then
                write(out_unit,*) 'Reading naked helium star tracks'
                call get_metallcity_file_from_Z(metallicity_file_list_he,Z_He,initial_Z,ierr)
                if (ierr/=0) then
                    write(out_unit,'(a,1p1e13.5)')" No matching Z_files found with Z_accuracy_limit =",Z_accuracy_limit
                    write(out_unit,*)"Switching to SSE formulae for helium stars "
                    ierr = 0
                    use_sse_NHe = .true.
                    cycle
                endif
                write(out_unit,'(a,1p1e13.5)')" Found matching Z_files ",initial_Z
              
                USE_DIR = METALLICITY_DIR_HE
                temp_filename = '.Mfilenames_He.txt'
    
            else
                write(out_unit,*) 'Reading main (hydrogen star) tracks'
                call get_metallcity_file_from_Z(metallicity_file_list,Z_H,initial_Z,ierr)
                if (ierr/=0) then
                    write(*,'(a,1p1e13.5)')" No matching Z_files found with Z_accuracy_limit =",Z_accuracy_limit
                    write(*,*)"If needed, Z_accuracy_limit can be increased to match one of the available Z_files "
                    write(*,'(1p100e13.5)') pack(Z_H, mask = Z_H>0)
                    return
                endif
    
                write(out_unit,'(a,1p1e13.5)')" Found matching Z_files ",initial_Z
                USE_DIR = METALLICITY_DIR
                temp_filename = '.Mfilenames_H.txt'
            endif
            
            !read file-format
            call read_format(USE_DIR,format_file,ierr); if (ierr/=0) return
                
            !get filenames from eep_tracks_dir
            if (read_eep_files) file_extension = '.eep'
            call get_files_from_path(eep_tracks_dir,file_extension,temp_filename,track_list,ierr)
            
            if (ierr/=0) then
                eep_tracks_dir = trim(USE_DIR)//'/'//trim(eep_tracks_dir)
                call get_files_from_path(eep_tracks_dir,file_extension,temp_filename,track_list,ierr)
            endif
            
            if (ierr/=0 ) then
                print*,'METISSE error: failed to read input files.'
                print*,'Check if eep_tracks_dir is correct.'
                return
            endif
            
            num_tracks = size(track_list)
            
            write(out_unit,*)"Found ", num_tracks, " tracks."
            allocate(xa(num_tracks))
            xa% filename = track_list
            get_cols = .true.
            
            ! set_eeps
            if (i == 2) then
                xa% is_he_track = .true.
                call read_key_eeps_he()
                if (debug) print*, "key eeps for he stars", key_eeps_he
            else
                xa% is_he_track = .false.
                call read_key_eeps()
                if (debug) print*, "key eeps", key_eeps
            endif
            
            ! read input files
            if (read_eep_files) then
                if (debug) print*,"reading eep files"
                do j=1,num_tracks
                    call read_eep(xa(j))
                    if (code_error) return
!                    if(debug) write(*,'(a100,f8.2,99i8)') trim(xa(j)% filename), xa(j)% initial_mass, xa(j)% ncol
                end do
            else
                !read and store column names in temp_cols from the the file if header location is not provided
                if (header_location<=0) then
                    if (debug) print*,"Reading column names from file"
                    call process_columns(column_name_file,temp_cols,ierr)
                    
                    if(ierr/=0) then
                        print*,"Failed while trying to read column_name_file"
                        print*,"Check if header location and column_name_file are correct "
                        return
                    endif
    
                    if (size(temp_cols) /= total_cols) then
                        print*,'Number of columns in the column_name_file does not matches with the total_cols'
                        print*,'Check if column_name_file and total_cols are correct'
                        return
                    endif
                end if
    
                do j=1,num_tracks
                    call read_input_file(xa(j))
                    if (code_error) return
!                    if(debug) write(*,'(a100,f8.2,99i8)') trim(xa(j)% filename), xa(j)% initial_mass, xa(j)% ncol
                end do
            endif
        end select
            
        call check_tracks(num_tracks)

        ! Process the input tracks
        if (i==2) then
            !sort the array based on intial mass if not sorted already
            call sort_minitial()
            !reset z parameters where available
            !and determine cutoff masses
            call set_zparameters_he(num_tracks)
            call copy_and_deallocatex(num_tracks,sa_he)
            call get_minmax(sa_he(1)% is_he_track,Mmax_he_array,Mmin_he_array)
            if (allocated(core_cols_he)) deallocate(core_cols_he)

            allocate(core_cols_he(4))
            core_cols_he = -1
            core_cols_he(1) = i_he_age
            core_cols_he(2) = i_logL
            core_cols_he(3) = i_co_core
            if (i_he_RCO>0) core_cols_he(4) = i_he_RCO
        else
            !sort the array based on intial mass if not sorted already
            call sort_minitial()
            !reset z parameters where available
            !and determine cutoff masses
            call set_zparameters(num_tracks,zpars)
            call copy_and_deallocatex(num_tracks,sa)
            
            call get_minmax(sa(1)% is_he_track,Mmax_array,Mmin_array)
            if (allocated(core_cols)) deallocate(core_cols)

            allocate(core_cols(6))
            core_cols = -1
            
            core_cols(1) = i_age
            core_cols(2) = i_logL
            core_cols(3) = i_he_core
            core_cols(4) = i_co_core

            if (i_RHe_core>0) core_cols(5) = i_RHe_core
            if (i_RCO_core>0) core_cols(6) = i_RCO_core
        endif
        deallocate(track_list)
        
    end do
    
    
    ! for main, commons are assigned within the METISSE_main
    if (front_end > main) call assign_commons()
        
end subroutine METISSE_zcnsts

