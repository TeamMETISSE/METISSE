subroutine METISSE_zcnsts(z,zpars,path_to_tracks,path_to_he_tracks,ierr)
    use track_support
    use z_support

    real(dp), intent(in) :: z
    character(len=*), intent(in) :: path_to_tracks, path_to_he_tracks
    real(dp), intent(out) :: zpars(20)
    integer, intent(out) :: ierr
    
    character(LEN=strlen), allocatable :: track_list(:)
    character(LEN=strlen) :: USE_DIR, find_cmd,rnd
    integer :: i,j,nloop,jerr
    integer :: num_tracks
    
    
    logical :: load_tracks
    logical :: debug, res
    
    debug = .false.
    ierr = 0
    
    ! At this point in the code front_end might not be assigned
    ! So we return ierr and let zcnsts.f of overlying code
    ! decide how to deal with errors.

    code_error = .false.
    
    if (front_end <0) then
        print*, 'Fatal error: front_end is not initialized for METISSE'
        ierr = 1; return
    endif
    
    
    if (debug) print*, 'in METISSE_zcsnts',z, trim(path_to_tracks),trim(path_to_he_tracks)
    
    ! read one set of stellar tracks (of input Z)
    load_tracks = .false.
    
    if (initial_Z <0) then
        load_tracks = .true.
    else
        ! tracks have been loaded at least once, for initial_Z
        ! check if they need to be reloaded
        
        ! if input metallicity 'z' has changed
        if (relative_diff(initial_Z,z) .ge. Z_accuracy_limit) load_tracks = .true.

        ! maybe metallicity is same, but paths may have changed
        ! (if paths are same, then metallicity doesn't matter)
        ! TODO: currently only for cosmic, can be modified to include others as well
        if (front_end == COSMIC)then
            if((trim(path_to_tracks)/=trim(TRACKS_DIR)) .or. &
            (trim(path_to_he_tracks)/=trim(TRACKS_DIR_HE))) load_tracks = .true.
        endif
            
    endif
            
    if (load_tracks.eqv. .false.) then
        if (debug) print*, '**** No change in metallicity or paths, exiting METISSE_zcnsts ****',initial_Z,z
        return
    endif
    
    if (debug) print*, 'Initializing METISSE_zcnsts'
    
    if (allocated(metallicity_file_list)) deallocate(metallicity_file_list,metallicity_file_list_he)

        
    ! use input file/path to locate list of *metallicity.in files
    ! these file contain information about eep tracks, their metallicity and format
    
    !read defaults option first
    call read_defaults()

    !read user inputs from evolve_metisse.in or use inputs from code directly
    if (front_end == main .or. front_end == bse) then
        call read_metisse_input(ierr); if (ierr/=0) call stop_code
    elseif (front_end == COSMIC) then
        TRACKS_DIR = path_to_tracks
        TRACKS_DIR_HE = path_to_he_tracks
        call get_metisse_input(TRACKS_DIR,metallicity_file_list)
        if (TRACKS_DIR_HE/='') call get_metisse_input(TRACKS_DIR_HE,metallicity_file_list_he)
    else
        print*, "Error: reading inputs; unrecognized front_end_name for METISSE"
        ierr = 1; return
    endif
    
        nloop = 2
    use_sse_NHe = .true.
    
    if (allocated(core_cols)) deallocate(core_cols)
    if (allocated(core_cols_he)) deallocate(core_cols_he)
    
    if (allocated(Mmax_array)) deallocate(Mmax_array, Mmin_array)
    
    !Some unit numbers are reserved: 5 is standard input, 6 is standard output.
    if (verbose) then
        out_unit = 6   !will write to screen
    else
        out_unit = 1000    !will write to file
        open(1000,file='tracks_log.txt',action='write',status='unknown')
    endif
    
    if (write_error_to_file) then
        err_unit = 99   !will write to fort.99
    else
        err_unit = 6      !will write to screen
    endif
    write(out_unit,'(a,f10.6)') ' Input Z is :', Z

    
    metallicity_file_list = pack(metallicity_file_list,mask=len_trim(metallicity_file_list)>0)
    
    metallicity_file_list_he = pack(metallicity_file_list_he,mask=len_trim(metallicity_file_list_he)>0)
    
    if (size(metallicity_file_list)<1) then
        print*, "Error: metallicity_file_list is empty"
        ierr = 1; return
    endif
    
    if (size(metallicity_file_list_he)<1) then
        write(out_unit,*) "Error: metallicity_file_list_he is empty"
        write(out_unit,*) "Switching to SSE formulae for helium stars "
        nloop = 1
    endif
    
    if(debug) print*,'metallicity files: ',metallicity_file_list
    if(debug) print*,'metallicity files he : ', metallicity_file_list_he
    
    if (front_end /= main) initial_Z = z

    ! need to intialize these seperately as they may be
    ! used uninitialized if he tracks are not present
    i_he_RCO = -1
    i_he_mcenv = -1
    i_he_Rcenv = -1
    i_he_age = -1
!    i_he_MoI = -1
        
    do i = nloop,1, -1
        !read metallicity related variables
        
        if (i == 2) then
            ZAMS_HE_EEP = -1
            TAMS_HE_EEP = -1
            GB_HE_EEP = -1
            TPAGB_HE_EEP = -1

            cCBurn_HE_EEP = -1
            post_AGB_HE_EEP = -1
            Initial_EEP_HE = -1
    
            Final_EEP_HE = -1
    
            write(out_unit,*) 'Reading naked helium star tracks'
            call get_metallcity_file_from_Z(initial_Z,metallicity_file_list_he,ierr)
            if (ierr/=0) then
                print*, "Switching to SSE formulae for helium stars "
                cycle
            endif
            USE_DIR = TRACKS_DIR_HE
        else
            write(out_unit,*) 'Reading main (hydrogen star) tracks'
            call get_metallcity_file_from_Z(initial_Z,metallicity_file_list,ierr)
            if (ierr/=0) return
            USE_DIR = TRACKS_DIR
        endif
        
        ! check if the format file exists
        inquire(file=trim(format_file), exist=res)
        
        if ((res .eqv. .False.) .and. (front_end == COSMIC)) then
            if (debug) print*, trim(format_file), 'not found; appending ',trim(USE_DIR)
            format_file = trim(USE_DIR)//'/'//trim(format_file)
        endif
    
        !read file-format
        call read_format(format_file,ierr); if (ierr/=0) return
            
        !get filenames from input_files_dir
        
        if (trim(INPUT_FILES_DIR) == '' )then
            print*,"Error: INPUT_FILES_DIR is not defined for Z= ", initial_Z
            ierr = 1; return
        endif
        
        if (front_end == COSMIC) then
!            find_cmd = 'find '//trim(INPUT_FILES_DIR)//'/*'//trim(file_extension)//' -maxdepth 1 > .file_name.txt'
!
!            call execute_command_line(find_cmd,exitstat=ierr,cmdstat=jerr)
!            if (ierr/=0) then
!            write(out_unit,*) trim(INPUT_FILES_DIR), 'not found; appending ',trim(USE_DIR), ierr, jerr
            INPUT_FILES_DIR = trim(USE_DIR)//'/'//trim(INPUT_FILES_DIR)
!            ierr = 0
!            endif
        endif
        
        write(out_unit,*)"Reading input files from: ", trim(INPUT_FILES_DIR)

        call get_files_from_path(INPUT_FILES_DIR,file_extension,track_list,ierr)
        
        if (ierr/=0) then
            print*,'Error: failed to read input files.'
            print*,'Check if INPUT_FILES_DIR is correct.'
            return
        endif

        num_tracks = size(track_list)
        
        write(out_unit,*)"Found: ", num_tracks, " tracks"
        allocate(xa(num_tracks))
        xa% filename = track_list
        set_cols = .true.
        
        if (i == 2) then
            xa% is_he_track = .true.
            call read_key_eeps_he()
            if (debug) print*, "key he eeps", key_eeps_he
        else
            xa% is_he_track = .false.
            call read_key_eeps()
            if (debug) print*, "key eeps", key_eeps
        endif
        
        if (read_eep_files) then
            if (debug) print*,"reading eep files"
            do j=1,num_tracks
                call read_eep(xa(j),ierr)
                if (ierr/=0) return
                if(debug) write(*,'(a100,f8.2,99i8)') trim(xa(j)% filename), xa(j)% initial_mass, xa(j)% ncol
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
                call read_input_file(xa(j),ierr)
                if (ierr/=0) return
                if(debug) write(*,'(a100,f8.2,99i8)') trim(xa(j)% filename), xa(j)% initial_mass, xa(j)% ncol
            end do
        endif
        
        call count_tracks(num_tracks)

        ! Processing the input tracks
        if (i==2) then
            call set_zparameters_he(num_tracks)
            call copy_and_deallocatex(num_tracks,sa_he)
            call get_minmax(sa_he(1)% is_he_track,Mmax_he_array,Mmin_he_array)

            use_sse_NHe = .false.
            allocate(core_cols_he(4))
            core_cols_he = -1
            core_cols_he(1) = i_he_age
            core_cols_he(2) = i_logL
            core_cols_he(3) = i_co_core
            if (i_he_RCO>0) core_cols_he(4) = i_he_RCO
        else
            !reset z parameters where available
            !and determine cutoff masses
            call set_zparameters(num_tracks,zpars)
            call copy_and_deallocatex(num_tracks,sa)
            
            call get_minmax(sa(1)% is_he_track,Mmax_array,Mmin_array)

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

    !TODO: 1. check BGB phase
    if (debug) print*,sa% initial_mass
    
    if (front_end /= main) call assign_commons()
    ! for main, commons are assigned within the METISSE_main 

        
        
end subroutine METISSE_zcnsts

