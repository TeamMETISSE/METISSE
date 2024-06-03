module z_support
    use track_support
    implicit none

    character(LEN=strlen) :: INPUT_FILES_DIR
    logical :: read_eep_files, read_all_columns, set_cols

    integer :: max_files = 50
    character(LEN=strlen) :: format_file, extra_columns_file
    character(LEN=strlen), allocatable :: metallicity_file_list(:),metallicity_file_list_he(:)
    character(LEN=col_width) :: extra_columns(100)

    real(dp) :: Z_files, Z_accuracy_limit

    !format_specifications
    character(LEN=5):: file_extension
    integer :: header_location, eep_location
    character(LEN=strlen) :: column_name_file
    type(column), allocatable :: key_cols(:), temp_cols(:)
    integer :: total_cols

    character:: extra_char

    type(eep_track), allocatable :: xa(:)

    real(dp) :: mass, max_age, min_mass, max_mass
    !used to set star_type_from_history
    ! central limits for high- / intermediate-mass stars, set these from input eep_controls nml
    real(dp) :: center_gamma_limit = 1d2
    real(dp) :: center_carbon_limit = 1d-4
    real(dp) :: log_center_T_limit = 9d0
    real(dp) :: high_mass_limit = 1d1 !Msun
    real(dp) :: he_core_mass_limit = 2.2

!    real(dp) ::  Mup_core,Mec_core


    namelist /SSE_input_controls/ initial_Z, max_age,read_mass_from_file,&
                        input_mass_file, number_of_tracks, max_mass, min_mass, &
                        WD_mass_scheme,use_initial_final_mass_relation, allow_electron_capture, &
                        BHNS_mass_scheme, max_NS_mass,pts_1, pts_2, pts_3, write_track_to_file

    namelist /METISSE_input_controls/ metallicity_file_list, Z_accuracy_limit, &
                        mass_accuracy_limit, construct_wd_track, verbose, &
                        write_eep_file,write_error_to_file, metallicity_file_list_he
            
    namelist /metallicity_controls/ INPUT_FILES_DIR, Z_files,format_file, extra_columns_file, &
                        read_all_columns, extra_columns,Mhook, Mhef, Mfgb, Mup, Mec, Mextra, Z_H, Z_He
                        
    namelist /format_controls/ file_extension, read_eep_files,total_cols,&
                        extra_char, header_location, column_name_file, &
                        PreMS_EEP, ZAMS_EEP, IAMS_EEP, TAMS_EEP, BGB_EEP, cHeIgnition_EEP, &
                        cHeBurn_EEP, TA_cHeB_EEP, TPAGB_EEP, cCBurn_EEP, post_AGB_EEP, &
                        Initial_EEP, Final_EEP, Extra_EEP1 ,Extra_EEP2, Extra_EEP3, &
                        fix_track, low_mass_final_eep, high_mass_final_eep, lookup_index, &
                        age_colname, mass_colname, log_L_colname ,log_T_colname, log_R_colname, &
                        Lum_colname, Teff_colname, Radius_colname, &
                        he_core_mass, co_core_mass, he_core_radius, co_core_radius, &
                        log_Tc, c12_mass_frac, o16_mass_frac,he4_mass_frac, &
                        mass_conv_envelope, radius_conv_envelope, moment_of_inertia, &
                        ZAMS_HE_EEP, TAMS_HE_EEP, GB_HE_EEP, cCBurn_HE_EEP, TPAGB_HE_EEP, &
                        post_AGB_HE_EEP, Initial_EEP_HE, Final_EEP_HE
            
    contains

    subroutine read_defaults()
        
        allocate(metallicity_file_list(max_files))
        allocate(metallicity_file_list_he(max_files))
        
        include 'defaults/evolve_metisse_defaults.inc'
    
    end subroutine read_defaults
    
    
    subroutine read_metisse_input(ierr)
        integer :: io
        integer, intent(out) :: ierr
        character(len=strlen) :: infile

        ierr = 0
        io = alloc_iounit(ierr)
        !reading user input
        
        infile = trim(METISSE_DIR)// '/evolve_metisse.in'
        open(io,FILE=infile,action="read",iostat=ierr)
            if (ierr /= 0) then
               print*, 'Error: Failed to open', trim(infile)
               call free_iounit(io)
               return
            end if
            if (front_end == main) then
                read(unit = io, nml = SSE_input_controls)
                if (.not. defined(initial_Z ))then
                    print*,"Error: initial_Z is not defined"
                    ierr = 1
                    return
                endif
            endif
            read(unit = io, nml = METISSE_input_controls)
        close(io)
        call free_iounit(io)
        
    end subroutine read_metisse_input
    
    
    subroutine get_metisse_input(path,file_list)

    character(LEN=strlen) :: path
    character(LEN=strlen), allocatable :: temp_list(:),file_list(:)

    integer :: ierr, n

        ierr = 0
        ! use inputs from COSMIC
        call get_files_from_path(path,'_metallicity.in',temp_list,ierr)
        
        if (.not. allocated(temp_list)) then
            print*, 'Could not find metallicity file(s) in ',trim(path)
            ierr = 1
            return
        endif
        
        n = size(temp_list)
        if (n > max_files) then
            print*, 'Too many metallicity files in ',trim(path)
            print*, 'Only using first ',max_files
            n = max_files
        endif
        file_list(1:n) = temp_list(1:n)
        deallocate(temp_list)
        
    end subroutine get_metisse_input
    
    subroutine get_metallcity_file_from_Z(Z_req,file_list,ierr)
        real(dp), intent(in) :: Z_req
        integer, intent(out) :: ierr
        character(LEN=strlen), allocatable :: file_list(:)

        integer :: i
        logical:: found_z,debug
        
        debug = .false.

        if (verbose) write(*,'(a,f7.3)') ' Input Z is', Z_req
        ierr = 0
        found_z = .false.
        
        do i = 1, size(file_list)
            if (len_trim(file_list(i))>0) then
                if (debug) print*, 'Reading : ', trim(file_list(i))
                call read_metallicity_file(file_list(i),ierr)
                if (ierr/=0) cycle
                if (.not. defined (Z_files)) then
                    print*, 'Warning: Z_files not defined in "'//trim(file_list(i))//'"'
                else
                    if (debug) print*, 'Z_files is', Z_files
                    if (relative_diff(Z_files,Z_req) < Z_accuracy_limit) then
                        if (debug) print*, 'Z_files matches with input Z'
                        found_z = .true.
                        exit
                    endif
                endif

            endif
        end do
        if ((found_z .eqv. .false.) .and. (ierr==0)) then
            print*, 'Error: metallicity value =', Z_req, 'not found amongst given Z_files'
            print*, 'Check metallicity_file_list and value of Z_files for each file'
            print*, 'If needed, Z_accuracy_limit can be relaxed (set to a greater value).'
            ierr = 1
            return
        endif
     
    end subroutine get_metallcity_file_from_Z
    
    subroutine read_metallicity_file(filename,ierr)
        character(LEN=strlen), intent(in) :: filename

        integer :: io
        integer, intent(out) :: ierr

        ierr = 0
        
        !initialize the defaults (even if already set)
        include 'defaults/metallicity_defaults.inc'
        
        io = alloc_iounit(ierr)
        open(io, file=filename, action='read', iostat=ierr)
            if (ierr /= 0) then
               print*, 'Error: failed to open metallicity_file: "'//trim(filename)//'"'
               call free_iounit(io)
               return
            end if
            read(io, nml = metallicity_controls)
        close(io)
        call free_iounit(io)
        
    end subroutine read_metallicity_file
    
    subroutine read_format(filename,ierr)
        character(LEN=strlen), intent(in) :: filename
        integer :: io
        integer, intent(out) :: ierr
        
        ierr = 0
        
        !initialize file format specs
        include 'defaults/format_defaults.inc'
        
        
        io = alloc_iounit(ierr)
        !read file format specs
        open(unit=io,file=trim(filename),action='read',iostat=ierr)
            if(ierr/=0)then
                print*,'Erorr: failed to open format file: "'//trim(filename)//'"'
                print*,'check if format file is correct'
                return
            endif
            read(unit = io, nml = format_controls)
        close(io)
        call free_iounit(io)

    end subroutine read_format

    
    subroutine get_files_from_path(path,extension,file_list,ierr)
    character(LEN=strlen), intent(in) :: path
    character(LEN=*), intent(in) :: extension
    character(LEN=strlen), allocatable :: file_list(:)

    integer, intent (out) ::  ierr

    character(LEN=strlen) :: str,find_cmd
    integer :: n,i, io
    
        ierr = 0
        
        find_cmd = 'find '//trim(path)//'/*'//trim(extension)//' -maxdepth 1 > .file_name.txt'
        call system(find_cmd,ierr)
        
        if (ierr/=0) return

        io = alloc_iounit(ierr)
        open(io,FILE='.file_name.txt',action="read")

        !count the number of tracks
        n = 0
        do while(.true.)
            read(io,*,iostat=ierr)
            if(ierr/=0) exit
            n = n+1
        end do

        allocate(file_list(n))
        rewind(io)
        ierr = 0
        do i = 1,n
            read(io,'(a)',iostat=ierr)str
            if (ierr/=0) exit
            file_list(i) = trim(str)
        end do
        
        close(io)
        call free_iounit(io)
        
    end subroutine get_files_from_path

    subroutine read_eep(x)      !from iso/iso_eep_support.f90
    type(eep_track), intent(inout) :: x
    integer :: ierr, io, j

    logical :: read_phase
    character(LEN=8) :: phase_info
    character(LEN=strlen) :: eepfile
    character(LEN=10) :: type_label
    
    logical :: debug

    ierr = 0
    debug = .false.
    read_phase = .false.

    eepfile = trim(x% filename)

    io = alloc_iounit(ierr)
    open(io,file=trim(eepfile),status='old',action='read',iostat=ierr)

    !check if the file was opened successfully; if not, then fail
    if(ierr/=0) then
       write(*,*) 'PROBLEM OPENING EEP FILE: ', trim(eepfile)
       close(io)
       call free_iounit(io)
       return
    endif

    read(io,'(25x,a8)') !x% version_string
    read(io,'(25x,i8)') !x% MESA_revision_number
    read(io,*) !comment line
    read(io,*) !comment line
    read(io,'(2x,f6.4,1p1e13.5,0p3f9.2)') x% initial_Y, x% initial_Z, x% Fe_div_H, x% alpha_div_Fe, x% v_div_vcrit
    read(io,*) !comment line
    read(io,*) !comment line
    read(io,'(2x,1p1e16.10,3i8,a8,2x,a10)') x% initial_mass, x% ntrack, x% neep, x% ncol, phase_info, type_label

    if (debug) print*,'reading',eepfile,x% neep
    
    if(index(phase_info,'YES')/=0) then
       x% has_phase = .true.
       allocate(x% phase(x% ntrack))
       x% ncol = x% ncol - 1
    else
       x% has_phase = .false.
       x% ncol = x% ncol
    endif
    
    allocate(x% tr(x% ncol, x% ntrack),x% cols(x% ncol))
    allocate(x% eep(x% neep))

    read(io,'(8x,299i8)') x% eep
    read(io,*) ! comment line
    read(io,*) ! column numbers

    !to exclude pms -read from eep(2)
    read(io,'(1x,299a32)') x% cols% name
    do j = 1, x% ntrack
        read(io,'(1x,299(1pes32.16e3))') x% tr(:,j)
    enddo

    !determine column of mass, age etc.
    if (set_cols) call locate_column_numbers(x% cols, x% ncol,x% is_he_track)
    x% initial_mass = max(x% initial_mass,maxval(x% tr(i_mass,:)))
!    call set_star_type_from_label(type_label,x)
    call set_star_type_from_history(x)

    close(io)
    call free_iounit(io)
    
  end subroutine read_eep
    
    !adapted from read_history_file of iso_eep_support
    subroutine read_input_file(x)
        type(eep_track), intent(inout) :: x
        character(LEN=8192) :: line
        integer :: i, io, j,ierr
        logical :: debug

        ierr = 0
        debug = .false.

        if (debug) print*,"in read_input_file",x% filename
        
        io = alloc_iounit(ierr)
        open(unit=io,file=trim(x% filename),status='old',action='read')
        !read lines of header as comments

        if (header_location >0)then
            do i = 1,header_location-1
                read(io,*) !header
            end do
            allocate(temp_cols(total_cols))
            !get column names
            read(io,'(a)') line
            do i =1, total_cols
                j = scan(line," ")
                temp_cols(i)% name = line(1:j)
                if (trim(temp_cols(i)% name)==extra_char) then
                    line = adjustl(line(j:))
                    j = scan(line," ")
                    temp_cols(i)% name = line(1:j)
                endif
                line = adjustl(line(j:))
            end do
        endif
        !figure out how many data lines
        j=0
        do while(.true.)
            read(io,*,iostat=ierr)
            if(ierr/=0) exit
            j=j+1
        enddo

        x% ntrack = j
        x% ncol = total_cols
        x% cols% name = temp_cols% name
        
        rewind(io)

        if (header_location >0)then
        !ignore file header, already read it once
            do i=1,header_location
               read(io,*) !header
            enddo
        endif

        allocate(x% tr(total_cols, x% ntrack))

        do j=1, x% ntrack
            read(io,'(a)') line
            call split(line, x% tr(:,j), x% ncol)
        enddo

        close(io)
        call free_iounit(io)

        !determine column of mass, age etc.
        if (set_cols) call locate_column_numbers(temp_cols, total_cols,x% is_he_track)
        
        if(header_location >0) deallocate(temp_cols)

        if (x% is_he_track) then
            x% neep = count(key_eeps_he .le. x% ntrack,1)
            allocate(x% eep(x% neep))
            x% eep = pack(key_eeps_he,mask = key_eeps_he .le. x% ntrack)
        else
            x% neep = count(key_eeps .le. x% ntrack,1)
            allocate(x% eep(x% neep))
            x% eep = pack(key_eeps,mask = key_eeps .le. x% ntrack)
        endif
        !print*,x% eep
        x% initial_mass = maxval(x% tr(i_mass,:))
        x% initial_mass = max(x% initial_mass,maxval(x% tr(i_mass,:)))
        call set_star_type_from_history(x)
        x% initial_Z = initial_Z
        x% initial_Y = Z_He

        if (debug) print*,x% initial_mass, x% initial_Z, x% ncol
    end subroutine read_input_file

    !from C.Flynn's driver routine

    subroutine split(line,values,ncol)
    character(LEN=*) :: line
    real(dp) :: values(:)
    integer:: i,ncol, iblankpos
    line = adjustl(line)
        do i =1, ncol
            !print*,i,trim(line)
            iblankpos = scan(line," ")
            if (trim(line)/= '') read(line(1:iblankpos),*) values(i)
            !print*, values(i)
            line = adjustl(line(iblankpos:))
        end do
    end subroutine split

    subroutine copy_and_deallocatex(y)
        type(eep_track), allocatable :: y(:)
        integer :: i, n, start
        logical :: debug

        debug = .false.
    
        if (allocated(y)) deallocate(y)
        allocate(y(size(xa)))
    
        !copy header
        y% filename = xa% filename
        y% initial_mass = xa% initial_mass
        y% ntrack = xa% ntrack
        y% initial_Y = xa% initial_Y
        y% initial_Z = xa% initial_Z
        y% Fe_div_H = xa% Fe_div_H
        y% alpha_div_Fe = xa% alpha_div_Fe
        y% v_div_vcrit = xa% v_div_vcrit
        y% star_type = xa% star_type
        y% is_he_track = xa% is_he_track
        
        !determine key columns
        call get_key_columns(xa(1)% cols, xa(1)% ncol, xa(1)% is_he_track)
        if (debug) print*,'key_cols', size(key_cols),xa(1)% ncol

        !copy columns and the track
        do n = 1,size(xa)
            if (read_all_columns) then
                if (debug) print*, 'using all columns'
                y(n)% ncol = xa(n)% ncol
                allocate(y(n)% tr(y(n)% ncol, y(n)% ntrack), y(n)% cols(y(n)% ncol))
                y(n)% cols% name = xa(n)% cols% name
                y(n)% tr = xa(n)% tr
            else
                y(n)% ncol = size(key_cols)
                allocate(y(n)% tr(y(n)% ncol, y(n)% ntrack), y(n)% cols(y(n)% ncol))
                do i = 1, y(n)% ncol
                    if (debug .and. n==1) print*, 'key column ',i,':',key_cols(i)% name,key_cols(i)% loc
                    y(n)% cols(i)% name = key_cols(i)% name
                    y(n)% tr(i,:) = xa(n)% tr(key_cols(i)% loc,:)
                end do
            endif
            ! copy eep information
            y(n)% neep = xa(n)% neep
            allocate(y(n)% eep(y(n)% neep))
            y(n)% eep = xa(n)% eep
            ! copy FSPS phase information if there
            y(n)% has_phase = xa(n)% has_phase
            if (y(n)% has_phase) then
                allocate(y(n)% phase(y(n)% ntrack))
                y(n)% phase = xa(n)% phase
            endif
            ! check for mass loss
            y(n)% has_mass_loss = check_mass_loss(y(n))
            
!             recalibrate age from ZAMS
            start = ZAMS_EEP
            if (y(n)% is_he_track)start = ZAMS_HE_EEP
            y(n)% tr(i_age2,:) = y(n)% tr(i_age2,:)- y(n)% tr(i_age2,start)
!
!            print*,'test input', y(n)% initial_mass,y(n)% ntrack
            if (.false.) then
            if (y(n)% ntrack< get_min_ntrack(y(n)% star_type, y(n)% is_he_track)) then
            print*, 'input track is incomplete', y(n)% initial_mass,y(n)% ntrack
            endif
            
            if (y(n)% star_type == star_high_mass .and. (y(n)% tr(i_co_core,y(n)% ntrack)<tiny)) then
            print*,'track has zero co_core',y(n)% initial_mass,y(n)% tr(i_co_core,y(n)% ntrack)
            endif
            endif
        end do
        
        !Now deallocate xa
        deallocate(xa)
        deallocate(key_cols)
    end subroutine

    subroutine get_minmax(x,Mmax,Mmin)
    
        type(eep_track) :: x(:)
        real(dp), allocatable :: Mmax(:), Mmin(:)
        integer :: i,j,nmax
        
        nmax = maxval(x% ntrack)
        allocate(Mmax(nmax), Mmin(nmax))
        Mmax = 0.d0
        Mmin = huge(0.0d0)    !largest float
        
        do i = 1,size(x)
            !Find maximum and minimum mass at each eep
            do j = 1, nmax
                if (x(i)% ntrack>=j) then
                    Mmax(j) = max(Mmax(j),x(i)% tr(i_mass,j))
                    Mmin(j) = min(Mmin(j),x(i)% tr(i_mass,j))
                endif
            end do
        end do
        
    end subroutine
    
    !locating essential columns here
    subroutine locate_column_numbers(cols,ncol,is_he_track)
        type(column), intent(in) :: cols(:)
        integer, intent(in) :: ncol
        logical, intent(in) :: is_he_track
        logical :: essential
        essential = .true.
        
        ! i_age is the extra age column for recording age values of new tracks
        ! it is used for interpolating in surface quantities after any explicit mass gain/loss
        ! It is the same as i_age2 if the input tracks already include mass loss due to winds/no mass loss
        
        i_age = ncol+1
        
        i_age2 = locate_column(cols, age_colname, essential)
        i_mass = locate_column(cols, mass_colname, essential)
        
        if (log_L_colname /= '') then
            !find the log luminosity column
            i_logL = locate_column(cols, log_L_colname, essential)

        else
            !find the luminosity column and convert it into log
            i_logL = locate_column(cols, Lum_colname, essential)
            call make_logcolumn(xa, i_logL)
        endif

        if (log_R_colname/= '') then
            i_logR = locate_column(cols, log_R_colname, essential)
        else
            i_logR = locate_column(cols, Radius_colname, essential)
            call make_logcolumn(xa, i_logR)
        endif

        !Teff gets calculated in the code, but the col used at several other places

        if (log_T_colname/= '') then
            i_logTe = locate_column(cols, log_T_colname, essential)
        else
            i_logTe = locate_column(cols, Teff_colname, essential)
            call make_logcolumn(xa, i_logTe)
        endif
            
        i_he_core = locate_column(cols, he_core_mass, essential)
        i_co_core = locate_column(cols, co_core_mass, essential)
        
        !optional columns
        if (is_he_track) then
            if (co_core_radius/= '') i_he_RCO = locate_column(cols, co_core_radius)
            if (mass_conv_envelope/= '') i_he_mcenv = locate_column(cols, mass_conv_envelope)
            if (radius_conv_envelope/= '') i_he_rcenv = locate_column(cols, radius_conv_envelope)
            if (moment_of_inertia/= '') i_he_MoI = locate_column(cols, moment_of_inertia)
            i_he_age = ncol+1
        else
            i_RHe_core = -1
            i_RCO_core = -1
            if (he_core_radius/= '') i_RHe_core = locate_column(cols, he_core_radius)
            if (co_core_radius/= '') i_RCO_core = locate_column(cols, co_core_radius)
            i_mcenv = -1
            if (mass_conv_envelope/= '') i_mcenv = locate_column(cols, mass_conv_envelope)
            i_Rcenv = -1
            if (radius_conv_envelope/= '') i_rcenv = locate_column(cols, radius_conv_envelope)
            i_MoI = -1
            if (moment_of_inertia/= '') i_MoI = locate_column(cols, moment_of_inertia)
        endif
        
        i_he4 = locate_column(cols, he4_mass_frac)
        i_c12 = locate_column(cols, c12_mass_frac)
        i_o16 = locate_column(cols, o16_mass_frac)
        i_Tc = locate_column(cols, log_Tc)
        
        set_cols = .false.
    end subroutine locate_column_numbers


    integer function locate_column(cols,colname,essential)
        character(LEN=col_width), intent(in) :: colname
        type(column), intent(in) :: cols(:)

        logical, intent(in),optional :: essential
        logical :: essential1

        integer :: i

        !unless explicitly specified
        !assume that the column is not essential
        essential1 = .false.
        if (present(essential)) essential1 = essential

        !now find the column
        locate_column = -1
        if (trim(colname)=='') return
        do i=1,size(cols)
           if(adjustl(adjustr(cols(i)% name))==trim(colname)) then
              locate_column = i
              return
           endif
        enddo

        !check whether the column has been successfully located
        if(locate_column<0) then
            write(0,*) 'Could not find column: ', trim(colname)
            !STOP the code if cannot locate one of the essential columns
            if(essential1) STOP
        endif
        
    end function locate_column
      
    subroutine make_logcolumn(x, itemp)
        type(eep_track) :: x(:)
        integer :: itemp,k
        do k = 1, size(x)
            x(k)% tr(itemp,:) = log10(x(k)% tr(itemp,:))
            x(k)% cols(itemp)% name = "log("//trim(x(k)% cols(itemp)% name)//")"
        end do
    end subroutine make_logcolumn
    
    subroutine make_pow10column(x, itemp,newname)
        type(eep_track) :: x(:)
        integer :: itemp,k
        character(LEN=col_width), intent(in), optional :: newname
        do k = 1, size(x)
            x(k)% tr(itemp,:) = 10.d0**(x(k)% tr(itemp,:))
            if (present(newname)) x(k)% cols(itemp)% name = trim(newname)
        end do
    end subroutine make_pow10column

    subroutine get_key_columns(cols,ncol,is_he_track)
        type(column), intent(in) :: cols(:)
        integer, intent(in) :: ncol
        logical, intent(in) :: is_he_track
        
        type(column) :: temp(ncol)
        type(column), allocatable :: temp_extra_columns(:)
        integer :: i,j,n,c,ierr
     
!        if (debug) print*, 'assigning key columns'

        ! Essential columns get reassigned here to match to reduced array format
        temp% loc = -1
        temp% name  =  ''

        n = 1
        
        call assign_sgl_col(temp, i_age2, age_colname,n)
        call assign_sgl_col(temp, i_mass, mass_colname,n)
        call assign_sgl_col(temp, i_logL, log_L_colname,n)
        call assign_sgl_col(temp, i_logR, log_R_colname,n)
        call assign_sgl_col(temp, i_he_core, he_core_mass,n)
        call assign_sgl_col(temp, i_co_core, co_core_mass,n)
        call assign_sgl_col(temp, i_logTe, log_T_colname,n)

        if(is_he_track) then
            if (i_he_RCO >0) call assign_sgl_col(temp, i_he_RCO, co_core_radius,n)
            if (i_he_mcenv>0) call assign_sgl_col(temp, i_he_mcenv, mass_conv_envelope,n)
            if (i_he_Rcenv>0) call assign_sgl_col(temp, i_he_Rcenv, radius_conv_envelope,n)
            if (i_he_MoI>0) call assign_sgl_col(temp, i_he_MoI, moment_of_inertia,n)
         else
        
            if (i_RHe_core >0) call assign_sgl_col(temp, i_RHe_core, he_core_radius,n)
            if (i_RCO_core >0) call assign_sgl_col(temp, i_RCO_core, co_core_radius,n)
            if (i_mcenv>0) call assign_sgl_col(temp, i_mcenv, mass_conv_envelope,n)
            if (i_Rcenv>0) call assign_sgl_col(temp, i_Rcenv, radius_conv_envelope,n)
            if (i_MoI>0) call assign_sgl_col(temp, i_MoI, moment_of_inertia,n)
                        
!            if (i_Tc >0) call assign_sgl_col(temp, i_Tc, log_Tc,n)
!            if (i_he4 >0) call assign_sgl_col(temp, i_he4, he4_mass_frac,n)
!            if (i_c12 >0) call assign_sgl_col(temp, i_c12,c12_mass_frac,n)
!            if (i_o16>0) call assign_sgl_col(temp, i_o16, o16_mass_frac,n)
        
        endif
        
        c = count(len_trim(extra_columns)>0)
        if (c>0) then
            do i = 1, size(extra_columns)
                if (len_trim(extra_columns(i))<1) exit
                j = locate_column(cols,extra_columns(i))
                if (j>0 )then
                    temp(n)% loc = j
                    temp(n)% name = extra_columns(i)
                    n = n+1
                endif
            end do
        else
            ierr = 0
            if (extra_columns_file /= '') call process_columns(extra_columns_file,temp_extra_columns,ierr)
            if(ierr/=0) print*, "Failed while trying to read extra_columns_file"
            
            if (allocated(temp_extra_columns)) then
                do i = 1, size(temp_extra_columns)
                    if (len_trim(temp_extra_columns(i)% name)<1) exit
                    temp(n)% loc = locate_column(cols,temp_extra_columns(i)% name)
                    temp(n)% name = temp_extra_columns(i)% name
                    n = n+1
                end do
                deallocate(temp_extra_columns)
            endif
        endif
    
!        if (allocated(key_cols)) deallocate(key_cols)
        allocate(key_cols(n-1))
        key_cols% name = temp(1:n-1)% name
        key_cols% loc = temp(1:n-1)% loc

        i_age = n
        if(is_he_track) i_he_age = n
        
    end subroutine get_key_columns

    subroutine assign_sgl_col(temp, col, colname,n)
        type(column) :: temp(:)
        character(len=col_width), intent(in) :: colname
        integer, intent(inout):: n,col
        
        temp(n)% loc = col
        temp(n)% name = colname
        col = n
        n = n+1
    end subroutine

    !reading column names from file - from iso_eep_support.f90
    subroutine process_columns(filename,cols,ierr)
        character(LEN=strlen), intent(in) :: filename
        integer, intent(out) :: ierr
        integer :: i, ncols(2), nchar, column_length, pass
        character(LEN=strlen) :: line, column_name
        logical :: is_int,debug
        type(column), allocatable, intent(out) :: cols(:)
        integer :: ncol,io

        debug =.false.
        ierr = 0
        io = alloc_iounit(ierr)
        open(io,file=trim(filename),action='read',status='old',iostat=ierr)
        if(ierr/=0) then
           write(*,*) 'failed to open columns list file: ', trim(filename)
           call free_iounit(io)
           return
        endif
        ncols=0
        do pass=1,2
           if(pass==2) allocate(cols(ncols(1)))
           inner_loop: do while(.true.)
              is_int = .false.
              read(io,'(a)',iostat=ierr) line
              if(ierr/=0) exit inner_loop

              !remove any nasty tabs
              do while(index(line,char(9))>0)
                 i=index(line,char(9))
                 line(i:i)=' '
              enddo

              nchar=len_trim(line)
              if(nchar==0) cycle inner_loop ! ignore blank line

              line=adjustl(line)
              i=index(line,'!')-1
              if(i<0) i=len_trim(line)

              if(i==0) then       !comment line
                 if(debug) write(*,*) ' comment: ', trim(line)
                 cycle inner_loop
              else if(index(line(1:i),'number')>0) then
                 if(debug) write(*,*) '****** ', trim(line)
                 if(debug) write(*,*) '****** index of number > 0 => integer'
                 is_int = .true.
              else if(index(line(1:i),'num_')==1)then
                 if(debug) write(*,*) '****** ', trim(line)
                 if(debug) write(*,*) '****** index of num_ == 1 => integer'
                 is_int = .true.
              endif

              column_name = line
              ncols(pass)=ncols(pass)+1
              if(i==0) then
                 column_length=len_trim(column_name)
              else
                 column_length=len_trim(column_name(1:i))
              endif
              do i=1,column_length
                 if(column_name(i:i)==' ') column_name(i:i)='_'
              enddo
              !if(debug) write(*,'(2i5,a32,i5)') pass, ncols(pass),trim(column_name(1:column_length)), column_length
              if(pass==2) then
                 cols(ncols(pass))% name = trim(column_name(1:column_length))
                 if(is_int) then
                    cols(ncols(pass))% type = column_int
                 else
                    cols(ncols(pass))% type = column_dbl
                 endif
                 cols(ncols(pass))% loc = ncols(pass)
              endif
           end do inner_loop
           if(pass==1) rewind(io)
           if(pass==2) close(io)
        end do
        if(ncols(1)==ncols(2)) then
           ierr=0
           ncol=ncols(1)
        endif
        call free_iounit(io)
        if(debug) write(*,*) 'process_columns: ncol = ', ncol

      end subroutine process_columns


    subroutine read_key_eeps()
    integer :: temp(15), neep,ieep
    
        temp = -1
        if (allocated(key_eeps)) deallocate(key_eeps)
        ieep = 1
        temp(ieep) = PreMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = ZAMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = IAMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = TAMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = BGB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = cHeIgnition_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = cHeBurn_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = TA_cHeB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = cCBurn_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = TPAGB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = post_AGB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = Extra_EEP1
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = Extra_EEP2
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = Extra_EEP3
        if(.not. add_eep(temp,ieep)) temp(ieep) = -1

        neep = count(temp > 0,1)
        allocate(key_eeps(neep))
        key_eeps = pack(temp,temp > 0)
        
        !define initial and final eep if not already defined
        if(Initial_EEP <0 .or. Initial_EEP< minval(key_eeps)) Initial_EEP = ZAMS_EEP
        if(Final_EEP < 0 .or. Final_EEP > maxval(key_eeps)) Final_EEP = maxval(key_eeps)
        
        if(low_mass_final_eep<0 .or. low_mass_final_eep>final_eep) low_mass_final_eep = final_eep
        if(high_mass_final_eep<0 .or. high_mass_final_eep>final_eep) high_mass_final_eep = final_eep
!        print*, 'eep main', Initial_EEP, final_eep, low_mass_final_eep, high_mass_final_eep
    end subroutine

    subroutine read_key_eeps_he()
    integer :: temp(15), neep,ieep
        temp = -1
        if (allocated(key_eeps_he)) deallocate(key_eeps_he)

        ieep = 1
        temp(ieep) = ZAMS_HE_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = TAMS_HE_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1
        
        temp(ieep) = GB_HE_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1
        
        temp(ieep) = TPAGB_HE_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1
        
        temp(ieep) = cCBurn_HE_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1
        
        temp(ieep) = post_AGB_HE_EEP
!        if(add_eep(temp,ieep)) ieep=ieep+1
        neep = count(temp > 0,1)
        allocate(key_eeps_he(neep))
        key_eeps_he = pack(temp,temp > 0)
        
        !define initial and final eep if not already defined
        if (Initial_EEP_HE < 0 .or. Initial_EEP_HE< minval(key_eeps_he)) Initial_EEP_HE = ZAMS_HE_EEP
        if (Final_EEP_HE < 0 .or. Final_EEP_HE > maxval(key_eeps_he)) Final_EEP_HE = maxval(key_eeps_he)
    
        if(low_mass_final_eep<0 .or. low_mass_final_eep>final_eep_he) low_mass_eep_he = Final_EEP_HE
        if(high_mass_final_eep<0 .or. low_mass_final_eep>final_eep_he) high_mass_eep_he = Final_EEP_HE
        
!        print*, 'eep he', Initial_EEP_he, final_eep_he, low_mass_eep_he, high_mass_eep_he
    end subroutine read_key_eeps_he
     
    logical function add_eep(temp, i)
    integer :: last, i, temp(:)
        last = 0
        add_eep = .false.
        if (i>1 ) last = temp(i-1)
        if (temp(i) >0 .and. temp(i)/= last) then
            add_eep = .true.
            if (Final_EEP>0 .and. temp(i)> Final_EEP) add_eep = .false.
        endif
    end function

    
  subroutine set_star_type_from_history(x)
    type(eep_track), intent(inout) :: x
    integer :: n

    !set the WDCS primary EEP if center_gamma < center_gamma_limit
    !center_gamma_limit = 19

    !set the CarbonBurn primary EEP if center_c12 < center_carbon_limit
    center_carbon_limit = 0.05

    !set star_type to high_mass_star if max(log_center_T) > this
    log_center_T_limit = 8.7 !changing from 8.5

    !set star_type to high mass star if M_init >= this
    high_mass_limit = 10.0 !Msun

    !from Pols et al. 1998- set star_type to high mass star if He core mass>= this
    he_core_mass_limit = 2.2 !Msun

    n = x% ntrack
    
    !- can be based on co mass in the end<1.4 maybe
    
    !only reach center_gamma_limit if the star evolves to a WD
    !if( x% tr(i_gamma,n) > center_gamma_limit) then
       !x% star_type = star_low_mass
       !return
    !endif

    x% star_type = star_low_mass

    !last gasp test for high-mass stars is the initial mass...
    if(x% is_he_track) then
        if(x% initial_mass >= he_core_mass_limit) then
           x% star_type = star_high_mass
           return
        endif
    else
        if(x% initial_mass >= high_mass_limit) then
           x% star_type = star_high_mass
           return
        endif
        
        !i_he_core
        if (x% tr(i_he_core,n)>= he_core_mass_limit) then
            x% star_type = star_high_mass
            return
        endif
        
    endif
    
    
    !simple test for high-mass stars is that central C is depleted
    if (i_c12 >0) then
        if(maxval(x% tr(i_c12,:)) > 0.4d0 .and. x% tr(i_c12,n) < center_carbon_limit)then
           x% star_type = star_high_mass
           return
        endif
    endif
    
    !alternative test for high-mass stars is that they reach a
    !central temperature threshhold
    if (i_Tc >0) then
        if(x% tr(i_Tc,n) > log_center_T_limit)then
            x% star_type = star_high_mass
            return
        endif
    endif

  end subroutine set_star_type_from_history

    logical function check_mass_loss(x)
    type(eep_track), intent(in) :: x
    real(dp) :: dm

    dm = (x% tr(i_mass,1)-x% tr(i_mass,x% ntrack))/x% initial_mass
    if (dm<tiny) then
        check_mass_loss = .false.
    else
        check_mass_loss = .true.
    endif
    end function

    subroutine set_star_type_from_label(label,x)
        character(LEN=10), intent(in) :: label
        type(eep_track), intent(inout) :: x
        integer :: n,i
        n = size(star_label)
        do i=1,n
            if(label==star_label(i)) x% star_type = i
        enddo
    end subroutine set_star_type_from_label

    subroutine set_zparameters(zpars)
        real(dp), intent(out) :: zpars(20)
        real(dp) :: old_co_frac,co_fraction,change_frac
        real(dp) :: smass,Teff,last_val,he_diff, mup_max
        real(dp), allocatable :: T_centre(:)
        integer :: len_track, i, min_index
        integer:: j_bagb, j_tagb, start
        real(dp), allocatable :: mass_list(:)
        integer :: num_tracks

        logical:: debug

        debug = .true.
        
        ! default is SSE
        Mup_core = 1.6d0
        Mec_core = 2.2d0
        !first calculate zpars the SSE way for use as backup
        call calculate_sse_zpars(initial_z,zpars)

        num_tracks = size(xa)
        Mcrit% mass= -1.d0
        Mcrit% loc = 0

        Mcrit(1)% mass = xa(1)% initial_mass
        Mcrit(1)% loc = 1

        Mcrit(2)% mass = very_low_mass_limit
        Mcrit(3)% mass = Mhook
        Mcrit(4)% mass = Mhef
        Mcrit(5)% mass = Mfgb
        Mcrit(6)% mass = Mup
        Mcrit(7)% mass = Mec
        Mcrit(8)% mass = Mextra

        Mcrit(9)% mass = xa(num_tracks)% initial_mass
        Mcrit(9)% loc = num_tracks+1 !TODO: explain why+1?

        if (verbose) write(*,'(a,f7.1)') ' Minimum initial mass', Mcrit(1)% mass
        if (verbose) write(*,'(a,f7.1)') ' Maximum initial mass', Mcrit(9)% mass

        allocate(mass_list(num_tracks),source=xa% initial_mass)
        old_co_frac = 0.d0

        !if already defined, do index search here otherwise search below
        do i = 2, size(Mcrit)-1
            if (.not. defined(Mcrit(i)% mass)) cycle
            call index_search (num_tracks, mass_list, Mcrit(i)% mass, min_index)
            Mcrit(i)% mass = xa(min_index)% initial_mass
            Mcrit(i)% loc = min_index
!            if (debug) print*, i, Mcrit(i)% mass
        end do

        start = max(Mcrit(1)% loc, Mcrit(2)% loc)
        
        do i = start, num_tracks
            smass = xa(i)% initial_mass
            !print*,smass, xa(i)% star_type
            len_track = xa(i)% ntrack
            
            if (smass<=3.0 .and. i_Tc>0) then
                !determining Mhook
                !where the maximum of the central temperature (Tc) between
                !the IAMS and the TAMS EEPs is greater than the Tc at
                !the TAMS EEP i.e., Tc,max>Tc,TAMS

                if ((Mcrit(3)% mass< 0.d0) .and. len_track >= TAMS_EEP) then
                    allocate(T_centre,source=xa(i)% tr(i_Tc,IAMS_EEP:TAMS_EEP))
                    last_val = T_centre(size(T_centre))
                    if (maxval(T_centre)>last_val) then
                        Mcrit(3)% mass = smass
                        Mcrit(3)% loc = i
                        if (debug) print*,"Mhook",smass,i
                    endif
                    deallocate(T_centre)
                endif

                !determining Mhef
                !the minimum temperature for core helium burning is about 100 million Kelvin.
                !In stars that undergo the helium flash, a slight expansion of the core
                ! following the flash causes the central temperature to decrease
                !a little before increasing again with stable helium burning.
                
                if ((Mcrit(4)% mass< 0.d0) .and. len_track>=TA_cHeB_EEP) then
                    allocate(T_centre,source=xa(i)% tr(i_Tc,cHeIgnition_EEP:TA_cHeB_EEP))
                    if (minval(T_centre)>7.4) then
                        Mcrit(4)% mass = smass
                        Mcrit(4)% loc = i
                        if (debug) print*,"Mhef",smass,i
                    endif
                    deallocate(T_centre)
                endif

            elseif(smass>=3.0) then
            
                !determining Mfgb
                if ((Mcrit(5)% mass< 0.d0) .and. len_track>=cHeIgnition_EEP) then
                    Teff = xa(i)% tr(i_logTe,cHeIgnition_EEP)       !temp at the end of HG/FGB
                    if (i_he4>0) then
                        he_diff = abs(xa(i)% tr(i_he4,cHeIgnition_EEP)-xa(i)% tr(i_he4,TAMS_EEP))
                    else
                        he_diff = 0.d0
                    endif
!                    print*,"bgb",smass,Teff, he_diff
                    if (Teff> T_bgb_limit .or. he_diff >0.01) then
                        Mcrit(5)% mass = smass
                        Mcrit(5)% loc = i
                        if (debug) print*,"Mfgb",smass,i
                    endif
                endif
            
                !determining Mup
                !where the absolute fractional change in the
                !central carbon-oxygen mass fraction exceeds
                !0.01 at the end of the AGB (TPAGB EEP)
                mup_max = 10.d0
                if (Mcrit(7)% mass> 0.d0) mup_max = Mcrit(7)% mass
                if ((Mcrit(6)% mass< 0.d0) .and. smass<mup_max) then
                    if (i_c12>0 .and. i_o16 >0) then
                        j_tagb = min(cCBurn_EEP,TPAGB_EEP)      !end of agb
                        j_tagb = min(len_track,j_tagb)
                        co_fraction = xa(i)% tr(i_c12,j_tagb)+xa(i)% tr(i_o16,j_tagb)
                        if (old_co_frac>0.0) then
                            change_frac = abs(co_fraction-old_co_frac)
                            change_frac = change_frac/old_co_frac
                            if (change_frac>0.01) then
                                ! this is the mass at which C/O ignition occur
                                ! we need the mass preceeding it
                                ! (Mup = M below which C/O ignition doesn't occur)
                                Mcrit(6)% loc = i-1
                                Mcrit(6)% mass = xa(Mcrit(6)% loc)% initial_mass
                                if (debug) print*,"Mup",Mcrit(6)% mass, Mcrit(6)% loc
                            endif
                        endif
                        old_co_frac = co_fraction
                    endif
                endif
            endif

            !determining Mec
            if (Mcrit(7)% mass< 0.d0)then
                if (xa(i)% star_type == star_high_mass) then
                    Mcrit(7)% mass = smass
                    Mcrit(7)% loc = i
                    if (debug) print*,"Mec",smass,i
                endif
            endif

        end do

        !If the tracks are beyond the zpars limits, above procedure
        !picks up the first or second track, which can lead to errors later,
        !hence those values need to be reverted

        do i = 2,size(Mcrit)-1
            if (Mcrit(i)% mass <= Mcrit(1)% mass) then
                Mcrit(i)% mass= -1.d0
                Mcrit(i)% loc = 0
            endif
        end do

        if (Mcrit(7)% loc < 1) then
            Mcrit(7)% mass = zpars(5)
        else
            j_bagb = min(xa(Mcrit(7)% loc)% ntrack,TA_cHeB_EEP)
            Mec_core = xa(Mcrit(7)% loc)% tr(i_he_core,j_bagb)
        endif
        
        !if cannot locate Mup or located it beyond Mec (which is incorrect),

        if (Mcrit(6)% loc < 1 .or. Mcrit(6)% loc >=  Mcrit(7)% loc) then
            if (debug) print*, 'Mcrit(6)/Mup not found'
            if (debug) print*, 'Using value closest to Mec-1.8 (the SSE way)'
            !modify Mup by SSE's way
            Mcrit(6)% mass = Mcrit(7)% mass - 1.8d0
            call index_search (num_tracks, mass_list, Mcrit(6)% mass, Mcrit(6)% loc)
            !make sure the new location for Mup does not exceed Mec
            if (Mcrit(7)% loc>1) Mcrit(6)% loc = min(Mcrit(6)% loc,Mcrit(7)% loc-1)
        endif

        if (Mcrit(6)% loc >= 1 .and. Mcrit(6)% loc <=num_tracks) then
            Mcrit(6)% mass = xa(Mcrit(6)% loc)% initial_mass
            if (debug) print*,"new Mup",Mcrit(6)% mass, Mcrit(6)% loc
            j_bagb = min(xa(Mcrit(6)% loc)% ntrack,TA_cHeB_EEP)
            Mup_core = xa(Mcrit(6)% loc)% tr(i_he_core,j_bagb)
        endif
        
        if (debug) print*,"Mup_core =", Mup_core
        if (debug) print*,"Mec_core =", Mec_core

        allocate (m_cutoff(size(Mcrit)))
        m_cutoff = Mcrit% loc
        call sort_mcutoff(m_cutoff)
        if (debug) print*, "m_cutoffs: ", m_cutoff
    
        !now redefine zpars where applicable
        do i = 3,7
            if(defined(Mcrit(i)% mass)) then
                zpars(i-2) = Mcrit(i)% mass
            else
                Mcrit(i)% mass = zpars(i-2)
            endif
        end do
    
        if (defined(Z_He)) then
            zpars(12) = Z_He
        elseif (xa(1)% initial_Y >0.d0)then
            zpars(12) = xa(1)% initial_Y
        endif

        if (defined(Z_H)) then
            zpars(11) = Z_H
        else
            zpars(11) = 1-initial_Z-Z_He
        endif
        
        Z04 = initial_Z**0.4

        
        !Redefine these for use later in the code
        Mhook = zpars(1)
        Mhef = zpars(2)
        Mfgb = zpars(3)
        Mup = zpars(4)
        Mec = zpars(5)
    
    
        if (debug) print*, 'zpars:',  zpars(1:5)
        deallocate(mass_list)
    end subroutine set_zparameters

    subroutine set_zparameters_he()
        real(dp) :: smass,frac_mcenv
        integer :: len_track, i, min_index, start!,j,jstart,jend
        real(dp), allocatable :: mass_list(:)
        integer :: num_tracks

        logical:: debug

        debug = .false.

        Mcrit_he% mass= -1.d0
        Mcrit_he% loc = 0

        Mcrit_he(1)% mass = xa(1)% initial_mass
        Mcrit_he(1)% loc = 1

!        Mcrit_he(2)% mass =!very_low_mass_limit
        !keep it empty for now

        Mcrit_he(3)% mass = Mhook
!        Mcrit_he(4)% mass = Mhef
        Mcrit_he(5)% mass = Mfgb
        Mcrit_he(6)% mass = Mup
        Mcrit_he(7)% mass = Mec
        Mcrit_he(8)% mass = Mextra
        
        num_tracks = size(xa)

        Mcrit_he(9)% mass = xa(num_tracks)% initial_mass
        Mcrit_he(9)% loc = num_tracks+1

        if (verbose) write(*,'(a,f7.1)') ' Minimum initial mass', Mcrit_he(1)% mass
        if (verbose) write(*,'(a,f7.1)') ' Maximum initial mass', Mcrit_he(9)% mass

        allocate(mass_list(num_tracks))
        mass_list = xa% initial_mass

        !if already defined, do index search here otherwise search below
        do i = 2, size(Mcrit_he)-1
            if (.not. defined(Mcrit_he(i)% mass)) cycle
            call index_search (num_tracks, mass_list, Mcrit_he(i)% mass, min_index)
            Mcrit_he(i)% mass = xa(min_index)% initial_mass
            Mcrit_he(i)% loc = min_index
!            if (debug) print*, i, Mcrit_he(i)% mass
        end do

        start = max(Mcrit_he(1)% loc, Mcrit_he(2)% loc)
        
        do i = start, num_tracks
            smass = xa(i)% initial_mass
            !print*,smass, xa(i)% star_type
            len_track = xa(i)% ntrack
             
            IF (.not. defined(Mcrit_he(3)% mass) .and. i_logTe>0)THEN
                if (len_track>= TPAGB_HE_EEP) then
                    !temp at the end of HG/FGB <= temp at TAMS
                    if (xa(i)% tr(i_logTe,TPAGB_HE_EEP) .le. xa(i)% tr(i_logTe,TAMS_HE_EEP)) then
                        Mcrit_he(3)% mass = smass
                        Mcrit_he(3)% loc = i
                        if (debug) print*,"Mhook",smass,i
                    endif
                endif
            ENDIF
            
            IF (.not. defined(Mcrit_he(5)% mass))THEN
                if (i_he_mcenv>0 .and. GB_HE_EEP>0 .and. len_track>=GB_HE_EEP) then
                    frac_mcenv = xa(i)% tr(i_he_mcenv,GB_HE_EEP)/xa(i)% tr(i_mass,GB_HE_EEP)
                    if (.not. defined(Mcrit_he(4)% mass) .and. frac_mcenv.ge.0.12d0) then
                        Mcrit_he(4)% mass = smass
                        if (debug) print*,"Mfgb1",smass,i

                    elseif(defined(Mcrit_he(4)% mass) .and. frac_mcenv.lt.0.12d0) then
                        Mcrit_he(5)% mass = smass
                        Mcrit_he(5)% loc = i
                        if (debug) print*,"Mfgb2",smass,i
                    endif
!                elseif (len_track>= Final_EEP_HE)) then
!                    jstart = TAMS_HE_EEP
!                    jend = Final_EEP_HE
!                    do j = jstart,jend
!                        if (xa(i)% tr(i_mcenv,j)/xa(i)% tr(i_mass,j).ge.0.12d0) then
!                            Mcrit_he(5)% mass = smass
!                            Mcrit_he(5)% loc = i
!                            if (debug) print*,"Mfgb",smass,i
!                            exit
!                        endif
!                    enddo
                endif
            ENDIF

            !determining Mec
            if (.not. defined(Mcrit_he(7)% mass))then
                if (xa(i)% star_type == star_high_mass) then
                Mcrit_he(7)% mass = smass
                Mcrit_he(7)% loc = i
                if (debug) print*,"Mec",smass,i
                endif
            endif
        end do

        !If the tracks are beyond the zpars limits, above procedure
        !picks up the first or second track, which can lead to errors later,
        !hence those values need to be reverted

        do i = 2,size(Mcrit_he)-1
            if (Mcrit_he(i)% mass <= Mcrit_he(1)% mass) then
                Mcrit_he(i)% mass= -1.d0
                Mcrit_he(i)% loc = 0
            endif
        end do

        !Mec
        Mcrit_he(7)% loc = max(Mcrit_he(7)% loc,1)

        allocate (m_cutoff_he(size(Mcrit_he)))
        m_cutoff_he = Mcrit_he% loc
        call sort_mcutoff(m_cutoff_he)
        if (debug) print*, "m_cutoffs he : ", m_cutoff_he

        deallocate(mass_list)
    end subroutine set_zparameters_he
    
    subroutine sort_mcutoff(m_cutoff)
     !subroutine to sort Mcutoffs, removing the ones who are at less than 2 distance from the last one
    !making sure there are at least 2 tracks between subsequent mcutoff
    
        integer, allocatable :: m_cutoff(:)
        integer, allocatable :: mloc(:)
        integer :: val,n
        integer :: i, k, loc,a
        
        n = size(m_cutoff)
        allocate (mloc(n))
        mloc = m_cutoff
        m_cutoff = 0
        m_cutoff(1) = 1
        k=1
        
        do i = 1, n
            val = minval(mloc(i:n))
            a = minloc(mloc(i:n),dim=1)
            loc = (i - 1) + a
            mloc(loc) = mloc(i)
            mloc(i) = val
            if ((val - m_cutoff(k))>=2) then
                k=k+1
                m_cutoff(k) = val
            endif
        end do
        
        m_cutoff = pack(m_cutoff,mask = m_cutoff .ne. 0)
        deallocate(mloc)
    end subroutine sort_mcutoff

    subroutine sort(mloc, list)
        integer, intent(in) :: mloc(:)
        integer, allocatable, intent(out) :: list(:)
        integer :: val,n
        integer :: i, a, loc

        n = size(mloc)
        allocate(list(n))
        list=pack(mloc,mask = mloc .ne. 0)
        
        !sort array
        do i = 0, n-1
            val = minval(list(i:n-1))
            a = minloc(list(i:n-1),dim=1)
            loc = (i - 1) + a
            list(loc) = list(i)
            list(i) = val
        end do
    end subroutine sort

    !ZPARS
    !finds critical masses and their locations
    !1; M below which hook doesn't appear on MS, Mhook. 3
    !2; M above which He ignition occurs non-degenerately, Mhef. 4
    !3; M above which He ignition occurs on the HG, Mfgb. 5
    !4; M below which C/O ignition doesn't occur, Mup. 6
    !5; M above which C ignites in the centre, Mec. 7

    subroutine calculate_sse_zpars(z,zpars)

    real(dp),intent(in) :: z
    real(dp),intent(out) :: zpars(14)
    real(dp) :: lzs,dlzs,lz,lzd

        lzs = log10(z/0.02d0)
        dlzs = 1.d0/(z*log(10.d0))
        lz = log10(z)
        lzd = lzs + 1.d0

        zpars(1) = 1.0185d0 + lzs*(0.16015d0 + lzs*0.0892d0)
        zpars(2) = 1.995d0 + lzs*(0.25d0 + lzs*0.087d0)
        zpars(3) = 16.5d0*z**0.06d0/(1.d0 + (1.0d-04/z)**1.27d0)
        zpars(4) = MAX(6.11044d0 + 1.02167d0*lzs, 5.d0)
        zpars(5) = zpars(4) + 1.8d0
        zpars(6) = 5.37d0 + lzs*0.135d0
        !* set the hydrogen and helium abundances
        zpars(11) = 0.76d0 - 3.d0*z
        zpars(12) = 0.24d0 + 2.d0*z
        !* set constant for low-mass CHeB stars
        zpars(14) = z**0.4d0

    end subroutine
    
    elemental function relative_diff(z1,z2) result(y)
        real(dp),intent(in) :: z1,z2
        real(dp) :: y
        !the formula can catch the difference in small values of metallicity
        !while allowing for differences due to precision errors
        y = abs(z1-z2)/MIN(z1,z2)
    end function

end module z_support
