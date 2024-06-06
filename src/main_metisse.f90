program metisse_main

    !This is the main program to use METISSE in standlaone mode
    !Evolves stars through evolv_metisse using inputs from SSE_input_controls namelist in evolve_metisse.in
    !For more details see Agrawal et al. 2020

    use track_support
    use z_support

    implicit none
    integer:: ierr,i
    real(dp):: zpars(20)
    real(dp), allocatable :: mass_array(:)

    !set the front end for METISSE
    call initialize_front_end('main')
    
    initial_Z = -1.d0
    
    ! read input metallicity and load the crresponding EEP tracks
    ! path for tracks are read from inlists
    call METISSE_zcnsts(initial_Z,zpars,'','',ierr)
    if (ierr/=0) STOP 'Fatal error: terminating METISSE'
    
    ! sets remnant schmeme from SSE_input_controls
    call assign_commons_main()
        
    allocate(mass_array(number_of_tracks))
    mass_array = 0.0

    if (read_mass_from_file) then
        !reads mass and age values from path for mass_file
        open(101, FILE= trim(input_mass_file), action="read",iostat =ierr)
        
        if (ierr/=0) then
            print*,'Erorr reading input masses from', trim(input_mass_file)
            print*,'check if input_mass_file is correct'
            STOP 'Fatal error: terminating METISSE'
        endif
        
        do i=1,number_of_tracks
            read(101,*) mass_array(i)
        end do
        close(101)
    else
        if (number_of_tracks>1) then
            call uniform_distribution(number_of_tracks,min_mass,max_mass,mass_array)
        else
            mass_array = min_mass
        endif
    endif
    
    allocate(t_incomplete(number_of_tracks), t_notfound(number_of_tracks))
    t_notfound = 0.d0
    t_incomplete = 0.d0
    
    
!    call allocate_track(number_of_tracks,mass_array)

    !evolve stars
    do i=1,number_of_tracks
        mass = mass_array(i)
        if (mass > Mcrit(9)% mass .or. mass< Mcrit(1)% mass) then
            t_notfound(i) = mass
            cycle
        endif
        if (verbose) write(*,'(a6, i9, a15,f7.3)') "count", i, "input mass = ", mass
        call allocate_track(1,mass)
        call evolv_metisse(mass,max_age,ierr,1)
        call dealloc_track()
        if (ierr/=0) t_incomplete(i) = mass
    end do
    
    if (verbose) print*,"Reached end of program"

    t_notfound = pack(t_notfound, mask = t_notfound >0)

    if (size(t_notfound)>0) then
        !write(*,'(a6, i9, a10,f7.3, a20)') "count", i, "ERROR: ", mass_array(i), "is out of bounds."
        write(*,*) "Stellar tracks of following initial masses were not interpolated."
        write(*,'(10f7.3)') t_notfound
        write(*,'(a25,f7.3, a5,f7.3)') "Reason: out of bounds for " , Mcrit(1)% mass," and ", Mcrit(9)% mass
    endif

    t_incomplete = pack(t_incomplete, mask = t_incomplete>0)
    if (size(t_incomplete)>0) then
        write(*,*) "Stellar tracks of following initial masses were rendered incomplete."
        write(*,'(5f7.3)') t_incomplete
    endif

end program

