program test_metisse

    !This is the main program to use METISSE in standlaone mode
    !Evolves stars through evolv_metisse using inputs from SSE_input_controls namelist in evolve_metisse.in
    !For more details see Agrawal et al. 2020

    use track_support
    use z_support
    use imf_support

    implicit none
    integer:: ierr,i, io, n
    real(dp):: zpars(20)
    real(dp), allocatable :: mass_array(:)

    !set the front end for METISSE
    call initialize_front_end('test')
    
    initial_Z = -1.d0
    
    ! read input metallicity and load the crresponding EEP tracks
    ! path for tracks are read from inlists
    call METISSE_zcnsts(initial_Z,zpars,ierr)
    if (ierr/=0 .or. code_error) STOP 1
    
    ! sets remnant schmeme from SSE_input_controls
    call assign_commons_main()
        
    ! number_of_tracks is read as real for user's ease, convert it to integer
    n = nint(number_of_tracks)
    allocate(mass_array(n))
    call sample_uniform(n,min_mass,max_mass,mass_array)
    

    !evolve stars
    do i = 1,n
        mass = mass_array(i)
        if (verbose) write(*,'(a6, i9, a15,f7.3)') "count", i, "input mass = ", mass
        call allocate_track(1,mass)
        call evolv_metisse(mass,max_age,ierr,1)
        call dealloc_track()
        if (ierr/=0) STOP 1
    end do
    
    if (verbose) print*,"Reached the end of the program"
    STOP 0
end program

