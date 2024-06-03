subroutine evolv_metisse(mass,max_age,ierr,id)

    ! evolve subroutine to use metisse in standlaone mode
    ! evolves one star at a time and writes output to file

    use track_support
    use sse_support, only:time_He_MS

    !variable declaration
    real(dp), intent(in):: mass,max_age
    integer, intent(out):: ierr
    integer, intent(in), optional :: id


    !type(track), target :: t
    integer :: str,lines,old_phase,idd
    real(dp):: tphys,timestep,dt,dtr
    character(len=strlen) :: output_file,eep_filename
    logical :: output, t_end
    real(dp) :: dms, M_env,dml,x
    type(track), pointer :: t


    ! dummy variables for bse/ sse
    real(dp) :: mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
    real(dp) :: mc,rc,menv,renv,k2,mcx,r,lum,epoch,age
    integer :: kw
    
    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)
   
!    initialize variables

    ierr = 0
    timestep = 0.d0
    mt = mass
    kw = 1
    tscls = 0.d0
    lums = 0.d0
    GB = 0.d0
    tm = 0.d0
    tn = 0.d0
    str = int(mass*100)
    call METISSE_star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,timestep,id)

    if (t% complete .eqv..false.) ierr = 1

  !write the mass interpolated track if write_eep_file is true
    if (write_eep_file) then
        write(eep_filename,"(a,a,i5.5,a)") trim(METISSE_DIR),"/output/",str,"M_full.track.eep"
        call write_eep_track(t,mass,eep_filename)
    end if

    lines = 0
    tphys = 0.d0
    timestep = 0.d0
    old_phase = -1
    t_end  = .false.
    epoch = 0.d0
    
    !SSE like output file if write_track_to_file is true
    output = write_track_to_file
    if (output) write (output_file,"(a,a,i5.5,a)") trim(METISSE_DIR), "/output/evolve_", str, "M.dat"
    if (output) open (120,FILE=trim(output_file),action="write")

    if (output) write(120,'(9a15,2a10)') "time", "age", "mass","core_mass","He_core" &
                ,"CO_core","log_L","log_Teff","log_radius", "phase","e"

    do while(.true.)
        
!       advance the time
        tphys = tphys+timestep

        !evolve the star- calculate stellar parameters at tphys
        t% pars% age = tphys - epoch
        age = t% pars% age
        
        call METISSE_hrdiag(mass,age,mt,tm,tn,tscls,&
            lums,GB,zpars,r,lum,kw,mc,rc,menv,renv,k2,mcx,id)
        
         if (t% pars% phase>6 .or. t% post_agb) then
            t% pars% log_L = log10(t% pars% luminosity)
            t% pars% Teff = 1000*((1130.d0*t% pars% luminosity/(t% pars% radius**2))**0.25)
            t% pars% log_Teff = log10(t% pars% Teff)
            t% pars% log_R =  log10(t% pars% radius)
        endif
       
       if (tphys>=max_age) then
            tphys = max_age
            t_end  = .true. ! to print output before exiting
            !it is different to end_of_file defined in hrdiag
        end if
        
        !write output if flag is true
        if (old_phase /=t% pars% phase .or. t_end) then
            if (verbose) write(*,'(a10,f10.1,3a10,f7.3)') "Time ", tphys, "Phase ", phase_label(t% pars% phase+1), &
                                                                                "Mass ", t% pars% mass
            if (output) call write_dat_track(tphys, t% pars)
        else if (lines<3000 .and. output) then
            call write_dat_track(tphys,t% pars)
        endif
    
        if (t_end) exit
         
        lines = lines+1
        old_phase = t% pars% phase
       
        !calculate next time step
        call METISSE_deltat(id,age,dt,dtr)
        timestep = min(dt,dtr)
        
        ! only for SSE_he stars
        ! Calculate mass loss and modify timestep if need be
        if (t% pars% phase>=He_MS .and. t% pars% phase<=He_GB .and. use_sse_NHe) then
            !LBV-like mass loss beyond the Humphreys-Davidson limit
            x = 1.0d-5*r*SQRT(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0) dms = 1.5d0*1.0d-04
            ! Mass loss of Hamann & Koesterke (1998, A&A, 335, 1003) for WR (naked helium) stars.
            dml = 1.0d-13*lum**(3.d0/2.d0)
            ! Add metallicity factor from Vink & de Koter (2005, A&A, 442, 587).
            dml = dml*(initial_z/0.02d0)**0.86d0
            ! Or use mass loss of Nugis & Lamers (2000, A&A, 360, 227).
    !       dml = 1.0d-11*(lum**1.29d0)*(z**0.5d0)

            dms = MAX(dms,dml)
            dms = dms *1.0d+06*timestep
            M_env = mt - t% pars% core_mass
            if(dms.ge.M_env)then
                timestep = (M_env/dms)*timestep
                dms = M_env
            endif

            !Limit to 1% mass loss.
            if(dms.gt.0.01d0*mt)then
                timestep = 0.01d0*timestep*mt/dms
                dms = 0.01d0*mt
            endif
            mt = mt-dms

            if (t% pars% phase == He_MS) then
                t% zams_mass = mt
                epoch = tphys-(t% pars% age*time_He_MS(t% zams_mass)/t% ms_time)
            endif
            t% pars% mass = mt
        endif
        
    end do
    
    if (output) close(120)
   
    nullify(t)
    if (verbose) write(*,*) "-------------------------------------------------------------------------"
    return
end subroutine evolv_metisse



