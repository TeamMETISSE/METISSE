 subroutine METISSE_hrdiag(mass,aj,mt,tm,tn,&
                            tscls,&
                            lums,&
                            GB,&
                            zpars,&
                                r,lum,kw,mc,rc,menv,renv,k2,mcx,id)
    use track_support
    use interp_support, only: interpolate_age
    use sse_support
    use remnant_support

    implicit none
    integer, intent(in), optional :: id
    real(dp) :: mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
    real(dp) :: r,lum,mc,rc,menv,renv,k2,mcx
    
    integer :: kw,i,idd,j_bagb,old_phase
    real(dp) :: rg,rzams,rtms
    real(dp) :: Mcbagb, mc_max,HeI_time,dt_hold
    real(dp) :: bhspin ! only for cosmic
    type(star_parameters) :: old_pars

    logical :: debug, post_agb
    logical :: has_become_remnant
    type(track), pointer :: t

    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)
    
    debug = .false.
!    if ((id == 1).and. kw>=1 )debug = .true.
!if(id ==1 .and. t% is_he_track)debug = .true.
    if (debug) print*, '-----------HRDIAG-------------'
    if (debug) print*,"started hrdiag",mt,mc,aj,tn,kw,id,t% post_agb

    end_of_file = .false. !this is just the end of eep track
    has_become_remnant = .false.
    mc_max= 0.d0
        
    t% pars% mass = mt
    t% pars% phase = kw
    t% pars% core_mass = mc
    dt_hold = aj - t% pars% age
    if (aj/=aj) aj = t% pars% age
    t% pars% age = aj
!    if (aj<0.d0) stop
    
    
    IF (t% pars% phase<=TPAGB) THEN
        if (t% post_agb) then
            ! contruct an artficial track until WD cooling phase is reached
            call evolve_after_agb(t)
        else
            !check if phase/type/kw of the star has changed

            old_phase = t% pars% phase
            do i = 6,1,-1
                if (.not. defined(t% times(i))) cycle
                if (t% pars% age .lt. t% times(i)) then
                    t% pars% phase = i
                    if (debug)print*,"phase",t% pars% age,t% times(t% pars% phase),t% pars% phase
                endif
            enddo
            if (debug .and. t% pars% phase /= old_phase .and.t% pars% phase>0) &
                    print*,"phase change",t% pars% age,kw,t% pars% phase

            if (t% initial_mass<very_low_mass_limit .and. t% pars% phase==1) t% pars% phase =0

             !interpolate in age
            call interpolate_age(t,t% pars% age)
            if (debug)print*, "mt difference",t% pars% mass, mt, mt-t% pars% mass,kw
            if (front_end==main) then
                mt = t% pars% mass
            else
                t% pars% mass = mt
            endif
            IF (check_ge(t% pars% age,t% times(11))) THEN
                !check if have reached the end of the eep track
                if (debug) print*,"end of file:aj,tn ",t% pars% age,t% times(11),t% times(kw)
                if (kw<5 .and. (t% initial_mass.gt.very_low_mass_limit)) then
                if (t% ierr==0 .and. (dt_hold.le.t% pars% dt)) then
                    write(UNIT=err_unit,fmt=*) 'WARNING: Early end of file at phase, mass and id',&
                    kw,t% initial_mass,id
                    t% ierr = -1
!                    call stop_code
                endif
                endif
                end_of_file = .true.
                
                j_bagb = min(t% ntrack, TA_cHeB_EEP)
                Mcbagb = t% tr(i_he_core, j_bagb)
                mc_max = MAX(M_ch,0.773* Mcbagb-0.35)

                if (check_remnant_phase(t% pars, mc_max)) has_become_remnant = .true.
            
            ELSEIF (check_ge(t% pars% core_mass,t% pars% mass)) THEN
                !check if envelope has been lost
    
                if (debug) print*, "envelope lost at",t% pars% mass, t% pars% age,t% pars% phase
                
                if (t% pars% phase == TPAGB) then
                    ! TPAGB star becomes a CO-WD/ONe WD upon losing envelope
                    j_bagb = min(t% ntrack, TA_cHeB_EEP)
                    Mcbagb = t% tr(i_he_core, j_bagb)
                    has_become_remnant = .true.
                    !TODO: add a check if it's not a white dwarf
                else
                    call assign_stripped_star_phase(t, HeI_time)
                    if(t% pars% phase == HeWD) then
                        has_become_remnant = .true.
                    elseif (use_sse_NHe) then
                        t% star_type = sse_he_star
                        t% zams_mass = t% pars% mass

                        call initialize_SSE_helium_star(t,HeI_time)
                        call calculate_SSE_He_star(t,tscls,lums,GB,tm,tn)
                    else
                        t% is_he_track = .true.
                        t% star_type = switch
                        
                        if (t% pars% phase >= He_HG) then
                            j_bagb = min(t% ntrack, TA_cHeB_EEP)
                            mass = t% tr(i_mass, j_bagb)
                        else
                            mass = t% pars% mass
                        endif
                        ! zams_mass is assigned in the star
                        call METISSE_star(t% pars% phase, mass,t% pars% mass,tm,tn,tscls,lums,GB,zpars,0.d0,id)
                        t% pars% age_old = t% pars% age

                        if (debug) print*, 'after env loss', t% pars% phase, t% pars% age,t% MS_time,t% nuc_time,id
                    endif
                endif
            ELSE
                ! Calculate mass and radius of convective envelope, and envelope gyration radius.
                if (t% pars% core_radius<0) CALL calculate_rc(t,tscls,zpars,t% pars% core_radius)
                rc = t% pars% core_radius
                call get_mcrenv_from_cols(t,lums,menv,renv,k2)
            ENDIF
        endif
    ENDIF
    
    ! Stripped/ naked helium stars phase 7:9
    IF(t% pars% phase >= He_MS .and. t% pars% phase <=He_GB) THEN
        if (use_sse_NHe) then
            call evolve_after_envelope_loss(t,zpars(10))
            Mcbagb = t% zams_mass
            mc_max = max_core_mass_he(t% pars% mass, Mcbagb)
            
            if(t% pars% phase >He_MS) then
                if(check_remnant_phase(t% pars, mc_max)) has_become_remnant = .true.
            endif
            
            if (has_become_remnant .eqv..false.) then
                ! Calculate mass and radius of convective envelope, and envelope gyration radius.
                rzams = t% He_pars% Rzams
                CALL calculate_rc(t,tscls,zpars,rc)
                t% pars% core_radius = rc
                CALL calculate_rg(t,rg)
                CALL mrenv(t% pars% phase,t% zams_mass,t% pars% mass,t% pars% core_mass, &
                    t% pars% luminosity,t% pars% radius,rc,t% pars% age,t% MS_time,lums(2),lums(3),&
                    lums(4),rzams,rtms,rg,menv,renv,k2)
            endif
        else
            !check if phase/type/kw of the star has changed
            old_phase = t% pars% phase
            do i =9,7,-1
                if (.not. defined(t% times(i))) cycle
                if (t% pars% age .lt. t% times(i)) then
                    t% pars% phase = i
                    if (debug)print*,"phase",t% pars% age,t% times(t% pars% phase),t% pars% phase

                endif
            enddo
            if (debug .and. t% pars% phase /= old_phase) &
                    print*,"phase change",t% pars% age,t% times(t% pars% phase),t% pars% phase
            !interpolate in age
            call interpolate_age(t,t% pars% age)
            if(debug)print*,"mt difference",t% pars% mass,mt,mt-t% pars% mass,t% pars% phase
            t% pars% mass = mt
            if (check_ge(t% pars% age,t% times(11)) .or. (t% pars% core_mass.ge.t% pars% mass)) then
                !have reached the end of the eep track; self explanatory
                if (debug) print*,"end of file:aj,tn ",t% pars% age,t% tr(i_he_age,t% ntrack),t% times(kw)
                end_of_file = .true.
                
                j_bagb = min(t% ntrack, TAMS_HE_EEP)
                Mcbagb = t% tr(i_mass, j_bagb)
                mc_max = MAX(M_ch,0.773* Mcbagb-0.35)
                if (check_remnant_phase(t% pars, mc_max)) has_become_remnant = .true.
            else
                ! Calculate mass and radius of convective envelope, and envelope gyration radius.
                if (t% pars% core_radius<0) CALL calculate_rc(t,tscls,zpars,t% pars% core_radius)
                rc = t% pars% core_radius
                call get_mcrenv_from_cols(t,lums,menv,renv,k2)
            endif
        endif
    ENDIF
      
    kw = t% pars% phase
    ! remnants phases 10:15
    IF(has_become_remnant) THEN
!        print*, 'star',id,'is remnant',t% pars% mass,mcbagb,t% pars% core_mass
        t% star_type = remnant
        if (front_end == main .or. front_end == BSE) then
            if(t% pars% phase /= HeWD) then
                call assign_remnant_METISSE(t% pars, mcbagb)
                ! kw at this point contains old phase of the star,
                ! before the star became a remnant or lost its envelope
                call post_agb_parameters(t,kw)
            endif
        elseif (front_end == COSMIC) then
            ! storing mass that remnant would be in mt
            mt = t% pars% mass
            if(t% pars% phase /= HeWD) then
                call assign_remnant(zpars,t% pars% core_mass,&
                                mcbagb,t% zams_mass,mt,t% pars% phase,bhspin,id)
                t% pars% bhspin = bhspin
                !kw at this point contains old phase of the star
                call post_agb_parameters(t,kw)
                ! if t% post_agb is .false., assign correct remnant mass
            endif
            if(t% pars% phase >=10) t% pars% mass = mt
        endif
        has_become_remnant = .false.
        !has_become_remnant is only for assigning remnants, setting it to false now
    ENDIF

    IF(t% pars% phase >= HeWD) THEN
        if (front_end == main .or. front_end == BSE) then
            call evolve_remnants_METISSE(t% pars)
        elseif (front_end == COSMIC) then
            call hrdiag_remnant(zpars,t% pars% mass,t% pars% core_mass,t% pars% luminosity,&
                            t% pars% radius,t% pars% age,t% pars% phase)
        endif
        rc = t% pars% radius
        menv = 1.0d-10
        renv = 1.0d-10
        k2 = 0.21d0
    ENDIF
    
    kw = t% pars% phase
    aj = t% pars% age
    lum = t% pars% luminosity
    r =  t% pars% radius
    mc = t% pars% core_mass
    mcx = t% pars% McCO

    mt = t% pars% mass
    mass= t% zams_mass

    t% pars% mcenv = menv
    t% pars% rcenv = renv
!    t% pars% env_frac = (mt-t% pars% McHe)/mt

    if (debug) print*,"finished hrdiag",mt,mc,aj,kw,id,tm,tn
!    if (id==1)print*,"finished hrdiag",t% pars% mass, t% pars% core_mass,t% pars% age,t% pars% radius,id

    nullify(t)
    end subroutine METISSE_hrdiag
