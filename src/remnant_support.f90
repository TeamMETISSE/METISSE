 module remnant_support
    use track_support
    use sse_support
    implicit none
    
    logical :: end_of_file, debug_rem

    !flags to be used while making decision for which method to use
    !for calculating properties of remnnants: neutron stars and black holes
    integer, parameter :: original_SSE = 0
    integer, parameter :: Belczynski2002 = 1
    integer, parameter :: Belczynski2008 = 2
    integer, parameter :: Eldridge_Tout2004 = 3
    integer, parameter :: Fryer2012 = 4

    !for white dwarfs
    integer, parameter :: Mestel = 0
    integer, parameter :: Modified_mestel = 1

    real(dp), parameter :: A_He = 4.d0
    real(dp), parameter :: A_CO = 16.d0
    real(dp), parameter :: tmin= 0.1d0

    !flags from SSE
    integer :: ns_flag = 0
    integer :: wd_flag = 0
    integer :: ec_flag = 0
    integer :: if_flag = 0
    contains
    

    logical function check_remnant_phase(pars,mc_max)
        type(star_parameters) :: pars
        real(dp) :: mc_max,mc_threshold

        debug_rem = .false.
        check_remnant_phase = .false.
        
        mc_threshold = pars% core_mass
        
        if (pars% phase <= TPAGB) then       !without envelope loss
            !mc = MAX(mc_max,mc_threshold)
            !mc_max = MIN(pars% mass,mc_max)
            mc_threshold = pars% McCO
        elseif (pars% phase >= He_MS) then
            mc_threshold = pars% core_mass
        else
            return
        endif

        if (mc_max<=0.d0 .or. mc_threshold<=0.d0) then
            write(UNIT=err_unit,fmt=*)"Fatal error: non-positive core mass",mc_max, mc_threshold
            code_error = .true.
            !assigning an ad-hoc non-zero core mass so the code doesn't break
            !TODO: fix very low-mass stars that form hewd and may get caught in this
            mc_threshold = 0.7*pars% McHe
            mc_max = mc_threshold
        endif
        
        if(mc_threshold>=mc_max .or. abs(mc_max-mc_threshold)<tiny .or. end_of_file)then
            !mc = MIN(mc_max,mc_threshold)
            pars% core_mass = mc_threshold
            pars% age_old = pars% age
            check_remnant_phase = .true.
            if (debug_rem) then
                print*, "check_remnant_phase is true"
                print*, "mass, core_mass, McCO, mc_max"
                print*, pars% mass, pars% core_mass, pars% McCO, mc_max
            end if
        endif
        
        end function check_remnant_phase
        
        
        subroutine assign_remnant_METISSE(pars, mcbagb)
            implicit none
            type(star_parameters) :: pars
            real(dp) :: Mcbagb

            !Mup_core, Mec_core are calculated in set_zparmeters routine of zfuncs
            if(pars% core_mass < M_ch)then
                if(mcbagb< Mup_core)then
                    pars% phase = CO_WD        !Zero-age Carbon/Oxygen White Dwarf
                else
                    if (ec_flag>0 .and. pars% McCO >= 1.372)then
                     !electron-capture collapse of an ONe core
                     !1.372<= mc< 1.44
                        pars% phase = NS
                        call initialize_ECSNe(pars)
                        if (debug_rem) print*,"ECSNe I: Mc< Mch, Mcbagb>Mup"
                    else
                        pars% phase = ONeWD    !Zero-age Oxygen/Neon White Dwarf
                    endif
                endif
            else !(mc>mch)
                !supernova 
                if(Mcbagb < Mup_core)then
                    ! Star is not massive enough to ignite C burning.
                    ! so no remnant is left after the SN
                    pars% phase = Massless_REM
                    call initialize_massless_rem(pars)

                else if(Mcbagb>= Mup_core .and. Mcbagb<= Mec_core)then
                    !Check for an electron-capture collapse of an ONe core.
                    if(ec_flag>0) then
                        pars% phase = NS
                        call initialize_ECSNe(pars)
                        if (debug_rem) print*,"ECSNe II: Mc> Mch, Mup< Mbagb< Mec"
                    else
                        call check_ns_bh(pars)
                    endif
                else
                    call check_ns_bh(pars)
                endif
            endif
        if (debug_rem .and. pars% phase>9) print*,"In remnant phase", phase_label(pars% phase+1)," , mass", pars% mass
        
    end subroutine assign_remnant_METISSE

    subroutine post_agb_parameters(t,old_phase)
        integer, intent(in) :: old_phase
        type(track), pointer :: t

        ! If a star becomes a WD, and construct_wd_track is true,
        ! then we use post_agb_parameters and evolve_after_agb
        ! to mimic post-agb evolution of star on HRD until WD cooling phase is reached.
        ! However naked helium stars don't go through this process
        ! and should directly jump to WD cooling track.
        ! So we use kw/old_phase as a check.
        
        !first check if the remnant is a WD
        if (t% pars% phase>HeWD .and. t% pars% phase<=ONeWD) then
            if (old_phase <=TPAGB .and. construct_wd_track) then
                ! contruct the track
                t% agb% phase_wd = t% pars% phase
                t% pars% phase = TPAGB
                t% post_agb = .true.
                t% pars% age_old = 0.0

                t% agb% tini = t% pars% age
                t% agb% lum = t% pars% luminosity
                t% agb% radius = t% pars% radius
                t% agb% mass = t% pars% mass
                call evolve_after_agb(t)
                if (debug_rem) print*, "In post-agb phase, mass = ", t% pars% mass
    !            print*,t% pars% luminosity, t% pars% radius
            else
                !jump to the WD cooling track
                t% zams_mass = t% pars% mass
                call initialize_white_dwarf(t% pars)
            endif
        endif
        
    end subroutine post_agb_parameters

    subroutine evolve_after_agb(t)
        type(track),pointer, intent(inout) :: t

        real(dp) :: alfa, beta,dt,r3,m0
        real(dp) :: radius_wd,lum_wd,mass_wd
        real(dp) :: t1,t2,t_post_agb
        
        t1 = fit_Z_t1(initial_Z)*t% MS_time
        t2 = 10*t1
        t_post_agb = t1+t2
        t% agb% tfinal = t% nuc_time + t_post_agb
        mass_wd = t% pars% core_mass
        radius_wd = calculate_wd_radius(mass_wd)
        lum_wd = calculate_wd_lum(mass_wd, 0.d0, A_CO)  !xx = a_co
        dt= t% pars% age- t% agb% tini
        m0 = t% pars% mass
!        print*,"aj",t% pars% age,t% agb% age, dt,t1
!            print*, t% nuc_time,t% agb% tfinal,dt
        alfa = 0d0; beta = 0d0
        r3 = 0.3*(t% agb% radius+radius_wd)
        t% pars% core_radius = radius_wd
        t% pars% core_mass = mass_wd
        if (dt<t1) then
            alfa = dt/t1
            beta = 1d0-alfa
            t% pars% radius = alfa* r3 + beta* t% agb% radius
            t% pars% luminosity = alfa* 0.9*t% agb% lum + beta*t% agb% lum
            t% pars% extra = 1
!            t% pars% mass = alfa* mass_wd + beta*t% agb% mass
!            t% pars% dms = (t% pars% mass-m0)/(dt*1.0d+06)
        else
            alfa = (dt-t1)/t2
            beta = 1d0-alfa
            t% pars% radius = (radius_wd*r3)/(alfa*(r3-radius_wd)+radius_wd)
            t% pars% luminosity = alfa* lum_wd + beta* 0.9*t% agb% lum
            t% pars% extra = 2
            t% pars% dms = 0.d0
        endif
        if (check_ge(t% pars% age,t% agb% tfinal)) then
            t% post_agb = .false.
            t% pars% extra = 0
            t% pars% age_old = t% pars% age
            t% pars% phase = t% agb% phase_wd
            t% zams_mass = t% pars% mass
!            call initialize_white_dwarf(t% pars)
        endif
!        print*, 'mass',t% pars% core_mass,t% pars% mass
    end subroutine

    subroutine initialize_white_dwarf(pars)
    type(star_parameters),intent(inout) :: pars

        if(if_flag>=1) call check_IFMR(pars% mass, pars% core_mass)
        pars% mass = pars% core_mass
        !pars% McHe = 0.0
        !pars% McCO = 0.0
        pars% age = 0.0

        call evolve_white_dwarf(pars)
        if (debug_rem) print*, "Remnant phase = ", phase_label(pars% phase+1), ", mass =", pars% mass
    end subroutine

    subroutine check_IFMR(mass, mc)
        real(dp), intent(in) :: mass
        real(dp) :: mc

        ! Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800.
          if(Z04>= 1E-08)then
             mc = MIN(0.36+0.104*mass,0.58+0.061*mass)
             mc = MAX(0.54+0.042*mass,mc)
             if(mass<1.0) mc= 0.46
          else
             mc= MIN(0.29+0.178*mass,0.65+0.062*mass)
             mc= MAX(0.54+0.073*mass,mc)
          endif
          mc= MIN(M_ch,mc)
    end subroutine

    subroutine check_ns_bh(pars)
        type(star_parameters) :: pars
        real(dp):: Mrem
        
        if (debug_rem) print*,"In CCSNe section",pars% core_mass,pars% mass
        pars% age = 0.0
        Mrem = calculate_NSBH_mass(pars% core_mass,pars% mass)
        if(debug_rem) print*, "Mrem= ", Mrem

        if(Mrem <= Max_NS_mass)then
            pars% phase = NS       !Zero-age Neutron star
            pars% mass = calculate_gravitational_mass(Mrem, pars% phase)
            call evolve_neutron_star(pars)
        else
            pars% phase = BH       !Zero-age Black hole
            pars% mass = calculate_gravitational_mass(Mrem, pars% phase)
            call evolve_black_hole(pars)

        endif
        if (debug_rem) print*, "NS/BH mass from", trim(BHNS_mass_scheme),"scheme = ",pars% mass
    end subroutine

    real(dp) function calculate_NSBH_mass(Mc,Mt) result(Mrem)
        real(dp), intent(in):: Mc,Mt
        real(dp) :: Mc_FeNi

        Mrem = Mc
        select case(ns_flag)
            case(Belczynski2002)
                !Use FeNi core mass given by Belczynski et al. 2002, ApJ, 572, 407.
                if(Mc< 2.5d0)then
                    Mc_FeNi = 0.161767d0*Mc+ 1.067055d0
                else
                    Mc_FeNi = 0.314154d0*Mc+ 0.686088d0
                endif
                Mrem = calculate_remnant_mass(Mc,Mc_FeNi,Mt)

            case(Belczynski2008)
                !Use FeNi core mass given by Belczynski et al. 2008, ApJSS, 174, 223.
                if(Mc< 4.82d0)then
                    Mc_FeNi = 1.5d0
                elseif(Mc< 6.31d0)then
                    Mc_FeNi = 2.11d0
                elseif(Mc< 6.75d0)then
                    Mc_FeNi = 0.69d0*Mc- 2.26d0
                else
                    Mc_FeNi = 0.37d0*Mc- 0.0828d0
                endif
                Mrem = calculate_remnant_mass(Mc,Mc_FeNi,Mt)

            case(Eldridge_Tout2004)
                !Use remnant masses based on Eldridge & Tout 2004, MNRAS, 353, 87.
                !NOT checked
                if(Mc<6.d0)then
                    Mc_FeNi = 1.44d0
                else
                    Mc_FeNi = (1.4512017d0*Mc)-(6.5913737d-03*Mc*Mc)-6.1073371d0
                endif
                Mrem = Mc_FeNi

            case(original_SSE)
                !Use the original SSE NS/BH mass.
                Mrem = 1.17d0 + 0.09d0*Mc
        end select
    end function

    real(dp) function calculate_remnant_mass(mc,mcfeni,Mt) result (Mrem)
        real(dp), intent(in) :: mc,mcfeni,mt
        real(dp) :: mc1 = 5.d0, mc2 = 7.6d0 !mass cutoffs for Belczynski methods

        !For Belczynski methods calculate the remnant mass from the FeNi core.
        if(mc<= mc1)then
            Mrem = mcfeni
        elseif(mc< mc2)then
            Mrem = mcfeni + (mc- mc1)*(mt - mcfeni)/(mc2-mc1)
        else
            Mrem = Mt
        endif
    end function calculate_remnant_mass

    real(dp) function calculate_gravitational_mass(mass,phase)
        real(dp) :: mass, mt
        integer,intent(in) :: phase
        mt = 0.0
        if(ns_flag>=2)then
            if (phase == NS) mt = quadratic(0.075d0, 1.d0, -mass)
            if (phase == BH) mt = 0.9d0* mass
        else
            mt = mass
        endif
        calculate_gravitational_mass = mt
    end function
    
    subroutine evolve_remnants_METISSE(pars)
    type(star_parameters):: pars
!        if (debug_rem)print*, 'evolving remnant in METISSE'
        select case(pars% phase)
        case(HeWD:ONeWD)
            call evolve_white_dwarf(pars)
        case(NS)
            call evolve_neutron_star(pars)
        case(BH)
            call evolve_black_hole(pars)
        case(Massless_Rem)
            call initialize_massless_rem(pars)
        end select
    end subroutine
            
    subroutine evolve_white_dwarf(pars)
    type(star_parameters):: pars
    real(dp) :: xx
    logical :: debug

    debug = .false.
    
    if (debug) print*,"In evolve_white_dwarf"
        pars% core_mass = pars% mass             !to account for accretion
        if(pars% phase == HeWD)then
            xx = a_he
        else
            xx = a_co
        endif

        if(pars% core_mass >= M_ch)then    !Mch or 1.372 !check for ec_flag
            if(pars% phase == ONeWD)then
                if (debug_rem) print*, "ECSNe:NS from accn onto WD", pars% mass
                pars% phase = NS
                call initialize_ECSNe(pars)
            else
                !Accretion induced supernova with no remnant
                if (debug_rem) print*, "ECSNe: massless remnant from accn onto WD", pars% mass
                pars% phase = Massless_REM
                call initialize_massless_rem(pars)
            endif
        else
            pars% luminosity = calculate_wd_lum(pars% mass, pars% age, xx)
            pars% radius = calculate_wd_radius(pars% mass)
            if(pars% mass < 0.0005) pars% radius= MIN(pars% radius,0.01d0)
            if (debug) print*, "Evolving WD", pars% mass, pars% luminosity, pars% radius
        endif
    end subroutine

    real(dp) function calculate_wd_lum(mass,age,xx) result(lum)
    real(dp), intent(in) :: mass, age, xx
    real(dp) :: fac
        lum = 0.d0
        select case(wd_flag)
        case(Mestel)            ! Mestel cooling
            Lum = 635.d0* mass* Z04/(xx*(age+tmin))**1.4
        case(Modified_mestel)       ! modified-Mestel cooling
            if(age < 9000.0)then
                Lum = 300* mass* Z04/(xx*(age+tmin))**1.18
            else
                fac = (9000.1*xx)**5.3
                Lum = 300*fac* mass* Z04/(xx*(age+tmin))**6.48
            endif
        end select

        !tmin= ((635.d0*pars% mass*Z04/lum)**(1.0/1.4))/xx
        !lum = (635.d0*mass*(Z04))/((xx*tmin)**1.4)
        !if (t% post_agb) lum = min(lum,t% agb% lum)
    end function

    real(dp) function calculate_wd_radius(mass) result(radius)
    real(dp), intent(in) :: mass
        radius = sqrt((M_ch/mass)**pow-(mass/M_ch)**pow)
        radius = max(1.4d-5,0.0115*radius)
        !radius  = 0.0115*SQRT(MAX(1.48204d-06,(M_ch/mass)**pow-(mass/M_ch)**pow))
        radius = MIN(0.1d0,radius )
    end function

    subroutine evolve_neutron_star(pars)
    type(star_parameters) :: pars
        !Neutron Star
        pars% core_mass = pars% mass
        !pars% McCO = 0.0
        !pars% McHe = 0.0
        if(pars% core_mass > Max_NS_mass)then
            pars% phase = BH  !Accretion induced Black Hole
            pars% age = 0.d0
            pars% luminosity = 1.0d-10
            pars% radius = 4.24d-06*pars% mass
        else
            pars% luminosity = 0.02*(pars% mass**0.67)/(MAX(pars% age,0.1d0))**2
            pars% radius= 1.4d-05
        endif
        
        !print*,"I am in evolve_neutron_star"
    end subroutine

    subroutine evolve_black_hole(pars)
    type(star_parameters) :: pars
        !Black hole
        !pars% mass has been calculated during remnant check
        pars% core_mass = pars% mass 
        pars% luminosity = 1.0d-10
        pars% radius = 4.24d-06*pars% mass
        !pars% McCO = 0.0
        !pars% McHe = 0.0

    end subroutine

    subroutine initialize_massless_rem(pars)
        type(star_parameters) :: pars
        pars% age = 0.d0
        pars% mass = 0.d0
        pars% luminosity = 1d-10
        pars% radius = 1d-10
        pars% Teff = 1d-10
        pars% core_mass = 0.d0
        pars% McCO = 0.d0
        pars% McHe = 0.d0
    end subroutine

    subroutine initialize_ECSNe(pars)
    type(star_parameters) :: pars
        !pars% age_old = pars% age
        pars% age = 0.d0
        pars% mass = 1.26d0
        pars% core_mass= pars% mass
        pars% luminosity = 0.02*(pars% mass**0.67)/(MAX(pars% age,0.1d0))**2
        pars% radius= 1.4d-05
    end subroutine
    
    subroutine assign_stripped_star_phase(t,HeI_time)
    
        type(track), pointer :: t
        real(dp) :: HeI_time
        logical :: debug

        debug = .false.
        
        HeI_time = 0.d0
        if (debug) print*,"Lost envelope at phase", t% pars% phase
        if (debug) print*,"age, core mass, mass", t% pars% age, t% pars% core_mass, t% pars% mass

        t% pars% age_old = t% pars% age

        select case(t% pars% phase)
            case(MS:RGB)   !MS,HG or RGB
            !ideally MS or Hg shouldn't directly jump to He_MS
            !There should be something like a He_PreMS
                if(t% zams_mass< Mhef)then
                    t% pars% phase = HeWD      !Zero-age helium white dwarf
                    t% pars% core_mass = t% pars% mass
!                    print*, 'hewd',t% pars% mass,Mhef
                else
                    t% pars% phase = He_MS       !Zero-age helium star
                    t% pars% core_mass = 0.d0
                    t% pars% McCO = 0.d0
                    t% pars% McHe = 0.d0
                    HeI_time = t% pars% age
                endif

            case(HeBurn)  !core he Burning
                t% pars% phase = He_MS
                t% pars% core_mass = 0.d0
                t% pars% McCO = 0.d0
                t% pars% McHe = 0.d0
                HeI_time = t% times(3)

            case(EAGB) !eAGB
                t% pars% phase = He_HG       !Evolved naked He star
                t% pars% mass = t% pars% core_mass
                t% pars% McHe = t% pars% mass
                t% pars% core_mass = t% pars% McCO
        end select
        
    end subroutine
    
    
    subroutine initialize_SSE_helium_star(t,HeI_time)
        type(track), pointer, intent(inout) :: t
        real(dp) :: HeI_time, HeB_time

        call calculate_SSE_He_timescales(t)
        
        if (t% pars% phase == He_MS) then
            HeB_time = t% times(4)-t% times(3)
            t% pars% age = t% MS_time*((t% pars% age- HeI_time)/HeB_time)
        else
            t% pars% age = He_GB_age(t% pars% core_mass,t% times(8), &
                            t% times(9),t% He_pars% D, t% He_pars% Mx)
            t% pars% age = MAX(t% pars% age,t% MS_time)
        endif
        
    end subroutine

    subroutine evolve_after_envelope_loss(t,McHeI)
    type(track),pointer, intent(inout) :: t

    real(dp) :: rg,tau,McHeI
    logical :: debug

    ! This is to prevent re-assigning of HeWD to HeMS phase
    ! if this function is called immediately after assign_stripped_star_phase
    if(t% pars% phase == HeWD) return

    debug = .false.
    if (debug) print*,"In evolve_after_envelope_loss: phase age",t% pars% phase, t% pars% age,t% ms_time,t% zams_mass
            
    t% He_pars% Rzams = radius_He_ZAMS(t% pars% mass)
    !if (check_le(t% pars% age ,t% MS_time)) then
    if(t% pars% age <t% MS_time .and. abs(t% pars% age-t% MS_time)>tiny)then
        ! Helium Main Sequence
        ! From SSE: "Star has no core mass and hence no memory of its past
        ! which is why we subject mass and mt to mass loss for this phase."

        t% pars% phase = He_MS
        tau = t% pars% age/t% MS_time
        t% He_pars% Lzams = lum_He_ZAMS(t% zams_mass)
        t% pars% luminosity = lum_He_MS(t% zams_mass, t% He_pars% Lzams,tau)
        t% pars% radius = radius_He_MS(t% pars% mass,t% He_pars% Rzams,tau)
        rg = t% He_pars% Rzams
        t% pars% core_mass = 0.0
        t% pars% age_old = t% pars% age
        !McHeI= He_McDu = core_mass_He_GB(lums(4), D,lx)   !core mass at He Ignition
!        print*, 'MCHEI',t% pars% mass,McHeI
        if(t% pars% mass < McHeI) t% pars% phase = HeWD
     else
        !Helium Shell Burning
        t% pars% phase = He_HG
        t% pars% luminosity = lum_He_GB(t% pars% age,t% times(8),t% times(9),&
                                t% times(10),t% He_pars% D)
        t% pars% radius = radius_He_HG(t% pars% mass,t% pars% luminosity,&
                                t% He_pars% Rzams,t% He_pars% LtMS)
        rg = radius_He_GB(t% pars% luminosity)  !radius on he giant branch
        if(t% pars% radius >= rg)then
           t% pars% phase = He_GB
           t% pars% radius = rg
        endif
        t% pars% core_mass = core_mass_He_GB(t% pars% luminosity,t% He_pars% Lx,t% He_pars% D)
        t% pars% McHe = t% pars% mass
        t% pars% McCO = t% pars% core_mass
    endif

    if (debug) print*,"End: Phase",t% pars% phase ," core mass",t% pars% core_mass
!        print*,"lum", t% pars% luminosity, "rad", t% pars% radius
    end subroutine evolve_after_envelope_loss


    !for constructing post-main sequence phase if construct_wd_track is true
    real(dp) function fit_Z_t1(Z)
    real(dp):: Z,x
        x = log10(Z)
        !y = 10**(-0.1*x*x-0.46*x-5.1)
        fit_Z_t1 = 10**(-0.136*x*x-1.117*x-6.256)
        return
    end function

    subroutine calculate_rc(t,tscls,zpars,rc)
     ! Calculate the core radius
     implicit none
     type(track), pointer, intent(in) :: t
     real(dp) :: tscls(20), zpars(20)
     real(dp) :: tau, lx,rx, rc,am,mt,mc,aj
     real(dp) :: tbagb,mass,lums1,lums2,tn
     real(dp) :: D, mx

        mass = t% zams_mass
        mt = t% pars% mass
        mc = t% pars% core_mass
        aj = t% pars% age
        tn = t% nuc_time
        tau = 0.d0
       select case(t% pars% phase)
         case(low_mass_MS:MS)
             rc = 0.d0
         case(HG:RGB)
             if(mass.gt.zpars(2))then
                 rc = radius_He_ZAMS(mc)
             else
                 rx = calculate_wd_radius(mc)
                 rc = 5.d0*rx
             endif
          case(HeBurn)
              tau = (aj - tscls(2))/tscls(3)     !tau = (aj - t_HeI)/t_He
              rx = radius_He_ZAMS(mc)
              !Following is akin to rzhef of SSE
              !TODO: following two lines may cause dt error
              ! might need interpolation
              am = MAX(0.d0,0.4d0-0.22d0*LOG10(mc))
              rc = rx*(1.d0+am*(tau-tau**6))
         case(EAGB)
    !            kwp = 9
             mc = t% pars% McHe
             tbagb = t% times(HeBurn)
             if(tn.gt.tbagb) tau = 3.d0*(aj-tbagb)/(tn-tbagb)
             D = 5.5d+04/(1.d0+0.4d0* t% zams_mass**4)
             Mx = (B/D)**(1.d0/(p-q))
             lx = lmcgbf(t% pars% McCO,D, Mx)
             lums1 = lum_He_ZAMS(mc)
             lums2 = lum_He_MS(mc,lums1,1.d0)
             if(tau.lt.1.d0) lx = lums2*(lx/lums2)**tau
             rx = radius_He_HG(mc,lx,radius_He_ZAMS(mc),lums2)
             rc = MIN(rx,radius_He_GB(lx))
         case(TPAGB: He_GB)
             if (t% pars% phase == He_MS) then
                 rc = 0.d0
             else
                 rx = calculate_wd_radius(mc)
                 rc = 5.d0*rx
             endif
         end select
        rc = MIN(rc,t% pars% radius)
    end subroutine calculate_rc

    subroutine calculate_rg(t,rg)
    !  rg = giant branch or Hayashi track radius, appropiate for the type.
    !       For kw=1 or 2 this is radius at BGB, and for kw=4 either GB or
    !       AGB radius at present luminosity.
    implicit none
    type(track), pointer, intent(in) :: t
    real(dp), intent(out):: rg
    real(dp) :: Rbgb, Rbagb, Lbgb, Lbagb, L, alfa
    integer :: j
    logical :: debug
    
        debug = .false.
        select case(t% pars% phase)
            case(low_mass_MS: HG)
                if (identified(BGB_EEP)) then
                    j = min(BGB_EEP,t% ntrack)
                    rg = t% tr(i_logR,j)   !TODO: check BGB  (or lums(3))
                else
                    j = min(cHeIgnition_EEP,t% ntrack)
                    rg = t% tr(i_logR,j)
                endif
                rg = 10.d0**rg
            case(RGB)
                rg = t% pars% radius
            case(HeBurn)
                !Linear interpolation between r(bgb) and r(bagb)
                !wrt luminosity l(bgb) and l(bagb) and L

                j = min(cHeIgnition_EEP,t% ntrack)
                Rbgb = t% tr(i_logR,j)
                Lbgb = t% tr(i_logL,j)
                j = min(TA_cHeB_EEP,t% ntrack)
                Rbagb = t% tr(i_logR,j)
                Lbagb = t% tr(i_logL,j)
                L = log10(t% pars% luminosity)
                alfa = (L - Lbgb)/(Lbagb-Lbgb)
                rg = (alfa*Rbagb)+((1d0 - alfa)*Rbgb)
                rg = 10**rg
                
                rg = max(rg, t% pars% radius )
!                if (rg .lt. t% pars% radius) then
!                    write(UNIT=err_unit,fmt=*)'Error in calculate_rg: Rg, R',rg,t% pars% radius,Rbgb, Rbagb
!                    write(UNIT=err_unit,fmt=*) t% pars% age, alfa, L, Lbgb, Lbagb
!                endif
                    
            case(EAGB:TPAGB)
                rg = t% pars% radius
            case(He_MS:He_GB)
                IF (use_sse_NHe) THEN
                    if (t% pars% phase==He_MS) then
                        rg = radius_He_ZAMS(t% pars% mass)
                    else
                        rg = radius_He_GB(t% pars% luminosity)
                    endif
                ELSE
                    if (t% pars% phase==He_GB) then
                        rg = t% pars% radius
                    elseif (t% initial_mass > Mcrit_he(4)% mass .and. &
                            t% initial_mass< Mcrit_he(5)% mass .and. &
                            identified(GB_HE_EEP)) then
                        j = min(GB_HE_EEP,t% ntrack)
                        rg = t% tr(i_logR,j)
                        rg = 10.d0**rg
                    else
                        rg = maxval(t% tr(i_logR,:))
                        rg = 10.d0**rg
                    endif
                    
                ENDIF
        end select
    end subroutine calculate_rg
    
    subroutine get_mcrenv_from_cols(t,lums,menv,renv,k2)
    
        type(track),pointer, intent(in) :: t
        real(dp), intent (in) :: lums(10)
        real(dp), intent(out) :: menv,renv,k2
        integer :: rcenv_col, mcenv_col, moi_col
        real(dp) :: rc, rg,rzams,rtms
    
        if (t% is_he_track) then
            mcenv_col = i_he_mcenv
            rcenv_col = i_he_rcenv
            moi_col = i_he_MoI

            rzams = 10.d0**t% tr(i_logR, ZAMS_HE_EEP)
            rtms = 10.d0**t% tr(i_logR, TAMS_HE_EEP)
        else
            mcenv_col = i_mcenv
            rcenv_col = i_rcenv
            moi_col = i_MoI
            rzams = 10.d0**t% tr(i_logR, ZAMS_EEP)
            rtms = 10.d0**t% tr(i_logR, TAMS_EEP)
        endif

        !rc, menv, renv and moi (moment of inertia) are calculated during age interpolation
        !if neccessary columns are present,
        !revert to SSE method if those columns are not present
        
        rc = t% pars% core_radius  ! it's calculated in hrdiag
        
        if ((.not. identified(mcenv_col)) .or. (.not. identified(rcenv_col)) .or. (.not. identified(moi_col)) .or. t% post_agb) then
            CALL calculate_rg(t,rg)
            CALL mrenv(t% pars% phase,t% zams_mass,t% pars% mass,t% pars% core_mass, &
            t% pars% luminosity,t% pars% radius,rc,t% pars% age,t% MS_time,lums(2),lums(3),&
            lums(4),rzams,rtms,rg,menv,renv,k2)
        endif
        
        if (mcenv_col>0) then
            !mass of convective envelope
            menv = t% pars% mcenv
            menv = min(menv,t% pars% mass-t% pars% core_mass)  ! limit it to the total envelope mass
            menv = MAX(menv,1.0d-10)
        endif
        
        if (rcenv_col>0) then
            renv = t% pars% rcenv
            renv = min(renv,t% pars% radius-rc)! limit it to the total envelope radius
        else
            if((t% pars% mass - t% pars% core_mass)>0) then
                renv = (t% pars% radius - rc)*menv/(t% pars% mass - t% pars% core_mass)
            else
                renv = 0.d0
            endif
        endif
        renv = MAX(renv,1.0d-10)
        ! radius of gyration, k2 given by sqrt(I/M*R*R)
!        if (moi_col>0) k2 = sqrt((t% pars% moi)/(t% pars% mass*t% pars% radius*t% pars% radius))
!        k2 = 0.21d0
    
    end subroutine


!    subroutine cutoffs_for_Belzynski_methods(ns_flag,mc1,mc2)
!        use track_support
!        use z_support, only : Mcrit
!        implicit none
!
!        integer, intent(in) :: ns_flag
!        real(dp), intent(out)  :: mc1, mc2
!        type(track) :: x
!        integer :: j_ntrack
!
!        if (ns_flag == Belczynski2002 .or. ns_flag == Belczynski2008) then
!            if(Mcrit(9)% mass>=20d0) then
!                call star(20d0,x)
!                j_ntrack = min(final_eep,x% ntrack)
!                mc1 = x% tr(i_co_core,j_ntrack)
!                call dealloc_track (x)
!            endif
!
!            if(Mcrit(9)% mass>=42d0) then
!                call star(42d0,x)
!                j_ntrack = min(final_eep,x% ntrack)
!                mc2 = x% tr(i_co_core,j_ntrack)
!                call dealloc_track (x)
!            endif
!        endif
!    end subroutine
end module remnant_support



