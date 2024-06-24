    subroutine METISSE_deltat(id,age,dt,dtr)
    
    !calculates timestep for evolution
    use track_support
    integer, intent(in), optional :: id

    INTEGER :: idd!,kw
    REAL(dp) :: age!,tm,tn,tscls(20)
    REAL(dp) :: pts1,pts2,pts3
    COMMON /POINTS/ pts1,pts2,pts3
    
    real(dp) :: dtr,dt, tMS_hook
    type(track), pointer :: t
    
        idd = 1
        if(present(id)) idd = id
        t => tarr(idd)
    
        dtr = -1.d0
        dt = -1.d0
        !Base new time scale for changes in radius & mass on stellar type.
        select case (t% pars% phase)
            case(low_mass_MS:MS)
                tMS_hook = (0.95d0* t% times(1))
                if (age .le. tMS_hook) then
                    dt = pts1* t% times(1)
                else
                    dt = pts2*(t% times(1)-tMS_hook)
                endif
                dtr = t% times(1)-age
            case(HG)
                dt = pts2*(t% times(2)-t% times(1))
                dtr = t% times(2)-age

            case(RGB)
                dt = pts2*(t% times(3)- t% times(2))
                dtr = MIN(t% nuc_time,t% times(3))-age
!                print*,"dt", t% times(2),t% times(3),age

            case(HeBurn)
                dt = pts2*(t% times(4)- t% times(3))
                dtr = MIN(t% nuc_time,t% times(4))-age
!                print*, t% times(3),t% times(4),t% nuc_time
!                print*, dt,dtr,age
            case(EAGB)
                dt = pts3*0.1*(t% times(5)- t% times(4))
                dtr = MIN(t% nuc_time,t% times(5))-age
!                print*, t% nuc_time, t% times(5), age
            case(TPAGB)
                if (t% post_agb) then
                    dt = pts3*(t% agb% tfinal - t% agb% tini)
                    dtr = t% agb% tfinal - age
!                    print*,"aj,tn",age,t% nuc_time,t% agb% tfinal,t% agb% age
                else
                    dt = pts3*(t% times(6)-t% times(5))
                    dt = MIN(dt,0.001d0)
                    dtr = MIN(t% nuc_time,t% times(6))- age
                endif
            case(He_MS)
                ! this is same for he stars with sse formulae as well as interpolated tracks

                dt = pts1* t% times(7)
                dtr = t% times(7) - age
!                if (id==2)print*, 'deltat',t% times(7), age,dt,dtr

            case(He_HG:He_GB)
                if (t% star_type == sse_he_star) then
                    if(age < t% times(10))then
                        dt = pts2*(t% times(8) - age)
                    else
                        dt = pts2*(t% times(9) - age)
                    endif
                     dtr = t% nuc_time -age
                    !Note that this is done to avoid negative timesteps
                    ! that result from more massive cores than what sse formulae can handle
                    dtr = max(dtr,1d-10)
                elseif(t% pars% phase== HE_HG) then
                    dt = pts2*(t% times(8)-t% times(7))
                    dtr = MIN(t% nuc_time,t% times(8))-age
                    
                elseif(t% pars% phase== HE_GB) then
                    dt = pts2*(t% times(9)-t% times(8))
                    dtr = t% nuc_time-age
                endif
!                    print*, 'deltat',t% nuc_time,t% times(8),t% times(9),age,dt,dtr
            case(HeWD:NS)
                dt = MAX(0.1d0,age*10.0)
                dt = MAX(0.1d0,dt*10.0)
                !dt = MIN(dt,5.0E+02)
                dt = MIN(dt,1.d0)
                dtr = dt
            case(BH:Massless_REM)
                !dt = MIN(dt,1.d0)
                dt = 10.d0
                dtr = dt
        end select

            t% pars% dt = min(dt,dtr)
!            if (id ==1)print*,'timestep',dt,dtr,t% pars% dt,t% times(4),age, t% pars% phase,id

            if (t% pars% dt<=0.0 .and. t% ierr==0) then
                write(UNIT=err_unit,fmt=*)"METISSE error: invalid timestep", t% pars% dt ,"for phase and id", t% pars% phase,id
                write(UNIT=err_unit,fmt=*)"zams_mass, pars% mass, nuc_time, age, dt,dtr"
                write(UNIT=err_unit,fmt=*)t% zams_mass, t% pars% mass,t% nuc_time, age, dt,dtr
                t% ierr = -1
                t% pars% dt = 1e+10
                ! forcing a large timestep to avoid crashing of code outside this function
!               call stop_code
            endif
            
            if(t% ierr<0) code_error = .true.

            nullify(t)
      end subroutine METISSE_deltat
