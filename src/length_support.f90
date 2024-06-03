module length_support

    use track_support
    use z_support, only: Mcrit, m_cutoff
    implicit none
    contains

 subroutine check_length_input_track(xa)
        type(track), allocatable :: xa(:)
        type(track), pointer :: a(:)
        integer :: bounds(4)
        real(dp) :: f(3), dx, x(4), y(4)
        integer :: i, j, k, mlo, mhi,min_ntrack,n,iend
        logical:: debug
        
        do n = 1,size(xa)
            
            min_ntrack = get_min_ntrack(xa(n)% star_type, xa(n)% is_he_track)
            !check length
            if ((xa(n)% ntrack >= min_ntrack) ) then
                if (debug) print*,"length ok", xa(n)% initial_mass, xa(n)% ntrack
                return
            else
               if (debug) print*,"not complete", xa(n)% initial_mass, xa(n)% ntrack, min_ntrack

                if(n>3 .and. n<size(xa)-2) then
                    bounds = [n-2,n-1,n+1,n+2]
                endif
                       
                mlo = 1
                mhi = 4
                
                a => xa(bounds(mlo):bounds(mhi))
                k = minloc(a(mlo:mhi)% ntrack,dim=1)

                iend = min(a(k)% ntrack, min_ntrack)
                
                x = a% initial_mass
                dx = xa(n)% initial_mass - x(2)
                do i=xa(n)% ntrack,iend
                    do j=1,xa(n)% ncol
                        do k=1,4
                            y(k) = a(k)% tr(j,i)
                        enddo
                        call interp_4pt_pm(x, y, f)
                        xa(n)% tr(j,i) = y(2) + dx*(f(1) + dx*(f(2) + dx*f(3)))
                    enddo
                enddo

                nullify(a)
            endif
        end do
    end subroutine check_length_input_track

    
    
    subroutine rev()
    
    
    if (t% pars% phase>= He_MS .and. t% pars% phase<=He_GB .and. kw<=TPAGB)then
                if (debug)print*,'rev to old initial_mass for phase',t% pars% phase,kw
                t% pars% phase = kw
                if (abs(t% zams_mass_old-t% zams_mass)>1.0d-12) then
                    ! need new track core properties
                    if (debug) print*, 'diff in mass', t% zams_mass_old,t% zams_mass,t% zams_mass_old-t% zams_mass

                    t% is_he_track = .false.
                    t% star_type = switch
                    mass_hold = t% initial_mass
                    mnew = t% zams_mass_old
                    eep_m = TAMS_HE_EEP
                    call get_initial_mass_for_new_track(t,idd,mnew,eep_m)
                    call interpolate_mass(t,exclude_core)
                    t% zams_mass_old = t% zams_mass
                    t% initial_mass = mass_hold
                    call calculate_timescales(t)
                endif
                t% initial_mass = t% initial_mass_old
                mass_check = .true.
            endif
            
            
            
            IF(kw>MS .and. kw/=He_MS) THEN
                exclude_core = .true.
            else
                t% zams_mass_old = t% zams_mass
            ENDIF
            
             t% pars% delta = 0.d0
                t% times_new = -1.d0
                
                if (exclude_core) then
                    !store core properties for post-main sequence evolution
                    if (debug) print*,'post-main-sequence star'
                    times_old = t% times

                    if ((t% pars% mcenv/t% pars% mass).ge.0.2d0) then
!                     .and.(t% pars% env_frac.ge.0.2)
                        allocate(rlist(t% ntrack))
                        rlist = t% tr(i_logR,:)
                        consvR = .true.
                    endif
    !
                    call interpolate_mass(t,exclude_core)
                   !write the mass interpolated track if write_eep_file is true
                    if (kw>=1 .and. kw<=4 .and. .false.) call write_eep_track(t,mt)
                    
                    call calculate_timescales(t)
                    t% times_new = t% times
                    t% times = times_old
                    t% nuc_time = t% times(11)
                    t% ms_time = t% times(MS)
                    if (t% is_he_track) t% ms_time = t% times(He_MS)

                    !TEST: R remains unchanged if Mcenv is significant
                    if (consvR) then
                        t% tr(i_logR,:) = rlist(1:t% ntrack)
                        deallocate(rlist)
                    endif
            endif
            t% pars% mass = mt
    
  end module length_support

