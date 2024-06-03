module interp_support
    use track_support
    implicit none

    integer, parameter :: no_interpolation = 0
    integer, parameter :: linear = 1
    integer, parameter :: Steffen1990 = 2
    !integer, parameter :: extrapolation = 3
    logical :: debug_mass

    contains

    !   interpolates a new track given initial mass
    subroutine interpolate_mass(t, exclude_core)
        implicit none

        logical, intent(in) :: exclude_core

        real(dp) :: mass
        type(track), pointer :: t

        integer :: iseg, keyword,min_index
        type(track), pointer :: a(:), s(:)
        real(dp) :: f(3), dx, x(4), y(4), alfa, beta
        integer :: i, j, k, mlo, mhi,nt,age_col,start
        integer, allocatable :: eeps(:), excl_cols(:)
        
        debug_mass = .false.
!        if(t% is_he_track)debug_mass = .true.

        if (debug_mass) print*, 'in interpolate_mass',t% initial_mass,t% pars% phase
        mass = t% initial_mass
        
        if (mass/=mass .or. mass<=0.d0) then
            write(UNIT=err_unit,fmt=*)"Fatal Error: Invalid mass",mass,t% pars% phase
            t% ierr = -1
!            call stop_code
            return
        endif
        
        dx=0d0; alfa=0d0; beta=0d0; x=0d0; y=0d0
        
        ! this line is to avoid array length problem with multiple calls to fix-track
        if (allocated(t% tr) .and. (.not.exclude_core)) call deallocate_arrays(t)
        if (allocated(t% bounds)) deallocate(t% bounds)
        
        ! takes a set of EEP-tracks and find tracks for interpolation (a)
        call findtracks_for_interpolation(mass,t% is_he_track,t% bounds,min_index,keyword,iseg)
        
        mlo = 1
        mhi = size(t% bounds)
        
        if(t% is_he_track) then
            a => sa_he(t% bounds(mlo):t% bounds(mhi))
            s => sa_he
            excl_cols = core_cols_he
            age_col = i_he_age
            start = ZAMS_HE_EEP
            eeps = key_eeps_he
        else
            a => sa(t% bounds(mlo):t% bounds(mhi))
            s => sa
            excl_cols = core_cols
            age_col = i_age
            start = ZAMS_EEP
            eeps = key_eeps
        endif
        
        k = minloc(a(mlo:mhi)% ntrack,dim=1)
        t% min_index = min_index
        
        if (debug_mass) print*,"mass, keyword", mass,keyword,t% is_he_track
        if (debug_mass) print*,"interpolate mass" , a% initial_mass
        if (debug_mass) print*,"interpolate ntrack" , a% ntrack
        
        if (exclude_core) then
            nt = t% ntrack
            t% ntrack = min(nt,a(k)% ntrack)
        else
            call write_header(t,s(min_index))
            !t% ntrack = min(a(k)% ntrack,get_min_ntrack(t% star_type, t% is_he_track))
            t% ntrack = a(k)% ntrack
            allocate(t% tr(t% ncol+1, t% ntrack))
            t% tr = 0d0
        endif
        
        
        ! interpolate the new track for given initial mass
        ! based on keyword

        select case(keyword)
        case(no_interpolation)
            do j=1,t% ncol
                if (exclude_core .and. any(j .eq. excl_cols,1)) cycle
                t% tr(j,start:t% ntrack) = a(1)% tr(j,start:t% ntrack)
            end do
            
        case(linear)
            alfa = (t% initial_mass - a(mlo)% initial_mass)/(a(mhi)% initial_mass - a(mlo)% initial_mass)
            beta = 1d0 - alfa
            do i = start,t% ntrack
                do j=1,t% ncol
                    if (exclude_core .and. any(j .eq. excl_cols,1)) cycle
                    t% tr(j,i) = alfa*a(mhi)% tr(j,i) + beta*a(mlo)% tr(j,i)
                enddo
            enddo

        case(Steffen1990)
            x = a(mlo:mhi)% initial_mass
            dx = t% initial_mass - x(2)
            do i = start,t% ntrack
                do j=1,t% ncol
                    if (exclude_core .and. any(j .eq. excl_cols,1)) cycle
                    do k=1,4
                        y(k) = a(k)% tr(j,i)
                    enddo
                    call interp_4pt_pm(x, y, f)
                    t% tr(j,i) = y(2) + dx*(f(1) + dx*(f(2) + dx*f(3)))
                enddo
            enddo
        end select
        
        if (fix_track) call check_length(iseg,t,min_index,exclude_core)

        t% tr(i_age2,:) = t% tr(i_age2,:)*1E-6          !Myrs
        if (start>1) t% tr(:,1:start-1) = -1.d0

    
        ! check if age is monotonically increasing
        call mod_PAV(t% tr(i_age2,start:t% ntrack))
        
        ! check if mass is same or monotonically decreasing
        call smooth_track(t,start)
        
        ! recalibrate age from ZAMS
        
!        t% tr(i_age2,:) = t% tr(i_age2,:)- t% tr(i_age2,start)
        
        t% j_bgb = -1
        if (t% is_he_track .eqv. .false.) then
            !determine the base of the giant branch, if present
            if (t% initial_mass > Mcrit(2)% mass .and. t% initial_mass< Mcrit(5)% mass) then
                if (BGB_EEP>0 .and. t% ntrack >= BGB_EEP) then
                    !Red giant Branch
                    t% j_bgb = BGB_EEP
                else
                    if(i_mcenv>0) then
                        t% j_bgb = bgb_mcenv(t,TAMS_EEP,cHeIgnition_EEP)
                    elseif(t% ntrack >= cHeIgnition_EEP) then
                        t% j_bgb = base_GB(t)
                    endif
                    if (t% j_bgb<0 .and. debug_mass)print*, "Unable to locate BGB for",t% initial_mass
                endif
            endif
            if (exclude_core .eqv. .false.) t% j_bgb0 = t% j_bgb
        endif
        
        if (exclude_core) then
!            if (t% ntrack/=nt) write(UNIT=err_unit,fmt=*)'WARNING: track length changed',t% initial_mass,nt,t% ntrack
            t% ntrack = nt
        else
            ! get eeps
            t% neep = count(eeps<= t% ntrack,1)
            allocate(t% eep(t% neep))
            t% eep = pack(eeps, eeps<= t% ntrack)
            if (debug_mass) print*, 'eeps',t% neep,t% eep(t% neep), eeps(size(eeps))
        endif
        
        nullify(a,s)
    end subroutine interpolate_mass

    subroutine findtracks_for_interpolation(mass,is_he_track,bounds,min_index,keyword,iseg)
        ! takes a set of EEP-defined tracks and find tracks for interpolation
        real(dp), intent(in) :: mass
        logical, intent(in) :: is_he_track
        integer, intent(out) :: iseg,min_index,keyword
        
        type(track), pointer :: s(:)
        real(dp), allocatable :: mass_list(:)
        integer, allocatable :: bounds(:), cutoff(:)
        integer :: m_low,m_high,num_list,min_index1,j

        m_low = 0
        num_list = 0
        min_index = -1
        
        if (is_he_track) then
            cutoff = m_cutoff_he
            s => sa_he
        else
            cutoff = m_cutoff
            s => sa
        endif
        
        if (mass > s(size(s))% initial_mass*1.01) then
            min_index = size(s)
            keyword = no_interpolation
            allocate(bounds(1))
            bounds = min_index
            if (debug_mass) print*,"No interpolation: mass EXCEEDS highest initial mass, ",mass
            iseg = size(cutoff)-1
            deallocate(cutoff); nullify(s)
            return
        endif
            
        ! we don't want to search the whole list, only between the mass cutoffs
        ! therefore we create smaller list of initial_masses
        ! mcutoff arrays are defined in the z_support module
        do iseg = 1,size(cutoff)-1
            m_low = cutoff(iseg)
            m_high = cutoff(iseg+1)-1
!            print*,'check:',m_low, m_high, iseg,is_he_track

            if (m_high == m_low) m_high = m_low+1
            if (mass > s(m_high)% initial_mass*1.01) cycle
            num_list = m_high-m_low+1      !+1 to include count for m_low point
            allocate(mass_list(num_list))
            mass_list = s(m_low:m_high)% initial_mass
            exit
        end do
        
        if (debug_mass) print*, 'finding nbrs', mass, num_list,m_high,m_low,is_he_track
        !search the smaller list for the nearest mass
        call index_search(num_list,mass_list,mass,min_index1)
        if(min_index1>num_list) min_index1 = num_list
        
        ! rescale min_index for the bigger list
        min_index = m_low+min_index1-1

        if (debug_mass) print*,"min_index for mass", min_index,min_index1
        if (debug_mass) print*,"mass(min_index)", s(min_index)% initial_mass
        
        if(abs(mass-mass_list(min_index1))< mass_accuracy_limit) then
            ! track is already in the database
            keyword = no_interpolation
            allocate(bounds(1))
            bounds = min_index
            if (debug_mass) print*,"Interpolation NOT required during interpolate mass",mass
        elseif(mass < mass_list(2)) then
            !index adjustments for mass at boundaries
            keyword = linear
            allocate(bounds(2))
            bounds = [m_low,m_low+1]
            if (debug_mass) print*, "Close to lower cutoff"

        elseif(mass > mass_list(num_list-1)) then
            keyword = linear
            allocate(bounds(2))
            bounds = [m_high-1,m_high]
            if (debug_mass) print*, "Close to uppper cutoff"
        else
            keyword = Steffen1990
            allocate(bounds(4))
            if (mass_list(min_index1)> mass) then
                bounds = (/(j, j=min_index-2,min_index+1)/)
            else
                bounds = (/(j, j=min_index-1,min_index+2)/)
            end if
        endif
        
        if (debug_mass .and. mass < mass_list(1)) print*, "doing extrapolation",mass , mass_list(1)

		deallocate(mass_list,cutoff)
        nullify(s)
    end subroutine findtracks_for_interpolation

    subroutine write_header(b,a)
        implicit none
        type(track):: a
        type(track), pointer :: b
        
        b% star_type = a% star_type
        !b% version_string = a% version_string
        b% initial_Z = a% initial_Z
        b% initial_Y = a% initial_Y
        b% Fe_div_H = a% Fe_div_H
        b% alpha_div_Fe = a% alpha_div_Fe
        b% v_div_vcrit = a% v_div_vcrit
        b% ncol = a% ncol
        allocate(b% cols(b% ncol+1))
        b% cols(1: b% ncol) = a% cols
        b% cols(i_age2)% name = a% cols(i_age2)% name
         
        b% cols(b% ncol+1)% name = 'age_old'
        b% complete = .true.
        b% has_mass_loss = a% has_mass_loss
        b% is_he_track = a% is_he_track
    end subroutine

    subroutine check_length(iseg,t,min_index,exclude_core)
        type(track), pointer :: t
        integer :: min_index,iseg
        logical, intent(in) :: exclude_core
        
        type(track), pointer :: s(:)
        real(dp), allocatable :: mass_list(:)
        integer :: m_low,m_high,num_list
        integer :: i, m1, up_count,low_count,temp(4)
        integer :: min_ntrack, low_lim, upp_lim
        real(dp) :: upper_tol, lower_tol
        type(track) :: a(2)
        
        min_ntrack = get_min_ntrack(t% star_type, t% is_he_track)
        !check length
        if (debug_mass) print*,"checking length now"
        if ((t% ntrack >= min_ntrack) ) then
            if (debug_mass) print*,"length ok", t% initial_mass, t% ntrack,min_ntrack,t% is_he_track
            return
        else
            if (debug_mass) print*,"not complete", t% initial_mass, t% ntrack, min_ntrack,t% is_he_track
            
            if (t% is_he_track) then
                m_low = m_cutoff_he(iseg)
                m_high = m_cutoff_he(iseg+1)-1
                num_list = m_high-m_low+1      !to include count for m_low point
                s => sa_he(m_low:m_high)
                allocate(mass_list(num_list))
                mass_list = sa_he(m_low:m_high)% initial_mass
            else
                m_low = m_cutoff(iseg)
                m_high = m_cutoff(iseg+1)-1
                num_list = m_high-m_low+1      !to include count for m_low point
                s => sa(m_low:m_high)
                allocate(mass_list(num_list))
                mass_list = sa(m_low:m_high)% initial_mass
            endif
            ! rescale min_index for the smaller list
            min_index = min_index-m_low+1

            temp = 0
            up_count = 0; low_count = 0

            upper_tol = t% initial_mass + t% initial_mass* lookup_index
            call index_search(num_list,mass_list,upper_tol,upp_lim)
            upp_lim = min(num_list,upp_lim)

            lower_tol = t% initial_mass - t% initial_mass* lookup_index
            call index_search(num_list,mass_list,lower_tol,low_lim)
            low_lim = max(1,low_lim)

            if (s(min_index)% initial_mass> t% initial_mass) then
                m1 = min_index
            else
                m1 = min_index+1
            end if

            do i= m1-1, low_lim,-1
                if (s(i)% ntrack >= min_ntrack) then
                    low_count = low_count+1
                    temp(low_count)=i        !new_a(count) = s(i)
                if (low_count == 2) exit
                endif
            end do
            do i= m1, upp_lim
                if (s(i)% ntrack >= min_ntrack) then
                    up_count = up_count+1
                    temp(low_count+up_count)=i        !new_a(count) = s(i)
                    if (up_count == 2) exit
                endif
            end do
            
            select case(low_count)
            case(0)
                if (up_count==2) then  !extrapolate
                    !m =1
                    a = s(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"0,2"
                else
                    t% complete = .false.
                endif

            case(1)
                if(up_count>0) then !interpolate
                    temp(3) = 0       !if already not so
                    a = s(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"1,1"
                else
                    t% complete = .false.
                endif

            case(2)
                if(up_count>0) then !interpolate
                    temp(2) = 0  !low_count =2 and we don't want to use more distant track
                    temp(4) = 0 !if already not so
                    a = s(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"2,1"
                else  !extrapolate
                    !m=1
                    temp(3) = temp(1) !linearly increasing order
                    temp(1) = 0
                    a = s(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"2,0"
                endif
            end select
            
            call fix_incomplete_tracks(a,t,min_ntrack,exclude_core)
            if (debug_mass) print*, "new length", t% ntrack

            deallocate(mass_list)
            nullify(s)
            
        endif
        
    end subroutine check_length
    
    subroutine fix_incomplete_tracks(a,t,min_ntrack,exclude_core)
    !this has been modified for use with fix_incomplete_tracks only
        implicit none
        type(track), intent(in) :: a(:)
        type(track), pointer :: t
        integer, intent(in) :: min_ntrack
        logical, intent(in) :: exclude_core
        real(dp) :: alfa,beta,bprime
        integer :: i,j,n
        integer, allocatable :: excl_cols(:)

        real(dp), allocatable :: c(:,:)
        
        n = t% ntrack
        t% ntrack = min_ntrack
        
        if(t% is_he_track) then
            excl_cols = core_cols_he
        else
            excl_cols = core_cols
        endif
        
        if (.not. exclude_core) then
            !store orginal track tr in c
            allocate(c(t% ncol, n))
            c(:,:) = t% tr(1: t% ncol,:)

            !reallocate tr for rewriting with new length
            deallocate(t% tr)

            allocate(t% tr(t% ncol+1, t% ntrack))

            t% tr = 0d0
            t% tr(1: t% ncol,1:n) = c(:,1:n)
            deallocate (c)
        endif
            
        
        do i= n+1, t% ntrack
            do j=1,t% ncol
                if (exclude_core .and. any(j .eq. excl_cols,1)) cycle
                t% tr(j,i) = t% tr(j,n)
            end do
        end do
        
        
        if (t% complete) then
            if (debug_mass) print*, "new masses for interpolate", a% initial_mass
            !complete the track between temp_ntrack and t% ntrack
!                    call linear_interp(a,t,temp_ntrack)
            alfa=0d0; beta=0d0
            bprime = 0.d0
            alfa = (t% initial_mass - a(1)% initial_mass)/(a(2)% initial_mass - a(1)% initial_mass)
            beta = 1d0 - alfa

            !bprime: the offest from previously calculated value

            do i = n, t% ntrack
                do j=1,t% ncol
                    if (exclude_core .and. any(j .eq. excl_cols,1)) cycle
                    bprime = t% tr(j,n)-(alfa*a(2)% tr(j,n) + beta*a(1)% tr(j,n))
                    t% tr(j,i) = alfa*a(2)% tr(j,i) + beta*a(1)% tr(j,i) +bprime
                enddo
            enddo
        end if
        deallocate(excl_cols)
                
    end subroutine fix_incomplete_tracks

    subroutine smooth_track(t,start)
    implicit none
    type(track), pointer :: t

    integer :: i, start
    
    real(dp), pointer :: mass_list(:)

        mass_list => t% tr(i_mass,start:t% ntrack)

        do i = 2, size(mass_list)
            if (mass_list(i).le.0.d0) then
            ! although rare, sometime extrapolation can cause negative mass values
                mass_list(i) = mass_list(i-1)
            else
                mass_list(i) = min(mass_list(i), mass_list(i-1))
            endif
        end do
        nullify(mass_list)
        
    end subroutine smooth_track
    
    subroutine mod_PAV(y)
        !PAV from ISO, modified for steps in time instead of smoothing over with average values
        real(dp), intent(inout) :: y(:)
        integer :: i,j, n, start, last, m,old_start
        real(dp), allocatable :: d(:)
        real(dp) :: diff, h
        logical :: debug

        debug = .false.
        n = size(y)
        allocate(d(n-1))
        old_start = 0

        do while(.true.)
           d = y(2:n)-y(1:n-1)
           if(all(d>=0)) exit !test for monotonicity
           i = locate(d) !finds the first point in d that is < 0
           start = i
            last = n
            if (debug) print*, 'in mod_PAV', i, d(i),y(i),y(n)
            !if start is the greatest value not in the end, take one step up and redo
            if (start == old_start) start=old_start-1
            do while(.true.)
                if (y(n)- y(start)>0) exit
!                print*,y(n)- y(start),y(start),y(n)
                start = start-1
                if (start <1)  then
                    write(UNIT=err_unit,fmt=*)"Error in mod_PAV, start<1",i,n
!                    call stop_code
                endif
            end do
            if (debug) print*,'start, old_start', start, old_start
            do j = start+1,n
                diff= y(j)-y(start)
                if (debug) print*,"i and diff",i,diff
                if (diff>0) then
                    last = j
                    if (debug) print*,j,y(j),y(i)
                    if (y(last)>y(n)) last = n
                    exit
                endif
            enddo

           m = last - start
           h = (y(last)-y(start))/real(m)
            do j= 1,m
                y(j+start) = y(start)+ j*h
            end do
            old_start = start
            if (debug) print*,last,start,y(last),y(start),h,m
        end do
        
        deallocate(d)
    end subroutine mod_PAV
  
    integer function locate(y)
    real(dp), intent(in) :: y(:)
    integer :: i, n
    n=size(y)
    do i=1,n
       if(y(i)<0.0)then
          locate=i
          return
       endif
    end do
    locate=0
    end function locate
    
    real(dp) function get_secondary_age(phase,input_age,times,times_new)
        real(dp), intent(in) :: input_age,times(:),times_new(:)
        integer :: phase
        real(dp) :: age2,frac
    
        frac = 1.d0
        select case (phase)
        case(low_mass_MS:MS)
            age2 = input_age
!            age2= input_age*(t% MS_time/t% ms_old)
        case(HG)
            frac = (times_new(phase)-times_new(phase-1))/(times(phase)-times(phase-1))
            age2 = times_new(phase-1)+((input_age-times(phase-1))*frac)
        case(RGB:HeBurn)
            ! since all stars don't have a GB, we skip the phase here
            frac = (times_new(HeBurn)-times_new(HG))/(times(HeBurn)-times(HG))
            age2 = times_new(HG)+((input_age-times(HG))*frac)
        case(EAGB:TPAGB)
            frac = (times_new(phase)-times_new(phase-1))/(times(phase)-times(phase-1))
            age2 = times_new(phase-1)+((input_age-times(phase-1))*frac)
        case(He_MS)
            age2 = input_age
!            age2= input_age*(t% MS_time/t% ms_old)
        case(HE_HG:HE_GB)
            frac = (times_new(HE_GB)-times_new(He_MS))/(times(HE_GB)-times(He_MS))
            age2 = times_new(He_MS)+((input_age-times(He_MS))*frac)
        end select
!        print*, "in interp2", input_age,age2,frac,phase
        get_secondary_age = age2
    end function
                    
    subroutine interpolate_age(t, input_age, icolumn, val)
        implicit none
        
        type(track), pointer :: t

        real(dp), intent(in) :: input_age
        integer, intent(in), optional :: icolumn
        real(dp), intent(out), optional :: val
        integer :: jstart, jend, age_col
        real(dp) :: age, age2

        real(dp) :: f(3), dx, x(4), y(4), alfa, beta
        integer :: j, k, mlo, mhi, pass, n_pass
        
        real(dp), allocatable :: new_line(:,:)
        integer, allocatable :: min_eeps(:)

        logical :: debug

        debug = .false.
!         if (t% is_he_track) debug = .true.
!        if (t% pars% phase>=4 .and. (present(icolumn).eqv..false.)) debug = .true.

        if (debug) print*,"in interpolate age",t% pars% phase
        dx=0d0; alfa=0d0; beta=0d0; x=0d0; y=0d0
        jstart = 1
        jend = t% ncol
        if (present(icolumn)) then
            jstart = icolumn
            jend = icolumn
            if (debug) print*,"only interpolating in column number",icolumn
        endif
        
        allocate (new_line(t% ncol,1))
        new_line = -1.d0
        n_pass = 1!2
        
        if (t% pars% phase<=1 .or. t% pars% phase ==7) n_pass = 1
            
        age2 = get_secondary_age(t% pars% phase,input_age,t% times,t% times_new)
        t% pars% age2 = age2

        do pass = 1, n_pass
            if (pass == 2) then
                !non core values post MS
                age = age2
                age_col = i_age2     !i_age2 = new age, age col in the main array
            else
                !core values post MS, all values during MS
                age = input_age
                age_col = i_age     !i_age = old age, stored at t% ncol+1
                if (t% is_he_track) age_col = i_he_age
            endif
            
            call find_nearest_eeps(t,min_eeps, age, age_col)

            mlo = minval(min_eeps)
            mhi = maxval(min_eeps)

            if (debug) print*, 'pass', n_pass,pass,age_col,age, t% pars% phase
            if (debug) print*,"neighbouring_eeps", min_eeps
!            if (debug) print*,"ages", t% tr(age_col,mlo:mhi)

            if (mhi == mlo) then
                if (debug) print*, "no interp in age needed"
                do j=jstart,jend
                    if (pass==2) then
                        ! check if it is core-related quantity
                        if (check_core_quant(j,t% is_he_track)) cycle
                    endif
                    new_line(j,1) = t% tr(j,mlo)
                end do

            elseif ((mhi-mlo)<4) then
                !linear interpolation
                alfa = (age - t% tr(age_col,mlo))/(t% tr(age_col,mhi) - t% tr(age_col,mlo))
                beta = 1d0 - alfa
                if (debug) print*, "doing linear interp in age",alfa,age,t% tr(age_col,mlo),t% tr(age_col,mhi)

                do j = jstart,jend
                    if (pass==2) then
                        ! check if it is core-related quantity
                        if (check_core_quant(j,t% is_he_track)) cycle
                    endif

                    new_line(j,1) = alfa*t% tr(j,mhi) + beta*t% tr(j,mlo)
                    if (new_line(j,1)/= new_line(j,1) .and. t% ierr==0) then
                        write(UNIT=err_unit,fmt=*) 'Warning: NaN encountered during interpolation age',&
                            t% initial_mass,input_age,j,mhi,mlo,t% pars% phase
                        t% ierr = -1
            !            call stop_code
                        return
                    endif
                end do
                if (debug) print*, "ending linear interp in age"

            else
                ! currently only linear interpolation is used
                if (debug) print*, "doing cubic interp in age"

                x = t% tr(age_col,mlo:mhi)
                dx = new_line(age_col,1) - x(2)

                if (age< x(2) .or. age> x(3)) then
                   write(UNIT=err_unit,fmt=*)"Error in cubic interpolation in interp_support"
!                   call stop_code
                endif

                do j = jstart,jend
                    if (pass==2) then
                        ! check if it is core-related quantity
                        if (check_core_quant(j,t% is_he_track)) cycle
                    endif
                    do k=1,4
                        y(k) = t% tr(j,mlo-1+k)
                    enddo
                    call interp_4pt_pm(x, y, f)
                    new_line(j,1) = y(2) + dx*(f(1) + dx*(f(2) + dx*f(3)))
                enddo
            endif
            if (allocated(min_eeps)) deallocate(min_eeps)
        end do
                    
        if (present(icolumn)) then
            val = new_line(icolumn,1)
        else
            call save_values(new_line,t% pars)
        endif        
        deallocate(new_line)

        if (t% pars% mass <0.0) then
            write(UNIT=err_unit,fmt=*)"Fatal Error: mass <0 in interpolate age",input_age,t% pars% phase
!            call stop_code
        endif
        if (debug) print*, 'exiting interpolate_age'
        
    end subroutine interpolate_age
    
    logical function check_core_quant(j,is_he_track)
        integer, intent(in) :: j
        logical, intent(in) :: is_he_track
        integer, allocatable :: excl_cols(:)

        if (is_he_track) then
            excl_cols = core_cols_he
        else
            excl_cols = core_cols
        endif

        if (any(j .eq. excl_cols,1) .or. j == i_age2) then
            check_core_quant = .true.
        else
            check_core_quant = .false.
        endif
        
        deallocate(excl_cols)
!        print*, 'core quant',j, pass,check_core_quant
    end function check_core_quant


    subroutine find_nearest_eeps(t,min_eeps,age,age_col)
        implicit none
        type(track), pointer :: t
        integer, allocatable, intent(out) :: min_eeps(:)
        real(dp), intent(in) :: age
        integer :: i, j, len_eep, min_index,age_col,initial_eep
        real(dp), allocatable :: age_list(:)
        real(dp) :: last_age
        logical :: debug
        
        debug =  .false.
        !Todo: min_eeps-> nbr_eeps
!         if (t% is_he_track) debug = .true.
        
        initial_eep = ZAMS_EEP
        if (t% is_he_track) initial_eep = ZAMS_HE_EEP
        
        if (age .lt. t% tr(age_col,initial_eep)) then
        ! check for lower boundary
            allocate(min_eeps(1))
            min_eeps = initial_eep
            if (debug)write(*,*)"age<initial_eep",age, t% tr(age_col,initial_eep)
            return
        elseif (age .gt. t% tr(age_col,t% eep(t% neep))) then
        ! check for upper boundary
            allocate(min_eeps(1))
            min_eeps = t% eep(t% neep)
            if (debug)write(*,*)"age>t%neep",age,t% ntrack,t%neep,t% eep(t% neep)
        endif
                
        last_age = 0.d0
        
        do i = 1,t% neep-1
            if (debug) print*,"loc_low", t% eep(i), t% eep(i+1),t%neep

            last_age = t% tr(age_col,t% eep(i+1))
            if (debug) print*,"ages", age,last_age, t% eep(i)< initial_eep

            if ((t% eep(i)< initial_eep) .or. (age .gt. last_age)) cycle
            len_eep = t% eep(i+1)-t% eep(i)+1
            allocate(age_list(len_eep))
            age_list = t% tr(age_col,t% eep(i):t% eep(i+1))

            if (debug) print*,"len_eep:",len_eep,"bounds:",age_list(1),age_list(len_eep)

            if (len_eep<5) then
            call index_search(len_eep,age_list,age,min_index)
        else
            min_index = binary_search(len_eep,age_list,age)
        endif
            if(abs(age_list(min_index)-age)< tiny) then        !less than a year
                if (debug) print*,"no interpolation, min_index", age_list(min_index)
                allocate(min_eeps(1))
                min_eeps = min_index
!                            a => b(:,min_index: min_index)
            elseif(age< age_list(2)) then
                if (debug) print*,"age< age_list(2)", age_list(1:2)
                allocate(min_eeps(2))
                min_eeps = [1,2]
!                            a => b(:,1:2)
            elseif(age> age_list(len_eep-1)) then
                if (debug) print*,"age> age_list(len_eep-1)", age_list(len_eep-1:len_eep)
                allocate(min_eeps(2))
                min_eeps = [len_eep-1,len_eep]
!                            a => b(:,len_eep-1:len_eep)
            elseif(age < age_list(min_index)) then
                if (debug) print*,"age< min_index", age_list(min_index-2:min_index+1)

!                            a => b(:,min_index-2:min_index+1)
                allocate(min_eeps(2))
                min_eeps = (/(j, j=min_index-1,min_index)/)
!                            a => b(:,min_index-1:min_index) !forcing linear

            else
                if (debug) print*,"age> min_index",age_list(min_index-1:min_index+2)
!                            a => b(:,min_index-1:min_index+2)
                allocate(min_eeps(2))
                min_eeps = (/(j, j=min_index,min_index+1)/)
!                            a => b(:,min_index:min_index+1)
            endif
            ! original min_eeps were only a given primary eep,
            ! scale them back to full track
            min_eeps = min_eeps +t% eep(i)-1
            deallocate(age_list)
            exit
        end do
    
        if(.not.allocated(min_eeps)) then
            write(UNIT=err_unit,fmt=*)'Error finding nearest eeps for age:',age,age_col
        endif
        
    end subroutine find_nearest_eeps
    
    subroutine save_values(new_line,pars)
        type(star_parameters) :: pars
        real(dp),intent (in) :: new_line(:,:)
        real(dp) :: lim_R
        
        pars% mass = new_line(i_mass,1)
        pars% McHe = new_line(i_he_core,1)
        pars% McCO = new_line(i_co_core,1)
        pars% log_L = new_line(i_logL,1)
        pars% luminosity = 10**pars% log_L
        pars% log_Teff = new_line(i_logTe,1)
        pars% Teff = 10**(pars% log_Teff)
        pars% log_R = new_line(i_logR,1)
    !        pars% log_R = 2*(3.762+(0.25*pars% log_L)-pars% log_Teff )
    
        !restrict radius from going beyond the hayashi
        !limit during extrapolation
        lim_R = 2*(3.762+(0.25*pars% log_L)-3.555 )
        if (pars% log_R>lim_R) pars% log_R = lim_R

        pars% radius = 10**pars% log_R
        pars% core_radius = -1.0
        
        if (pars% phase <= TPAGB) then
            if (pars% phase == TPAGB ) then
                pars% core_mass = pars% McCO
                if (i_RCO_core>0) pars% core_radius = new_line(i_RCO_core,1)
            else
                pars% core_mass = pars% McHe
                if (i_RHe_core>0) pars% core_radius = new_line(i_RHe_core,1)
            endif
            if (i_mcenv>0) pars% mcenv = new_line(i_mcenv,1)
            if (i_rcenv>0) pars% rcenv = new_line(i_rcenv,1)
            if (i_MoI>0) pars% moi = new_line(i_MoI,1)
        elseif(pars% phase >= He_MS) then
            pars% core_mass = pars% McCO
            if (i_he_RCO>0) pars% core_radius = new_line(i_he_RCO,1)
            if (i_he_mcenv>0) pars% mcenv = new_line(i_he_mcenv,1)
            if (i_he_rcenv>0) pars% rcenv = new_line(i_he_rcenv,1)
            if (i_he_MoI>0) pars% moi = new_line(i_he_MoI,1)
        endif
            
    end subroutine
                    
                    
    subroutine calculate_timescales(t)
        !calculate timescales associated with different phases (0-6)
        implicit none
        
        type(track), pointer :: t
        integer :: i
        real(dp), pointer :: age(:)=> NULL()


        if (t% is_he_track) then
            call calculate_he_timescales(t)
            return
        endif
        
        age => t% tr(i_age2,:)
        t% times = undefined

        do i = 1, t% neep
            if (t% eep(i) == TAMS_EEP) then    !MS
                t% times(MS) = age(TAMS_EEP)

            elseif (t% eep(i) == cHeIgnition_EEP) then
                !Herztsprung gap
                t% times(HG) = age(cHeIgnition_EEP)
                !t% times(HG) gets modified to t(BGB)
                ! if RGB phase is present
                t% times(RGB) = t% times(HG)

            elseif (t% eep(i) == TA_cHeB_EEP) then
                !red_HB_clump /core He Burning
                t% times(HeBurn) = age(TA_cHeB_EEP)

            elseif (t% eep(i) == cCBurn_EEP) then
                !EAGB/ core C burning
                t% times(EAGB) = age(cCBurn_EEP)
                
            elseif (t% eep(i) == TPAGB_EEP) then
                !AGB
                t% times(EAGB) = age(TPAGB_EEP)
                
            elseif (t% eep(i) == post_AGB_EEP) then
                !TP-AGB :only for low_inter mass stars
                t% times(TPAGB) = age(post_AGB_EEP)
            endif
        enddo
        !print*,"bgb",BGB_EEP,identified(BGB_EEP)

        t% times(11) = age(min(Final_EEP,t% ntrack))
        t% MS_time = t% times(MS)
        !Todo: nuc_time should be for WR phase
        t% nuc_time = t% times(11)
        
        !Red giant Branch
        if (t% j_bgb > 1) t% times(HG) = age(t% j_bgb)
                    

        nullify(age)
    end subroutine calculate_timescales
    
    subroutine calculate_he_timescales(t)
        !calculate timescales associated with different he star phases (7,8,9)
        implicit none
        
        type(track), pointer :: t
        integer :: i
        real(dp), pointer :: age(:)=> NULL()

        age => t% tr(i_age2,:)
        
        t% times(He_MS:) = undefined

        do i = 1, t% neep
            if (t% eep(i) == TAMS_HE_EEP) then    !MS
                t% times(He_MS) = age(TAMS_HE_EEP)

            elseif (t% eep(i) == TPAGB_HE_EEP) then
                !helium Herztsprung gap
                t% times(HE_HG) = age(TPAGB_HE_EEP)

            elseif (t% eep(i) == cCBurn_HE_EEP) then
                !EAGB/ core C burning
                t% times(HE_HG) = age(cCBurn_HE_EEP)
                
!            elseif (t% eep(i) == post_AGB_HE_EEP) then
!                !TP-AGB :only for low_inter mass stars
!                t% times(10) = age(post_AGB_HE_EEP)
            endif
        enddo
        t% times(HE_GB) = t% times(HE_HG)

        t% times(11) = age(min(Final_EEP_HE,t% ntrack))
        
        t% nuc_time = t% times(11)
        t% MS_time = age(TAMS_HE_EEP)

        !Red giant Branch
        if (t% j_bgb > 1) t% times(HG) = age(t% j_bgb)
        
        
        if (t% initial_mass > Mcrit_he(4)% mass .and. t% initial_mass< Mcrit_he(5)% mass) then
            if (identified(GB_HE_EEP)) then
                t% times(HE_HG) = age(GB_HE_EEP)
            ! TODO: this needs to be written
!            elseif (t% ntrack >TAMS_EEP) then
!                j_bgb =  base_GB(t)
!                if (j_bgb>0) then
!                    j_bgb = j_bgb+TAMS_EEP-1
!                    !Red giant Branch
!                    t% times(HG) = age(j_bgb)
!                elseif (debug) then
!                    print*, "Unable to locate BGB ", j_bgb
!                end if
            endif
        endif
        
        nullify(age)
    end subroutine calculate_he_timescales
    
    !from MESA-r7503/1d_interp/
    subroutine interp_4pt_pm(x, y, a)
        ! returns coefficients for monotonic cubic interpolation from x(2) to x(3)
        real(dp), intent(in)    :: x(4)    ! junction points, strictly monotonic
        real(dp), intent(in)    :: y(4)    ! data values at x's
        real(dp), intent(inout)   :: a(3)    ! coefficients
        real(dp) :: h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3

        !integer, parameter :: pm_work_size = 3  !from mesa/interp_1d_def.f90

        ! for x(2) <= x <= x(3) and dx = x-x(2),
        ! y(x) = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
        h1 = x(2)-x(1)
        h2 = x(3)-x(2)
        h3 = x(4)-x(3)
        s1 = (y(2)-y(1))/h1
        s2 = (y(3)-y(2))/h2
        s3 = (y(4)-y(3))/h3
        p2 = (s1*h2+s2*h1)/(h1+h2)
        p3 = (s2*h3+s3*h2)/(h2+h3)
        as2 = abs(s2)
        ss2 = sign(1d0, s2)
        yp2 = (sign(1d0, s1)+ss2)*min(abs(s1), as2, 0.5d0*abs(p2))
        yp3 = (ss2+sign(1d0, s3))*min(as2, abs(s3), 0.5d0*abs(p3))
        a(1) = yp2
        a(2) = (3*s2-2*yp2-yp3)/h2
        a(3) = (yp2+yp3-2*s2)/(h2*h2)
    end subroutine interp_4pt_pm
    
    integer function eqv_eep(EEP2,EEP1,EEP_OLD2,EEP_OLD1,m)
        integer, intent(in) :: EEP2,EEP1,EEP_OLD2,EEP_OLD1,m
        real(dp) ::  frac

        frac = (M-EEP1)* 1.d0 /(EEP2-EEP1)
        eqv_eep = nint(EEP_OLD1+(frac* 1.d0 *(EEP_OLD2-EEP_OLD1)))
        
    end function
    
    subroutine get_initial_mass_for_new_track(t,id,mnew,eep_m)

        real(dp), intent(in) :: Mnew

        type(track), pointer :: t
        integer :: eep_m,id

        type(track), pointer :: s(:)
        real(dp), pointer :: age_list(:)
        real(dp), allocatable:: mlist(:), Mmax(:), Mmin(:)
        real(dp) :: alfa,beta,age
        integer :: i,j,k,nt,eep_n,num_list,Mupp,Mlow

        logical :: debug

        debug = .false.
!        if(id ==2 .and. t% is_he_track) debug = .true.
!        if (id ==1) debug = .true.
!        if (t% star_type==rejuvenated .or. t%star_type== switch) debug = .true.
                    
        nt = t% ntrack
        
        eep_n = -1
        Mlow = -1
        Mupp = -1
    
        IF (eep_m<0) THEN
            !using the original age of the star to keep core properties comparable
            age = t% pars% age
            
            ! create appropiate age pointers
            if (t% is_he_track .and. t% star_type /= switch) then
                age_list => t% tr(i_he_age,1:nt)
                initial_eep = ZAMS_HE_EEP
            elseif (t% star_type == switch .and. (.not. t% is_he_track)) then
                if (debug) print*, 'switching from he to h star'
                age_list => t% tr(i_he_age,1:nt)
                initial_eep = ZAMS_HE_EEP
            else
                if (debug.and.(t% star_type == switch)) print*, 'switching from h to he star'
                age_list => t% tr(i_age,1:nt)
                initial_eep = ZAMS_EEP
            endif
            if (debug) print*,"getting new initial mass mnew at age and phase: ",mnew,age,t% pars% phase,id,t% is_he_track
            
!            call index_search(nt-initial_eep+1,age_list(initial_eep:),age,eep_m)
            eep_m = binary_search(nt-initial_eep+1,age_list(initial_eep:),age)

            eep_m = eep_m+initial_eep-1

            if (eep_m > nt) eep_m = nt
            if (debug) print*,"nearest index eep_m, ntrack : ",eep_m,nt, age, age_list(eep_m)
            
            nullify(age_list)
        
        ENDIF
        
        if (t% star_type == switch) then
            if (t% is_he_track)then
                if (eep_m <= cHeIgnition_EEP) then
                    eep_m = ZAMS_HE_EEP
                elseif (eep_m <= TA_cHeB_EEP) then
                    eep_m = eqv_eep(TA_cHeB_EEP,cHeIgnition_EEP,TAMS_HE_EEP,ZAMS_HE_EEP,eep_m)
                elseif (eep_m<=TPAGB_EEP) then
                    eep_m = eqv_eep(TPAGB_EEP,TA_cHeB_EEP,TPAGB_HE_EEP,TAMS_HE_EEP,eep_m)
                elseif (eep_m<=cCBurn_EEP) then
                    eep_m = eqv_eep(cCBurn_EEP,TA_cHeB_EEP,cCBurn_HE_EEP,TAMS_HE_EEP,eep_m)
                elseif (eep_m >cCBurn_EEP) then
                    eep_m  = cCBurn_HE_EEP
                 elseif (eep_m >TPAGB_EEP) then
                    eep_m  = TPAGB_HE_EEP
                endif
            else
                if (eep_m<=TAMS_HE_EEP)then
                    eep_m = eqv_eep(TAMS_HE_EEP,ZAMS_HE_EEP,TA_cHeB_EEP,cHeIgnition_EEP,eep_m)
                elseif (eep_m<=TPAGB_HE_EEP) then
                    eep_m = eqv_eep(TPAGB_HE_EEP,TAMS_HE_EEP,TPAGB_EEP,TA_cHeB_EEP,eep_m)
                elseif (eep_m<=cCBurn_HE_EEP) then
                    eep_m = eqv_eep(cCBurn_HE_EEP,TAMS_HE_EEP,cCBurn_EEP,TA_cHeB_EEP,eep_m)
                endif
                
             endif
             ! get correct t% min_index
            t% min_index = 1
            if (debug) print*, 'new eep_m for switch', eep_m
        endif
        
        if (t% is_he_track) then
            s => sa_he
            Mmin = Mmin_he_array
            Mmax = Mmax_he_array
        else
            s => sa
            Mmin = Mmin_array
            Mmax = Mmax_array
        endif
        eep_m = min(eep_m, size(Mmin))
        ! It's crucial to avoid extrapolation here as it results in serious issues
        ! Mmax_array is the maximum mass at given eep amongst all input tracks, similarily Mmin_array has minimum
        ! first check if mass bounds exist at eep_m,
        ! if not check higher eeps (older age) in case of mass loss, and lower eeps for mass gain
        
        if (debug)print*, 'check',Mmin(eep_m),Mmax(eep_m)
        if (t% pars% delta.lt.0.d0) then
        ! mass loss
        if (debug) print*, 'mass loss',t% pars% delta
            do j = 0,nt
                if (eep_m+j<= nt) then
                    if (Mnew>= Mmin(eep_m+j).and.Mnew<=Mmax(eep_m+j)) then
                        eep_n = eep_m+j
                        exit
                    endif
                endif
            end do
        else
        ! mass gain
        if (debug) print*, 'mass gain',t% pars% delta
            do j = 0,nt
                if (eep_m-j>=1) then
                    if (Mnew>= Mmin(eep_m-j).and.Mnew<=Mmax(eep_m-j)) then
                        eep_n = eep_m-j
                        exit
                    endif
                endif
            end do
        endif
        if (debug .and. eep_n/=eep_m) print*,"modified eep_n : ",eep_n

        if (t% star_type == switch .and. eep_n<0) eep_n = eep_m
        
        if (eep_n >0) then
            ! get mass bounds for Mnew at eep_n
            num_list = size(s)
            allocate(mlist(size(s)))
            
            do i = 1, num_list
                if (s(i)% ntrack >= eep_n) then
                    mlist(i) = s(i)% tr(i_mass,eep_n)
                else
                    mlist(i) = -1.d0
                end if
            end do
            
            !find lower bound, start at Min_index
            do i = 0, num_list
                k = t% min_index+i

                if (k <= num_list) then
                    if (mlist(k) >0 .and. mlist(k).le.Mnew) then
                        Mlow = k
                        exit
                    endif
                endif
                ! search the other side now
                k = t% min_index-i
                if (k >=1 .and. i>0) then
                    if (mlist(k) >0 .and. mlist(k).le.Mnew) then
                        Mlow = k
                        exit
                    endif
                endif
            end do
            
!            find upper bound, start at Min_index
            do i = 0, num_list
                k = t% min_index+i
                if (k <= num_list .and. k>=1) then
                    if (mlist(k) >0 .and. mlist(k).gt.Mnew) then
                        Mupp = k
                        exit
                    endif
                endif

                ! search the other side now
                k = t% min_index-i
                if (k >=1 .and. k<=num_list) then
                    if (mlist(k) >0 .and. mlist(k).gt.Mnew) then
                        Mupp = k
                        exit
                    endif
                endif
            end do
        endif
        
        if (debug)  print*, "Mup", Mupp, "mlow",Mlow,"min_index",t% min_index
        ! if no solution is found, we keep using the old tracks for interpolation

        if (t% star_type == switch) then
            if (Mlow<1) then
                Mlow = 1
                Mupp = 2
            elseif (Mupp> num_list) then
                Mupp = num_list
                mlow = num_list -1
            endif
        endif
            
        if(Mlow < 0 .or. Mupp <0) then
            if (debug) print*,"Error: beyond the bounds for interpolation"
            if (debug) print*, "Mlow,Mupp,num_list,mnew,eep_n", &
                    Mlow,Mupp,num_list,mnew,eep_n
        else
            if (debug) print*,"ini_old",t% initial_mass,"mnew =",mnew,"masses at Mup =",mlist(Mupp),"mlow = ",mlist(Mlow)
            
            !intrepolate the bounds and their initial masses to get the initial mass for the new track
            alfa = (Mnew - mlist(Mlow))/(mlist(Mupp) - mlist(Mlow))
            beta = 1d0 - alfa
            t% initial_mass = alfa*s(Mupp)% initial_mass + beta*s(Mlow)% initial_mass
            
            if (debug) print*, "new ini mass",t% initial_mass, s(Mupp)% initial_mass, s(Mlow)% initial_mass
        endif
        
        !eep_m not needed except for switch cases
        ! sending out a dummy value to avoid segmentation faults
        nullify(s)
        deallocate(Mmax,Mmin)
        if (allocated(mlist)) deallocate(mlist)
    end subroutine get_initial_mass_for_new_track
    
    
  end module interp_support
