! This file contains miscellaneous subroutines needed by METISSE
! to work in stand alone or otherwise
! Ideally these should be packed in a module
! But Fortran 77 does not know how to use modules
! So here we are

subroutine initialize_front_end(front_end_name)
    use track_support
    character (len=*), intent(in) :: front_end_name

!    print*, 'in front end',trim(front_end_name)

    if (ANY((/'MAIN','main'/)== trim(front_end_name))) then
        ! METISSE's main code as described in Agrawal et al. 2020
        ! Can be used to evolve single stars and/or debugging purposes.
        front_end = main

    elseif (ANY((/'SSE','sse','BSE','bse'/)== trim(front_end_name))) then
        ! SSE (Single Star Evolution) from Hurley et al. 2000
        ! BSE (Binary Star Evolution) from Hurley et al. 2002
        front_end = BSE
        
    elseif (ANY((/'COSMIC','cosmic'/)== trim(front_end_name))) then
        ! COSMIC (Compact Object Synthesis and Monte Carlo Investigation Code)
        ! Binary evolution code from Breivik et al. 2020
        front_end = COSMIC
        
    else
        print*, "Error: Unrecongnized front_end_name for METISSE"
        print*, "Choose from 'MAIN', 'SSE', 'BSE', 'COSMIC' "
        STOP
    endif
    
end subroutine initialize_front_end


subroutine allocate_track(n,mass)
    use track_support
    implicit none

    integer, intent(in):: n
    real(dp), intent(in), optional :: mass(:)
        
!    n= 1
!    n = size(mass)
!    print*,"I am in alloc_track with ", n,mass

    allocate(tarr(n))
    tarr% star_type = unknown
    tarr% pars% age = 0.d0
    tarr% pars% extra = 0
    tarr% pars% bhspin = 0.d0
    tarr% ierr = 0
    tarr% pars% dms = 0.d0
    tarr% pars% delta = 0.d0
    
end subroutine allocate_track


subroutine dealloc_track()
    use track_support
    implicit none
    
    integer:: n,i

    n = size(tarr)
!    print*," deallocating", n
!    n = 1
    do i = 1,n
!        deallocate(tarr(i)% times_new)
!        deallocate(tarr(i)% times)
        deallocate(tarr(i)% eep)
        deallocate(tarr(i)% tr)
        deallocate(tarr(i)% cols)
        deallocate(tarr(i)% bounds)
        if((tarr(i)% ierr/=0).and.verbose) write(UNIT=err_unit,fmt=*)'Error in evolving the system',i

    end do
    deallocate(tarr)
end subroutine dealloc_track


subroutine set_star_type(id)

! set star type to rejuvenated before calling star
    use track_support
        implicit none
        integer, intent(in) :: id
!        print*, 'setting star to reju',tarr(id)% pars% age,id

        tarr(id)% star_type = rejuvenated
end subroutine set_star_type

