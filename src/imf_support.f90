module imf_support
  use track_support, only: dp
  implicit none
  contains

  subroutine sample_uniform(n, mmin, mmax, marray)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: mmin, mmax
    real(dp), intent(out) :: marray(n)
    real(dp) :: h
    integer :: i

    marray = 0.d0
    !linearly spaced
    h = abs(mmax- mmin)/(n-1)
    do i= 1,n
        marray(i) = mmin + (i-1)*h
    end do
  end subroutine sample_uniform
    
!  subroutine sample_uniform_log()
!  !for log spaced
!  end subroutine

  subroutine sample_kroupa_imf(marray, n, mmin, mmax)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: mmin,  mmax
    real(dp), intent(out) :: marray(n)
    real(dp) :: r, x, norm1, norm2, total, mbreak
    real(dp), parameter :: alpha1 = 1.3d0, alpha2 = 2.3d0
    integer :: i

    marray = 0.d0
    mbreak = max(mmin,0.5d0)

    ! Compute cumulative normalization for piecewise power-law
    norm1 = (mbreak**(1.d0 - alpha1) - mmin**(1.d0 - alpha1)) / (1.d0 - alpha1)
    norm2 = (mmax**(1.d0 - alpha2) - mbreak**(1.d0 - alpha2)) / (1.d0 - alpha2)
    total = norm1 + norm2

    call random_seed()  ! Initialize RNG

    i=1
    do while(i<=n)
       call random_number(r)
       if ((norm1>0.d0) .and. (r < norm1 / total)) then
        print*, 'i m here'
          call random_number(x)
          marray(i) = ((x * (mbreak**(1.d0 - alpha1) - mmin**(1.d0 - alpha1)) + mmin**(1.d0 - alpha1)))**(1.d0 / (1.d0 - alpha1))
          i=i+1
       else
          call random_number(x)
          marray(i) = ((x * (mmax**(1.d0 - alpha2) - mbreak**(1.d0 - alpha2)) + mbreak**(1.d0 - alpha2)))**(1.d0 / (1.d0 - alpha2))
          i=i+1
       end if
    end do
    
  end subroutine sample_kroupa_imf

end module imf_support
