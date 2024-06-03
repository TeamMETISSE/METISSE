subroutine METISSE_gntage(mc,mt,kw,zpars,m0,aj,id)

    use track_support
    use interp_support
    
    implicit none
    integer, intent(in), optional :: id
    real(dp) :: m0,aj,mt,mc
    real(dp):: tscls(20),lums(10),GB(10),zpars(20)
    real(dp):: dtm,mcy,tm,tn,mt0
    integer :: kw,idd,j
    type(track), pointer :: t
    logical :: debug
    
    debug = .false.
    if (debug)write(UNIT=err_unit,fmt=*) 'in gntage',kw,m0,mt,mc,aj,id
!    write(UNIT=*,fmt=*) 'in gntage',kw,m0,mt,mc,aj,id

    if (kw>7 .and. use_sse_NHe)  then
        CALL SSE_gntage(mc,mt,kw,zpars,m0,aj,id)
        return
    endif

    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)

    dtm = 0.d0
    mcy = 0.d0
    mt0 = mt
    ! for very low mass stars MS extends beyond the age of the universe
    ! so higher phases are often missing
    ! below simplification avoids error in length
    if(mt<Mmin_array(TA_cHeB_EEP) .and. kw>1 .and. kw<7) then != very_low_mass_limit
        mt = Mmin_array(TA_cHeB_EEP)
!        print*,mt,mt0
    endif
    
    !this is just to signal star that gnatge is calling it
    !pars% phase will get updated to its correct value in star
!    if (t% pars% phase<=kw) t% pars% phase = kw+1
    t% star_type = rejuvenated

    !TODO: provide a backup in case one of the mcrits are not defined

    if(kw.eq.4)then
    ! Set the minimum CHeB core mass (at BGB or He ignition)
    ! using M = Mflash/Mhef
         if(Mcrit(4)% loc >0) then
             j = min(sa(Mcrit(4)% loc)% ntrack,cHeIgnition_EEP)
             mcy = sa(Mcrit(4)% loc)% tr(i_he_core,j)
         endif
         if(mc.le.mcy) then
            kw = 3
            if (debug) WRITE(*,*)' GNTAGE4: changed to 3'
         endif
      endif

    !Next we check that we don't have a GB star for M => Mfgb
      if(kw.eq.3)then
        ! Set the maximum GB core mass using M = Mfgb
        if(Mcrit(5)% loc >0) then
            j = min(sa(Mcrit(5)% loc)% ntrack,cHeIgnition_EEP)
            mcy = sa(Mcrit(5)% loc)% tr(i_he_core,j)
        endif
        if(mc.ge.mcy)then
            kw = 4
            aj = 0.d0
           if (debug) WRITE(*,*)' GNTAGE3: changed to 4'
         endif
      endif

    select case(kw)
    case(8:9)
        if (debug) WRITE(*,*)'he stars'

        if (t% is_he_track) then
            t% pars% age = t% times(He_MS)
        else
            t% pars% age = t% times(HeBurn)
            t% star_type = switch
        endif
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tm + 1.0d-10*tm
        
    case (7)
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        if (debug)write(UNIT=err_unit,fmt=*) '7 in gntage',kw,m0,mt,mc,aj,id
        aj = t% pars% age
    case(6)
        ! We try to start the star from the start of the SAGB.
        if (debug) WRITE(*,*)'TPAGB'

        if (t% is_he_track) then
            t% pars% age = t% times(He_HG)
            t% star_type = switch
        else
            t% pars% age = t% times(EAGB)
            if (t% pars% age<0.d0) t% pars% age = maxval(t% times(1:5))
            
        endif
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(13)
    case(5)
        ! We fit a Helium core mass at the base of the AGB.
        if (debug) WRITE(*,*)'AGB'

        if (t% is_he_track) then
            t% pars% age = t% times(He_MS)
            t% star_type = switch
        else
            t% pars% age = t% times(HeBurn)
            if (t% pars% age<0.d0) t% pars% age = maxval(t% times(1:4))

        endif

        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(2) + tscls(3)
    case(4)
    
        if (debug) WRITE(*,*)'chb'

        ! The supplied age is actually the fractional age, fage, of CHeB lifetime
        ! that has been completed, ie. 0 <= aj <= 1.
        if(aj.lt.0.d0) aj = 0.d0
        if(aj.gt.1.d0) aj = 1.d0

        if (t% is_he_track) then
            t% pars% age = t% times(He_MS)+ aj*(t% times(He_MS))
            t% star_type = switch
        else
            t% pars% age = t% times(RGB) + aj*(t% times(HeBurn)-t% times(RGB))
        endif
        
        !for stars that don't have RGB phase, times(RGB) corresponds to times(HG)
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(2) + aj*tscls(3)
                
        
    case(3)
        if (debug) WRITE(*,*)'rgb'

        !Place the star at the BGB
        if (t% is_he_track) then
            t% pars% age = 0.d0
            t% star_type = switch
        else
            t% pars% age = t% times(HG)
        endif
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(1) + 1.0d-06*(tscls(2) - tscls(1))
    end select
    
    t% pars% age = aj
    t% pars% phase = kw
    mt = mt0
    m0 = t% zams_mass

    nullify(t)
    
    if (debug)write(UNIT=err_unit,fmt=*)'exit gntage',kw,m0,mt,mc,aj
!    write(UNIT=*,fmt=*)'exit gntage',kw,m0,mt,mc,aj

end subroutine METISSE_gntage

