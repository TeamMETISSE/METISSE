subroutine assign_commons_main()
    use track_support
    use remnant_support
    implicit none
    ! assign values for variables used in METISSE from SSE_input_controls
    
    REAL(dp) :: pts1,pts2,pts3
    COMMON /POINTS/ pts1,pts2,pts3
        
    if (front_end == main) then
       if (WD_mass_scheme == 'Mestel') then
            wd_flag = Mestel
        elseif (WD_mass_scheme == 'Modified_mestel') then
            wd_flag = Modified_mestel
        else
            print*,"METISSE error: Invalid Option for WD_mass_scheme."
            print*,"Choose from 'Mestel' and 'Modified_mestel'. "
            call stop_code
        endif

        if (BHNS_mass_scheme == 'original_SSE') then
            ns_flag = original_SSE
        elseif (BHNS_mass_scheme == 'Belczynski2002') then
            ns_flag = Belczynski2002
        elseif (BHNS_mass_scheme == 'Belczynski2008') then
            ns_flag = Belczynski2008
        elseif (BHNS_mass_scheme == 'Eldridge_Tout2004') then
            ns_flag = Eldridge_Tout2004
        else
            print*,"METISSE error: Invalid Option for BHNS_mass_scheme."
            print*,"Choose from 'original_SSE','Belczynski2002','Belczynski2008','Eldridge_Tout2004'. "
            call stop_code
        endif
        !call cutoffs_for_Belzynski_methods(ns_flag,mc1,mc2)

        !if (verbose) print*,"For Belzynski methods", mc1,mc2
        if (allow_electron_capture) ec_flag = 1
        if (Use_Initial_final_mass_relation) if_flag = 1
        
        !timesteps
        pts1 = pts_1
        pts2 = pts_2
        pts3 = pts_3
        
    else
        print*,'METISSE error: Front end mismtach in assign_commons_main'
        print*,'expected 0 (main); got ', front_end
    endif

    end subroutine

