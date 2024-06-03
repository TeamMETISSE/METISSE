module track_support

    !some subroutines of this module have been adpated from iso_eep_support module of ISO package (Dotter 2016).

    implicit none
    integer, parameter :: min_io_unit = 29
    integer, parameter :: max_io_unit = 99
    logical :: assigned(max_io_unit) = .false.
    
    !----from mesa const_def.f90
    ! real number precision options: single, double
    integer, parameter :: sp = selected_real_kind(p=5)
    integer, parameter :: dp = selected_real_kind(p=15)

    integer, parameter :: strlen = 256 ! for character (len=strlen)
    integer, parameter :: col_width = 32

    real(dp), parameter :: ln10 = log(1.0d1)
    real(sp), parameter :: ln10_sp = log(10.0)
    real(dp), parameter :: tiny = 1.0d-6
    real(dp), parameter :: undefined  =  -1.d0
    integer, parameter :: undefined_i = -1
    
    logical :: verbose, use_sse_NHe
    logical :: write_track_to_file, write_eep_file, write_error_to_file
    integer :: err_unit

    integer :: front_end = -1
    integer, parameter :: main = 0
    integer, parameter :: BSE = 1
    integer, parameter :: COSMIC = 2

    character(len=strlen) :: METISSE_DIR,TRACKS_DIR,TRACKS_DIR_HE

    ! for use when constructing EEP distance
    logical :: weight_center_rho_T_by_Xc
    real(dp) :: Teff_scale=2d0
    real(dp) :: logL_scale=0.125d0
    real(dp) :: age_scale=0.05d0
    real(dp) :: Rhoc_scale=1d0
    real(dp) :: Tc_scale=1d0

    !stellar types for handling primary eeps
    integer, parameter :: unknown           =  1 !for initialization only
    integer, parameter :: sub_stellar       =  2 !no fusion = brown dwarf
    integer, parameter :: star_low_mass     =  3 !ends as a WD  (applies to both H and He tracks)
    integer, parameter :: star_high_mass    =  4 !does not end as a WD (applies to both H and He tracks)
!    integer, parameter :: post_agb          =  5
    integer, parameter :: remnant           =  6
    integer, parameter :: sse_he_star       =  7    ! special type for He star evolved using sse formulae
    integer, parameter :: rejuvenated       =  8
    integer, parameter :: switch            =  9     ! temporary type when switching from hydrogen to helium stars
    
    
    character(len=10) :: star_label(4) = ['   unknown', 'substellar', '  low-mass', ' high-mass']
    character(len=5) :: phase_label(16) = ['lm_MS','   MS','   HG','  FGB',' CHeB',' EAGB',&
    'TPAGB','He_MS', 'He_HG','He_GB','He_WD','CO_WD','ONeWD','   NS','   BH','   MR']
    
    real(dp) :: T_bgb_limit = 3.85
    real(dp) :: very_low_mass_limit = 0.75d0 !Msun

    !sse phases

    integer, parameter :: low_mass_MS = 0
    integer, parameter :: MS = 1
    integer, parameter :: HG = 2
    integer, parameter :: RGB = 3
    integer, parameter :: HeBurn = 4
    integer, parameter :: EAGB = 5
    integer, parameter :: TPAGB = 6
    integer, parameter :: He_MS = 7
    integer, parameter :: He_HG = 8
    integer, parameter :: He_GB = 9
    integer, parameter :: HeWD = 10
    integer, parameter :: CO_WD = 11
    integer, parameter :: ONeWD = 12
    integer, parameter :: NS = 13
    integer, parameter :: BH = 14
    integer, parameter :: Massless_REM = 15

    !EEPs

    integer :: PreMS_EEP
    integer :: ZAMS_EEP
    integer :: IAMS_EEP
    integer :: TAMS_EEP
    integer :: BGB_EEP
    integer :: cHeIgnition_EEP
    integer :: cHeBurn_EEP
    integer :: TA_cHeB_EEP

    integer :: cCBurn_EEP
    integer :: TPAGB_EEP
    integer :: post_AGB_EEP

    
    integer :: Extra_EEP1
    integer :: Extra_EEP2
    integer :: Extra_EEP3

    ! for he stars
    
    integer :: ZAMS_HE_EEP
    integer :: TAMS_HE_EEP
    integer :: GB_HE_EEP
    integer :: TPAGB_HE_EEP
    integer :: cCBurn_HE_EEP
    integer :: post_AGB_HE_EEP
    
    integer :: Initial_EEP, Initial_EEP_HE
    integer :: Final_EEP, Final_EEP_HE
    integer :: low_mass_final_eep, high_mass_final_eep
    integer :: low_mass_eep_he, high_mass_eep_he

    integer, allocatable :: key_eeps(:),key_eeps_he(:)
    
    !quantities from history file that are needed directly in the code

    character(len=col_width) :: age_colname, mass_colname, log_L_colname,log_T_colname, &
                                log_R_colname, he_core_mass,co_core_mass, &
                                log_Tc,c12_mass_frac,o16_mass_frac, he4_mass_frac, &
                                Lum_colname,Teff_colname,Radius_colname, &
                                he_core_radius, co_core_radius, mass_conv_envelope, &
                                radius_conv_envelope, moment_of_inertia

    integer :: i_age, i_age2, i_mass, i_logTe, i_logL, i_logR, i_he_core, i_co_core
    integer :: i_RHe_core,i_RCO_core,i_mcenv, i_Rcenv,i_MoI
    integer :: i_he_RCO,i_he_mcenv, i_he_Rcenv,i_he_MoI,i_he_age

    integer :: i_Tc, i_he4, i_c12,i_o16
    integer :: i_Xc, i_Yc, i_Cc,i_Rhoc, i_gamma, i_surfH

    integer :: number_of_core_columns
    integer, allocatable :: core_cols(:), core_cols_he(:)
    !for columns
    integer, parameter :: max_col = 180
    integer, parameter :: column_int=0
    integer, parameter :: column_dbl=1
!    character(len=strlen) :: extra_core_columns_names        !TODO: make it flexible

    type column
     character(len=col_width) :: name
     integer :: type, loc
    end type column

    !EEP arrays
    integer, parameter :: primary = 10 ! number of primary EEPs !TODO: --change this
    ! as set by primary_eep
    integer :: eep_interval(primary-1) ! number of secondary EEPs
    ! between the primaries

    real(dp), allocatable :: t_incomplete(:), t_notfound(:)
    real(dp), allocatable :: Mmax_array(:), Mmin_array(:)
    real(dp), allocatable :: Mmax_he_array(:), Mmin_he_array(:)

  !holds an evolutionary track for input, use an array of these for multiple tracks

    type eep_track
        character(len=strlen) :: filename
        type(column), allocatable :: cols(:)

        logical :: has_phase = .false., ignore=.false.
        logical :: has_mass_loss, is_he_track
        integer :: ncol, ntrack, neep
        integer :: star_type = unknown

        integer, allocatable :: eep(:), phase(:)
        real(dp) :: initial_mass, initial_Z, initial_Y, Fe_div_H,  v_div_vcrit, alpha_div_Fe
        real(dp), allocatable :: tr(:,:)

    end type eep_track

    !holds current parameters of star-- used by track
    type star_parameters
        integer :: phase,extra
        real(dp) :: mass,core_mass,core_radius, McHe, McCO
        real(dp) :: luminosity,Teff,radius
        real(dp) :: log_L,log_Teff,log_R                !log values
        real(dp) :: epoch, age, age_old,age2
        real(dp) :: delta, dt, dms, mcenv, rcenv,moi,bhspin
    end type star_parameters
    

    !holds values of agb parameters for constructing AGB to WD track
    type agb_parameters
        real(dp) :: mass,radius,lum
        real(dp) :: tfinal,tini
        integer :: phase_wd
    end type agb_parameters

    !holds values of SSE related parameters
    type sse_parameters
        real(dp) :: D,Mx,Lx,LtMS
        real(dp) :: Rzams, Lzams      !zams values
    end type

    !holds interpolated track
    type track
        logical :: complete = .true., post_agb = .false.
        logical :: has_mass_loss = .false., is_he_track = .false.

        integer :: ncol, ntrack, neep,min_index,j_bgb,j_bgb0
        integer :: star_type = unknown, ierr = 0
        integer, allocatable :: eep(:), bounds(:)
        type(column), allocatable :: cols(:)

        real(dp) :: initial_mass, initial_Z, initial_Y, Fe_div_H,  v_div_vcrit, alpha_div_Fe
        real(dp) :: initial_mass_old,initial_mass0,zams_mass
        
        ! initial_mass is the initial mass of the track
        ! zams_mass is the effective initial mass (M0/mass0 of SSE)
        ! pars% mass contains current total mass of the star (mt of sse)
        ! initial_mass0 is the initial_mass of the track for which mass at tams is zams_mass
        ! initial_mass_old is the last initial mass of the track used for interpolation
        
        ! Initial_mass/Initial_mass_old is needed for surface parameters,
        ! while Initial_mass0 is for core parameters(except MS when initial_mass0 is enough)
        
        real(dp), allocatable :: tr(:,:)
        real(dp) :: times(11), times_new(11)           !timescales
        real(dp) :: MS_time, nuc_time, MS_old
        type(star_parameters) :: pars    ! parameters at any instant
        
        type(agb_parameters) :: agb   ! parameters of stars used when constructing post-agb to WD track
        
        type(sse_parameters) :: He_pars   !parameters of naked helium stars, used when using SSE formulae
    end type track
    
    !defining array for input tracks
    type(eep_track), allocatable, target :: sa(:), sa_he(:)
    type(track), allocatable, target :: tarr(:)
    real(dp) :: initial_Z

    !variable declaration-- for main
    integer :: number_of_tracks
    character(len=strlen) :: input_mass_file
    logical :: read_mass_from_file
    
    !for z_support
    real(dp) :: Mhook, Mhef,Mfgb, Mup, Mec, Mextra,Mup_core,Mec_core
    real(dp) :: Z04, Z_H, Z_He
    integer, allocatable :: m_cutoff(:), m_cutoff_he(:)

    type critical_mass
        integer :: loc
        real(dp) :: mass
    end type critical_mass

    type(critical_mass) :: Mcrit(9), Mcrit_he(9)
    
    !for interp_support
    logical :: fix_track
    real(dp) :: lookup_index, mass_accuracy_limit
    
    !for remnant support
    real(dp), parameter :: M_ch = 1.44d0

    !in case of direct call
    real(dp) :: max_NS_mass         !maximum NS mass
    logical :: construct_wd_track, allow_electron_capture, use_Initial_final_mass_relation
    character (len=strlen) :: BHNS_mass_scheme, WD_mass_scheme
!    real(dp) :: mc1, mc2 !mass cutoffs for Belczynski methods
    real(dp) :: pts_1,pts_2,pts_3
    
    contains

    !linear search alogorithm
    subroutine index_search(size_list,list,value,min_index,debug)
        integer, intent(in) :: size_list
        real(dp), intent(in) :: list(:)
        real(dp), intent(in) :: value
        integer, intent(out) :: min_index
        logical, optional :: debug

        if (present(debug)) then
            if (debug)  then
    !            print*, "in index search",value,size_list
            print*, "in index search; list(1),list(size_list-1),list(size_list)"
            print*, list(1),list(size_list-1),list(size_list)
            endif
        endif
        
        if (size(list)<1) print*, 'error in list size',size(list),size_list
        if (value < list(1)) then             !from num_binary_search.inc
            min_index = 1; return
        elseif (check_equal(value, list(size_list)))then
            min_index = size_list; return
        elseif (value > list(size_list)) then
            min_index = size_list+1; return
        end if

        min_index = minloc(abs(list-value), dim=1)
    end subroutine index_search


    ! from ISO (Dotter et al. 2016), adapted from MESA; modified by PA to avoid 0 loc

    ! if vec contains decreasing values,
    ! returns 1 if val > vec(1); returns n if val <= vec(n)
    ! else returns k between 1 and n-1 such that
    !     vec(k) >= val > vec(k+1)

    ! if vec contains increasing values,
    ! returns 0 if val < vec(1); returns n if val >= vec(n)
    ! else returns k between 1 and n-1 such that
    !     vec(k) <= val < vec(k+1)

    integer function binary_search(n, vec, val) result(loc)
     integer, intent(in) :: n
     real(dp), intent(in) :: val
     real(dp), intent(in) :: vec(:) ! (n)
    !     real(dp), parameter :: tiny = 1.0d-13
     integer :: first, last, mid
     real(dp) :: list(3)
     if (n <= 1) then
        loc = n; return
     end if

     if (vec(n) < vec(1)) then ! decreasing values

        if (val > vec(1)) then
           loc = 1; return
        else if (abs(val - vec(1)) < tiny ) then
           loc = 1; return
        else if (val <= vec(n)) then
           loc = n; return
        end if


        first = 1
        last = n-1
        loc = -1
        do while (first <= last)
           mid = (first + last)/2
           if (vec(mid) >= val) then
              if (val > vec(mid+1)) then
                 loc = mid
                 exit
              end if
              first = mid + 1
           else
              last = mid - 1
           end if
        end do

     else ! increasing values
     
!        print*, 'test', vec(1),vec(n),val
        if (val < vec(1)) then
           loc = 1; return
        else if (abs(val - vec(1)) < tiny) then
           loc = 1; return
        else if (check_equal(val,vec(n))) then
            loc = n; return
        else if (val >= vec(n)) then
           loc = n; return
        end if

        first = 1
        last = n-1
        loc = -1
        do while (first <= last)
           mid = (first + last)/2
           if (vec(mid) <= val) then
              if (val < vec(mid+1)) then
                 loc = mid
                 exit
              end if
              first = mid + 1
           else
              last = mid - 1
           end if
        end do
     end if

    
        if (loc>1 .and. loc<n) then
            list(1) = vec(loc-1)
            list(2) = vec(loc)
            list(3) = vec(loc+1)
!            print*, 'list',list,val
            loc = loc+minloc(abs(list-val), dim=1)-2
!            print*,'loc',minloc(abs(list-val), dim=1)
        endif
    end function binary_search
    
    subroutine stop_code()
        print*, 'Error encountered; stopping METISSE'
        print*, 'See error file (fort.99) or terminal for details'
        STOP
    end subroutine stop_code
    
    subroutine write_eep_track(x,mt,filename)
    !modified from ISO subroutine of same name

    type(track), pointer, intent(in) :: x
    character(len=strlen), intent(in), optional :: filename
    real(dp), intent(in), optional :: mt
    character(len=strlen) :: eep_filename
    integer :: str, io, ierr, j! ncol
    real(dp) :: real_phase
    character(len=8) :: have_phase
    real(dp), allocatable :: phase(:)
    io = alloc_iounit(ierr)

    if (present(filename)) then
        eep_filename = trim(filename)
    elseif (present(mt)) then
        str = int(mt*100)
        write(eep_filename,"(a,a,i5.5,a)") trim(METISSE_DIR),"/output_eep/",str,"M.track.eep"
    else
        print*, 'ERROR: NO EEP FILE WRITTEN, either provide FILENAME or mass of the star'
        return
    ENDIF
    
    if (verbose) print*,'writing',str,'M.track.eep'
    
    call calculate_sse_phases(x,phase)
    
    open(io,file=trim(eep_filename),action='write',status='unknown')
    have_phase = 'YES'

    !write(io,'(a25,a8)') '# MIST version number  = ', x% version_string
    !write(io,'(a25,i8)') '# MESA revision number = ', x% MESA_revision_number
    !                     123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    ! write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a88)') '#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit                                        '
    write(io,'(a2,f6.4,1p1e13.5,0p3f9.2)') '# ', x% initial_Y, x% initial_Z, x% Fe_div_H, x% alpha_div_Fe, x% v_div_vcrit
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a1,1x,a16,4a8,2x,a10)') '#','initial_mass', 'N_pts', 'N_EEP', 'N_col', 'phase', 'type'
    write(io,'(a1,1x,1p1e16.10,3i8,a8,2x,a10)') '#', x% initial_mass, x% ntrack, x% neep, x% ncol+2, have_phase, &
         star_label(x% star_type)
    write(io,'(a8,20i8)') '# EEPs: ', x% eep
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a1,299(27x,i5))') '#', (j,j=1,x% ncol+2)
    write(io,'(a1,299a32)') '#', adjustr(x% cols(:)% name), 'phase'
    
    do j=x% eep(1),x% ntrack
      real_phase = real(phase(j))
      write(io,'(1x,299(1pes32.16e3))') x% tr(:,j), real_phase
    enddo
       
       
    deallocate(phase)
    close(io)
    call free_iounit(io)
    end subroutine write_eep_track

    subroutine calculate_sse_phases(t,phase)
        !subroutine to assign sse phases to the interpolated track
        implicit none
        
        real(dp), allocatable :: phase(:)

        type(track), pointer,intent(in) :: t
        integer :: i,j_bgb
        real(dp) :: mass
        
        logical :: debug
        debug = .false.


        mass = t% initial_mass
        allocate(phase(t% ntrack))

        do i = 1, t% neep

        if (t% eep(i) == PreMS_EEP) then   !pre_MS
            phase(PreMS_EEP:ZAMS_EEP-1) = -1

        elseif (t% eep(i) == TAMS_EEP) then    !MS
            if (mass< Mhook-0.3) then
                phase(ZAMS_EEP:TAMS_EEP) = low_mass_MS
            else
                phase(ZAMS_EEP:TAMS_EEP) = MS
            endif

        elseif (t% eep(i) == cHeIgnition_EEP) then
            phase(TAMS_EEP+1: cHeIgnition_EEP) = HG          !Herztsrung gap

        elseif (t% eep(i) == TA_cHeB_EEP) then
            phase(cHeIgnition_EEP+1:TA_cHeB_EEP) = HeBurn                !red_HB_clump /core He Burning

        elseif (t% eep(i) == TPAGB_EEP) then
            phase(TA_cHeB_EEP+1:TPAGB_EEP) = EAGB                !AGB :massive stars' evolution ends here

        elseif (t% eep(i) == post_AGB_EEP) then
            phase(TPAGB_EEP+1:post_AGB_EEP) = TPAGB              !TP-AGB :only for low_inter mass stars

        elseif (t% eep(i) == cCBurn_EEP) then
            phase(TA_cHeB_EEP+1: cCBurn_EEP) = EAGB
        endif
        enddo

        !determine the base of the giant branch times, if present
        if (mass > very_low_mass_limit .and. mass< Mfgb) then
            if (identified(BGB_EEP)) then
                phase(BGB_EEP: cHeIgnition_EEP) = RGB           !Red giant Branch
            elseif (t% ntrack >TAMS_EEP) then
                j_bgb =  base_GB(t)
                if (j_bgb>0) then
                    phase(j_bgb: cHeIgnition_EEP) = RGB           !Red giant Branch
                else
                    write(UNIT=err_unit,fmt=*) "Unable to locate BGB ", j_bgb, t% initial_mass
                end if
            endif
        endif
    end subroutine
    
    integer function bgb_mcenv(t,jstart,jend) result(j_bgb)
        type(track), intent(in) :: t
        integer :: j, jstart,jend

        j_bgb = -1
!        print*, 'teff',s% tr(i_logTe,jend)
        if ((t% is_he_track .eqv. .false. ).and.t% tr(i_logTe,jend)> T_bgb_limit) return

        do j = jstart,jend
            if (t% tr(i_mcenv,j)/t% tr(i_mass,j).ge.0.12d0) then
                j_bgb = j+ jstart-1
                exit
            endif
        enddo
    end function
    
    integer function base_GB(t) result(j_bgb)
        type(track), pointer,intent(in) :: t
        integer :: peak, jfinal, jini, j_diff,k

        real(dp) :: mass,l_calc,diff
        real(dp), allocatable ::Lum(:),Teff(:),core_mass(:)
        real(dp), allocatable ::diff_L(:),diff_Te(:),dLdTe(:)

        j_bgb = -1
        
        if (t% is_he_track) then
            jini = TAMS_HE_EEP
            jfinal = Final_EEP_HE
        else
            jini = TAMS_EEP
            jfinal = min(cHeIgnition_EEP,cHeBurn_EEP)
        endif
        
        j_diff = jfinal - jini

        if ((t% tr(i_logTe,jfinal)> T_bgb_limit) .or. j_diff<=0) return

        allocate(Lum(j_diff),Teff(j_diff),core_mass(j_diff))
        allocate (diff_L(j_diff -1),diff_Te(j_diff -1),dLdTe(j_diff -1))

        
        mass = t% initial_mass
        Lum =  t% tr(i_logL,jini:jfinal)
        Teff = t% tr(i_logTe,jini:jfinal)            !TODO: make this wrt radius
        core_mass = t% tr(i_he_core,jini:jfinal)
        
        
        diff_L = Lum(2:j_diff)-Lum(1:j_diff-1)
        diff_Te = Teff(2:j_diff)-Teff(1:j_diff-1)
        dLdTe = diff_L/diff_Te

        peak = maxloc(dLdTe,dim = 1)
!        print*,"peak = ",mass,peak, Teff(peak),j_diff

        if(peak>=j_diff) then
            deallocate(diff_L,diff_Te,dLdTe)
            deallocate(Lum,Teff,core_mass)
            return
        endif
        
        if (dLdTe(peak)>0.d0 .and. dLdTe(peak+1)<0.d0) then         !checking for oscillations
            if (peak>1) peak = maxloc(dLdTe(:peak-1),dim = 1)
        endif
        
        if (dLdTe(peak)>0.d0) then          !convective core
            do k = peak,j_diff-1
                if (dLdTe(k)<1E-12 .and. Teff(k)<=T_bgb_limit) then
                    j_bgb = k+TAMS_EEP-1
                    exit
                endif
            end do
        elseif (mass<Mhook) then            !radiative core
            do k = 1,j_diff-1
                l_calc = log10(2.3d+5*(core_mass(k)**6))
                diff = l_calc-Lum(k)
                if (abs(diff)<0.12) then
                    j_bgb = k+TAMS_EEP-1
                    exit
                endif
            end do
        end if

        deallocate(diff_L,diff_Te,dLdTe)
        deallocate(Lum,Teff,core_mass)
        
    end function
    
    subroutine write_dat_track(tphys, pars)
        real(dp), intent(in) :: tphys
        type(star_parameters), intent(in) :: pars
        character(LEN=*), PARAMETER  :: FMT= '(1p9e15.6,2i10)'
        write(120,FMT) tphys,pars% age,pars% mass,pars% core_mass,pars% McHe, pars% McCO &
                    ,pars% log_L,pars% log_Teff,pars% log_R,pars% phase ,pars% extra
    end subroutine write_dat_track

!    subroutine alloc_track(filename,x)
!        character(len=strlen), intent(in) :: filename
!        type(eep_track), pointer :: x
!        allocate(x)
!        x% neep = primary
!        x% filename = trim(filename)
!        allocate(x% eep(x% neep))
!      end subroutine alloc_track

    subroutine distance_along_track(t)
      type(track), intent(inout) :: t
      real(dp) :: tmp_dist, weight, max_center_h1
      integer :: j

      if(weight_center_rho_T_by_Xc)then
         max_center_h1 = maxval(t% tr(i_Xc,:))
         if(max_center_h1 <= 0d0) max_center_h1 = 1d0
      else
         max_center_h1 = 1d0
         weight = 1d0
      endif

!      t% dist(1) = 0d0
      if(t% ntrack > 3)then
         do j = 2, t% ntrack
            
            if(weight_center_rho_T_by_Xc)then
               weight = max(0d0, t% tr(i_Xc,j)/max_center_h1)
            endif
            
            !build up the distance between EEPs piece by piece
            tmp_dist =            Teff_scale*sqdiff(t% tr(i_logTe,j) , t% tr(i_logTe,j-1))
            tmp_dist = tmp_dist + logL_scale*sqdiff(t% tr(i_logL, j) , t% tr(i_logL, j-1))
            tmp_dist = tmp_dist + weight * Rhoc_scale * sqdiff(t% tr(i_Rhoc, j) , t% tr(i_Rhoc, j-1))
            tmp_dist = tmp_dist + weight * Tc_scale*  sqdiff(t% tr(i_Tc,   j) , t% tr(i_Tc,   j-1))
            tmp_dist = tmp_dist + age_scale* sqdiff(log10(t% tr(i_age,j)) , log10(t% tr(i_age,j-1)))

!            t% dist(j) = t% dist(j-1) + sqrt(tmp_dist)
         enddo
      endif
    end subroutine distance_along_track

    elemental function sqdiff(x0,x1) result(y) !square of x, y=x*x
      real(dp), intent(in) :: x0, x1
      real(dp) :: y, dx
      dx=x0-x1
      y = dx*dx
    end function sqdiff
    
    elemental function pow10_sg(x) result(y)
        real(sp), intent(in) :: x
        real(sp) :: y
        y = exp(ln10_sp*x)
    end function pow10_sg

    elemental function pow10(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = exp(ln10*x)
    end function pow10
    
    !from mesa/utils_lib.f
    integer function alloc_iounit(ierr)
        !use utils_def
        integer, intent(out) :: ierr
        integer :: i
        ierr = 0
        alloc_iounit = -1
        do i = min_io_unit, max_io_unit
            if (.not. assigned(i)) then
                assigned(i) = .true.
                alloc_iounit = i
                exit
            end if
        end do
        if (alloc_iounit == -1) then
            ierr = -1
        end if
    end function alloc_iounit

    subroutine free_iounit(iounit)
        !use utils_def
        integer, intent(in) :: iounit
        logical :: bad_iounit
        bad_iounit = .false.
        !$omp critical (utils_alloc_io_unit)
        if (iounit >= min_io_unit .and. iounit <= max_io_unit) then
            assigned(iounit) = .false.
        else
            bad_iounit = .true.
        end if
        !$omp end critical (utils_alloc_io_unit)
        if (bad_iounit) then
            write(*,*) 'called free_iounit with invalid arg', iounit
            stop 'free_iounit'
        end if
    end subroutine free_iounit
    
    
    logical function check_ge(x,y) result(z)
    !TODO: needs to be checked before use
    real(dp), intent(in) :: x,y
        if (x.ge.y .or. abs(x-y)<tiny) then
            z = .true.
        else
            z = .false.
        endif
        return
    end function check_ge
    
    logical function check_le(x,y) result(z)
    !TODO: needs to be checked before use
    real(dp), intent(in) :: x,y
        if (x.le.y .or. abs(x-y)<tiny) then
            z = .true.
        else
            z = .false.
        endif
        return
    end function check_le
    
    logical function check_equal(x,y,limit) result(z)
    real(dp), intent(in) :: x,y
    real(dp), intent(in), optional :: limit
    real(dp) :: diff, threshold

    if (present(limit)) then
        threshold = limit
        else
    threshold = 1d-4
    endif
        diff = abs(x-y)
        if (diff.le.threshold) then
            z = .true.
        else
            z = .false.
        endif
        return
    end function check_equal

    logical function defined(x) result(y)
        real(dp), intent(in) :: x
        if (abs(x-undefined).le.tiny) then
            y = .false.
        else
            y = .true.
        endif
    return
    end function defined

    !same as the function 'defined' above but for integers
    ! used for checking if locations are identified
    logical function identified(x) result(y)
        integer, intent(in) :: x
        if (abs(x-undefined_i).le.0.d0) then
            y = .false.
        else
            y = .true.
        endif
        return
    end function identified
    
    subroutine uniform_distribution(n, minval, maxval, marray)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: minval, maxval
        real(dp), intent(out) :: marray(n)
        real(dp) :: h
        integer :: i

        marray = 0.d0
        !linearly spaced
        h = abs(maxval- minval)/(n-1)
        do i= 1,n
            marray(i) = minval + (i-1)*h
        end do
        
        !for log spaced
    end subroutine uniform_distribution

    !from COMPAS (Team Compas 2020)
    real(dp) function quadratic(a,b,c) result(x)
      real(dp),intent(in) :: a,b,c
      real(dp) :: D,sqrtD,x1,x2

        D = (B*B)-(4.d0*A*C)
        x = 0.0

        if (D< 0.0) then
            write(UNIT=err_unit,fmt=*)"fatal error: Non-real roots"
            call stop_code
        else if (D > 0.0) then
            sqrtD = sqrt(D)
            x1 = (-B + sqrtD)/(2*A)
            x2 = (-B - sqrtD)/(2*A)
            x= max(x1, x2)
        else
            x = -B/(2*A)
        endif
    end function quadratic
    
    subroutine deallocate_arrays(t)
    type(track), pointer :: t
        deallocate(t% eep)
        deallocate(t% tr)
!        deallocate(t% times)
        deallocate(t% cols)

    end subroutine deallocate_arrays
    
    integer function get_min_ntrack(star_type,is_he_track)
    integer, intent(in):: star_type
    logical, intent(in):: is_he_track
    
        get_min_ntrack = 0
        !calculating min required length to the new track
        !use Mec here?
        
        if (is_he_track) then
            if (star_type == star_high_mass) then
                get_min_ntrack = high_mass_eep_he
            else
                get_min_ntrack = low_mass_eep_he
            endif
        else
            if (star_type == star_high_mass) then
                get_min_ntrack = high_mass_final_eep
            else
                get_min_ntrack = low_mass_final_eep
            endif
        endif
        
    end function
end module track_support
