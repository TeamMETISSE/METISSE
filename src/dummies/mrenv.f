      SUBROUTINE mrenv(kw,mass,mt,mc,lum,rad,rc,aj,tm,ltms,lbgb,lhei,
     & rzams,rtms,rg,menv,renv,k2)
        use track_support, only: dp
        !dummy subroutine

        integer :: kw
        real(dp):: mass,mt,mc,lum,rad,rc,aj,tm
        real(dp):: ltms,lbgb,lhei,rzams,rtms,rg,menv,renv,k2

        menv = (mt - mc)
        renv = (rad - rc)
        k2 = 0.21
          
      END SUBROUTINE mrenv
