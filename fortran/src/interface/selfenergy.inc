      interface
      subroutine selfenergy(job1,ikpt,eval,n,coul,Win(ctrafo),Win(screen), head, wing,
     &                                                        Win(screenc),headc,wingc,plasma)

      use global !inc
      use arrays !inc

      implicit none
      type(jobtype), intent(in)           :: job1
      integer,       intent(in)           :: n,ikpt
      logical,       intent(in)           :: eval(2)
      real_dp,       intent(in), optional :: coul(n)
      complex_dp,    intent(in), optional :: headc(3,3,nfreqc),wingc(3,n,2,nfreqc)
      MCOMPLEX_dp,   intent(in), optional :: head(6,nfreq),wing(3,n,nfreq)
      MCOMPLEX_dp,               optional :: ctrafo(nbasm(ikpt),n)   ! can be modified but leaves subroutine unchanged
      MCOMPLEX_dp,               optional :: screen(n*(n+1)/2,nfreq) ! modified by Gamma
      complex_dp,                optional :: screenc(n,n,nfreqc)     ! modified by Gamma
      real_dp,                   optional :: plasma                  ! changed if plasma=-1d0 (PLASMA METAL)
# ifdef MPI
      integer,       intent(in), optional :: win_screen,win_screenc,win_ctrafo
# endif
      end subroutine selfenergy
      end interface
