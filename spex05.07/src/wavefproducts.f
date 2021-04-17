c WAVE-FUNCTION PRODUCTS  <M phi | phi >
c
c wavefproducts1_mt : Muffin-tin part, for susceptibility
c wavefproducts2_mt : Muffin-tin part, for selfenergy
c wavefproducts3_mt : Muffin-tin part, for Wannier
c wavefproducts1_pw : Plane-wave part, for susceptibility
c wavefproducts2_pw : Plane-wave part, for selfenergy
c wavefproducts3_pw : Plane-wave part, for Wannier
c
c --------------------------------------------------------------------------------

c The parallelization assumes ik1,b1 (respectively k1,b1low,b1up) to be identical over the Mcomm processes (but not b2 etc.).
c It requires a relatively large communication overhead.
c Comment the following to avoid parallelization.
# ifdef MPI
#   define PARA_MT
#   define PARA_PW
# endif

# ifndef SUBTYPE
#   define SUBTYPE 1
#   include "wavefproducts.f"
#   undef SUBTYPE
#   define SUBTYPE 2
#   include "wavefproducts.f"
#   define SUBWANNIER
# endif

# ifdef SUBWANNIER
#   ifdef PARA_MT
#     undef PARA_MT
#   endif
#   ifdef PARA_PW
#     undef PARA_PW
#   endif
# endif

# ifndef noWARN
# if SUBTYPE == 1 && defined(MPI)
#   ifdef PARA_MT
#     warning MT products will be parallelized.
#   else
#     warning MT products won't be parallelized.
#   endif
#   ifdef PARA_PW
#     warning PW products will be parallelized.
#   else
#     warning PW products won't be parallelized.
#   endif
# endif
# endif

# if SUBTYPE != 1 && SUBTYPE != 2
#   error "Unknown subtype"
# endif

c LOAD: cmt and cpw are loaded (different array structure)
c (0) Standard definitions
# define CMT1 CMT
# define CPW1 CPW
# define IKX ,ikx
# define ADJUST_KINDX_IF_NEEDED
# define UNDO_ADJUST_KINDX
# ifdef LOAD
c   (1) In case of SUBTYPE==1 (selfenergy), ik1,b1 wave functions are in cmtq,cpwq
#   if SUBTYPE == 1
#     undef CMT1
#     undef CPW1
#     undef IKX
#     define CMT1 cmtq
#     define CPW1 cpwq
#     define IKX
#   endif
c   (2) cmt,cpw spin index = 1
#   ifndef SUBWANNIER
#     define ISPIN1 1
#     define ISPIN2 1
#   endif
c   (3) macros for redefining kindx to adjust to different kindx ordering of ik2
#   if SUBTYPE == 2 && !defined(SUBWANNIER)
#     undef ADJUST_KINDX_IF_NEEDED
#     undef UNDO_ADJUST_KINDX
#     define ADJUST_KINDX_IF_NEEDED if(.not.storeibz) then ; ikx = kindx(ik2) ; kindx(ik2) = kindx(ik1) ; endif
#     define UNDO_ADJUST_KINDX      if(.not.storeibz) kindx(ik2) = ikx
#   endif
# endif

c Wannier Bloch functions
c (1) Replace cmt->cmtu and cpw->cpwu
# ifdef SUBWANNIER
#   define LWAN .true.
#   define WSUB _w
#   define CMT cmtu
#   define CPW cpwu
# else
#   define LWAN .false.
#   define WSUB
#   define CMT cmt
#   define CPW cpw
# endif
c (2) PW coefficients cannot be defined real for Wannier Bloch functions
# if defined SUBWANNIER && defined INV
#   define WSUFFIX _w
#   define INVW
#   undef INV
# else
#   define WSUFFIX
# endif

# include "cppmacro.h"

# ifdef PARA_MT
#   define ifMOD_ ifMOD
#   define endMOD_ endMOD
#   define Mpi_ Mpi
# else
#   define ifMOD_(arg)
#   define endMOD_
#   define Mpi_(arg)
# endif

c --------------------------------------------------------------------------------
c
c Muffin-tin part of wave-function products
c
c cprod = < M(k0) phi(b1,k1,s1) | phi(b2,k2,s2) >            ! k2 = k0+k1
c
c SUBTYPE = 1 : k1 fixed, for selfenergy
c SUBTYPE = 2 : k0 fixed, for susceptibility
c
c SUBTYPE     1                 2
c b1          b1(:nb1)          b1low:b1up
c b2          (b2up-b2low)/M+1  b2low:b2up     (striped storage M=Msize if MPI and not LOAD, otherwise M=1)
c k0          k0(:nk)           ik0
c k1          ik1               k1(:nk)
c s1          ispin             ispin1
c s2          ispin             ispin2
c

# if SUBTYPE == 1
      subroutine wavefproducts1_mt(cprod,dim,b1,nb1,ik1,ispin,k0,nk,b2low,b2up,b2inc)
# else
#   ifndef SUBWANNIER
      subroutine wavefproducts2_mt(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   else
      subroutine wavefproducts3_mt(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   endif
# endif
      use global
      use wrapper, only: macmat
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )
      implicit none
# if SUBTYPE == 1
      integer,     intent(in)    :: nb1,b1(nb1),ik1,ispin,nk,k0(nk),b2low,b2up,b2inc,dim
      MCOMPLEX_dp, intent(inout) :: cprod(dim, (b2up-b2low) / b2inc + 1 ,nb1,nk)
      integer                    :: ispin1,ispin2
      integer                    :: ik0,ib,ic1,isym
      logical                    :: trafo,trs
      complex_dp                 :: cdum
      complex_dp,  allocatable   :: dwgn1(:,:)
# else
      integer,     intent(in)    :: nk,k1(nk),ik0,ispin1,ispin2,b1low,b1up,b2low,b2up
      MCOMPLEX_dp, intent(inout) :: cprod(b2low:b2up,b1low:b1up,nk,nbasp)
      integer                    :: ik1,ib2
      complex_dp,  allocatable   :: cmt4(:,:,:)
# endif
      complex_dp,  allocatable   :: mat(:,:),mat1(:,:,:),mat2(:,:,:),mat0(:,:),cmt1(:),cmt2(:,:),cmt3(:,:)
      complex_dp                 :: cexp
      integer                    :: lmstart(-maxlcut:maxlcut,0:maxlcut),pnt(maxlmindx)
      integer                    :: ik,ik2,k2(nk),ib1,s,ikx
      integer                    :: itype,ieq,ic
      integer                    :: l,l1,ll,m,m1,mm,n,n1,nn,lln,ln,ln1,lm,lm1,llm,llm0,nx,nx1,maxlm InvC(maxlmm)
      logical                    :: def
      real_dp,     allocatable   :: integral(:,:,:)
      real_dp                    :: gaunt1(0:maxlcut,0:maxlcut,0:maxlcutm,-maxlcut:maxlcut,-maxlcutm:maxlcutm),gnt
      real_dp                    :: gaunt,intgrf
      integer                    :: kptsum
c      logical trf1,trf2
c      trf1 = .false.
c      trf2 = .false.

      if(nk==0.or.b2low>b2up) return

# if SUBTYPE == 1
      ispin1 = ispin
      ispin2 = ispin
      if(storeibz) then
        l = maxval(lcutp)
        allocate ( dwgn1(-l:l,-l:l) )
      endif
# else
      if(l_soc.and.(ispin1/=1.or.ispin2/=1)) Bug('Wrong spin indices (SOC).')
c      if(ispin1/=ispin2) call teststop('ispin1/=ispin2')
# endif

      ! Precalculate Gaunt coefficients
      gaunt1 = 0
      do ll = 0,maxlcutm
        do mm = -ll,ll
          do l = 0,maxlcut
            do m = -l,l
              do l1 = abs(l-ll),min(l+ll,maxlcut),2
                m1                   = mm + m
                gaunt1(l1,l,ll,m,mm) = gaunt(ll,l1,l,mm,m1,m)
              enddo
            enddo
          enddo
        enddo
      enddo

      ! k2 = q+k
# if SUBTYPE == 1
      trafo = .false.
# endif
      do ik = 1,nk
# if SUBTYPE == 1
        ik0    = k0(ik)
        k2(ik) = kptsum(ik0,ik1) ; if(storeibz.and.kptp(k2(ik))/=k2(ik)) trafo = .true.
# else
        ik1    = k1(ik)
        k2(ik) = kptsum(ik0,ik1)
# endif
      enddo

# if SUBTYPE == 1
      if(dim<nbasp) Bug('Dimension dim too small.')
      cprod(:nbasp,:,:,:) = 0
# else
      cprod = 0
# endif
      llm0  = 0
      ic    = 0
      do itype = 1,ntype
             maxlm  = sum ( [ ((2*l+1)*nindx(l,itype), l=0,lcutp(itype)) ] )
        Inv( maxlmm = sum ( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] ) )

        ! precalculate radial integrals (->integral1)
        def = .false.
 1      lln = 0 ; llm = 0 ; lm = 0 ; lm1 = 0
        do ll = 0,lcutm(itype)                         ; do nn = 1,nindxm(ll,itype) ; lln = lln + 1 ; ln  = 0
          do l = 0,lcutp(itype)                        ; do n  = 1,nindx(l,itype)   ; ln  = ln  + 1 ; ln1 = 0
            do l1 = abs(ll-l),min(ll+l,lcutp(itype)),2 ; do n1 = 1,nindx(l1,itype)  ; ln1 = ln1 + 1
              if(def) then
                integral(ln,ln1,lln) = intgrf( basm(:,nn,ll,itype) / rgrid(:,itype) *
     &                               ( bas1(:,n,l,itype,ispin1)*bas1(:,n1,l1,itype,ispin2) +
     &                                 bas2(:,n,l,itype,ispin1)*bas2(:,n1,l1,itype,ispin2) ) , itype )
              else
                llm = max(llm,lln)
                lm  = max(lm,ln)
                lm1 = max(lm1,ln1)
              endif
            enddo ; enddo
          enddo ; enddo
        enddo ; enddo
        if(.not.def) then
          def = .true.
          allocate ( integral(lm,lm1,llm) )
          goto 1
        endif

        ! start index
        lm = 0
        do l = 0,lcutp(itype)
          do m = -l,l
            lmstart(m,l) = lm
            lm           = lm + nindx(l,itype)
          enddo
        enddo

                         allocate ( cmt1(maxlm)   )
        if(l_soc)        allocate ( cmt3(maxlm,2) )
# if SUBTYPE == 2
        if(l_soc) then ; allocate ( cmt4(maxlm,b2low:b2up,2) )
        else           ; allocate ( cmt4(maxlm,b2low:b2up,ispin2:ispin2) )
        endif
# endif

        do ieq = 1,neq(itype)
          ic = ic + 1

# if SUBTYPE == 2
          cexp = exp ( -img * 2*pi * dot_product(kpt(:,ik0),cent(:,ic)) )
          do ik = 1,nk
            ik1 = k1(ik)
            ik2 = k2(ik)
#   ifndef old_trafo
            ADJUST_KINDX_IF_NEEDED
            do ib2 = b2low,b2up
              call wavefunction_mt WSUB (cmt4(:,ib2,:),maxlm,ic,ib2,ik2,ISPIN2)
            enddo
            UNDO_ADJUST_KINDX
#   else
            do ib2 = b2low,b2up
              if(storeibz.and.kptp(ik2)/=ik2) then
                if(l_soc) then
                  call waveftrafo_mt(cmt4(:,ib2,1),maxlm,ik2,ib2,1,ic,LWAN)
                  call waveftrafo_mt(cmt4(:,ib2,2),maxlm,ik2,ib2,2,ic,LWAN)
                  call waveftrafo_soc(cmt4(:,ib2,:),maxlm,symkpt(ik2))
                else
                  call waveftrafo_mt(cmt4(:,ib2,ispin2),maxlm,ik2,ib2,ISPIN2,ic,LWAN)
                endif
              else
                ikx = kindx(ik2) ; Load( if(.not.storeibz.and..not.LWAN) ikx = kindx(ik1) )
                if(l_soc) then ; cmt4(:,ib2,:)      = CMT(:maxlm,ic,ib2,ikx,:)
                else           ; cmt4(:,ib2,ispin2) = CMT(:maxlm,ic,ib2,ikx,ISPIN2)
                endif
              endif
            enddo
#   endif
            do ib1 = b1low,b1up
# else
          do ib = 1,nb1
            ib1 = b1(ib)
# endif

# ifndef old_trafo
#   if defined(LOAD) && SUBTYPE == 1
            if(storeibz.and.kptp(ik1)/=ik1) Error('Not implemented: LOAD & ik1 not in IBZ.')
            if(l_soc) then ; cmt3 = cmtq(:maxlm,ic,ib1,:)
            else           ; cmt1 = cmtq(:maxlm,ic,ib1,ISPIN1)
            endif
#   else
            if(l_soc) then ; call wavefunction_mt WSUB (cmt3,maxlm,ic,ib1,ik1,ISPIN1)
            else           ; call wavefunction_mt WSUB (cmt1,maxlm,ic,ib1,ik1,ISPIN1)
            endif
#   endif
# else
            if(storeibz.and.kptp(ik1)/=ik1) then
              if(l_soc) then
                call waveftrafo_mt(cmt3(:,1),maxlm,ik1,ib1,1,ic,LWAN)
                call waveftrafo_mt(cmt3(:,2),maxlm,ik1,ib1,2,ic,LWAN)
                call waveftrafo_soc(cmt3,maxlm,symkpt(ik1))
              else
                call waveftrafo_mt(cmt1,maxlm,ik1,ib1,ISPIN1,ic,LWAN)
              endif
            else
              ikx = kindx(ik1)
              if(l_soc) then ; cmt3 = CMT1(:maxlm,ic,ib1 IKX ,:)
              else           ; cmt1 = CMT1(:maxlm,ic,ib1 IKX ,ISPIN1)
              endif
            endif
# endif

c
c           Loop over MPB function

            llm = llm0
            lln = 0
            do ll = 0,lcutm(itype)
              do mm = -ll,ll

c
c               For each MPB function, calculate  <M phi(ib1) | basis >  ( -> mat, mat1 for SOC)

                allocate ( mat(maxlm,nindxm(ll,itype)) )
                if(l_soc) allocate ( mat1(maxlm,nindxm(ll,itype),2) )
                if(storeibz) then
                  if(l_soc) then ; allocate ( mat2(maxlm,nindxm(ll,itype),2) )
                  else           ; allocate ( mat0(maxlm,nindxm(ll,itype))   )
                  endif
                endif

                do s = 1,2 ! SOC spin loop
                if     (l_soc)     then ; cmt1 = cmt3(:,s)
                else if(s/=ispin1) then ; cycle
                endif

                mat = 0
                lm  = 0
                ln  = 0
                do l = 0,lcutp(itype)
                  nx = nindx(l,itype)
                  do m = -l,l

                    ifMOD_(l*(l+1)+m)

                    m1 = m + mm

                    ln1 = 0
                    do l1 = abs(l-ll),min(l+ll,lcutp(itype)),2
                      nx1 = nindx(l1,itype)

                      if(l1>=abs(m1)) then
                        gnt = gaunt1(l1,l,ll,m,mm)
                        if(gnt/=0) then
                          lm1 = lmstart(m1,l1)
                          do nn = 1,nindxm(ll,itype)
                            mat(lm1+1:lm1+nx1,nn) = mat(lm1+1:lm1+nx1,nn) +
     &                                              matmul(cmt1(lm+1:lm+nx),integral(ln+1:ln+nx,ln1+1:ln1+nx1,lln+nn)) * gnt
                          enddo
                        endif
                      endif

                      ln1 = ln1 + nx1
                    enddo ! l1

                    endMOD_

                    lm = lm + nx
                  enddo ! m
                  ln = ln + nx
                enddo ! l

                if(l_soc) mat1(:,:,s) = mat
                enddo ! SOC spin loop

                Mpi_( if(l_soc) then ; call Msum(mat1) ; else ; call Msum(mat) ; endif )

c
c               Multiply  cprod = <M phi(ib1) | phi(ib2)> = <M phi(ib1) | basis> <basis | phi(ib2)>
c               (Different cases: STOREIBZ, SOC)

# if SUBTYPE == 1
                if(trafo) then

                  if(.not.l_soc) mat0 = mat

                  if(l_soc) then ; do nn = 1,nindxm(ll,itype) ; call sort(mat1(1,nn,1)) ; call sort(mat1(1,nn,2)) ; enddo
                  else           ; do nn = 1,nindxm(ll,itype) ; call sort(mat0(1,nn))                             ; enddo
                  endif

                  do ik = 1,nk
                    ik0  = k0(ik)
                    ik2  = k2(ik)
                    ikx  = kindx(ik2)
                    cexp = exp ( -img * 2*pi * dot_product(kpt(:,ik0),cent(:,ic)) ) / phase(ik2)

                    ! rotate mat (no SOC: mat0->mat, SOC: mat1->mat2) if needed (otherwise copy)
                    ! - note that we use macmat later; so, all complex matrices have to be conjugated: transp(dwgn)->conjg(transp(dwgn)), esoc->conjg(esoc)
                    if(kptp(ik2)/=ik2) then
c                      trf2  = .true.
                      isym  = symkpt(ik2)
                      trs   = isym>nsymt !; if(trs) call teststop('time reversal symmetry')
                      ic1   = pcent(ic,sym(isym)%inv)
                      cdum  = exp( img * (2*pi) * dot_product(kpt(:,ik2),tcent(:,ic1,isym)) )
                      lm    = 1
                      do l = 0,lcutp(itype)
                        nx               = nindx(l,itype)
                        dwgn1(-l:l,-l:l) = conjg(transpose(dwgn(-l:l,-l:l,l,isym)))
                        do n = 1,nx
                          lm1 = lm + 2*l
                          if(l_soc) then
                            mat2(lm:lm1,:,1) = matmul ( dwgn1(-l:l,-l:l) , mat1(lm:lm1,:,1) ) * cdum
                            mat2(lm:lm1,:,2) = matmul ( dwgn1(-l:l,-l:l) , mat1(lm:lm1,:,2) ) * cdum
                          else
                            mat(lm:lm1,:) = matmul ( dwgn1(-l:l,-l:l) , mat0(lm:lm1,:) ) * cdum
                          endif
                          lm = lm + 2*l+1
                        enddo
                      enddo
                      if(l_soc) then
                        mat         = mat2(:,:,1) * conjg(sym(isym)%esoc(1,1)) + mat2(:,:,2) * conjg(sym(isym)%esoc(1,2))
                        mat2(:,:,2) = mat2(:,:,1) * conjg(sym(isym)%esoc(2,1)) + mat2(:,:,2) * conjg(sym(isym)%esoc(2,2))
                        mat2(:,:,1) = mat
                      endif
                    else
                      if(l_soc) then ; mat2 = mat1
                      else           ; mat  = mat0
                      endif
                    endif

                    do s = 1,2 ! SOC spin loop
                    if     (l_soc)     then ; mat = mat2(:,:,s)
                    else if(s/=ISPIN2) then ; cycle
                    endif

                    do nn = 1,nindxm(ll,itype)
                      call resort(mat(1,nn))
                    enddo

                    lm1 = maxlm
                    if(mtthr>0) then
                      lm1 = 0
                      do lm = 1,maxlm
                        if(sum(abs(real(mat(lm,:)))+abs(imag(mat(lm,:))))>mtthr) then
                          lm1        = lm1 + 1
                          pnt(lm1)   = lm
                          mat(lm1,:) = mat(lm,:)
                        endif
                      enddo
                    endif

                    allocate(cmt2(lm1,(b2up-b2low)/b2inc+1))
                    nn = nindxm(ll,itype)
                    if(lm1<maxlm) then
                      if(kptp(ik2)/=ik2) then ; cmt2 = CMT(pnt(:lm1),ic1,b2low:b2up:b2inc,ikx,s) ; if(trs) cmt2 = conjg(cmt2)
                      else                    ; cmt2 = CMT(pnt(:lm1),ic ,b2low:b2up:b2inc,ikx,s)
                      endif
                      NoInv( cprod(llm+1:llm+nn,:,ib,ik) = cprod(llm+1:llm+nn,:,ib,ik) + cexp * macmat(mat(:lm1,:),cmt2)   )
                      Inv(   call symmetrize_to_cprod                                  ( cexp * macmat(mat(:lm1,:),cmt2) ) )
                    else
                      if(kptp(ik2)/=ik2) then ; cmt2 = CMT(    :lm1 ,ic1,b2low:b2up:b2inc,ikx,s) ; if(trs) cmt2 = conjg(cmt2)
                      else                    ; cmt2 = CMT(    :lm1 ,ic ,b2low:b2up:b2inc,ikx,s)
                      endif
                      NoInv( cprod(llm+1:llm+nn,:,ib,ik) = cprod(llm+1:llm+nn,:,ib,ik) + cexp * macmat(mat,cmt2)   )
                      Inv(   call symmetrize_to_cprod                                  ( cexp * macmat(mat,cmt2) ) )
                    endif
                    deallocate(cmt2)

                    enddo ! SOC spin loop

                  enddo ! ik

                else
# endif

                  do s = 1,2 ! SOC spin loop
                  if     (l_soc)     then ; mat = mat1(:,:,s)
# if SUBTYPE == 1
                  else if(s/=ISPIN2) then ; cycle ! CMT is accessed (therefore ISPIN2)
# else
                  else if(s/=ispin2) then ; cycle ! cmt4 is accessed (therefore ispin2)
# endif
                  endif

                  lm1 = maxlm
                  if(mtthr>0) then
                    lm1 = 0
                    do lm = 1,maxlm
                      if(sum(abs(real(mat(lm,:)))+abs(imag(mat(lm,:))))>mtthr) then
                        lm1        = lm1 + 1
                        pnt(lm1)   = lm
                        mat(lm1,:) = mat(lm,:)
                      endif
                    enddo
                  endif

                  nn = nindxm(ll,itype)
# if SUBTYPE == 1
                  allocate(cmt2(lm1,(b2up-b2low)/b2inc+1))
                  do ik = 1,nk
                    ik0  = k0(ik)
                    ik2  = k2(ik)
                    ikx  = kindx(ik2)
                    cexp = exp ( -img * 2*pi * dot_product(kpt(:,ik0),cent(:,ic)) )
                    if(lm1<maxlm) then
                      cmt2 = CMT(pnt(:lm1),ic,b2low:b2up:b2inc,ikx,s)
                      NoInv( cprod(llm+1:llm+nn,:,ib,ik) = cprod(llm+1:llm+nn,:,ib,ik) + cexp * macmat(mat(:lm1,:),cmt2)   )
                      Inv(   call symmetrize_to_cprod                                  ( cexp * macmat(mat(:lm1,:),cmt2) ) )
                    else
                      cmt2 = CMT(    :lm1 ,ic,b2low:b2up:b2inc,ikx,s)
                      NoInv( cprod(llm+1:llm+nn,:,ib,ik) = cprod(llm+1:llm+nn,:,ib,ik) + cexp * macmat(mat,cmt2)   )
                      Inv(   call symmetrize_to_cprod                                  ( cexp * macmat(mat,cmt2) ) )
                    endif
                  enddo
                  deallocate(cmt2)
# else
                  if(lm1<maxlm) then
                    allocate(cmt2(lm1,b2up-b2low+1))
                    cmt2 = cmt4(pnt(:lm1),b2low:b2up,s)
                    NoInv( cprod(:,ib1,ik,llm+1:llm+nn) = cprod(:,ib1,ik,llm+1:llm+nn) + cexp * conjg(macmat(cmt2,mat(:lm1,:)))   )
                    Inv(   call symmetrize_to_cprod                                    ( cexp * conjg(macmat(cmt2,mat(:lm1,:))) ) )
                    deallocate(cmt2)
                  else
                    NoInv( cprod(:,ib1,ik,llm+1:llm+nn) = cprod(:,ib1,ik,llm+1:llm+nn) + cexp * conjg(macmat(cmt4(:,:,s),mat))   )
                    Inv(   call symmetrize_to_cprod                                    ( cexp * conjg(macmat(cmt4(:,:,s),mat)) ) )
                  endif
# endif

                  enddo ! SOC spin loop

# if SUBTYPE == 1
                endif
# endif

                deallocate(mat)
                if(l_soc) deallocate(mat1)
                if(storeibz) then
                  if(l_soc) then ; deallocate(mat2)
                  else           ; deallocate(mat0)
                  endif
                endif

                llm = llm + nindxm(ll,itype)
              enddo ! mm
              lln = lln + nindxm(ll,itype)
            enddo   ! ll
          enddo     ! ib(1)
# if SUBTYPE == 2
          enddo     ! ik(1)
# endif
          llm0 = llm
        enddo ! ieq
        deallocate ( integral )

        deallocate(cmt1)
        if(l_soc) deallocate(cmt3)
# if SUBTYPE == 2
        deallocate(cmt4)
# endif

      enddo ! itype

# if SUBTYPE == 1
      if(storeibz) deallocate ( dwgn1 )
# endif

# ifdef SUBWANNIER
      if(storeibz) call transform_wan(cprod,k1,k2,nk,ispin1,ispin2,b2up-b2low+1,nbasp)
# endif

c      write(*,*) trf1,trf2

      contains
      subroutine sort(mat)
      implicit none
      complex_dp, intent(inout) :: mat(maxlm)
      complex_dp                :: hlp(maxindx*(2*maxlcut+1))
      integer                   :: l,n,nx,mm,lm,lm1
      lm = 0
      do l = 0,lcutp(itype)
        nx  = nindx(l,itype)
        mm  = 2*l+1
        lm1 = 0
        do n = 1,nx
          hlp(lm1+1:lm1+mm) = mat(lm+n:lm+mm*nx:nx)
          lm1               = lm1 + mm
        enddo
        mat(lm+1:lm+mm*nx) = hlp(:lm1)
        lm                 = lm + mm*nx
      enddo
      end subroutine sort

      subroutine resort(mat)
      implicit none
      complex_dp, intent(inout) :: mat(maxlm)
      complex_dp                :: hlp(maxindx*(2*maxlcut+1))
      integer                   :: l,m,mm,nx,lm,lm1
      lm = 0
      do l = 0,lcutp(itype)
        nx  = nindx(l,itype)
        mm  = 2*l+1
        lm1 = 0
        do m = 1,mm
          hlp(lm1+1:lm1+nx) = mat(lm+m:lm+mm*nx:mm)
          lm1               = lm1 + nx
        enddo
        mat(lm+1:lm+mm*nx) = hlp(:lm1)
        lm                 = lm + mm*nx
      enddo
      end subroutine resort

# ifdef INV
      subroutine symmetrize_to_cprod(mat)
      implicit none
      complex_dp, intent(in) :: mat(:,:)
      integer                :: llm1,ic1,sgn,nn
      real_dp,    parameter  :: sq = sqrt(0.5d0)
      nn   = nindxm(ll,itype)
      ic1  = pcent(ic,invsym)
      llm1 = llm - 2*mm*nindxm(ll,itype) + (ic1-ic)*maxlmm
      if(llm1/=llm) then
        if(llm1>llm) then
#   if SUBTYPE == 1
          cprod(llm +1:llm +nn,:,ib,ik)  = cprod(llm +1:llm +nn,:,ib,ik)  +       mat * sq
          cprod(llm1+1:llm1+nn,:,ib,ik)  = cprod(llm1+1:llm1+nn,:,ib,ik)  - img * mat * sq
#   else
          cprod(:,ib1,ik,llm +1:llm +nn) = cprod(:,ib1,ik,llm +1:llm +nn) +       mat * sq
          cprod(:,ib1,ik,llm1+1:llm1+nn) = cprod(:,ib1,ik,llm1+1:llm1+nn) - img * mat * sq
#   endif
        else
          sgn = (-1)**(mm+ll)
#   if SUBTYPE == 1
          cprod(llm1+1:llm1+nn,:,ib,ik)  = cprod(llm1+1:llm1+nn,:,ib,ik)  +       mat * (sq*sgn)
          cprod(llm +1:llm +nn,:,ib,ik)  = cprod(llm +1:llm +nn,:,ib,ik)  + img * mat * (sq*sgn)
#   else
          cprod(:,ib1,ik,llm1+1:llm1+nn) = cprod(:,ib1,ik,llm1+1:llm1+nn) +       mat * (sq*sgn)
          cprod(:,ib1,ik,llm +1:llm +nn) = cprod(:,ib1,ik,llm +1:llm +nn) + img * mat * (sq*sgn)
#   endif
        endif
      else
        if(mod(ll,2)==0) then
#   if SUBTYPE == 1
          cprod(llm+1:llm+nn,:,ib,ik)  = cprod(llm+1:llm+nn,:,ib,ik)  + mat
#   else
          cprod(:,ib1,ik,llm+1:llm+nn) = cprod(:,ib1,ik,llm+1:llm+nn) + mat
#   endif
        else
#   if SUBTYPE == 1
          cprod(llm+1:llm+nn,:,ib,ik)  = cprod(llm+1:llm+nn,:,ib,ik)  - img * mat
#   else
          cprod(:,ib1,ik,llm+1:llm+nn) = cprod(:,ib1,ik,llm+1:llm+nn) - img * mat
#   endif
        endif
      endif
      end subroutine symmetrize_to_cprod
# endif

      end

# undef ifMOD_
# undef endMOD_
# undef Mpi_

# ifdef PARA_PW
#   define Mcol1_ Mcol1
#   define McolD1_ McolD1
#   define MrangeDistr_ MrangeDistr
#   define Mrange1_ Mrange1
# else
#   define Mcol1_(arg) 1:arg
#   define McolD1_(arg1,arg2)
#   define MrangeDistr_(arg1,arg2)
#   define Mrange1_(arg) 1,arg
# endif

c --------------------------------------------------------------------------------
c
c Plane-wave part of wave-function products
c
c cprod = < k0+G0 phi(b1,k1) | phi(b2,k2) >
c
c
c (1) Fast Fourier transformation (FFT)
c
c FULLPW: Multiply with step function phi(b1,k1)*theta -> phi(b1,k1)
c
c Transform phi(b1,k1,G) -> phi(b1,k1,R)
c           phi(b2,k2,G) -> phi(b2,k2,R)
c Multiply  phi(b1,k1,R)*phi(b2,k2,R) -> prod(b1,b2,k1,k2,R)
c Backtrafo prod(b1,b2,k1,k2,R) -> cprod(b1,b2,k1,k2,G)
c
c
c (2) Explicit Convolutions (default)
c
c FULLPW: Multiply with step function
c         phi~(b1,k1,G) = SUM(G') phi(b1,k1,G') * theta(G-G')
c
c Convolute   < k0+G0     phi~(b1,k1,G1) | phi(b2,k2,G2)      >  ! k2 = k1+k0+g (g=backfolding vector)
c           = < k0+G0     phi~(b1,k1,G1) | phi(b2,k0+k1+g,G2) >  ! k0+g+k1 on the left balances with k0+k1+g on the right
c           = < k0+g+G0-g phi~(b1,k1,G1) | phi(b2,k0+k1+g,G2) >  ! -> G1 = G2-G0+g
c           = SUM(G2) conjg[phi~(b1,k1,G2-G0+g)] * phi(b2,k2,G2) -> cprod(b1,b2,k1,k2,G0)
c

# if SUBTYPE == 1
      subroutine wavefproducts1_pw(cprod,dim,b1,nb1,ik1,ispin,k0,nk,b2low,b2up,b2inc)
# else
#   ifndef SUBWANNIER
      subroutine wavefproducts2_pw(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   else
      subroutine wavefproducts3_pw(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   endif
# endif
      use global
      use wrapper, only: dotprod
      use m_fft
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )
      implicit none
# if SUBTYPE == 1
      integer,     intent(in)    :: nb1,b1(nb1),ik1,ispin,nk,k0(nk),b2low,b2up,b2inc,dim
      MCOMPLEX_dp, intent(inout) :: cprod(dim, (b2up-b2low) / b2inc + 1 ,nb1,nk)
      integer                    :: ispin1,ispin2
      integer                    :: ik0,ib2_
# else
      integer,     intent(in)    :: nk,k1(nk),ik0,ispin1,ispin2,b1low,b1up,b2low,b2up
      MCOMPLEX_dp, intent(inout) :: cprod(b2low:b2up,b1low:b1up,nk,ngptm(ik0))
      integer                    :: ik1,nb1,nb2
# endif
      complex_dp,  allocatable   :: rfft(:,:,:,:)
      MCOMPLEX_dp, allocatable   :: cpw1(:),cpw1s(:,:),cpw_1(:,:,:),cpw0(:,:,:),cpw3(:,:),cpw2(:,:,:)
      real_dp                    :: sqrtvoli
      integer,     allocatable   :: point(:,:,:),gpt0(:,:),pgptinv(:)
      integer                    :: k2(nk),ik,ik2,s,ib1,ib2,ib,ikx
      integer                    :: gpt2(3,nk),g(3),g1(3),ig,ig1,ig2
      integer                    :: ngpt0,ig0
      integer                    :: ngpt1,ngpt2,pgpt1(maxgpt),pgpt2(maxgpt),pg(maxgpt)
      integer                    :: i,j,k
      integer                    :: def
      integer                    :: kptsum
# if defined INV || defined INVW
      real_dp,     allocatable   :: step(:)
# else
      complex_dp,  allocatable   :: step(:)
# endif

      real cputime
c      logical trf1,trf2
c      trf1 = .false.
c      trf2 = .false.

      if(nk==0.or.b2low>b2up) return

      call cpu_time(cputime)

      sqrtvoli = 1/sqrt(vol)

# if SUBTYPE == 1
      ispin1 = ispin
      ispin2 = ispin
# else
      if(l_soc.and.(ispin1/=1.or.ispin2/=1)) Bug('Wrong spin indices (SOC).')
c      if(ispin1/=ispin2) call teststop('ispin1/=ispin2')
      nb1 = b1up - b1low + 1
      nb2 = b2up - b2low + 1
# endif

      if(gcutf/=0) then
        allocate ( cpw1(maxgpt) )
        if(l_soc) allocate ( cpw1s(maxgpt,2) )
      endif

      ! k2 = q+k
      do ik = 1,nk
# if SUBTYPE == 1
        ik0 = k0(ik)
# else
        ik1 = k1(ik)
# endif
        k2(ik)     = kptsum(ik0,ik1)
        gpt2(:,ik) = nint ( kpt(:,k2(ik)) - kpt(:,ik0) - kpt(:,ik1) )
        if(any(abs(kpt(:,ik1)+kpt(:,ik0)+gpt2(:,ik)-kpt(:,k2(ik)))>1d-12)) Bug('error')
      enddo

# if SUBTYPE == 1
      if(dim<minval(nbasm(k0))) Bug('Dimension dim too small.')
      if(dim==nbasp) Error('No IPWs?')
      cprod(nbasp+1:,:,:,:) = 0
# else
      cprod = 0
# endif

      if(.not.fullpw) allocate ( pgptinv(ngptall) )

c
c (1) Fast Fourier Transformation
c

      if(gcutf/=0) then
        ! Calculate step function with Fourier components |G|<=2*gcut+gcutm
        if(fullpw) then
          cffts WSUFFIX = 0
          do k = 0,nfft(3)-1     ; g(3) = k ; if(g(3)>nfft(3)/2) g(3) = g(3) - nfft(3)
            do j = 0,nfft(2)-1   ; g(2) = j ; if(g(2)>nfft(2)/2) g(2) = g(2) - nfft(2)
              do i = 0,nfft(1)-1 ; g(1) = i ; if(g(1)>nfft(1)/2) g(1) = g(1) - nfft(1)
                g1                    = matmul(imat,g) ; if(sum(matmul(rlat,g1)**2)>(2*gcut+gcutm)**2) cycle
                cffts WSUFFIX (i,j,k) = cstep(g1(1),g1(2),g1(3))
              enddo
            enddo
          enddo
          call dfftw_execute_dft Inv(_r2c) (plan_ffts WSUFFIX,cffts WSUFFIX,rffts)
        endif

# if SUBTYPE == 2
c       Either ib1 or ib2 as the inner loop (two blocks of code)
c       --- First block ---
        if(b2up-b2low>b1up-b1low) then
          i = 1 ; if(l_soc) i = 2
          allocate ( rfft(size(rfft1,1),size(rfft1,2),size(rfft1,3),nb1*i) )
          rfft = 0
          do ik = 1,nk
            ik1 = k1(ik)
            ik2 = k2(ik)
            do ib = 1,nb1
              ib1 = b1low + ib - 1
# else
        i = 1 ; if(l_soc) i = 2
        allocate ( rfft(size(rfft1,1),size(rfft1,2),size(rfft1,3),nb1*i) )
        rfft = 0
        do ib = 1,nb1
          ib1 = b1(ib)
# endif
# ifndef old_trafo
#   if defined(LOAD) && SUBTYPE == 1
          if(storeibz.and.kptp(ik1)/=ik1) Error('Not implemented: LOAD & ik1 not in IBZ.')
          if(l_soc) then ; cpw1s = cpwq(:,ib1,:)
          else           ; cpw1  = cpwq(:,ib1,ISPIN1)
          endif
#   else
          if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib1,ik1,ISPIN1)
          else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib1,ik1,ISPIN1)
          endif
#   endif
# else
          if(storeibz.and.kptp(ik1)/=ik1) then
            if(l_soc) then
              call waveftrafo_pw(cpw1s(:,1),ik1,ib1,1,LWAN)
              call waveftrafo_pw(cpw1s(:,2),ik1,ib1,2,LWAN)
              call waveftrafo_soc(cpw1s,maxgpt,symkpt(ik1))
            else
              call waveftrafo_pw(cpw1,ik1,ib1,ISPIN1,LWAN)
            endif
          else
            ikx = kindx(ik1)
            if(l_soc) then ; cpw1s = CPW1(:,ib1 IKX ,:)
            else           ; cpw1  = CPW1(:,ib1 IKX ,ISPIN1)
            endif
          endif
# endif
          do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
          cfft1 WSUFFIX = 0
          do ig1 = 1,ngpt(ik1)
            g                              = matmul(imati,gpt(:,pgpt(ig1,ik1)))
            g                              = modulo(g,nfft)
            cfft1 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig1)
          enddo
          call dfftw_execute_dft Inv(_r2c) (plan_fft1 WSUFFIX,cfft1 WSUFFIX,rfft1)
          rfft1 = conjg ( rfft1 * sqrtvoli )
          if(fullpw) then
            rfft1 = rfft1 * rffts
          endif
          if(l_soc) then ; rfft(:,:,:,ib+(s-1)*nb1) = rfft1
          else           ; rfft(:,:,:,ib)           = rfft1
          endif
          enddo ! SOC spin loop
        enddo ! ib1
# if SUBTYPE == 1
        do ik = 1,nk
          ik0  = k0(ik)
          ik2  = k2(ik)
          ib2_ = 0
          do ib2 = b2low,b2up,b2inc
            ib2_ = ib2_ + 1
# else
          do ib2 = b2low,b2up
# endif
# ifndef old_trafo
          ADJUST_KINDX_IF_NEEDED
          if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib2,ik2,ISPIN2)
          else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib2,ik2,ISPIN2)
          endif
          UNDO_ADJUST_KINDX
# else
          if(storeibz.and.kptp(ik2)/=ik2) then
            if(l_soc) then
              call waveftrafo_pw(cpw1s(:,1),ik2,ib2,1,LWAN)
              call waveftrafo_pw(cpw1s(:,2),ik2,ib2,2,LWAN)
              call waveftrafo_soc(cpw1s,maxgpt,symkpt(ik2))
            else
              call waveftrafo_pw(cpw1,ik2,ib2,ISPIN2,LWAN)
            endif
          else
            ikx = kindx(ik2)
#   if SUBTYPE == 2 && defined(LOAD)
            if(.not.storeibz.and..not.LWAN) ikx = kindx(ik1)
#   endif
            if(l_soc) then ; cpw1s = CPW(:,ib2,ikx,:)
            else           ; cpw1  = CPW(:,ib2,ikx,ISPIN2)
            endif
          endif
# endif
          do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
          cfft2 WSUFFIX = 0
          do ig2 = 1,ngpt(ik2)
            g                              = matmul(imati,gpt(:,pgpt(ig2,ik2)))
            g                              = modulo(g,nfft)
            cfft2 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig2)
          enddo
          call dfftw_execute_dft Inv(_r2c) (plan_fft2 WSUFFIX,cfft2 WSUFFIX,rfft2)
          rfft2 = rfft2 / nnfft
          do ib = 1,nb1
# if SUBTYPE == 2
            ib1 = b1low + ib - 1
# endif
            if(l_soc) then ; rfftp = rfft(:,:,:,ib+(s-1)*nb1) * rfft2
            else           ; rfftp = rfft(:,:,:,ib)           * rfft2
            endif
            call dfftw_execute_dft Inv(_c2r) (plan_fftp WSUFFIX,rfftp,cfftp WSUFFIX)
            do ig = 1,ngptm(ik0)
              g                          = matmul(imati,gptm(:,pgptm(ig,ik0))-gpt2(:,ik))
              g                          = modulo(g,nfft)
# if SUBTYPE == 1
              cprod(nbasp+ig,ib2_,ib,ik) = cprod(nbasp+ig,ib2_,ib,ik) + cfftp WSUFFIX (g(1),g(2),g(3))
# else
              cprod(ib2,ib1,ik,ig)       = cprod(ib2,ib1,ik,ig)       + cfftp WSUFFIX (g(1),g(2),g(3))
# endif
            enddo
          enddo ! ib(1)
          enddo ! SOC spin loop
        enddo   ! ib2
        enddo   ! ik0 / ik

# if SUBTYPE == 2
c       --- Second block ---
        else
          i = 1 ; if(l_soc) i = 2
          allocate ( rfft(size(rfft2,1),size(rfft2,2),size(rfft2,3),nb2*i) )
          rfft = 0
          do ik = 1,nk
            ik1 = k1(ik)
            ik2 = k2(ik)
            do ib2 = b2low,b2up
              ib = ib2 - b2low + 1
# ifndef old_trafo
              ADJUST_KINDX_IF_NEEDED
              if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib2,ik2,ISPIN2)
              else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib2,ik2,ISPIN2)
              endif
              UNDO_ADJUST_KINDX
# else
              if(storeibz.and.kptp(ik2)/=ik2) then
                if(l_soc) then
                  call waveftrafo_pw(cpw1s(:,1),ik2,ib2,1,LWAN)
                  call waveftrafo_pw(cpw1s(:,2),ik2,ib2,2,LWAN)
                  call waveftrafo_soc(cpw1s,maxgpt,symkpt(ik2))
                else
                  call waveftrafo_pw(cpw1,ik2,ib2,ISPIN2,LWAN)
                endif
              else
                ikx = kindx(ik2) ; Load( if(.not.storeibz.and..not.LWAN) ikx = kindx(ik1) )
                if(l_soc) then ; cpw1s = CPW(:,ib2,ikx,:)
                else           ; cpw1  = CPW(:,ib2,ikx,ISPIN2)
                endif
              endif
# endif
              do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
              cfft2 WSUFFIX = 0
              do ig2 = 1,ngpt(ik2)
                g                              = matmul(imati,gpt(:,pgpt(ig2,ik2)))
                g                              = modulo(g,nfft)
                cfft2 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig2)
              enddo
              call dfftw_execute_dft Inv(_r2c) (plan_fft2 WSUFFIX,cfft2 WSUFFIX,rfft2)
              if(fullpw) then
                rfft2 = rfft2 * rffts
              endif
              if(l_soc) then ; rfft(:,:,:,ib+(s-1)*nb2) = rfft2
              else           ; rfft(:,:,:,ib)           = rfft2
              endif
              enddo ! SOC spin loop
            enddo ! ib2
            do ib1 = b1low,b1up
# ifndef old_trafo
              if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib1,ik1,ISPIN1)
              else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib1,ik1,ISPIN1)
              endif
# else
              if(storeibz.and.kptp(ik1)/=ik1) then
                if(l_soc) then
                  call waveftrafo_pw(cpw1s(:,1),ik1,ib1,1,LWAN)
                  call waveftrafo_pw(cpw1s(:,2),ik1,ib1,2,LWAN)
                  call waveftrafo_soc(cpw1s,maxgpt,symkpt(ik1))
                else
                  call waveftrafo_pw(cpw1,ik1,ib1,ISPIN1,LWAN)
                endif
              else
                ikx = kindx(ik1)
                if(l_soc) then ; cpw1s = CPW1(:,ib1 IKX ,:)
                else           ; cpw1  = CPW1(:,ib1 IKX ,ISPIN1)
                endif
              endif
# endif
              do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
              cfft1 WSUFFIX = 0
              do ig1 = 1,ngpt(ik1)
                g                              = matmul(imati,gpt(:,pgpt(ig1,ik1)))
                g                              = modulo(g,nfft)
                cfft1 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig1)
              enddo
              call dfftw_execute_dft Inv(_r2c) (plan_fft1 WSUFFIX,cfft1 WSUFFIX,rfft1)
              rfft1 = conjg ( rfft1 * sqrtvoli ) / nnfft
              do ib = 1,nb2
                ib2 = b2low + ib - 1
                if(l_soc) then ; rfftp = rfft(:,:,:,ib+(s-1)*nb2) * rfft1
                else           ; rfftp = rfft(:,:,:,ib)           * rfft1
                endif
                call dfftw_execute_dft Inv(_c2r) (plan_fftp WSUFFIX,rfftp,cfftp WSUFFIX)
                do ig = 1,ngptm(ik0)
                  g                    = matmul(imati,gptm(:,pgptm(ig,ik0))-gpt2(:,ik))
                  g                    = modulo(g,nfft)
                  cprod(ib2,ib1,ik,ig) = cprod(ib2,ib1,ik,ig) + cfftp WSUFFIX (g(1),g(2),g(3))
                enddo
              enddo ! ib(2)
              enddo ! SOC spin loop
            enddo   ! ib1
          enddo     ! ik
        endif
# endif
        deallocate ( rfft )
        deallocate ( cpw1 )
        if(allocated(cpw1s)) deallocate ( cpw1s )

      else

c
c (2) Explicit Convolutions
c

# if SUBTYPE == 2
        do ik = 1,nk
          ik1 = k1(ik) ; ngpt1 = ngpt(ik1) ; pgpt1 = pgpt(:,ik1)
          ik2 = k2(ik) ; ngpt2 = ngpt(ik2) ; pgpt2 = pgpt(:,ik2)
          allocate ( cpw2(ngpt2,b2low:b2up,nspin3),cpw3(ngpt2,maxgptm) )
# ifndef old_trafo
          ADJUST_KINDX_IF_NEEDED
          do ib2 = b2low,b2up
            call wavefunction_pw WSUB(cpw2(:,ib2,:),ngpt2,ib2,ik2,ISPIN2)
          enddo
          UNDO_ADJUST_KINDX
# else
          do ib2 = b2low,b2up
            if(storeibz.and.kptp(ik2)/=ik2) then
              if(l_soc) then
                call waveftrafo_pw(cpw2(:,ib2,1),ik2,ib2,1,LWAN)
                call waveftrafo_pw(cpw2(:,ib2,2),ik2,ib2,2,LWAN)
                call waveftrafo_soc(cpw2(:,ib2,:),ngpt2,symkpt(ik2))
              else
                call waveftrafo_pw(cpw2(:,ib2,1),ik2,ib2,ISPIN2,LWAN)
              endif
            else
              ikx = kindx(ik2) ; Load( if(.not.storeibz.and..not.LWAN) ikx = kindx(ik1) )
              cpw2(:,ib2,1) = CPW(:ngpt(ik2),ib2,ikx,ISPIN2)
              if(l_soc) cpw2(:,ib2,2) = CPW(:ngpt(ik2),ib2,ikx,2)
            endif
          enddo
# endif
# else
        ngpt1 = ngpt(ik1)
        pgpt1 = pgpt(:,ik1)
# endif

        if(.not.fullpw) then
          pgptinv                       = 0
          pgptinv(pgpt(:ngpt(ik1),ik1)) = [ (i,i=1,ngpt(ik1)) ]
        endif

c        write(*,'(A'NoA) '1';call cpu_done(cputime)

c
c       Prepare multiplication with step function (->point,gpw0)
        def   =  0
        g1    = -1
 1      ngpt0 =  0
# if SUBTYPE == 1
        do ik = 1,nk
          ik0 = k0(ik)
          ik2 = k2(ik) ; ngpt2 = ngpt(ik2) ; pgpt2 = pgpt(:,ik2)
# endif
        do ig2 = 1,ngpt2
          do ig = 1,ngptm(ik0)
            g = gpt(:,pgpt2(ig2)) - gptm(:,pgptm(ig,ik0)) + gpt2(:,ik)
            if(def==0) then
              g1(1) = max(g1(1),abs(g(1)))
              g1(2) = max(g1(2),abs(g(2)))
              g1(3) = max(g1(3),abs(g(3)))
            else
              if(point(g(1),g(2),g(3))==0) then
                if(fullpw) then
                  ngpt0                 = ngpt0 + 1
                  point(g(1),g(2),g(3)) = ngpt0
                else if(sum(matmul(rlat,kpt(:,ik1)+g)**2)<=gcut**2) then
                  ngpt0                 = ngpt0 + 1
                  point(g(1),g(2),g(3)) = pgptinv ( pntgpt(g(1),g(2),g(3)) )
                endif
                if(def==2) then
                  gpt0(:,ngpt0) = g
                endif
              endif
            endif
          enddo
        enddo
# if SUBTYPE == 1
        enddo
# endif
        if(def==0) then
          allocate ( point(-g1(1):g1(1),-g1(2):g1(2),-g1(3):g1(3)) )
          point = 0
          def   = 1
          goto 1
        else if(def==1.and.fullpw) then
          allocate ( gpt0(3,ngpt0) )
          point = 0
          def   = 2
          goto 1
        endif

        if(fullpw) allocate ( cpw0(ngpt0,nspin3,nb1) )
        allocate ( cpw_1(ngpt1,nspin3,nb1) )

c
c       Get coefficients for states b1 (->cpw_1)

        do ib = 1,nb1
# if SUBTYPE == 1
          ib1 = b1(ib)
# else
          ib1 = b1low + ib - 1
# endif
# ifndef old_trafo
#   if defined(LOAD) && SUBTYPE == 1
          if(storeibz.and.kptp(ik1)/=ik1) Error('Not implemented: LOAD & ik1 not in IBZ.')
          if(l_soc) then ; cpw_1(:,:,ib) = cpwq(:ngpt1,ib1,:)
          else           ; cpw_1(:,1,ib) = cpwq(:ngpt1,ib1,ISPIN1)
          endif
#   else
          call wavefunction_pw WSUB (cpw_1(:,:,ib),ngpt1,ib1,ik1,ISPIN1)
#   endif
# else
          if(storeibz.and.kptp(ik1)/=ik1) then
            if(l_soc) then
              call waveftrafo_pw(cpw_1(:,1,ib),ik1,ib1,1,LWAN)
              call waveftrafo_pw(cpw_1(:,2,ib),ik1,ib1,2,LWAN)
              call waveftrafo_soc(cpw_1(:,:,ib),ngpt1,symkpt(ik1))
            else
              call waveftrafo_pw(cpw_1(:,1,ib),ik1,ib1,ISPIN1,LWAN)
            endif
          else
            ikx = kindx(ik1)
            if(l_soc) then ; cpw_1(:,:,ib) = CPW1(:ngpt1,ib1 IKX ,:)
            else           ; cpw_1(:,1,ib) = CPW1(:ngpt1,ib1 IKX ,ISPIN1)
            endif
          endif
# endif
        enddo

c
c       fullpw: Multiply cpw_1 with step function (-> cpw0)

        if(fullpw) then
          allocate(step(ngpt1))
          do ig0 = Mrange1_(ngpt0)
            do ig1 = 1,ngpt1
              g         = gpt0(:,ig0) - gpt(:,pgpt1(ig1))
              step(ig1) = cstep(g(1),g(2),g(3))
            enddo
            cpw0(ig0,:,:) = reshape ( matmul(step,reshape(cpw_1,(/ngpt1,nspin3*nb1/))) , (/nspin3,nb1/) )
          enddo
          deallocate(step)
          MrangeDistr_(cpw0(McolD1_(ngpt0,i),:,:),i)
        endif
        
c
c       Loop over ib1

        do ib = 1,nb1
# if SUBTYPE == 1
          ib1 = b1(ib)
# else
          ib1 = b1low + ib - 1
# endif        

# if SUBTYPE == 1
          do ik = 1,nk
            ik0 = k0(ik)
            ik2 = k2(ik) ; ngpt2 = ngpt(ik2) ; pgpt2 = pgpt(:,ik2)
            allocate ( cpw2(ngpt2,size(cprod,2),nspin3),cpw3(ngpt2,maxgptm) )
#   ifndef old_trafo
            ib2_ = 0
            do ib2 = b2low,b2up,b2inc
              ib2_ = ib2_ + 1
              call wavefunction_pw WSUB(cpw2(:,ib2_,:),ngpt2,ib2,ik2,ISPIN2)
            enddo
#   else
            if(storeibz.and.kptp(ik2)/=ik2) then
              ib2_ = 0
              do ib2 = b2low,b2up,b2inc
                ib2_ = ib2_ + 1
                if(l_soc) then
                  call waveftrafo_pw(cpw2(:,ib2_,1),ik2,ib2,1,LWAN)
                  call waveftrafo_pw(cpw2(:,ib2_,2),ik2,ib2,2,LWAN)
                  call waveftrafo_soc(cpw2(:,ib2_,:),ngpt2,symkpt(ik2))
                else
                  call waveftrafo_pw(cpw2(:,ib2_,1),ik2,ib2,ISPIN2,LWAN)
                endif
              enddo
            else
              ikx = kindx(ik2)
              cpw2(:,:,1) = CPW(:ngpt2,b2low:b2up:b2inc,ikx,ISPIN2)
              if(l_soc) cpw2(:,:,2) = CPW(:ngpt2,b2low:b2up:b2inc,ikx,2)
            endif
#   endif
# endif         
          do s = 1,nspin3 ! SOC spin loop

          ! define cpw3 with which we multiply phi(ib2,ik2)
          if(fullpw) then
            do ig = Mrange1_(ngptm(ik0))
              do ig2 = 1,ngpt2
                g            = gpt(:,pgpt2(ig2)) - gptm(:,pgptm(ig,ik0)) + gpt2(:,ik)
                cpw3(ig2,ig) = cpw0(point(g(1),g(2),g(3)),s,ib)
              enddo
            enddo
            MrangeDistr_(cpw3(:,McolD1_(ngptm(ik0),i)),i)
          endif

c
c         Loop over mixed-basis IPWs and define cprod
          do ig = 1,ngptm(ik0)
            if(fullpw) then
c             (1) FULLPW
              ! multiply with phi(ib2,ik2)
# if SUBTYPE == 1
              do ib2 = 1,size(cprod,2)
                cprod(nbasp+ig,ib2,ib,ik) = cprod(nbasp+ig,ib2,ib,ik) + dotprod(cpw3(:,ig),cpw2(:,ib2,s)) * sqrtvoli
# else
              do ib2 = b2low,b2up
                cprod(ib2,ib1,ik,ig)      = cprod(ib2,ib1,ik,ig)      + dotprod(cpw3(:,ig),cpw2(:,ib2,s)) * sqrtvoli
# endif
              enddo
            else
c             (2) APPROXPW
              ! define cpw3 with which we multiply phi(ib2,ik2)
              k = 0
              do ig2 = 1,ngpt2
                g         = gpt(:,pgpt2(ig2)) - gptm(:,pgptm(ig,ik0)) + gpt2(:,ik)
                ig1       = point(g(1),g(2),g(3))
                if(ig1/=0) then
                  k         = k + 1
                  cpw3(k,1) = cpw_1(ig1,s,ib)
                  pg(k)     = ig2
                endif
              enddo
              ! multiply with phi(ib2,ik2)
# if SUBTYPE == 1
              do ib2 = 1,size(cprod,2)
                cprod(nbasp+ig,ib2,ib,ik) = cprod(nbasp+ig,ib2,ib,ik) + dotprod(cpw3(:k,1),cpw2(pg(:k),ib2,s)) * sqrtvoli
# else
              do ib2 = b2low,b2up
                cprod(ib2,ib1,ik,ig)      = cprod(ib2,ib1,ik,ig)      + dotprod(cpw3(:k,1),cpw2(pg(:k),ib2,s)) * sqrtvoli
# endif
              enddo
            endif
          enddo

          enddo ! SOC spin loop
# if SUBTYPE == 1
            deallocate ( cpw2,cpw3 )
          enddo ! ik(0)
# endif

        enddo ! ib(1)

c        write(*,'(A'NoA) '4';call cpu_done(cputime)

        deallocate ( cpw_1 )

        if(fullpw) deallocate ( cpw0,gpt0 )
        deallocate ( point )
# if SUBTYPE == 2
          deallocate ( cpw2,cpw3 )
        enddo ! ik(1)
# endif

      endif

      if(.not.fullpw) deallocate ( pgptinv )

# ifdef SUBWANNIER
      if(storeibz) call transform_wan(cprod,k1,k2,nk,ispin1,ispin2,b2up-b2low+1,ngptm(ik0))
# endif

      end

c --------------------------------------------------------------------------------

# ifdef SUBWANNIER

c --------------------------------------------------------------------------------

c If storeibz=.true., wavefproducts returns < M(k0) P(k1) phi(b1,k1p) | P(k2) phi(b2,k2p) >,
c where P(k1) [P(k2)] rotates from k1p (k2p), the parent of k1 (k2), to k1 (k2).
c
c In the case of Wannier functions, this is < M(k0) P(k1) wan(n1,k1p) | P(k2) wan(n2,k2p) >,
c but we need < M(k0) wan(n1,k1) | wan(n2,k2) >.
c
c Here, we transform both Wannier functions using:      
c | wan(n,k) > = sum(b) U(b,n,k) | phi(b,k) > = sum(b) U(b,n,k) | P(k) phi(b,kp) >
c              = sum(b) U(b,n,k) [ sum(n') U*(b,n',kp)   | P(k) wan(n',kp) > ]
c              = sum(n') [ sum(b) U(b,n,k) U*(b,n',kp) ] | P(k) wan(n',kp) >
c                \-----       = trafo(n,n')       -----/
c
      subroutine transform_wan(cprod,k1,k2,nk,s1,s2,dim1,dim4)
      use global
      use wrapper, only: matmat
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: nk,k1(nk),k2(nk),dim1,dim4,s1,s2
      complex_dp, intent(inout) :: cprod(dim1,nwan,nk,dim4)
      complex_dp                :: trafo(nwan,nwan)
      integer                   :: ik,ik1,ik2,ik1p,ik2p,ibas
      do ik = 1,nk
        ik1 = k1(ik) ; ik1p = kptp(ik1)
        ik2 = k2(ik) ; ik2p = kptp(ik2)
        if(ik1>nkpti) then
          if(symkpt(ik1)>nsymt) then
            trafo = matmat ( conjg(transpose(uwan(:,:,ik1,s1))) ,    conjg(uwan(:,:,ik1p,s1)) )
          else
            trafo = matmat ( conjg(transpose(uwan(:,:,ik1,s1))) ,          uwan(:,:,ik1p,s1)  ) * phase(ik1)
          endif
          do ibas = 1,dim4
            cprod(:,:,ik,ibas) = matmat ( cprod(:,:,ik,ibas) , transpose(trafo) )
          enddo
        endif
        if(ik2>nkpti) then
          if(symkpt(ik2)>nsymt) then
            trafo = matmat ( transpose(uwan(:,:,ik2,s2)) ,       uwan(:,:,ik2p,s2) )
          else
            trafo = matmat ( transpose(uwan(:,:,ik2,s2)) , conjg(uwan(:,:,ik2p,s2) ) ) * conjg(phase(ik2))
          endif
          do ibas = 1,dim4
            cprod(:,:,ik,ibas) = matmat ( trafo , cprod(:,:,ik,ibas) )
          enddo
        endif
      enddo
      end

c --------------------------------------------------------------------------------

c   Special wavefunction routines for Wannier Bloch functions (SUBWANNIER)

#   ifndef old_trafo

      subroutine wavefunction_mt_w(cmtout,dim,ic,iband,ikpt,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ic,ikpt,iband,ispin,dim
      complex_dp, intent(out) :: cmtout(dim,*)
      complex_dp, pointer_cnt :: cmtin(:,:)
      integer                 :: itype,ieq,icent,i,ic0

      if(ispin<1.or.ispin>nspin2) Bug('ispin out of range')
      if(l_soc) Error('Not implemented for SOC.')

      if(storeibz.and.kptp(ikpt)/=ikpt) then
        if     (ic==0) then ; ic0 = ic                              ; icent = 1
        else if(ic> 0) then ; ic0 = pcent(ic,sym(symkpt(ikpt))%inv) ; icent = ic0
        else                ; ic0 = ic                              ; icent = 1 + sum(neq(:-ic-1))
        endif
        cmtin => cmtu(:,:,iband,kindx(ikpt),ispin)
        call waveftrafo_mt_io(cmtout,cmtin(1,icent),dim,maxlmindx,ic,symkpt(ikpt),kptp(ikpt))
        nullify(cmtin)
        if     (ic==0) then ; cmtout(:,:ncent)    = phase(ikpt) * cmtout(:,:ncent) ! undo division by phase factor
        else if(ic> 0) then ; cmtout(:,1     )    = phase(ikpt) * cmtout(:,1)
        else                ; cmtout(:,:neq(-ic)) = phase(ikpt) * cmtout(:,:neq(-ic))
        endif
      else
        i     = 0
        icent = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            icent       = icent + 1 ; if((ic>0.and.icent/=ic).or.(ic<0.and.itype/=-ic)) cycle
            i           = i + 1
            cmtout(:,i) = cmtu(:dim,icent,iband,kindx(ikpt),ispin)
          enddo
        enddo
      endif

      end

c     ------------------

      subroutine wavefunction_pw_w(cpwout,dim,iband,ikpt,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ikpt,iband,ispin,dim
      complex_dp, intent(out) :: cpwout(dim)
      complex_dp, pointer_cnt :: cpwin(:)

      if(ispin<1.or.ispin>nspin2) Bug('ispin out of range')
      if(l_soc) Error('Not implemented for SOC.')

      if(storeibz.and.kptp(ikpt)/=ikpt) then
        cpwin => cpwu(:,iband,kindx(ikpt),ispin)
        call waveftrafo_pw_io WSUFFIX (cpwout,cpwin,symkpt(ikpt),kptp(ikpt),ispin)
        nullify(cpwin)
      else
        cpwout = cpwu(:dim,iband,kindx(ikpt),ispin)
      endif

      end

#   endif

c --------------------------------------------------------------------------------

# endif

c --------------------------------------------------------------------------------

# undef Mcol1_
# undef McolD1_
# undef MrangeDistr_
# undef Mrange1_

c --------------------------------------------------------------------------------

# undef CMT
# undef CPW
# undef WSUB
# undef WSUFFIX
# undef LWAN
# ifdef INVW
#   undef INVW
#   define INV
# endif
# undef CMT1
# undef CPW1
# undef IKX
# undef ADJUST_KINDX_IF_NEEDED
# undef UNDO_ADJUST_KINDX
# ifdef ISPIN1
#   undef ISPIN1
#   undef ISPIN2
# endif

# ifdef PARA_MT
#   undef PARA_MT
# endif
# ifdef PARA_PW
#   undef PARA_PW
# endif
