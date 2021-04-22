c Adds the contribution of the current k-point (ikpt) (and its symmetry-equivalent points) to the matrix elements
c of the correlation self-energy (-> selfc).
c
c job1%full = .true.  : full self-energy matrix including off-diagonal elements
c                       in block-matrix form according to block(:,:)
c job1%full = .false. : only diagonal elements
c
c (1) If the self-energy is later analytically continued from the imaginary to the real frequency axis, we evaluate
c     the matrix elements <n|SIGMA(iw)|m> by frequency convolution (along the imaginary frequency axis)
c
c          s         s          s        -OMEGA        3   u/o  +inf             s                  c               s          s  ~         ~       s       s
c     < phi   | SIGMA (iw) | phi   >  =  ------ INT   d k  SUM  INT dw' (iw+iw'-E      )^(-1) SUM  W  (k,iw')  < phi      | phi   M    >  < M    phi   | phi      >
c          nq        c          mq       16pi^4    BZ       n'  -inf             n'q+k         IJ   IJ              n'q+k      mq  k,I       k,J    nq      n'q+k
c
c (2) If the frequency integration is performed by contour integration, we evaluate
c     (a) the integrand of the integral along the imaginary frequency axis <n|SIGMA(w0+iw)|m>
c         using the same formulas as above with the substitution
c
c          s         s
c         E      -> E     - w0 ,  where the Fermi energy is "shifted" downwards by w0,
c          n'q+k     n'q+k
c
c     (b) the sum over the residues, i.e., the integrated values of W at the poles of the Green function
c
c             OMEGA        3  u/o             s                 c     s               s         s  ~        ~       s       s
c         +/- ----- INT   d k SUM  theta[+/-(E     +w0)]  SUM  W  (k,E     -w0)  < phi     | phi   M    > < M    phi   | phi      >
c             8pi^3    BZ      n'             n'q+k        IJ   IJ    n'q+k           n'q+k     nq  k,I      k,J    nq      n'q+k
c
c         for w0 < 0 / w0 > 0. Both contributions are calculated at the frequencies freqr1(:nfreqr) which are calculated
c         from freqr(:nfreqr) and the band energies (see below).
c         theta is the Heaviside function.
c         Contributions (a) and (b) are summed.
c
c Input:
c   screen = screened interaction W
c   head   = head of W
c   ctrafo = transformation coefficients of eigenvectors of v wrt the mixed basis (screen is defined in the basis of eigenvectors)
c   n      = rank of W
c   ikpt   = current summation k-point (from the irreducible zone)
c
c Imported from global:
c   block  = definition of the block matrices (taking into account irreps, see irrep.f)
c   nfreq  = number of mesh points along imaginary frequency axis
c   freq   = mesh points along imaginary frequency axis
c   nfreqr = number of (real) frequencies w0 at which the self-energy is evaluated
c   freqr  = list of those frequencies (relative to efermi or the KS energy, see correlation.f)
c   nfreqc = number of (real) frequencies which are used to interpolate W (i.e., screenc and headc)
c            use Pade extrapolation for nfreqc = 0
c   freqc  = list of those frequencies
c
c The divergence at the Gamma point is integrated as in exchange.f (but anisotropy is taken into account in divergence_h/c).
c
c
c
c USE OF SYMMETRY:
c
c
c We explain the calculation of <i|GW|j>(ikpt), which is the contribution of the current k point ikpt to <i|GW|j>.
c i(r) and j(r) are states at q. (Prefactors such as "img" are omitted for simplicity.)
c
c <i|GW|j>(ikpt) = SUM(n) SUM(k0) INT dw <i (nk) | W | (nk) j >
c                = SUM(n) SUM(k0) INT dw INT j(r) nk(r)* W(r,r') nk(r') i(r')* dr dr',
c
c where we have used the symmetry W(r,r')=W(r',r). INT dw is the frequency convolution. The k0 summation runs over all k0 that
c are symmetry equivalent to ikpt. The definition is such that k=q+k0 rather than k=q-k0.
c
c Using symmetry, we can restrict the summation to the subset {k1} of {k0}. The set {k1} [kpt1(:nkpt1)] includes only those
c points k1 that are NOT pairwise related by a symop of the little group, which is defined as the set of symops P that leave
c q invariant [sym1(:nsym1)]. In other words, we have removed all symmetry equivalents (wrt sym1) from the set {k0} and have
c retained only one representive of each group [size nkpts(k1)] of symmetry-equivalent k0 points. The contribution of those
c that have been removed is taken into account by symmetrization. (The k1 are elements of EIBZ and equivalent to ikpt.)
c The little group also contains the inverses P^(-1). We use that W is invariant under the spatial part R of P, i.e.,
c W=RWR. (If P does not include time reversal (tr), R=P.)
c
c <i|GW|j>(ikpt) = SUM(n) SUM(k0) INT dw                                < i (nk) | W | (nk) j >
c                = SUM(n) SUM(k1) INT dw nkpts(k1)/nsym1   SUM(P)  <i P^(-1)(nk) | W | P^(-1)(nk) j>
c                = SUM(n) SUM(k1) INT dw nkpts(k1)/nsym1 { SUM(P,no tr) <Pi (nk) | W | (nk) Pj>
c                                                        + SUM(P,tr)    <Pj (nk) | W | (nk) Pi> } ,
c
c where k stands for q+k0 and q+k1, respectively.
c
c The tr case follows from
c <i P^(-1)(nk) | W | P^(-1)(nk) j> = <Ri (nk)* | W | (nk)* Rj> = <Pi* (nk)* | W | (nk)* Pj* > = <Pj (nk) | W | (nk) Pi> .
c Note that, if W is Hermitian (e.g., W=W(iw) or W=v), one can also write
c <i P^(-1)(nk) | W | P^(-1)(nk) j> = <Pi (nk) | W | (nk) Pj>* .
c
c Introducing the irreps P(i'i) = <i'|Pi> (defined in the subspace including state i), we have
c <i|GW|j>(ikpt) = SUM(n) SUM(k1) INT dw nkpts(k)/nsym1 SUM(i'j') { [ SUM(P,no tr) P(i'i)* P(j'j) ]  <i' (nk) | W | (nk) j'>
c                                                                 + [ SUM(P,tr)    P(i'i)* P(j'j) ]* <j' (nk) | W | (nk) i'> }
c                = SUM(n) SUM(k1) INT dw nkpts(k)/nsym1 SUM(i'j') {       irrep_contr(ij,i'j')       <i' (nk) | W | (nk) j'>
c                                                                 +       irrep_contt(ij,j'i')       <j' (nk) | W | (nk) i'> } ,
c which defines the contracted irreps irrep_contr and irrep_contt.
c
c In the following, we discuss the evaluation of <i(nk)|W|(nk)j>. The operator W is given in an auxiliary Bloch basis,
c <a|W|b>. The auxiliary basis can be the mixed basis (W=v) or the basis of Coulomb eigenvectors (W=W(w)).
c The matrices of W are only known in the irreducible wedge. Therefore, we use a suitable transformation [symkpt(k1)]
c P: IBZ->EIBZ, ikpt->k1, PI=I', Pa=a', etc.. In the following, summing over I, J, a, b is implicit.
c
c First case: P does not include tr:
c   <i|(nk)W(nk)|j> = <nk|jI'> <I'|a'> <a'|W|b'> <b'|J'> <J'i|nk>
c                   = <nk|jPI> <PI|Pa> <Pa|W|Pb> <Pb|PJ> <PJi|nk>
c                   = <nk|jPI>  <I|a>   <a|W|b>   <b|J>  <PJi|nk> .
c
c Second case: P includes tr:
c   <i|(nk)W(nk)|j> = <nk|jI'> <I'|a'> <a'|W|b'> <b'|J'> <J'i|nk>
c                   = <nk|jPI> <PI|Pa> <Pa|W|Pb> <Pb|PJ> <PJi|nk>
c                   = <nk|jPI>  <a|I>   <b|W|a>   <J|b>  <PJi|nk>
c                   = <nk|iPJ>* <J|b>   <b|W|a>   <a|I>  <PIj|nk>* ,
c
c which differs from the first case only by a complex conjugation of the vectors and a transpose operator. If W is Hermitian
c (i.e, for W(iw) and v), we may also use
c
c   <i|(nk)W(nk)|j> =   <nk|jPI> <PI|Pa> <Pa|W|Pb> <Pb|PJ> <PJi|nk>
c                   =   <nk|jPI>  <I|a>*  <a|W|b>*  <b|J>* <PJi|nk>
c                   = [ <nk|jPI>* <I|a>   <a|W|b>   <b|J>  <PJi|nk>* ]* ,
c
c i.e., we may replace the transpose operator by a complex conjugate.
c
c When testing symmetry, note that Pade interpolation is not a linear operation. This means that CONTOUR with Pade and FREQINT PADE
c will produce slightly different results if symmetry is switched off or not (switch_off_symmetry macro in susceptibility.f).
c
c For the same reason, one can get different results in the serial and parallelized version: In the serial version, the Green-function
c band summation is (degenerate-)subspace-wise and the k-point summation can be over the full Brillouin zone instead of the EIBZ.
c In the parallelized version, the former is band-wise and the latter is always over EIBZ. The sum is identical, but the non-linear
c Pade operation on summands can lead to deviations.
c
c Note:
c - wings of W are used only to give the Gamma point contribution to the body of W
c - on exit, the arrays screen and screenc are destroyed.
c

c Only calculate one of the identical diagonal self-energy matrix elements of
c degenerate states and copy the others (only serial version and .not.full).
c This can be disabled here by commenting out the following.
# define COPYDEG

# ifndef COPYDEG
#   warning COPYDEG unset for testing
# endif

# ifdef MPI
#   define NFENCE Nfence(screen) ; if(nfreqc>0) then ; Nfence(screenc) ; endif
# else
#   define NFENCE continue
# endif

# if defined(MPI) && !defined(LOAD)
#   define INC Msize
# else
#   define INC 1
# endif

# include "cppmacro.h"
# include "jobtype.h"

c begin interface
      subroutine selfenergy(job1,ikpt,eval,n,coul,Win(ctrafo),Win(screen), head, wing,
     &                                                        Win(screenc),headc,wingc,plasma)

      use global !inc
      use arrays !inc
      use util
      use wrapper
      use freq_integral
      use file
      use, intrinsic :: iso_fortran_env
      Mpi ( use Mwrapper )
      Mpi ( use key )
      Mpi2( use, intrinsic :: iso_c_binding )
      Load( use readwrite ) LoadC( only: read_wavef0 ) LoadC( read_wavef2 )

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
c end interface
      complex_dp,    allocatable :: matc(:,:,:),cprod1w(:),cvec(:),cprod_ibc(:,:,:,:)
      MCOMPLEX_dp,   allocatable :: mat(:,:,:),matx(:,:),vec(:)
      MCOMPLEX_dp,   allocatable :: cprod(:,:,:,:),cprod1(:),irrep_contr(:,:,:,:),irrep_contt(:,:,:,:)
      MCOMPLEX_dp,   allocatable :: moment(:,:,:),moment_diag(:,:)
      MCOMPLEX_dp,   pointer_cnt :: eigv(:,:),olap(:,:)
      real_dp,       allocatable :: eig(:)
      real_dp,       allocatable :: wintgrc(:,:,:,:)
      integer                    :: nwblock
      integer,       allocatable :: band(:),pnt(:),wblock(:)
      complex_dp,    allocatable :: wfreqintegral(:,:,:),wfreqintegral1(:,:,:),wfreqintegral2(:,:,:)
      complex_dp,    allocatable :: pm(:,:,:,:),pv(:,:,:),zm(:,:,:,:),cpade_head(:,:,:)
      complex_dp,    allocatable :: selfc_ibc(:,:,:)
      complex_dp                 :: h0(3,3),h1(3,3),h2(3,3)
      complex_dp                 :: cpade(nfreq+1),pole(nfreq*smooth(2)),resid(nfreq*smooth(2)),wfreq(0:3,nfreq)
      complex_dp                 :: cdum,cdum1,cdum2
      complex_dp                 :: cmat(3,3),hw(3,3,max(nfreq,nfreqr),0:2),divergence_c(nfreqc),iheadc(3,3,nfreqc)
      MCOMPLEX_dp                :: aspline(0:3,nfreq)
      MCOMPLEX_dp                :: ihead(6,nfreq)
      real_dp                    :: dene(3),dene2(3,3),ddene(3,3),divfac_x,divergence_x,divergence_h(nfreq)
      real_dp                    :: eneq,enek,enediff,freqr1(nfreqr),weight,rdum,rdum1,coeff(nfreqc)
      real_dp                    :: frq(27),pfrq(0:3,26),emin,emax,f1,f2,mom(3),logdiv1,logdiv2
      integer                    :: nfrq
      integer                    :: nkpt1,kpt1(nkpt),sym1(nsym),nsym1,nkpts(nkpt),nkpt0
      integer                    :: minb(nfreqr,nspin1),maxb(nfreqr,nspin1),npole1
      integer                    :: nfreq1,ifreq,ifreqr,ifreqc,ifreqc1,ifreqc2
      integer                    :: iself,iselfc,iselfx,iblock,ib,ibandq,ikptq,ispin,s,iband,ideg,id
      integer                    :: ndeg,deg0,deg1,band1,band2,bandup,pack,maxpack,nbnd,ibnd
      integer                    :: isub,jsub,isym
      integer                    :: ikpt1,ikpt2,ikpt2_
      integer                    :: i,j,k,l,m,p,ka,kb,kptp2,nkpt_2,sgn
      integer                    :: ipro,npro
      logical                    :: drude_separate,newkptq,evalx,evalc,ldum
      real                       :: time_mt,time_pw,time_trafo,time_exch,time_gamma,time_mat,time_equiv,time_freq,time_ibc,time_tot
      real                       :: time_tetra MpiC(time_idle) MpiC(time_maxidle)
      real                       :: time1,time2
      integer                    :: kptsum
      complex_dp                 :: residues,residues_pade
      complex_dp                 :: pade_func
      real_dp                    :: logdiv
# ifdef ifort_vector_bug
      real_dp                    :: ggg1(100),ggg2(100),g1,g2
# endif
# ifdef LOAD
      integer                    :: kpt2(nkpt)
# endif
# ifdef MPI
      integer                    :: win_eigv,win_olap
      logical                    :: lsplit
      logical, save              :: first = .true.,mpisym = .false.
      real_dp, save              :: Mover
      integer                    :: Mcolor,Merr,Msub(size(block,2)),Mstat(mpi_status_size),rank
      MCOMPLEX_dp,   allocatable :: math(:,:,:)
      type(c_ptr)                :: ptr
      if(first) then
        Rbegin
        if(job1%type==J_GW) then ; Mover = 0  ! parallelize over n'     (default for GW)
        else                     ; Mover = 10 ! parallelize over blocks (default for HF, ...)
        endif
        call getkey(inp,'MPIBLK',Mover,section='SENERGY',status=i,mini=0d0)
        if(i==1) Mover = 10 ! default
        if(nfreqc==0) call getkey(inp,'MPISYM',mpisym,section='SENERGY',default=.false.)
        if(mpisym.and.nfreqc/=0.and.freqint==2) then
          Info('MPISYM disabled for present case (nfreqc/=0 and FREQINT PADE)')
          mpisym = .false.
        endif
        Rend
        first = .false.
        call Mcast(Mover)
        call Mcast(mpisym)
      endif
# endif
      
# ifdef CHECK_SENERGY
#   warning CHECK_SENERGY defined
                    Rwrite(6,'(''screen: '',F30.14)') sum(abs(screen)**2)
      Rif(ikpt==1)   write(6,'(''head:   '',F30.14)') sum(abs(head)**2)
      Rif(nfreqc>0) then
                     write(6,'(''screen:'',F30.14)') sum(abs(screenc)**2)
        Rif(ikpt==1) write(6,'(''head:  '',F30.14)') sum(abs(headc)**2)
      endif
      Rwrite(6,*)
# endif

      evalx = eval(1)
      evalc = eval(2)

      if(all(job1%type/=[J_HF,J_GW,J_RPA,J_HFE,J_SX,J_COSX,J_PBE0]))
     &                                  Bug('Wrong type of calculation (job1%type).')
      if(.not.(evalx.or.evalc))         Bug('evalx and evalc both false.')
      if(evalc.and.nfreqr==0)           Bug('nfreqr is zero.')
      if(nfreqc>0) then
        if(any(abs(imag(freqc))>1d-10)) Bug('Nonzero imaginary part in freqc.')
      endif

      Rif(ikpt==1.and.evalc.and.freqint>=2) then
        if(ozero)     Info('FREQINT SPLINE assumed for zero-order corrections.')
        if(job1%full) Warn('FREQINT NONLIN never tested for FULL calculations. Use with care!')
      endif

      if(oibc/=0.and.allocated(cblock)) then
        write(*,*) 'cblock deallocated for keyword IBC'
        deallocate ( cblock )
      endif

      maxpack    = 1
      time_mt    = 0
      time_pw    = 0
      time_trafo = 0
      time_exch  = 0
      time_gamma = 0
      time_mat   = 0
      time_equiv = 0
      time_freq  = 0
      time_ibc   = 0
      time_tetra = 0

c
c     Define inverse head (-> ihead/c)
c
      if(ikpt==1) then
        if(present(screen)) then
          Mpi( ihead = 0 )
          do i = Mrange(1,nfreq)
            h0         = unpackmat(head(:,i)) ; call invert_angular(h0)
            ihead(:,i) = 4*pi * packmat(h0)
          enddo
          Mpi( call Msum(ihead) )
          if(nfreqc>0) then
            Mpi( iheadc = 0 )
            do i = Mrange(1,nfreqc)
              h0            = headc(:,:,i) ; call invert_angular(h0)
              iheadc(:,:,i) = 4*pi * h0
            enddo
            Mpi( call Msum(iheadc) )
          endif
        endif
      endif

      call cpu_time(time_tot)

c
c     Define Gamma divergence (->divergence/_h/_c)
      if(ikpt==1) then
        call cpu_time(time1)
        if(divergence==0) then
          RWarn('Contribution of Coulomb divergence has not been calculated yet. This should not happen!')
          call gamma_divergence(.false.)
        endif
        if(present(screen)) then
          do i = 1,nfreq  ; rdum = trace(ihead(:,i))/3    ; call gamma_divergence_h(divergence_h(i),head(:,i),   rdum) ; enddo
          do i = 1,nfreqc ; cdum = trace(iheadc(:,:,i))/3 ; call gamma_divergence_c(divergence_c(i),headc(:,:,i),cdum) ; enddo
        endif
        call cpu_time(time2) ; time_gamma = time_gamma + time2 - time1
      endif

c
c     Prepare Coulomb matrix for exchange calculation (->coulomb0)
c     - case SX and COSX: calculate static screened interaction
c     - if needed, multiply with the inverse of the overlap matrix for current k point
# define NC Ncol(1,ngptm(ikpt))
      if(evalx) then
        call cpu_time(time1)
        Nfence(coulomb0)
        Obegin
        if(any(job1%type==[J_SX,J_COSX])) call unitarytrafo ( coulomb0 , screen(:,1) , ctrafo , 2 )
        if(ikpt==1) then
          if(any(job1%type==[J_SX,J_COSX])) then ; divergence_x =   divergence_h(1) / vol ; divfac_x = trace(ihead(:,1))/3 / vol
          else                                   ; divergence_x = 4*pi * divergence / vol ; divfac_x = 4*pi / vol
          endif
        endif
        Oend
        Nfence(coulomb0)
        Mpi( call Mcast(divergence_x) ; call Mcast(divfac_x) )
        if(      fullpw .and. all(job1%type/=[J_SX,J_COSX]) .or.
     &      .not.fullpw .and. any(job1%type==[J_SX,J_COSX]) ) then
          allocate  ( eig(ngptm(ikpt)) )
          Nallocate ( eigv,(S_ ngptm(ikpt),ngptm(ikpt) S_) )
          Nallocate ( olap,(S_ ngptm(ikpt),ngptm(ikpt) S_) )
          Ocall olap_gptm(olap,ngptm(ikpt),ikpt,ikpt)        ; Nfence(olap)
          ! Calculate pseudoinverse of olap -> olap**(-1)
          call Mdiagonalize(eigv,eig,olap)
          if(any(eig<-cutzero)) Error('Overlap matrix not positive definite.')
          where(eig<cutzero) eig = 0
          where(eig/=0)      eig = 1/eig
          Nfence(eigv)
          do i = Nrange(1,ngptm(ikpt))
            eigv(:,i) = sqrt(eig(i)) * eigv(:,i)
          enddo
          Nfence(eigv)
          Nfence(olap) ; if(size(eigv(NC,:))>0) olap(:,NC) = matmac(eigv,eigv(NC,:))
          Nfence(olap)
          Ndeallocate ( eigv )
          ! Multiply coulomb0 with olap**(-1) from both sides
          Nallocate ( eigv,(S_ nbasm(ikpt),nbasm(ikpt) S_) ) ! copy of coulomb0
          ifO call p_unpackmat(eigv,coulomb0)
# ifndef INV
          Rbegin
          rdum = 0
          do i = 1,nbasm(ikpt)
            rdum      = rdum + abs(imag(eigv(i,i)))
            eigv(i,i) = dble(eigv(i,i))
          enddo
          rdum = rdum / nbasm(ikpt)
          if(rdum>1d-6) Error('Very large imaginary part on Coulomb diagonal: '//Chf(rdum,'F20.12'))
          if(rdum>1d-12) Warn('Large imaginary part on Coulomb diagonal: '//Chf(rdum,'F20.12'))
          Rend
# endif
# undef NC
# define NC Ncol(1,nbasm(ikpt))
          Nfence(eigv) ; eigv(nbasp+1:,NC) = matmat( olap , eigv(nbasp+1:,NC) )
          Nfence(eigv) ; eigv(NC,nbasp+1:) = matmat( eigv(NC,nbasp+1:) , olap )
          Nfence(eigv)
          Rbegin
          rdum = 0
          do j = 1,nbasm(ikpt)
            do i = 1,j
              rdum = rdum + abs(eigv(i,j)-MCONJG(eigv(j,i)))
            enddo
          enddo
          rdum = rdum / (nbasm(ikpt)*(nbasm(ikpt)+1)/2)
          write(6,'(A,ES8.1)') 'Hermiticity of transformed Coulomb matrix:',rdum
          if(rdum>1d-8) then
            Warn('Large average deviation from Hermiticity in Coulomb: '//Chf(rdum,'F11.8')//'. CUTZERO might help.')
          endif
          Rend
          Nfence(coulomb0)
          ifO call p_packmat ( coulomb0, eigv, .false. )
          Nfence(coulomb0)
          deallocate (eig)
          Ndeallocate(eigv)
          Ndeallocate(olap)
        endif
        call cpu_time(time2) ; time_exch = time_exch + time2 - time1
      endif
# undef NC

c
c     Define  W^c = W - v  (->screen,ihead/w,divergence_h/c) in order to remove the divergence
      if(evalc.or.job1%type==J_COSX) then
        NFENCE
        if(nfreq>0) then
          do i = 1,n
            ifO screen(i*(i+1)/2,:) = screen(i*(i+1)/2,:) - coul(i)
          enddo
          if(ikpt==1) then
            ihead(1,:) = ihead(1,:) - 4*pi
            ihead(3,:) = ihead(3,:) - 4*pi
            ihead(6,:) = ihead(6,:) - 4*pi
          endif
        endif
        if(oselfc/=1.and.nfreqc>0) then
          do i = 1,n
            ifO screenc(i,i,:) = screenc(i,i,:) - coul(i)
          enddo
          if(ikpt==1) then
            iheadc(1,1,:) = iheadc(1,1,:) - 4*pi
            iheadc(2,2,:) = iheadc(2,2,:) - 4*pi
            iheadc(3,3,:) = iheadc(3,3,:) - 4*pi
          endif
        endif
        if(ikpt==1) then
          if(nfreq >0) divergence_h = divergence_h - 4*pi * divergence
          if(nfreqc>0) divergence_c = divergence_c - 4*pi * divergence
        endif
        NFENCE
      endif

      if(any(job1%type==[J_SX,J_COSX])) then
        if(nfreq/=1) Error('nfreq not set to 1. Frequency mesh read from spex.cor?')
        nfreq = 0 ! we only need the static screened interaction, which is stored in coulomb0; avoid loops over ifreq=1,..,nfreq
      endif

c
c     WBLOCK (see susceptibility): Decompose W into blocks (redefine screen,screenc,ctrafo)
      if(allocated(cblock).and.evalc) then
        NFENCE ; Nfence(ctrafo)
        nwblock = maxval(cblock) ; allocate ( wblock(nwblock),pnt(n) )
        k       = 0
        do i = 1,nwblock
          wblock(i) = count(cblock==i)
          do j = 1,n
            if(cblock(j)==i) then
              k      = k + 1
              pnt(k) = j
            endif
          enddo
        enddo
        Obegin
        ctrafo = ctrafo(:,pnt)
        do i = 1,nfreq
          call block_diag_pack(screen(:,i))
        enddo
        do i = 1,nfreqc
          do j = 1,n ; screenc(:,j,i) = screenc(pnt,j,i) ; enddo
          do j = 1,n ; screenc(j,:,i) = screenc(j,pnt,i) ; enddo
        enddo
        Oend  
        deallocate ( pnt )
        NFENCE ; Nfence(ctrafo)
      endif

      if(.not.evalc)     then ; nfreq1 = 0
      else if(oselfc/=1) then ; nfreq1 = 1     ! contour integration:   self-energy arguments  w + img*freq(:1)      (freq(1)=0)
      else                    ; nfreq1 = nfreq ! analytic continuation: self-energy arguments      img*freq(:nfreq)
      endif

      if(evalc.and.(freqint<=1.or.ozero)) then
        if(oselfc==1) then ; m = nfreq
        else               ; m = nfreqr
        endif
        allocate ( wfreqintegral(0:3,nfreq,m) )
      endif

c
c     Remove Drude term, which is treated analytically
      if(ikpt==1.and.evalc) then
        ! Does Drude term exist?
        drude_separate = metal.and.plasma/=0
        ! Write head of screened interaction to spex.head
        Rif(wrtinfo) then
          i = fopen('spex.head',status='unknown',numbered=.true.)
          write(i,'(A)')       '# Head element of screened interaction on imaginary axis'
          if(drude_separate)
     &      write(i,'(A,F10.5)') '# Compare with -4pi * plasma**2 / (plasma**2 + freq**2),  plasma =',plasma
          write(i,'(F10.5,F15.10)') (freq(ifreq),dble((trace(ihead(:,ifreq)))/3),ifreq=1,nfreq)
          call fclose(i)
        endif
        ! If Padé is used, we cannot separate off the Drude term.
        drude_separate = drude_separate.and.freqint<2.and.(oselfc==1.or.nfreqc/=0)
        ! Check if separating off makes the head smoother
        if(drude_separate) then
          rdum  = 0
          rdum1 = 0
          do ifreq = 1,nfreq-1
            enediff = 4*pi * plasma**2 * ( ( freq(ifreq+1)**2 + plasma**2 )**(-1) -
     &                                     ( freq(ifreq)  **2 + plasma**2 )**(-1) )
            rdum    = rdum  + abs( trace(ihead(:,ifreq+1)-ihead(:,ifreq))/3           )
            rdum1   = rdum1 + abs( trace(ihead(:,ifreq+1)-ihead(:,ifreq))/3 + enediff )
          enddo
          Rwrite(6,'(A,F5.1,A,F5.1)') 'Variation in W(head):',rdum,' ->',rdum1
          if(rdum1>rdum) drude_separate = .false.
        endif
        ! Remove Drude term from head and Gamma divergence
        if(drude_separate) then
          Rwrite(6,'(A)') 'Drude term is separated off and treated analytically.'
          do ifreq = 1,nfreq
            rdum                = -plasma**2 / ( plasma**2 + freq(ifreq)**2 )
            ihead(:,ifreq)      = ihead(:,ifreq) - packmat ( 4*pi * rdum * identity(3) )
            divergence_h(ifreq) = divergence_h(ifreq) - divergence * 4*pi * rdum
          enddo
        endif
      endif

      Rbegin
      write(6,'(/A'NoA) 'Add contribution to self-energy'
      if(.not.(evalc.and.job1%full)) write(6,'(A'NoA) '... '
      Rend

c
c     Pade coefficients from head (-> cpade_head)
      if(ikpt==1.and.ozero.and.evalc.and.oselfc/=1.and.nfreqc==0) then
        allocate ( cpade_head(nfreq+1,3,3) )
        cpade_head = 0
        k          = 0
        do j = 1,3
          do i = 1,j
            k = k + 1
            if(all(abs(ihead(k,:))>1d-10)) then
                       call pade_init(cpade_head(:,i,j),img*freq,       ihead(k,:) *(1d0,0d0),nfreq,-1)
              if(i/=j) call pade_init(cpade_head(:,j,i),img*freq,MCONJG(ihead(k,:))*(1d0,0d0),nfreq,-1)
            endif
          enddo
        enddo
      endif

c
c     Test if screenc is symmetric (case inversion symmetry)
# ifdef INV
      if(evalc) then
        do ifreqc = 1,nfreqc
          if(sum(abs(screenc(:,:,ifreqc)-transpose(screenc(:,:,ifreqc))))>1d-10)
     &      Bug('Screened interaction not symmetric.')
        enddo
      endif
# endif

c
c     Distribute blocks over processes (Mover/=0)
# ifdef MPI
      call Mdistribute
      i = Mrank
      if(Mcolor/=0) call begin_split(Mcolor)
# endif

c
c     Initialize progress bar
      if(evalc.and.job1%full Mpi(.and.Mrank==0.and.Mcolor<=1) ) then
        ipro = 0
        npro = 0
        do iblock = 1,size(block,2)
          Mpi( if(Msub(iblock)/=Mcolor) cycle )
          ikptq = job1%kpt(block(1,iblock))
          ispin = job1%spin(block(1,iblock))
          call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikptq,ikpt,.false.)
          do k = 1,nkpt1
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikpt)
            do i = 1,size(block,1)
              if(block(i,iblock)/=0) then
                do iband = 1,nband(ikpt2,ispin) ; Mcycle(iband-1)
                  npro = npro + 1
                enddo
              endif
            enddo
          enddo
        enddo
        npro = max(npro,1) ! avoid division by zero
      endif

c
c
c     LOOP OVER BLOCKS   ( job1%full = .true.  : block matrices of full self-energy matrix;
c                          job1%full = .false. : dummy definition: each state is a "block" of its own )

      ikptq = 0
      ispin = 0

      do iblock = 1,size(block,2)

        Mpi( if(Msub(iblock)/=Mcolor) cycle )

        iselfc = 0
        iselfx = 0
        do i = 1,iblock-1
          m      = count(block(:,i)/=0)
          iselfc = iselfc + m**2
          iselfx = iselfx + m*(m+1)/2
        enddo

c        write(*,'(5I5)') iblock,Mcolor,Mrank,iselfc,size(selfc,1)

        ib = block(1,iblock) ! first band of block

c       Just copy the self-energy if the state is degenerate with the previous one
# if !defined(MPI) && defined(COPYDEG)
        if(.not.job1%full.and.freqint/=2) then
          if(ib/=1) then
            if(ikptq==job1%kpt(ib).and.ispin==job1%spin(ib)) then
              if(same_eigenspace(ibandq,job1%band(ib),ikptq,ispin)) then
                iselfc = iselfc + 1
                iselfx = iselfx + 1
                Oif(evalc) selfc(iselfc,:) = selfc(iselfc-1,:)
                if (evalx) selfx(iselfx)   = selfx(iselfx-1)
                cycle
              endif
            endif
          endif
        endif
# endif

        newkptq = job1%kpt(ib)/=ikptq.or.job1%spin(ib)/=ispin
        ibandq  = job1%band(ib)
        ikptq   = job1%kpt(ib)
        ispin   = job1%spin(ib)

c
c       Get irreducible k points kpt1(:nkpt1) wrt current q point ikptq (EIBZ)
        if(.not.job1%full) then
          call getkpt1_fullBZ(kpt1,nkpt0,nkpts,sym1,nsym1,ikpt) ! determine nkpt0
        endif
        call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikptq,ikpt,.false.)

c
c       Determine bands that are contained in current block (->band)
        if(job1%full) then
          ! according to block (full self-energy calculations)
          i = size(block,1)
          do while(block(i,iblock)==0)
            i = i - 1
          enddo
          allocate ( band(i) )
          band = job1%band(block(:i,iblock))
        else
          ! Determine first and last state of degenerate subspace (only diagonal self-energy elements)
          deg0 = ibandq
          deg1 = deg(ibandq,ikptq,ispin)
          if(deg1<deg0) then
            deg0 = deg1
            deg1 = deg(deg0,ikptq,ispin)
          endif
          ndeg = deg1 - deg0 + 1
          ! double spin degeneracies (SOC+inv.sym.) are not due to spatial symmetry and don't have to be averaged over: leave out.
          if(l_soc.and.invsym/=0.and.ndeg==2) then
            ndeg = 1
            deg1 = deg0
          endif
          ! define array band accordingly
          if(ndeg*nkpt1>=nkpt0 Load(.and..false.) ) then ! in this case summing over all k points is faster (disabled for LOAD)
            allocate ( band(1) )
            band = ibandq
            deg0 = ibandq
            deg1 = ibandq
            call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt)
          else
            allocate ( band(ndeg) )
            band = [ (i,i=deg0,deg1) ]
          endif
        endif

c
c       Define frequency arguments w of self-energy SIGMA(w)  (-> freqr1)
        if(evalc) then
          if(oselfc==1.or.oselfc==4) then
            freqr1 = freqr
          else
            eneq = ene(ibandq,ikptq,ispin)
            if(job1%full)
     &        Error('Full self-energy calculations only implemented for CONTINUE and CONTOUR [{...}] <...>')
            if(eneq<=efermi) then ; freqr1 = eneq - efermi - freqr
            else                  ; freqr1 = eneq - efermi + freqr
            endif
          endif
        endif

c        if(option==1) goto 123 ! skip val/cond (testing)

c
c       Initialize addition of equivalent k points (define contracted irreps -> irrep_contr/t)
        if(use_sym.and.job1%full) then!.and.any(nkpts(:nkpt1)/=1)) then
          call cpu_time(time1)
          ndeg = deg(ibandq,ikptq,ispin) - ibandq + 1 ! dimension of irrep
          if(.not.trsoff) then
            allocate ( irrep_contt(ndeg**2,ndeg**2,size(band)/ndeg,size(band)/ndeg) )
            irrep_contt = 0
          endif
          allocate ( irrep_contr(ndeg**2,ndeg**2,size(band)/ndeg,size(band)/ndeg) )
          irrep_contr = 0
          do i = 1,size(band)/ndeg
            do j = 1,size(band)/ndeg
              isub = psub(block((i-1)*ndeg+1,iblock))
              jsub = psub(block((j-1)*ndeg+1,iblock))
              do k = 1,ndeg ; do l = 1,ndeg ; do m = 1,ndeg ; do p = 1,ndeg
                deg0 = k+(p-1)*ndeg
                deg1 = l+(m-1)*ndeg
                do isym = 1,nsym1
                  cdum = conjg ( irrep_sub(l,k,isub,sym1(isym)) ) * irrep_sub(m,p,jsub,sym1(isym))
                  if(sym1(isym)>nsymt) then
                    irrep_contt(m+(l-1)*ndeg,k+(p-1)*ndeg,i,j) = irrep_contt(m+(l-1)*ndeg,k+(p-1)*ndeg,i,j) + conjg(cdum)
                  else
                    irrep_contr(l+(m-1)*ndeg,k+(p-1)*ndeg,i,j) = irrep_contr(l+(m-1)*ndeg,k+(p-1)*ndeg,i,j) + cdum
                  endif
                enddo
              enddo ; enddo ; enddo ; enddo
            enddo
          enddo
          call cpu_time(time2) ; time_equiv = time_equiv + time2 - time1
        endif

c
c       Gamma point
        if(ikpt==1) then
          call cpu_time(time1)
          if(ozero) then
c           Calculate momentum matrix <nq|-i\nabla|n'q>
            Load( allocate ( cmtq(maxlmindx,ncent,nband(ikptq,ispin),nspin3) ) )
            Load( allocate ( cpwq(maxgpt,         nband(ikptq,ispin),nspin3) ) )
            Load( call read_wavef2([(i,i=1,nband(ikptq,ispin))],nband(ikptq,ispin),ikptq,ispin,cmtq,cpwq)  )
            allocate ( moment(nband(ikptq,ispin),size(band),3) )
            do i = Mrange1(size(band))
              call momentum_matrix(moment(:,i,:),[ikptq],1,ispin,ispin,band(i),band(i),1,nband(ikptq,ispin),.false.
     &                             MpiC(.false.) LoadC(.true.) )
            enddo
            MrangeDistr( moment (:, McolD1(size(band),i) ,:) ,i )
            if(evalc.and.newkptq) then ! all diagonal elements needed for correlation
              if(allocated(moment_diag)) deallocate ( moment_diag )
              allocate ( moment_diag(nband(ikptq,ispin),3) )
              do iband = Mrange1(nband(ikptq,ispin))
                call momentum_matrix(moment_diag(iband,:),[ikptq],1,ispin,ispin,iband,iband,iband,iband,.false.
     &                               MpiC(.false.) LoadC(.true.) )
              enddo
              MrangeDistr( moment_diag ( McolD1(nband(ikptq,ispin),i) ,:) ,i)
            endif
          endif
          Rbegin
c         (A) Pole contribution 1/k**2 at k->0 (head element)
          ! Correlation
          if(evalc) then
            if(freqint>=2.or.nfreqc==0) call pade_init(cpade,img*freq,divergence_h*(1d0,0d0),nfreq,-1)
            if(freqint>=2) call pade_poles(pole,resid,npole1,img*freq,divergence_h*(1d0,0d0),cpade,nfreq,smooth(2),-1,.false.)
          endif
          do ib = 1,size(band)
            eneq = ene(band(ib),ikptq,ispin) - efermi
            if(ologdiv) then
              if(allocated(moment)) then
                mom = moment(band(ib),band(ib),:)
              else
                call momentum_matrix(mom,[ikptq],1,ispin,ispin,band(ib),band(ib),band(ib),band(ib),.false. MpiC(.false.) )
              endif
              logdiv1 = logdiv(eneq,mom)
            endif
            ! Correlation
            if(evalc) then
              iself = iselfc + (ib-1)*size(band) + ib
              do ifreqr = 1,nfreqr
                enediff = eneq - freqr1(ifreqr)
                if(ologdiv) logdiv2 = logdiv(enediff,mom)
                ! Frequency convolution
                do ifreq = 1,nfreq1
                  if     (freqint<=1) then ; cdum = freqintegral(divergence_h,freq,nfreq,1,freq(ifreq)+img*enediff,0)
                  else if(freqint>=2) then ; cdum = freqintegral_poles(pole,resid,npole1,1,freq(ifreq)+img*enediff,0)
                  endif
                  if(drude_separate) cdum = cdum + img*plasma/2 * drudefunc(freq(ifreq)+img*enediff,1) * 4*pi * divergence
                  if(ologdiv)        cdum = real(cdum) * abs(2*logdiv2-1) + img * imag(cdum) ! 2*logdiv-1 switches from 1 (enediff<<0) to -1 (enediff>>0); abs() because, if 2*logdiv-1<0, the frequency integral has switched sign of the real part
                  if(oselfc==1) then ; Nacc1_c( selfc,(iself,ifreq),  -cdum/vol ) ! The minus sign (-cdum/vol) comes from the prefactor i 
                  else               ; Nacc1_c( selfc,(iself,ifreqr), -cdum/vol ) ! and the integration differential d(iw) = i*dw.
                  endif
                enddo
                ! Residues / eneq*enediff<=0: (un)occupied states require a positive (negative) argument to W, i.e., eneq and enediff must not have identical signs.
                ldum = ologdiv
                if(ldum) ldum = freqr1(ifreqr)<0  .and. logdiv1*(1-logdiv2)>1d-8 .or.
     &                          freqr1(ifreqr)>=0 .and. logdiv2*(1-logdiv1)>1d-8                                
                if(oselfc/=1.and.(eneq*enediff<=0.or.ldum)) then
                  if(nfreqc==0) then
                    cdum = pade_func((1d0,0d0)*abs(enediff),img*freq,cpade,nfreq)
                  else
                    call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
                    cdum = sum(coeff*divergence_c)
                  endif
                  if(ologdiv) then
                    if(freqr1(ifreqr)<0) then ; cdum = cdum * logdiv1 * (1-logdiv2)
                    else                      ; cdum = cdum * logdiv2 * (1-logdiv1)
                    endif
                  endif
                  if(enediff==0)       cdum =  cdum / 2 ! the other half comes from the integral over iw
                  if(freqr1(ifreqr)>0) cdum = -cdum     ! integration runs the perimeter of the rectangle clock-wise (factor -1)
                  Nacc1_c( selfc,(iself,ifreqr), -cdum/vol ) ! The minus sign comes from the prefactor i/(2*pi) times the factor 2*pi*i from the residues
                endif
              enddo
            endif
            ! Exchange (case NOSTORE)
            if(evalx) then
              if(lomit(1)) then ; if(any(omit==band(ib))) cycle ; endif
              iself = iselfx + ib*(ib+1)/2
              if(ologdiv)     then ; selfx(iself) = selfx(iself) - divergence_x * logdiv1
              else if(eneq<0) then ; selfx(iself) = selfx(iself) - divergence_x
              endif                
              if(job1%type==J_COSX)  selfx(iself) = selfx(iself) + (divergence_x-4*pi/vol*divergence)/2
            endif
            if(.not.job1%full) exit
          enddo
          Rend
c         (B) Zero-order contribution 1/k**2 * k**2 at k->0 (OZERO)
          if(ozero) then
            ! Loop over bands
            if(oselfc==1) then ; allocate ( wfreqintegral1(0:3,nfreq,nfreq), wfreqintegral2(0:3,nfreq,nfreq) )
            else               ; allocate ( wfreqintegral1(0:3,nfreq,nfreqr),wfreqintegral2(0:3,nfreq,nfreqr) )
            endif
            allocate ( pm(3,3,size(band),size(band)),pv(3,size(band),size(band)),zm(3,3,size(band),size(band)) )
            ! Calculate second order term of <E phi_i | phi_n> = O(1) + O(k) + k*zm*k (->zm)
            zm = 0
            do j = 1,size(band)   ! phi_n (zm needed only if n is one of band)
              do i = 1,size(band) ! phi_i
                ldum = same_eigenspace( band(i) , band(j) , ikptq,ispin ,edeg)
                do iband = 1,nband(ikptq,ispin)
                  cmat = reshape ( [ ((moment(iband,i,k)*MCONJG(moment(iband,j,l)),k=1,3),l=1,3) ] , [3,3] )
                  if(ldum) then ! band(i), band(j) degenerate
                    if(.not.same_eigenspace( iband , band(j) , ikptq,ispin ,edeg)) ! contributes only if iband not degenerate to band(i),band(j)
     &                zm(:,:,i,j) = zm(:,:,i,j) - cmat /   (ene(band(j),ikptq,ispin) - ene(iband, ikptq,ispin))**2 / 2
                  else          ! band(i), band(j) not degenerate                     
                    if(same_eigenspace( iband , band(j) , ikptq,ispin ,edeg)) then ! contribution if iband, band(j) degenerate
                      zm(:,:,i,j) = zm(:,:,i,j) - cmat /   (ene(band(j),ikptq,ispin) - ene(band(i),ikptq,ispin))**2
                    else                                                           ! contribution if iband, band(j) not degenerate
                      zm(:,:,i,j) = zm(:,:,i,j) + cmat / ( (ene(band(j),ikptq,ispin) - ene(band(i),ikptq,ispin)) *
     &                                                     (ene(band(j),ikptq,ispin) - ene(iband, ikptq,ispin)) )
                    endif
                  endif
                enddo
              enddo
            enddo
            if(evalc) then ; bandup = nband(ikptq,ispin)
            else           ; bandup = bando
            endif
            m    = 0
            deg0 = 1
            do while(deg0<=bandup)
              enek = ene(deg0,ikptq,ispin) - efermi
              deg1 = deg0
              do while(deg1<bandup)
                if(abs(ene(deg1+1,ikptq,ispin)-ene(deg0,ikptq,ispin))>edeg) exit
                deg1 = deg1 + 1
              enddo
              if(deg1<deg0)   Bug('Count error.')
              if(deg1>bandup) Bug('deg1 exceeded upper bound.')
              ifMODP(m)
              if(evalc) then
                ! Initialize frequency convolutions
                if(oselfc==1) then ; allocate ( cvec(nfreq)  ) ; cvec = freq + img*enek
                else               ; allocate ( cvec(nfreqr) ) ; cvec = img*(enek-freqr1)
                endif
                call freqintegral_init(wfreqintegral1,freq,nfreq,cvec-1d-4,size(cvec))
                call freqintegral_init(wfreqintegral2,freq,nfreq,cvec+1d-4,size(cvec))
                call freqintegral_init(wfreqintegral, freq,nfreq,cvec,     size(cvec))
                deallocate ( cvec )
                wfreqintegral2 =   wfreqintegral2 + wfreqintegral1
                wfreqintegral1 = ( wfreqintegral2 - wfreqintegral1*2 ) / 2d-4
                wfreqintegral2 = ( wfreqintegral2 - wfreqintegral *2 ) / 1d-8
                ! Frequency convolutions (->hw)
                do j = 1,3
                  do i = 1,3
                    k = min(i,j) + max(i,j)*(max(i,j)-1)/2
                    call spline_init(aspline, freq, ihead(k,:), nfreq)
                    if(i>j) aspline = MCONJG(aspline)
                    do ifreq = 1,size(selfc,2)
                      hw(i,j,ifreq,0) = sum(wfreqintegral (:,:,ifreq)*aspline) / ( 2*pi*img)
                      hw(i,j,ifreq,1) = sum(wfreqintegral1(:,:,ifreq)*aspline) / ( 2*pi)
                      hw(i,j,ifreq,2) = sum(wfreqintegral2(:,:,ifreq)*aspline) / (-2*pi*img)
                    enddo
                  enddo
                enddo
                if(drude_separate) then
                  if(oselfc==1) then ; allocate ( cvec(nfreq)  ) ; cvec = freq + img* enek
                  else               ; allocate ( cvec(nfreqr) ) ; cvec =        img*(enek-freqr1)
                  endif
                  do ifreq = 1,size(selfc,2)
                    hw(:,:,ifreq,0) = hw(:,:,ifreq,0) +     identity(3) * (2*pi*plasma) * drudefunc(cvec(ifreq),1) * img
                    hw(:,:,ifreq,1) = hw(:,:,ifreq,1) +     identity(3) * (2*pi*plasma) * drudefunc(cvec(ifreq),2)
                    hw(:,:,ifreq,2) = hw(:,:,ifreq,2) - 2 * identity(3) * (2*pi*plasma) * drudefunc(cvec(ifreq),3) * img
                  enddo
                  deallocate ( cvec )
                endif
                ! Calculate derivatives of KS energies (->dene, dene2=dene*dene^T, ddene)
# if defined(__GFORTRAN__) && __GNUC__ < 5
#   warning rewritten to avoid gfortran 4.9.1 bug
                dene  = [ ( sum ( [ (moment_diag(id,i),id=deg0,deg1) ] ) , i=1,3 ) ]
                dene  = dene / (deg1-deg0+1)
# else
                dene  = [ ( sum ( [ (moment_diag(id,i),id=deg0,deg1) ] ) , i=1,3 ) ] / (deg1-deg0+1) ! average over degenerate states
# endif
                dene2 = reshape ( [ ((dene(i)*dene(j),i=1,3),j=1,3) ] , [3,3] )
                do k = 1,size(band) ; band1 = band(k) ! ddene only needed if one of band(:) falls inside deg0..deg1
                  if(band1>=deg0.and.band1<=deg1) then
                    ddene = identity(3) + 2 * reshape ( (/ ((
     &                sum   ( moment(:nband(ikptq,ispin),k,i) *
     &                MCONJG (moment(:nband(ikptq,ispin),k,j)) / (ene(band1,ikptq,ispin)-ene(:nband(ikptq,ispin),ikptq,ispin)),
     &                abs(ene(band1,ikptq,ispin)-ene(:nband(ikptq,ispin),ikptq,ispin))>edeg ), i=1,3),j=1,3) /) , [3,3] )
                    exit
                  endif
                enddo
              endif
              ! Calculate k expansion of <E phi_i | phi_n> <phi_n | phi_j E> = O(1) + k*pv + k*pm*k (->pm,pv)
              do iband = deg0,deg1 ! phi_n
                allocate ( cvec(3) )
                do j = 1,size(band)
                  do i = 1,size(band)
                    cvec = 0
                    cmat = 0
                    if( (band(i)<deg0.or.band(i)>deg1) .and. (band(j)<deg0.or.band(j)>deg1) ) then
                      cmat = cmat + reshape ( [((moment(iband,i,k)*MCONJG(moment(iband,j,l)),k=1,3),l=1,3)], [3,3] )
     &                              / ( (ene(iband,ikptq,ispin) - ene(band(i),ikptq,ispin)) *
     &                                  (ene(iband,ikptq,ispin) - ene(band(j),ikptq,ispin)) )
                    endif
                    if(band(i)==iband) then
                      if(band(j)<deg0.or.band(j)>deg1)
     &                cvec = cvec + MCONJG(moment(iband,j,:)) / (ene(iband,ikptq,ispin) - ene(band(j),ikptq,ispin))
                      cmat = cmat + conjg(zm(:,:,j,i))
                    endif
                    if(band(j)==iband) then
                      if(band(i)<deg0.or.band(i)>deg1)
     &                cvec = cvec + moment(iband,i,:) / (ene(iband,ikptq,ispin) - ene(band(i),ikptq,ispin))
                      cmat = cmat + zm(:,:,i,j)
                    endif
                    pm(:,:,i,j) = cmat
                    pv(:,  i,j) = cvec
                  enddo
                enddo
                deallocate ( cvec )
                ! Exchange self-energy (NOSTORE)
                if(evalx.and.enek<0) then
                  iself = iselfx
                  do j = 1,size(band)
                    do i = 1,j
                      iself        = iself + 1
                      selfx(iself) = selfx(iself) - divfac_x/(3*nkpt) * ( pm(1,1,i,j)+pm(2,2,i,j)+pm(3,3,i,j) )
                    enddo
                    if(.not.job1%full) exit
                  enddo
                endif
                ! Correlation self-energy
                if(evalc) then
                  do ifreq = 1,size(selfc,2)
# ifdef CMP_0405adj
                    Error('ddene undefined due to change in selfenergy!')
                    h0     =        hw(:,:,ifreq,0)
                    h2     = matmul(hw(:,:,ifreq,1),ddene) / 2 + matmul(hw(:,:,ifreq,2),dene2) / 2
# else
                    h0     = hw(:,:,ifreq,0)
                    h1     = hw(:,:,ifreq,1)
                    h2     = hw(:,:,ifreq,2)
# endif
                    iself  = iselfc
                    if(oselfc/=1) then
                      enediff = enek - freqr1(ifreq)
                      if(enek*enediff<=0) then
                        p = 1
                        if(enediff==0) p = 2
                        if(enek>0)     p = -p
                        do i = 1,3 ; do j = 1,3
                          cdum = 0 ; cdum1 = 0 ; cdum2 = 0
                          if(nfreqc==0) then
                            if(any(cpade_head(:,i,j)/=0)) then
                              cdum  = pade_func((1d0,0d0)* abs(enediff),      img*freq,cpade_head(:,i,j),nfreq)
                              cdum1 = pade_func((1d0,0d0)*(abs(enediff)-1d-4),img*freq,cpade_head(:,i,j),nfreq)
                              cdum2 = pade_func((1d0,0d0)*(abs(enediff)+1d-4),img*freq,cpade_head(:,i,j),nfreq)
                            endif
                          else
                            call getcoeff(coeff,abs(enediff),     real(freqc),nfreqc) ; cdum  = sum(coeff*iheadc(i,j,:))
                            call getcoeff(coeff,abs(enediff-1d-4),real(freqc),nfreqc) ; cdum1 = sum(coeff*iheadc(i,j,:))
                            call getcoeff(coeff,abs(enediff+1d-4),real(freqc),nfreqc) ; cdum2 = sum(coeff*iheadc(i,j,:))
                          endif
                          cdum2   = cdum1 + cdum2
                          cdum1   = ( cdum2 - cdum1*2 ) / 2d-4
                          cdum2   = ( cdum2 - cdum *2 ) / 1d-8
                          h0(i,j) = h0(i,j) + cdum  / p
# ifndef CMP_0405adj
c                          h1(i,j) = h1(i,j) + cdum1 / p ! currently disabled because a Taylor expansion of W along the real axis is probably unreliable
c                          h2(i,j) = h2(i,j) + cdum2 / p
# endif
                        enddo ; enddo
                      endif
                    endif
                    do j = 1,size(band)
                      do i = 1,size(band)
                        iself = iself + 1
# ifdef CMP_0405adj
                        if(i==j.and.iband==band(i))
     &                  Nacc1_c( selfc,(iself,ifreq), - ( h2(1,1)+h2(2,2)+h2(3,3) ) / (3*vol*nkpt) )
                        Nacc1_c( selfc,(iself,ifreq), - ( sum(h0*pm(:,:,i,j))     ) / (3*vol*nkpt) )
# else
                        cmat = reshape ( [((dene(k)*pv(l,i,j),k=1,3),l=1,3)], [3,3] )
                        if(i==j.and.iband==band(i))
     &                  Nacc1_c( selfc,(iself,ifreq), - ( avg(h1,(1d0,0d0)*ddene) + avg(h2,(1d0,0d0)*dene2) ) / (vol*nkpt*2) )
                        Nacc1_c( selfc,(iself,ifreq), - ( avg(h1,cmat)            + avg(h0,pm(:,:,i,j))     ) / (vol*nkpt)   )
# endif
                        if(.not.job1%full) exit
                      enddo
                      if(.not.job1%full) exit
                    enddo
                  enddo
                endif
              enddo
              endMOD
              deg0 = deg1 + 1
            enddo
            deallocate ( pm,pv,zm,wfreqintegral1,wfreqintegral2 )
          endif
          call cpu_time(time2) ; time_gamma = time_gamma + time2 - time1
        endif

c
c       Define tetrahedron weights
        if(evalc.and.freqint/=3) then          
          if(oselfc==1) then
            if(.not.allocated(wintgrc)) then
              ! copy from wintgr
              minb(1,:) = 1
              maxb(1,:) = bando
              allocate ( wintgrc(nkpti+nkpti2,bando,1,nspin1) )
              wintgrc(:nkpti,:,1,:) = wintgr(:nkpti,:,:)
              if(lkptadd) wintgrc(nkpti+1:,:,1,:) = wintgr(nkpt+1:nkpt+nkpti2,:,:)
            endif
          else if(oselfc/=4.or.oselfc==4.and..not.allocated(wintgrc)) then
            call cpu_time(time1)
            ! determine minimum and maximum bands
            do ifreqr = 1,nfreqr
              do s = 1,nspin1 ; if(oselfc/=4.and.s/=ispin) cycle
                ! (1) minimum
                i = 1
                do while(all(ene(i,:,s)-freqr1(ifreqr)<=efermi).and.i/=maxeband) ! all states below efermi(+freqr1) will get weight 1/nkpt
                  i = i + 1
                enddo
                j = i
                do k = 1,size(ene,2) ! be sure to catch all degenerate states
                  if(j<=nband(k,s)) i = min(i,deg(j,k,s))
                enddo
                minb(ifreqr,s) = i
                ! (2) maximum
                do while(any(ene(i,:,s)-freqr1(ifreqr)<efermi).and.i/=maxeband) ! all states above efermi(+freqr1) will get weight 0
                  i = i + 1
                enddo
                i = i - 1
                if(i>0) then
                  j = i
                  do k = 1,size(ene,2) ! be sure to catch all degenerate states
                    if(j<=nband(k,s)) then
                      l = deg(j,k,s) ; if(l<j) l = deg(l,k,s)
                      i = max(i,l)
                    endif
                  enddo
                endif
                maxb(ifreqr,s) = i
              enddo
            enddo
            ! calculate weights
            i = maxval(maxb-minb)+1
            if(oselfc/=4) then ; allocate ( wintgrc(nkpti+nkpti2,i,nfreqr,ispin:ispin) )
            else               ; allocate ( wintgrc(nkpti+nkpti2,i,nfreqr,nspin1) )
            endif
            if(i>0) then              
              wintgrc = 0
              do s = 1,nspin1 ; if(oselfc/=4.and.s/=ispin) cycle
                do ifreqr = 1,nfreqr
                  call tetrahedron_init(wintgrc(:nkpti,:,ifreqr,s),nkpti,minb(ifreqr,s),maxb(ifreqr,s),
     &                                  efermi+freqr1(ifreqr),s,.false.)
                  if(lkptadd)
     &              call tetrahedron_init(wintgrc(nkpti+1:,:,ifreqr,s),nkpti2,minb(ifreqr,s),maxb(ifreqr,s),
     &                                    efermi+freqr1(ifreqr),s,.true.)
                enddo
              enddo
            endif
            call cpu_time(time2) ; time_tetra = time_tetra + time2 - time1
          endif
        endif

c       IBC: allocate auxiliary selfc_ibc array
        if(oibc/=0.and.evalc) then
          Mpi( Error('IBC&MPI not implemented') )
          allocate ( selfc_ibc(size(band),size(band),nfreq) )
          selfc_ibc = 0
        endif

c       MPI: read cmtq/cpwq
        nkpt_2 = 0
# ifdef LOAD
        if(allocated(cmtq)) then
          cmtq(:,:,:size(band),:) = cmtq(:,:,band,:)
          cpwq(:,  :size(band),:) = cpwq(:,  band,:)
          call reallocate ( cmtq,maxlmindx,ncent,size(band),nspin3 )
          call reallocate ( cpwq,maxgpt,         size(band),nspin3 )
        else
          allocate ( cmtq(maxlmindx,ncent,size(band),nspin3) )
          allocate ( cpwq(maxgpt,         size(band),nspin3) )
          call read_wavef2(band,size(band),ikptq,ispin,cmtq,cpwq)
        endif
        band = [ (i,i=1,size(band)) ]
# endif

c
c       MEM work packages
        if(evalc) then ; bandup = maxval ( [ (nband(kptsum(ikptq,kpt1(k)),ispin),k=1,nkpt1) ] ) ! bandup = maximal number of bands
        else           ; bandup = bando
        endif
        Rbegin
        rdum = ( maxmem - mem ) Mpi(/Nsize) Mpi(*Msize) ! available memory in bytes
        pack = min ( bandup , int ( rdum / (    MBYTES*nbasm(ikpt)*nkpt1*size(band)                     ! for cprod
     &                                       + (MBYTES*maxgpt+16d0*maxlmindx*ncent)*nkpt_2*nspin3 ) ) ) ! for cmt/cpw (only MPI)
        if(pack<1) Error('Insufficient memory; increase MEM'//Mpi(' or number of processes'//)'!')
        maxpack = max(maxpack,ceiling(1d0*bandup/pack-1d-12))
        Rend
        Mpi ( call Mcast(pack) )

c
c       LOAD: allocate cmt/cpw
# ifdef LOAD
        if(pack/=bandup.or.newkptq) then                    ! MEM separation or newkptq
          newkptq = .true.                                  ! force execution of read_wavef below
          nbnd    = ceiling( (1d0*pack-1d-12) Mpi(/Msize) ) ! maximal number of bands per process
          kindx   = 0                                       ! prepare kindx and kpt2
          i       = 0
          do k = 1,nkpt1
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikptq)
            if(storeibz) ikpt2 = kptp(ikpt2)
            if(kindx(ikpt2)==0) then ; i = i + 1 ; kindx(ikpt2) = i ; kpt2(i) = ikpt2 ; endif
            if(storeibz) kindx(kptsum(ikpt1,ikptq)) = kindx(ikpt2)
          enddo
          nkpt_2 = i
          if(associated(cmt)) deallocate(cmt,cpw)
          allocate ( cmt(maxlmindx,ncent,nbnd,nkpt_2,nspin3) )
          allocate ( cpw(maxgpt,         nbnd,nkpt_2,nspin3) )
        endif
# endif

c
c       LOAD&MPI: Split in case of unused processes: Msize > bandup (avoid deadlock)
# ifdef MPI
        lsplit = .false.
        if(Msize>bandup) then
          RInfo('More processes than bands: '//Chr(Msize-bandup)//' cores unused. Use less processes or try MPIBLK or MPIKPT.') 
          lsplit = .true.
          if(1+Mrank>bandup) then ; i = mpi_undefined
          else                    ; i = 0
          endif
          call begin_split(i)
        endif
# endif

c
c       Loop over MEM packets (only one packet if memory sufficient)
        do band1 = 1 Mpi(+Mrank),bandup,pack  ! band1 = lower bound of band package
        band2 = min(band1+pack-1,bandup)      ! band2 = upper bound of band package
        nbnd  = (band2-band1) Mpi(/Msize) + 1 ! nbnd  = number of bands actually treated by current rank (or just band2-band1+1): bands band1,... treated by ranks 0,...,Msize-1,0,...

        if(evalx.and.band1<=bando) allocate ( matx(size(band),size(band))       )
        if(evalc)                  allocate ( mat (size(band),size(band),nfreq) )

# ifdef LOAD
        if(newkptq) call read_wavef0([(i,i=band1,band2)],kpt2(:nkpt_2),ispin,cmt,cpw,.true.)
# endif

c
c       Wave-function products at kpt1
        call cpu_time(time1)
        allocate ( cprod(nbasm(ikpt),nbnd,size(band),nkpt1), cprod1(nbasm(ikpt)) )
        if(evalc.and.oibc/=0) then
          Error('IBC currently disabled')
          allocate ( cprod_ibc(nbasp,maxlmindx,size(band),nkpt1) )
          cprod_ibc = 0
        endif

        call cpu_time(time1)
        call wavefproducts1_mt(cprod,nbasm(ikpt),band,size(band),ikptq,ispin,kpt1,nkpt1, ifLoad(1,band1), ifLoad(nbnd,band2), INC )
        call cpu_time(time2) ; time_mt = time_mt + time2 - time1 ; time1 = time2
        call wavefproducts1_pw(cprod,nbasm(ikpt),band,size(band),ikptq,ispin,kpt1,nkpt1, ifLoad(1,band1), ifLoad(nbnd,band2), INC )
        call cpu_time(time2) ; time_pw = time_pw + time2 - time1 ; time1 = time2

c       Rotate to kptp(ikpt1), and transform to eigenbasis
        do k = 1,nkpt1
          do i = 1,size(band)
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikptq)
            do iband = band1,min(band2,nband(ikpt2,ispin)) MpiC(Msize) ; ibnd = (iband-band1) Mpi(/Msize) + 1
              if(lomit(3)) then ; if(any(omit==iband)) then ; cprod(:,ibnd,i,k) = 0 ; cycle ; endif ; endif
              call mtrafo_r(cprod1,cprod(:,ibnd,i,k),nbasm(ikpt),1,kptp(ikpt1),-symkpt(ikpt1),1,.false.)
              cprod(:,ibnd,i,k) = cprod1
            enddo
          enddo
        enddo
        call cpu_time(time2) ; time_trafo = time_trafo + time2 - time1 ; time1 = time2

c
c       Exchange self-energy (case NOSTORE)
        if(evalx.and.band1<=bando) then
          matx = 0
          do k = 1,nkpt1
            do i = 1,size(band)
              ikpt1 = kpt1(k)
              ikpt2 = kptsum(ikpt1,ikptq)
              do iband = band1,min(bando,band2) MpiC(Msize) ; ibnd = (iband-band1) Mpi(/Msize) + 1
                if(lomit(1)) then ; if(any(omit==iband)) cycle ; endif
                ! matrix-vector product: cprod1 = coulomb * cprod(i)
                cprod1 = matvec ( coulomb0 , cprod(:,ibnd,i,k) )
                do j = 1,i
                  ! scalar products: conjg(cprod(j)) * cprod1
                  cdum = dotprod ( cprod(:,ibnd,j,k) , cprod1 ) * (nkpts(k)*wintgr(kptp(ikpt2),iband,ispin))
                  if(i/=j.and.symkpt(ikpt1)>nsymt) cdum = conjg(cdum)
                  matx(i,j) = matx(i,j) + cdum
                  if(i/=j) then ; matx(j,i) = MCONJG(matx(i,j))
                  else          ; matx(i,i) = dble(matx(i,i))
                  endif
                enddo
              enddo
            enddo
          enddo
          call cpu_time(time2) ; time_exch = time_exch + time2 - time1 ; time1 = time2
          ! Add contribution of equivalent k points
          if(job1%full.and.any(nkpts(:nkpt1)/=1)) then
            matx = matx / nsym1
            allocate ( vec(ndeg**2),cvec(ndeg**2) )
            do j = 1,size(band)/ndeg
              do i = 1,size(band)/ndeg
                l = (i-1)*ndeg
                m = (j-1)*ndeg
                if(trsoff) then ! no time-reversal symmetry
                  vec                         = reshape ( matx(l+1:l+ndeg,m+1:m+ndeg)       , [ ndeg**2     ] )
                  matx(l+1:l+ndeg,m+1:m+ndeg) = reshape ( matmul(vec,irrep_contr(:,:,i,j))  , [ ndeg , ndeg ] )
                else            ! with time-reversal symmetry
                  if(i>j) cycle
                  vec                         = reshape ( matx(l+1:l+ndeg,m+1:m+ndeg)       , [ ndeg**2     ] )
                  cvec                        = reshape ( matx(m+1:m+ndeg,l+1:l+ndeg)       , [ ndeg**2     ] )
                  matx(l+1:l+ndeg,m+1:m+ndeg) = reshape ( matmul( vec,irrep_contr(:,:,i,j)) , [ ndeg , ndeg ] )
     &                                        + reshape ( matmul(cvec,irrep_contt(:,:,i,j)) , [ ndeg , ndeg ] )
                  if(i/=j)
     &            matx(m+1:m+ndeg,l+1:l+ndeg) = reshape ( matmul(cvec,irrep_contr(:,:,j,i)) , [ ndeg , ndeg ] )
     &                                        + reshape ( matmul( vec,irrep_contt(:,:,j,i)) , [ ndeg , ndeg ] )
                endif
              enddo
            enddo
            deallocate ( vec,cvec )
            call cpu_time(time2) ; time_equiv = time_equiv + time2 - time1 ; time1 = time2
          endif
          ! use great orthogonality theorem for diagonal elements
          if(.not.job1%full) matx(1,1) = sum( [ (matx(i,i),i=1,size(band)) ] ) / size(band)
          ! Add to selfx
          iself = iselfx
          do j = 1,size(band)
            do i = 1,j
              iself        = iself + 1
              selfx(iself) = selfx(iself) - matx(i,j)
            enddo
            if(.not.job1%full) exit
          enddo
          deallocate ( matx )
        endif
        deallocate ( cprod1 )
        if(.not.evalc) goto 2 ! skip the rest for HF only

c
c       Transform to eigenfunctions of Coulomb
        call cpu_time(time1)
        do k = 1,nkpt1
          do i = 1,size(band)
            cprod(:n,:,i,k) = macmat ( ctrafo , cprod(:,:,i,k) )
          enddo
        enddo
        call cpu_time(time2) ; time_trafo = time_trafo + time2 - time1 ; time1 = time2

c
c       Loop over the equivalent k points of EIBZ (Green function)
c       (We will sum up the matrices mat(c) for ikpt1=kpt1(ka:kb) because they all lead to symmetry-equiv. ikpt2, which makes the frequ. convolution identical.)
        kb = 0
        do while(kb<nkpt1)
          ka    = kb + 1
          kb    = ka
          kptp2 = kptp(kptsum(kpt1(ka),ikptq))
          do while(kb<nkpt1)
            if(kptp(kptsum(kpt1(kb+1),ikptq))==kptp2) then ; kb = kb + 1
            else                                           ; exit
            endif
          enddo
          ikpt2  = kptp2
          ikpt2_ = ikpt2 ; if(ikpt2>nkpt) ikpt2_ = ikpt2 - nkpt + nkpti

c
c         Loop over bands (Green function)

          iband = band1
          do while(iband<=min(band2,nband(ikpt2,ispin)))
            enek = ene(iband,ikpt2,ispin) - efermi
            deg1 = deg(iband,ikpt2,ispin) ; if(deg1<iband) deg1 = deg(deg1,ikpt2,ispin)
            deg0 = deg(deg1, ikpt2,ispin)
            deg1 = min(deg1,band2)
            deg0 = max(deg0,band1 Mpi(-Mrank) )
            ideg = ifMpi( iband , deg1 )

c           Progress bar
            if(evalc.and.job1%full Mpi(.and.Mrank==0.and.Mcolor<=1) ) then
              j = (ideg-iband+1) * size(band) * (kb-ka+1)
              do i = 49*ipro/npro+1,49*(ipro+j)/npro ; write(6,'(''.'''NoA) ; enddo
              ipro = ipro + j
            endif

c           FREQINT NONLIN: Get tetrahedron weight function
            if(freqint==3) call tetrahedron_nonlin(nfrq,frq,pfrq,enek,iband,ikpt2,ispin,.false.)

c           Determine for which real frequencies {w} we must calculate matrices matc(w) = <nn''|W(w)|n''n'>, n''=iband  (->ifreqc1/2)
            ifreqc1 = nfreqc+1
            ifreqc2 = 0
            if(nfreqc/=0) then
              if(freqint==3) then
                f1 = frq(1)    + enek + 1d-12
                f2 = frq(nfrq) + enek - 1d-12
                do ifreqr = 1,nfreqr
                  emin = min(0d0,freqr1(ifreqr))
                  emax = max(0d0,freqr1(ifreqr))
                  sgn  = nint( -sign(1d0,freqr1(ifreqr)) )
                  if(emax-emin>1d-8.and.f1<emax.and.f2>emin) then
                    do ifreqc = 1,nfreqc-1
c# define ifort_vector_bug
# ifdef ifort_vector_bug
                      g1 = real(freqc(ifreqc))   ! Strange ifort12 bug: vectorization (-O2) makes g2 undefined (except at loop unrolling)
                      g2 = real(freqc(ifreqc+1)) ! Even stranger: works with different variable names, e.g., rdum, rdum1 (see below).
                      ggg1(ifreqc) = g1
                      ggg2(ifreqc) = g2
# else
                      rdum  = sgn*real(freqc(ifreqc))   + freqr1(ifreqr)
                      rdum1 = sgn*real(freqc(ifreqc+1)) + freqr1(ifreqr)
                      if(min(rdum,rdum1)<min(emax,f2).and.max(rdum,rdum1)>max(emin,f1)) then
                        ifreqc1 = min(ifreqc1,ifreqc)
                        ifreqc2 = max(ifreqc2,ifreqc+1)
                      endif
# endif
                    enddo
                  endif
# ifdef ifort_vector_bug
                  write(201,'(I3,2F20.10)') (ifreqc,ggg1(ifreqc),ggg2(ifreqc),ifreqc=1,nfreqc)
# endif
                enddo
              else
                do ifreqr = 1,nfreqr
                  enediff = enek - freqr1(ifreqr)
                  weight  = 0
                  if     (iband<=bando)                    weight = -wintgr(ikpt2,iband,ispin)
                  if     (iband< minb(ifreqr,ispin)) then; weight = weight + 1d0/nkpt
                  else if(iband<=maxb(ifreqr,ispin)) then; weight = weight + wintgrc(ikpt2_,iband-minb(ifreqr,ispin)+1,ifreqr,ispin)
                  endif
                  if(abs(weight)>1d-10) then
                    call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
                    do ifreqc = 1,nfreqc
                      if(coeff(ifreqc)/=0) then
                        ifreqc1 = min(ifreqc1,ifreqc)
                        ifreqc2 = max(ifreqc2,ifreqc)
                      endif
                    enddo
                  endif
                enddo
              endif
# ifdef ifort_vector_bug              
              write(*,*) ifreqc1,ifreqc2  
              stop
# endif
              if(ifreqc1<=ifreqc2) allocate ( matc(size(band),size(band),ifreqc1:ifreqc2) )
            endif

c           Vector-matrix-vector products : cprod(i)*screen*cprod(j)
            call cpu_time(time1)
            ! (a) imaginary frequencies (->mat)
            Mpi( if(.not.mpisym.or.iband-Msize<deg0) ) mat = 0
            allocate ( cprod1(n) )
            do ifreq = 1,nfreq
              do k = ka,kb ; ikpt1 = kpt1(k)
                do id = iband,ideg ; ibnd = (id-band1) Mpi(/Msize) + 1
                  do i = 1,size(band)
                    ! matrix-vector product: cprod1 = screen * cprod(i)
                    if(allocated(cblock)) then
                      call matvec_sym( cprod1 , screen(:,ifreq) , cprod(:n,ibnd,i,k) )
                    else
                      cprod1 = matvec ( screen(:,ifreq) , cprod(:n,ibnd,i,k) )
                    endif
                    do j = 1,i
                      ! scalar products: conjg(cprod(j)) * cprod1
                      cdum = dotprod ( cprod(:n,ibnd,j,k) , cprod1 ) * nkpts(k)
                      if(i/=j.and.symkpt(ikpt1)>nsymt) cdum = conjg(cdum)
                      mat(i,j,ifreq) = mat(i,j,ifreq) + cdum
                      if(i/=j) then  ; mat(j,i,ifreq) = MCONJG(mat(i,j,ifreq))
                      else           ; mat(i,i,ifreq) = dble(mat(i,i,ifreq))
                      endif
                    enddo
                  enddo
                enddo
              enddo
            enddo
            deallocate ( cprod1 )
            ! (b) real frequencies (->matc)
            if(nfreqc/=0) then
              allocate ( cprod1w(n) )
              do ifreqc = ifreqc1,ifreqc2
                matc(:,:,ifreqc) = 0
                do k = ka,kb ; ikpt1 = kpt1(k)
                  do id = iband,ideg ; ibnd = (id-band1) Mpi(/Msize) + 1
                    do i = 1,size(band)
                      ! matrix-vector product: cprod1w = screenc * cprod(i)
                      if(allocated(cblock)) then
                        call matvec_gen( cprod1w , screenc(:,:,ifreqc) , cprod(:n,ibnd,i,k) )
                      else
                        cprod1w = matvec ( screenc(:,:,ifreqc) , cprod(:n,ibnd,i,k) )
                      endif
# ifdef INV
                      do j = 1,i
# else
                      do j = 1,size(band)
# endif
                        ! scalar products: conjg(cprod(j)) * cprod1
                        cdum = dotprod ( cprod(:n,ibnd,j,k) , cprod1w ) * nkpts(k)
                        if(i/=j.and.symkpt(ikpt1)>nsymt) then
                          matc(j,i,ifreqc) = matc(j,i,ifreqc) + cdum
                        else
                          matc(i,j,ifreqc) = matc(i,j,ifreqc) + cdum
# ifdef INV
                          matc(j,i,ifreqc) = matc(i,j,ifreqc)
# endif
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              deallocate ( cprod1w )
            endif
            call cpu_time(time2) ; time_mat = time_mat + time2 - time1

# ifdef MPI
c           MPISYM: In case of CONTOUR and nfreqc=0, symmetrize over degenerate subspace to avoid symmetry breaking.
            if(mpisym.and.deg0/=deg1) then
              if(iband==deg1) then ! Last degenerate receives contribution from others
                allocate ( math(size(band),size(band),nfreq) )
                do ideg = deg0,deg1-1
                  rank = modulo(Mrank+ideg-deg1,Msize)
                  if(rank/=Mrank.and.ideg+Msize>deg1) then
                    call mpi_recv(math,size(math),MMCOMPLEX,rank,ideg,Mcomm,Mstat,Merr)
                    mat = mat + math
                  endif
                enddo
                deallocate ( math )
              else                 ! Others send their contribution and cycle
                rank  = modulo(Mrank+deg1-iband,Msize)
                if(rank/=Mrank.and.iband+Msize>deg1) call mpi_send(mat,size(mat),MMCOMPLEX,rank,iband,Mcomm,Merr)
                iband = iband + Msize
                cycle
              endif
            endif
# endif

c           Add contribution of equivalent k points
            if(job1%full.and.any(nkpts(ka:kb)/=1)) then
              ! use contracted irreps for full matrix
              call cpu_time(time1)
              mat = mat / nsym1
              if(ifreqc1<=ifreqc2) matc = matc / nsym1
              allocate ( vec(ndeg**2),cvec(ndeg**2) )
              do ifreq = 1,nfreq + nfreqc
                ifreqc = ifreq - nfreq ; if(ifreqc>0.and.(ifreqc<ifreqc1.or.ifreqc>ifreqc2)) cycle
                do j = 1,size(band)/ndeg
                  do i = 1,size(band)/ndeg
                    l = (i-1)*ndeg
                    m = (j-1)*ndeg
                    if(trsoff) then ! no time-reversal symmetry
                      if(ifreq<=nfreq) then
                        vec                                = reshape ( mat (l+1:l+ndeg,m+1:m+ndeg,ifreq) , [ ndeg**2     ] )
                        mat(l+1:l+ndeg,m+1:m+ndeg,ifreq)   = reshape ( matmul(vec,irrep_contr(:,:,i,j))  , [ ndeg , ndeg ] )
                      else
                        cvec                               = reshape ( matc(l+1:l+ndeg,m+1:m+ndeg,ifreqc), [ ndeg**2     ] )
                        matc(l+1:l+ndeg,m+1:m+ndeg,ifreqc) = reshape ( matmul(cvec,irrep_contr(:,:,i,j)) , [ ndeg , ndeg ] )
                      endif
                    else            ! with time-reversal symmetry
                      if(i>j) cycle
                      if(ifreq<=nfreq) then
                        vec                                = reshape ( mat (l+1:l+ndeg,m+1:m+ndeg,ifreq) , [ ndeg**2     ] )
                        cvec                               = reshape ( mat (m+1:m+ndeg,l+1:l+ndeg,ifreq) , [ ndeg**2     ] )
                        mat(l+1:l+ndeg,m+1:m+ndeg,ifreq)   = reshape ( matmul( vec,irrep_contr(:,:,i,j)) , [ ndeg , ndeg ] )
     &                                                     + reshape ( matmul(cvec,irrep_contt(:,:,i,j)) , [ ndeg , ndeg ] )
                        if(i/=j)
     &                  mat(m+1:m+ndeg,l+1:l+ndeg,ifreq)   = reshape ( matmul(cvec,irrep_contr(:,:,j,i)) , [ ndeg , ndeg ] )
     &                                                     + reshape ( matmul( vec,irrep_contt(:,:,j,i)) , [ ndeg , ndeg ] )
                      else
                        vec                                = reshape ( matc(l+1:l+ndeg,m+1:m+ndeg,ifreqc), [ ndeg**2     ] )
                        cvec                               = reshape ( matc(m+1:m+ndeg,l+1:l+ndeg,ifreqc), [ ndeg**2     ] )
                        matc(l+1:l+ndeg,m+1:m+ndeg,ifreqc) = reshape ( matmul( vec,irrep_contr(:,:,i,j)) , [ ndeg , ndeg ] )
     &                                                     + reshape ( matmul(cvec,irrep_contt(:,:,i,j)) , [ ndeg , ndeg ] )
                        if(i/=j)
     &                  matc(m+1:m+ndeg,l+1:l+ndeg,ifreqc) = reshape ( matmul(cvec,irrep_contr(:,:,j,i)) , [ ndeg , ndeg ] )
     &                                                     + reshape ( matmul( vec,irrep_contt(:,:,j,i)) , [ ndeg , ndeg ] )
                      endif
                    endif
                  enddo
                enddo
              enddo
              deallocate ( vec,cvec )
              call cpu_time(time2) ; time_equiv = time_equiv + time2 - time1
            endif
            if(.not.job1%full.and.size(band)>1) then
              ! use great orthogonality theorem for diagonal elements
              do ifreq = 1,nfreq
                mat(1,1,ifreq)   = sum( [ (mat(i,i,ifreq),i=1,size(band)) ] )   / size(band)
              enddo
              do ifreqc = ifreqc1,ifreqc2
                matc(1,1,ifreqc) = sum( [ (matc(i,i,ifreqc),i=1,size(band)) ] ) / size(band)
              enddo
            endif

c
c           Multiply Green function and perform frequency convolution (iw integral and residues)

            call cpu_time(time1)

c           Initialize frequency convolutions
            if(freqint<=1) then
              if(oselfc==1) then ; allocate ( cvec(nfreq)  ) ; cvec = freq + img*enek
              else               ; allocate ( cvec(nfreqr) ) ; cvec = img*(enek-freqr1)
              endif
              call freqintegral_init(wfreqintegral,freq,nfreq,cvec,size(cvec))
              deallocate ( cvec )
            endif

c           Loop over self-energy matrix elements i,j
            call cpu_time(time1)
            do j = 1,size(band)
              do i = 1,j!size(band)
                if(maxval(abs(mat(i,j,1:)))>1d-10) then
                  if(freqint>=2.or.nfreqc==0) call pade_init(cpade,img*freq,(1d0,0d0)*mat(i,j,1:),nfreq,-1)
                  if(freqint>=2) then
                    call pade_poles(pole,resid,npole1,img*freq,(1d0,0d0)*mat(i,j,1:),cpade,nfreq,smooth(2),-1,.false.)
                  else
                    call spline_init(aspline,freq,mat(i,j,1:),nfreq)
                  endif

c                 (A) Frequency convolutions
                  
                  if(freqint<=2) then ! FREQINT SPLINE/PADE
                    do ifreq = 1,size(selfc,2)
                      if(oselfc==1) then ; ifreqr = 1     ; rdum = freq(ifreq) ! the latter only for Pade integration
                      else               ; ifreqr = ifreq ; rdum = 0           !              - " -
                      endif
                      if      (iband<minb(ifreqr,ispin)) then ; weight = 1d0/nkpt
                      else if (iband>maxb(ifreqr,ispin)) then ; weight = 0
                      else                                    ; weight = wintgrc(ikpt2_,iband-minb(ifreqr,ispin)+1,ifreqr,ispin)
                      endif
                      enediff = enek - freqr1(ifreqr)
                      ! occ. states
                      if(weight>1d-10) then
                        if(freqint<=1) then
                          if(enediff>0) then ; wfreq = conjg(wfreqintegral(:,:,ifreq)) ! protrudes above Fermi energy
                          else               ; wfreq =       wfreqintegral(:,:,ifreq)  ! normal case
                          endif
                          cdum = sum(wfreq*aspline) / (2*pi*img)
                        else
                          cdum = freqintegral_poles(pole,resid,npole1,1,rdum-img*abs(enediff),0)
                        endif
                        if(rdum==0.and.enediff==0) cdum = cdum + mat(i,j,1)/2 ! add delta function (due to negative infinitesimal imaginary part)
                        iself = iselfc + i + (j-1)*size(band)
c                        write(*,*) '2'
                        Nacc1_c( selfc,(iself,ifreq), - cdum * weight )
                        if(i/=j) then
# ifndef INV
                          if(freqint<=1)then; cdum = sum(wfreq*conjg(aspline)) / (2*pi*img)
                          else              ; cdum = freqintegral_poles(-conjg(pole),-conjg(resid),npole1,1,rdum-img*abs(enediff),0)
                          endif
# endif
                          iself = iselfc + j + (i-1)*size(band)
                          Nacc1_c( selfc,(iself,ifreq), - cdum * weight )
                        endif
                      endif
                      ! unocc. states
                      weight = 1d0/nkpt - weight
                      if(weight>1d-10) then
                        if(freqint<=1) then
                          if(enediff<0) then ; wfreq = conjg(wfreqintegral(:,:,ifreq)) ! protrudes below Fermi energy
                          else               ; wfreq =       wfreqintegral(:,:,ifreq)  ! normal case
                          endif
                          cdum = sum(wfreq*aspline) / (2*pi*img)
                        else
                          cdum = freqintegral_poles(pole,resid,npole1,1,rdum+img*abs(enediff),0)
                        endif
                        if(rdum==0.and.enediff==0) cdum = cdum - mat(i,j,1)/2 ! add delta function (due to positive infinitesimal imaginary part)
                        iself = iselfc + i + (j-1)*size(band)
c                        write(*,*) '3',iself,ifreq,-cdum*weight,iband
                        Nacc1_c( selfc,(iself,ifreq), - cdum * weight )
                        if(i/=j) then
# ifndef INV
                          if(freqint<=1)then; cdum = sum(wfreq*conjg(aspline)) / (2*pi*img)
                          else              ; cdum = freqintegral_poles(-conjg(pole),-conjg(resid),npole1,1,rdum+img*abs(enediff),0)
                          endif
# endif
                          iself = iselfc + j + (i-1)*size(band)
                          Nacc1_c( selfc,(iself,ifreq), - cdum * weight )
                        endif
                      endif
                    enddo
                    
                  else ! FREQINT NONLIN

                    if(npole1>0) then
                      do ifreq = 1,size(selfc,2)
                        if(oselfc==1) then ; cdum1 = freq(ifreq) * img
                        else               ; cdum1 = freqr1(ifreq)
                        endif
                        cdum  = freqintegral_poles_nonlin(pole,resid,npole1,1,cdum1,0,enek,nfrq,frq,pfrq)
                        iself = iselfc + i + (j-1)*size(band)
                        Nacc1_c( selfc,(iself,ifreq) , -cdum )
                        if(i/=j) then
                          cdum  = freqintegral_poles_nonlin(-conjg(pole),-conjg(resid),npole1,1,cdum1,0,enek,nfrq,frq,pfrq)
                          iself = iselfc + j + (i-1)*size(band)
                          Nacc1_c( selfc,(iself,ifreq) , -cdum )
                        endif
                      enddo
                    endif
                    
                  endif
                  
                endif

c               (B) Residues
                if(oselfc/=1) then
                  if(nfreqc==0.and.maxval(abs(mat(i,j,1:)))<=1d-10.or.nfreqc>0.and..not.allocated(matc)) cycle
                  
                  if(freqint<=2) then ! FREQINT SPLINE/PADE

                    do ifreqr = 1,nfreqr
                      weight = 0
                      if     (iband<=bando)                     weight = -wintgr(ikpt2,iband,ispin)
                      if     (iband< minb(ifreqr,ispin)) then ; weight = weight + 1d0/nkpt
                      else if(iband<=maxb(ifreqr,ispin)) then ; weight = weight +
     &                                                                   wintgrc(ikpt2_,iband-minb(ifreqr,ispin)+1,ifreqr,ispin)
                      endif
                      if(abs(weight)>1d-10) then
                        enediff = enek - freqr1(ifreqr)
                        if(nfreqc==0) then
                          cdum = pade_func((1d0,0d0)*abs(enediff),img*freq,cpade,nfreq)
                        else
                          call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
                          cdum = sum(coeff(ifreqc1:ifreqc2)*matc(i,j,:))
                        endif
                        iself = iselfc + i + (j-1)*size(band)
                        Nacc1_c( selfc,(iself,ifreqr), weight * cdum )
                        if(i/=j) then
# ifndef INV
                          if(nfreqc==0) then
                            cdum = pade_func(-(1d0,0d0)*abs(enediff),-img*freq,conjg(cpade),nfreq)
                          else
                            call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
                            cdum = sum(coeff(ifreqc1:ifreqc2)*matc(j,i,:))
                          endif
# endif
                          iself = iselfc + j + (i-1)*size(band)
                          Nacc1_c( selfc,(iself,ifreqr), weight * cdum )
                        endif
                      endif
                    enddo

                  else ! FREQINT NONLIN

                    allocate ( cvec(nfreqc) ) ; cvec = 0
                    f1 = frq(1)    + enek + 1d-10
                    f2 = frq(nfrq) + enek - 1d-10
                    do ifreqr = 1,nfreqr
                      emin = min(0d0,freqr1(ifreqr))
                      emax = max(0d0,freqr1(ifreqr))
                      sgn  = nint( -sign(1d0,freqr1(ifreqr)) )
                      if(emax-emin>1d-8.and.f1<emax.and.f2>emin) then
                        if(nfreqc>0) then
                          rdum  = sgn*real(freqc(ifreqc1)) + freqr1(ifreqr)
                          rdum1 = sgn*real(freqc(ifreqc2)) + freqr1(ifreqr)
                          if(max(emin,f1)<min(rdum,rdum1)) then
                            write(*,*) 'maxof',emin,f1
                            write(*,*) min(rdum,rdum1),sgn,'lower',ifreqc1,ifreqc2
                            write(*,*) 'test',sgn*real(freqc(ifreqc2)) + freqr1(ifreqr)
                            read(*,*)
                            Warn('Lower bound error (freqc).')
                          endif
                          if(min(emax,f2)>max(rdum,rdum1)) then
                            write(*,*) min(emax,f2),max(rdum,rdum1),sgn,'upper'
                            read(*,*)
                            Warn('Upper bound error (freqc).')
                          endif
                          cvec(ifreqc1:ifreqc2) = matc(i,j,:)
                          cdum = residues(emin,emax,enek,nfrq,frq,pfrq,nfreqc,real(freqc)+sgn*freqr1(ifreqr),cvec,sgn==-1)
                        else
                          cdum = residues_pade(emin,emax,enek,nfrq,frq,pfrq,npole1,sgn*pole+freqr1(ifreqr),sgn*resid)
                        endif
                        if(emax==0d0) cdum = -cdum
                        iself = iselfc + i + (j-1)*size(band)
                        Nacc1_c( selfc,(iself,ifreqr) , cdum )
                        if(i/=j) then
                          if(nfreqc>0) then
                            cvec(ifreqc1:ifreqc2) = matc(j,i,:)
                            cdum = residues(emin,emax,enek,nfrq,frq,pfrq,nfreqc,real(freqc)+sgn*freqr1(ifreqr),cvec,sgn==-1)
                          else
                            cdum = residues_pade(emin,emax,enek,nfrq,frq,pfrq,npole1,sgn*(-conjg(pole))+freqr1(ifreqr),
     &                                                                                                  sgn*(-conjg(resid)))
                          endif
                          if(emax==0d0) cdum = -cdum
                          iself = iselfc + j + (i-1)*size(band)
                          Nacc1_c( selfc,(iself,ifreqr) , cdum )
                        endif
                      endif
                    enddo
                    deallocate(cvec)

                  endif
                  
                endif

              enddo
              if(.not.job1%full) exit
            enddo
            call cpu_time(time2) ; time_freq = time_freq + time2 - time1

            if(allocated(matc)) deallocate ( matc )

# ifdef MPI
            iband = iband + Msize ! MPI:  strided loop over states (to share the work load evenly among the processes)
# else
            iband = ideg + 1      ! else: loop over degenerate subspaces
# endif

          enddo ! iband loop

        enddo ! ka..kb loop

 2      deallocate ( cprod )
        if(allocated(mat)) deallocate ( mat )

        enddo ! band1..band2 loop

        Mpi( if(lsplit) call end_split )

# if 0
        if(job1%full) then
          j = 0
          do i = 1,size(band)
            j = j + 1
            write(*,'(2F30.20)') selfc(iselfc+j,:min(7,size(selfc,2)))
            j = j + size(band)
          enddo
        else
          write(*,'(2F30.20)') selfc(iselfc+1,:min(7,size(selfc,2)))
        endif
        read(*,*)
# endif
c        Error(' ')

c
c       IBC: Add equivalent k points and add to ibc_selfc
        if(oibc/=0.and.evalc) then
          write(*,*) 'Add equivalent k points for IBC.'
          if(job1%full.and.any(nkpts(:nkpt1)/=1)) then
            Error('Not implemented job1%full and IBC.')
          endif
          if(.not.job1%full) then
            ! use great orthogonality theorem for diagonal elements
            do ifreq = 1,nfreq
              ibc_selfc(iselfc+1,ifreq) = sum( [ (selfc_ibc(i,i,ifreq),i=1,size(band)) ] ) / size(band)
            enddo
            write(*,'(3F20.10)') (freq(ifreq),-ibc_selfc(iselfc+1,ifreq),ifreq=1,nfreq)
            write(800,'(3F20.10)') (freq(ifreq),-ibc_selfc(iselfc+1,ifreq),ifreq=1,nfreq)
          endif
          deallocate ( selfc_ibc )
        endif

        if(allocated(wintgrc)) then
          if(oselfc==2.or.oselfc==3) deallocate ( wintgrc )
        endif
        if(allocated(irrep_contr)) deallocate ( irrep_contr )
        if(allocated(irrep_contt)) deallocate ( irrep_contt )

        Load(deallocate(cmtq,cpwq))

 123    deallocate ( band )
        if(allocated(moment)) deallocate ( moment )

      enddo ! loop blocks

# ifdef MPI
      call cpu_time(time1)
      if(Mcolor/=0) call end_split
      if(evalc) then
        Nfence(selfc)
        Ocall Msum(selfc,0,comm=Ocomm) ; MnoR(ifO selfc = 0 )
        Nfence(selfc)
      endif
      if(evalx) then
        call Msum(selfx,0) ; MnoR( selfx = 0 )
      endif
      call cpu_time(time2) ; time_idle = time2 - time1
      call mpi_reduce(time_idle,time_maxidle,1,mpi_real,mpi_max,0,Mcomm,Merr)
# endif

      if(allocated(moment_diag))    deallocate ( moment_diag )
      if(allocated(wfreqintegral))  deallocate ( wfreqintegral )
      if(allocated(cpade_head))     deallocate ( cpade_head )
      if(allocated(wintgrc))        deallocate ( wintgrc )
      Load( if(associated(cmt))     deallocate ( cmt,cpw ) )

      if(any(job1%type==[J_SX,J_COSX])) nfreq = size(freq) ! nfreq was temporarily set to 0

      if(allocated(wblock)) deallocate ( wblock )

      Rbegin
      call cpu_time(time2)
      time_tot = time2 - time_tot
      if(evalc.and.job1%full) then ; write(6,*)
      else                         ; write(6,'(A)') 'done'
      endif
      write(6,'(31X,A)')   'Timings'
      write(6,'(A,F13.5)') '  MT wave-func. products:',time_mt
      write(6,'(A,F13.5)') '  PW wave-func. products:',time_pw
      if(evalc)
     &write(6,'(A,F13.5)') '  Matrix-vector products:',time_mat
      write(6,'(A,F13.5)') '  Transformations:       ',time_trafo
      if(evalx)
     &write(6,'(A,F13.5)') '  Exchange self-energy:  ',time_exch
      if(job1%full)
     &write(6,'(A,F13.5)') '  Equivalent k points:   ',time_equiv
      if(evalc)
     &write(6,'(A,F13.5)') '  Frequency integration: ',time_freq
      if(oselfc/=0)
     &write(6,'(A,F13.5)') '  Tetrahedron weights:   ',time_tetra
      if(ikpt==1)
     &write(6,'(A,F13.5)') '  Gamma-point correction:',time_gamma
      if(evalc.and.oibc/=0)
     &write(6,'(A,F13.5)') '  IBC double-counting:   ',time_ibc
# ifdef MPI
      write(6,'(A,2F13.5)')'  MPI idle time:         ',time_idle,time_maxidle
# endif
      write(6,'(A,F13.5)') '  Other:                 ',time_tot-(time_mt+time_pw+time_mat+time_trafo+time_exch+time_equiv+time_freq+
     &                                                           time_tetra+time_gamma+time_ibc Mpi(+time_idle))
      write(6,'(A)')       '  ------------------------------------'
      write(6,'(A,F13.5)') '  Total:                 ',time_tot
      if(maxpack>1) write(6,'(/A,I15)') '  MEM separation:',maxpack
      if(job1%type==J_RPA) write(6,*)
      Rend

c      read(*,*)

      contains

c -----------------------------

      function drudefunc(c,n)
      implicit none
      complex_dp             :: drudefunc
      complex_dp, intent(in) :: c
      integer,    intent(in) :: n
      real_dp                :: c2
      c2 = imag(c)
      if (plasma==0) then ; drudefunc =   0
      else if(c2>0)  then ; drudefunc =   1/( img*plasma+c)**n
      else if(c2<0)  then ; drudefunc =   1/(-img*plasma+c)**n
      else                ; drudefunc = ( 1/( img*plasma+c)**n + 1/(-img*plasma+c)**n ) / 2
      endif
      end function drudefunc

c -----------------------------

      subroutine block_diag_pack(mat)
      implicit none
      MCOMPLEX_dp, intent(inout) :: mat(n*(n+1)/2)
      MCOMPLEX_dp                :: mat1(n,n)
      integer                    :: i,n0,n1,m0,m1
      call p_unpackmat(mat1,mat)
      do i = 1,n ! following two loops replace mat1 = mat1(pnt,pnt), which is stack-intensive.
        mat1(:,i) = mat1(pnt,i)
      enddo
      do i = 1,n
        mat1(i,:) = mat1(i,pnt)
      enddo      
      n0 = 1
      m0 = 1
      do i = 1,nwblock
        n1 = n0 + wblock(i)                 - 1
        m1 = m0 + wblock(i)*(wblock(i)+1)/2 - 1 ; call p_packmat(mat(m0:m1),mat1(n0:n1,n0:n1))
        n0 = n1 + 1
        m0 = m1 + 1
      enddo
      end subroutine block_diag_pack

      subroutine matvec_sym(res,mat,vec)
      implicit none
      MCOMPLEX_dp, intent(out) :: res(n)
      MCOMPLEX_dp, intent(in)  :: mat(n*(n+1)/2),vec(n)
      integer                  :: i,n0,n1,m0,m1
      res = 0
      n0  = 1
      m0  = 1
      do i = 1,nwblock
        n1         = n0 + wblock(i)                 - 1
        m1         = m0 + wblock(i)*(wblock(i)+1)/2 - 1
        res(n0:n1) = res(n0:n1) + matvec ( mat(m0:m1) , vec(n0:n1) )
        n0         = n1 + 1
        m0         = m1 + 1
      enddo
      end subroutine matvec_sym

      subroutine matvec_gen(res,mat,vec)
      implicit none
      complex_dp,  intent(out) :: res(n)
      complex_dp,  intent(in)  :: mat(n,n)
      MCOMPLEX_dp, intent(in)  :: vec(n)
      integer                  :: i,n0,n1
      res = 0
      n0  = 1
      do i = 1,nwblock
        n1         = n0 + wblock(i) - 1
        res(n0:n1) = res(n0:n1) + matvec ( mat(n0:n1,n0:n1) , vec(n0:n1) )
        n0         = n1 + 1
      enddo
      end subroutine matvec_gen

c ------------------------------

# ifdef MPI
c     Distribute job over processes. The current communicator (Mcomm) is split into subgroups s of P(s) processors.
c     ( sum(s) P(s) is the total number of processors. )
c
c     Assume computation time for each block b : T(b) = N(b) * [ M + K(b)*N'/P(s) ]
c     N(b) - number of bands in block b
c     K(b) - number of Green-function k-points for block b (nkpt1)
c     N'   - number of Green-function bands
c     M    - relative cost of N' independent term (Mover, keyword MPIBLK)
c
c     We determine
c     - the number of subgroups
c     - the number of processors of each subgroup P(s)
c     - the subgroup index S(b) for the block b
c     so as to minimize
c     R(S,P) = max(s) [ SUM[b,S(b)=s] T(b) ]
c
c     Output
c     - Mcolor  : Subgroup index of current rank
c     - Msub(i) : Subgroup index of ith block
c
      subroutine Mdistribute
      implicit none
      integer :: s,proc(Msize),sub(size(block,2)),i,j,iblock,nn,nblock,nk(size(block,2)),nb(size(block,2))
      integer :: pnt(size(block,2))
      real_dp :: time(size(block,2)),stime(Msize),mintime,maxtime
      if(Mover==0) then ! do not distribute jobs
        Mcolor = 0
        Msub   = 0
        return
      endif
      nblock  = size(block,2)
      mintime = huge(0d0)
      ! information from blocks
      do iblock = 1,nblock
        ikptq      = job1%kpt(block(1,iblock)) ; call getkpt1(kpt1,nk(iblock),nkpts,sym1,nsym1,ikptq,ikpt,.false.) ! -> nk(iblock)
        nb(iblock) = count(block(:,iblock)/=0)
      enddo
      ! number of Green-function bands
      if(any(job1%type==[J_GW,J_RPA])) then ; nn = maxband
      else                                  ; nn = bando
      endif
      ! loop over subgroups
      do s = 1,Msize
c        write(*,*) 's=',s
        ! distribute Msize processors over subgroups -> proc
        proc = Msize / s
        if(sum(proc(:s))/=Msize) then
          i        = Msize - Msize/s*s
          proc(:i) = proc(:i) + 1
        endif
        if(sum(proc(:s))/=Msize) Bug('Count error.')
c        write(*,*) 'proc=',proc(:s)
        ! computation time for each block -> time(:)
        time = nb * ( Mover + 1d0 * nk * nn / Msize * s )
        call rorderp(pnt,time,nblock)
        ! distribute blocks over subgroups
        stime = -proc*1d-15
        do j = nblock,1,-1
          iblock      = pnt(j)
          i           = minloc(stime(:s),1)
          stime(i)    = stime(i) + nb(iblock) * ( Mover + 1d0 * nk(iblock) * nn / proc(i) )
          sub(iblock) = i
c          ifR write(*,'(I5'NoA) iblock,nint(nb(iblock) * ( Mover + 1d0 * nk(iblock) * nn / proc(i) ))
c          ifR write(*,'(F10.3'NoA) stime(:s)
c          ifR write(*,*)
        enddo
        maxtime = maxval(stime(:s))
c        ifR write(*,*) 'TIME',maxtime
c        ifR write(*,*) 'proc',proc(:s)
c        ifR write(*,*) 'sub',sub
c        ifR write(*,*) 'time',nb * ( Mover + 1d0 * nk * nn / Msize * s )
c        ifR write(*,*)
c        ifR read(*,*)
        if(maxtime<mintime) then
c        if(s==3) then
          mintime = maxtime
          Msub    = sub
          Mcolor  = 0 ; do while(sum(proc(:Mcolor))<=Mrank) ; Mcolor = Mcolor + 1 ; enddo
        endif
      enddo
c      write(*,*) Mrank,Mcolor
c      call mpi_barrier(mpi_comm_world,i)
c      Error('done')
      end subroutine Mdistribute
# endif

c ------------------------------

c Returns  1/(4*pi) * INT (kAk) (kBk) d²k = [ tr(a)*tr(b) + 2*tr(ab) ] / 15
c - Simplifies to tr(ab) if a (or b) is a multiple of the identity matrix.
c - At least one of the matrices must be symmetric.
      function avg(a,b)
      implicit none
      complex_dp             :: avg
      complex_dp, intent(in) :: a(3,3),b(3,3)
      if(maxval(abs(a-transpose(a)))>1d-8.and.maxval(abs(b-transpose(b)))>1d-8) Error('both matrices asymmetric')
      avg = ( (a(1,1)+a(2,2)+a(3,3)) * (b(1,1)+b(2,2)+b(3,3)) + 2*sum(a*b) ) / 15
      end function avg

c ------------------------------

# if 0
c     Converts data from double precision to single precision (to regularize Pade extrapolation). Currently not used.
      subroutine regularize(val,n)
      implicit none
      integer,     intent(in)    :: n
      MCOMPLEX_dp, intent(inout) :: val(n)
      integer                    :: i
      do i = 1,n
        val(i) = ifInv( real , cmplx ) (val(i))
      enddo
      end subroutine regularize
# endif

c ------------------------------

      end
