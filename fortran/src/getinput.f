c     Get parameters from input file (spex.inp) and data from DFT run.
c
c
c     Program sequence:
c
c       Command line options
c       read_param
c       read_symmetry
c       checkkeys
c       kpt definition (BZ, KPT...)
c       Mcast
c
c       getinput.inc
c
c       read_radial
c       read_core
c       read_pot
c       Mcast
c
c       if(it_mode/=0) then
c         iterate
c       else
c         read_gpt
c         read_ene
c       endif
c       Mcast
c
c       read_wavef
c
c       MINCPW
c       ENERGY
c       DELTAEX
c       deg
c       waveftrafo
c       ALIGNBD
c       Mcast
c
c       wannier
c       efermi
c       BANDINFO
c       DIPOLE
c       checksum
c       Mcast

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"

      subroutine getinput

      use key
      use file
      use global
      use util
      use wrapper
      use readwrite
      use m_fft
      use arrays, only: selfc
      use, intrinsic :: iso_fortran_env
      Mpi ( use Mwrapper )
      Mpi2( use, intrinsic :: iso_c_binding )

      implicit none
      integer                     :: it_mode,it_stop
      integer                     :: iunit,idum,ispin,itype,ineq,ineq1,i,j,k,l,ll,m1,m2,s
      integer                     :: ikpt,ieq,icent,icent1,m,n,nn,ind,ijob
      integer                     :: iband1,iband2,iband,bandi,bandf,nsym1,mincpw
      integer,    allocatable     :: iarr(:),iarr1(:),iarr2(:,:)
      integer,    allocatable     :: band1(:),spin1(:),kpt1(:)
      real_dp                     :: rdum,rdum1,rdum2,kvec(3),f1,f2,f
      complex_dp                  :: cdum
      real_dp,    allocatable     :: rarr(:),rarr2(:,:),rarr3(:,:,:)
      real_dp,    allocatable     :: eig(:),evec(:,:),ham(:,:)
      complex_dp, allocatable     :: olapmt(:,:,:,:)
      MCOMPLEX_dp,allocatable     :: cpwhlp(:,:),olappw(:,:)
      complex_dp, allocatable     :: olap(:,:)
      integer                     :: isym
      logical                     :: wrtkpt,ldum
      logical,    allocatable     :: done(:)
      character(2)                :: zlabel(112) = (/
     &  'H ',                                                                                'He',
     &  'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne',
     &  'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar',
     &  'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     &  'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',
     &  'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &            'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     &  'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
     &            'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn' /)
      character(7)                :: jobarg(15) = (/'HF     ','GW     ','GT     ','GWT    ','RPA    ','HFE    ','SUSCEP ',
     &                                              'SUSCEPR','DIELEC ','SCREEN ','SCREENW','SX     ','COSX   ','PBE0   ',
     &                                              'KS     '/)
      character(80), allocatable  :: lines(:)
      character(80)               :: line,parse_out(9)
      character(256)              :: longline
      real_dp                     :: ch2r,intgrf
      complex_dp                  :: wfolap
      integer                     :: istat
      integer                     :: ch2i
      type(symtype), allocatable  :: sym1(:)
      real                        :: time1,time2
      MCOMPLEX_dp                 :: stepfunction
      logical                     :: isinteger,isreal

# ifndef noPID
      integer                     :: getpid
# endif
# ifdef MPI
      type(c_ptr)                 :: ptr
      integer                     :: Merr,nnode
      MCOMPLEX_dp, pointer_cnt    :: Ninit_cpw(:)
      complex_dp,  pointer_cnt    :: Ninit_cmt(:)
# endif
# include "interface/band_info.inc"

      if(storage_size(i)    /= 4*8) Error('Datatype integer is not 4 bytes.')
      if(storage_size(ldum) /= 4*8) Error('Datatype logical is not 4 bytes.')
      if(storage_size(time1)/= 4*8) Error('Datatype real is not 4 bytes.')
      if(storage_size(rdum) /= 8*8) Error('Datatype real(real64) is not 8 bytes.')
      if(storage_size(cdum) /=16*8) Error('Datatype complex(real64) is not 16 bytes.')

# ifdef MPI
      call mpi_type_size(mpi_integer,         i,Merr) ; if(i/=4)  Error('MPI datatype integer is not 4 bytes.')
      call mpi_type_size(mpi_logical,         i,Merr) ; if(i/=4)  Error('MPI datatype logical is not 4 bytes.')
      call mpi_type_size(mpi_real,            i,Merr) ; if(i/=4)  Error('MPI datatype real is not 4 bytes.')
      call mpi_type_size(mpi_double_precision,i,Merr) ; if(i/=8)  Error('MPI datatype double precision is not 8 bytes.')
      call mpi_type_size(mpi_double_complex,  i,Merr) ; if(i/=16) Error('MPI datatype double complex is not 16 bytes.')
# endif

      if( error_unit/=0) Error('Unit for standard error is different from 0.')
      if( input_unit/=5) Error('Unit for standard input is different from 5.')
      if(output_unit/=6) Error('Unit for standard output is different from 6.')

      maxmem     = 0
      nalloc_    = 0
      it_mode    = 0
      it_stop    = 0
      nkpt_path  = 0
      kptadd     = 0
      fac_kptadd = 1
      lkptadd    = .false.
      latpar     = 1
      divergence = 0
      wrtinfo    = .false.
      wrtext     = .false.
      wrtkpt     = .false.
      inp        = 1
      option     = 0
      kpt_order  = 0
      nwan       = 0
# ifdef MPI
      Ncomm      = mpi_comm_null
      Ocomm      = mpi_comm_null
      nnode      = 1
# endif

      nullify(selfc)
      nullify(coulomb)
      nullify(cmt)
      nullify(cpw)
      nullify(vmt)
      nullify(vmt_xc)
      nullify(vpw)
      nullify(vpw_xc)
      nullify(cmtu)
      nullify(cpwu)

c
c Command line options
      inpfile  = 'spex.inp'
      i        = 0
      do
        i = i + 1
        Rcall get_command_argument(i,longline)
        Mpi( call Mcast(longline) )
        if(longline(:6)=='--inp=') then
          inpfile = longline(7:)
        else if(longline=='-w') then
          wrtkpt  = .true.
        else if(longline(:3)=='-o=') then
          option  = ch2i(longline(4:),' ') ; Rwrite(6,'(A,I3)') 'Special option: ',option
        else if(longline/=' ') then
          Error('Unknown command line argument: '//trim(longline))
        else
          exit
        endif
      enddo

      Rbegin

      ! tmpdir currently not used
# if 0
      call get_environment_variable('TMPDIR',tmpdir)
      if(tmpdir(len(tmpdir):len(tmpdir))/=' ') then
        Warn('Environment variable TMPDIR too long. Will not be used.')
        tmpdir = ' '
      endif
      if(tmpdir==' ') then
        open(1,file='/tmp/.spextouch',iostat=istat)
        if(istat==0) then
          close(1,status='delete')
          tmpdir='/tmp'
        endif
      endif       
# endif

      inp = fopen(trim(inpfile),action='read')
# ifdef MPI
#   ifdef noSHARED
      Info('noSHARED version: Autoset MPISPLIT SHRD=1')
      line = 'SHRD=1'
#   else
      call getkey(inp,'MPISPLIT',line,default='NODE')
#   endif
      if(line=='NODE'.or.line=='SHARED') then
        Msplit = mpi_comm_type_shared  ; if(Msplit<=-huge(0)+Msize0) Error('MPISPLIT NODE/SHARED handle in SHRD range.')
      else if(line=='SOCKET') then
#   ifdef OPEN_MPI
        Msplit = ompi_comm_type_socket ; if(Msplit<=-huge(0)+Msize0) Error('MPISPLIT SOCKET handle in SHRD range.')        
#   else
        Error('MPISPLIT SOCKET only available for OpenMPI.')
#   endif
      else if(line(:5)=='TYPE=') then
        Msplit = ch2i(line(6:),'MPISPLIT TYPE=...') ; if(Msplit<=-huge(0)+Msize0) Error('MPISPLIT TYPE handle in SHRD range.')        
      else if(line(:5)=='SHRD=') then
        Msplit = ch2i(line(6:),'MPISPLIT SHRD=...') ; if(Msplit<1) Error('SHRD parameter out of range.')
        Msplit = min(Msize,Msplit)
      else
        Error('MPISPLIT argument unknown.')
      endif
      Rend
      call Mcast(line)
      call Mcast(Msplit)      
      call mpi_comm_split_type(Mcomm,mpi_comm_type_shared,0,mpi_info_null,Ncomm,Merr) ! Define physical node (->Nsize0,Nrank0,Orank0)
      call mpi_comm_size(Ncomm,Nsize,Merr)                                            ! Ncomm will be redefined as communicator in virtual node below
      call mpi_comm_rank(Ncomm,Nrank,Merr)
      Nsize0 = Nsize
      Nrank0 = Nrank
#   ifndef noPID
      allocate(Mpid(0:Nsize0-1))
      Mpid         = 0
      Mpid(Nrank0) = getpid() ! get PIDs of all processes in physical node
      call Msum(Mpid,comm=Ncomm)      
#   endif
#   ifdef WRTMEM
      Ncommon(mem0)     ! define common variable for total (private) memory used in physical node
      Ncommon(mem_prv0) !
      if(Nrank0==0) then
        mem0     = 0
        mem_prv0 = 0
      endif
      Nfence_(mem0)
      Nfence_(mem_prv0)
#   endif
      if(Nrank0==0) then
        call mpi_comm_split(Mcomm,0,0,Ocomm,Merr)
        call mpi_comm_rank(Ocomm,Orank0,Merr)
        call mpi_comm_free(Ocomm,Merr)
      else
        call mpi_comm_split(Mcomm,MPI_UNDEFINED,0,Ocomm,Merr)
      endif
      call Mcast(Orank0,comm=Ncomm)  ! Orank0: rank of physical node
      call mpi_comm_free(Ncomm,Merr) ! Ncomm freed, will be redefined in def_nodes below
      Ocomm = mpi_comm_null
      Ncomm = mpi_comm_null
      if(line(:5)=='SHRD=') then
        if(Nsize0<Msplit) Error('SHRD parameter larger than number of physical cores.')
        if(Mrank0/=Msize-1.and.Nrank0==Nsize-1) then
          if(mod(Nsize0,Msplit)/=0) Error('SHRD parameter incomensurate with node.')
        endif
        Msplit = Msplit - huge(0) ! handle for def_nodes
      endif
      call def_nodes
      Ncommon(mem) ! define common variable for total memory used in virtual node
      ifO mem = 0
      Nfence_(mem)
      nnode = (Nsize0-1) / Nsize + 1
      if(nnode>1) then ; Rwrite(6,'(A)') 'MPI: '//Chr(Msize)//','//Chr(Osize)//' ('//Chr(nnode)//')'
      else             ; Rwrite(6,'(A)') 'MPI: '//Chr(Msize)//','//Chr(Osize)
      endif
# endif

      ! unique identification number for current run (part of scratch file names to avoid conflict of concurrent runs)
      uin = 0
# ifndef noPID
      ifR uin = getpid()
      Mpi( call Mcast(uin) )
# endif

      CHKMEM0(MPInodes)

      beginSingle

      write(6,'(//A)') '### subroutine: getinput ###'


      ! Allocation of parameters whose dimension is known. (Routine getkey needs allocatable arrays.)
      allocate ( nkpt3(3),gauss(2),smooth(2) )

c
c     Check input file (spex.inp)
      write(6,'(/A)') '##################'
      write(6,'(A)')  '### Input file ###'
      write(6,'(A)')  '##################'
      call checkkeys(inp)

c
c Read basic LAPW parameters
      call read_param
      call vprod(rlat(:,1),lat(:,2),lat(:,3))
      call vprod(rlat(:,2),lat(:,3),lat(:,1))
      call vprod(rlat(:,3),lat(:,1),lat(:,2))
      vol     = dot_product(lat(:,1),rlat(:,1))
      rvol    = 8*pi**3 / vol
      rlat    = 2*pi * rlat / vol
      maxlcut = maxval(lcut) ; allocate(nindx(0:maxlcut,ntype))
      nindx   = 2
      do itype = 1,ntype
        do n = 1,nlo(itype)
          l              = llo(n,itype) ; if(l>lcut(itype)) Error('l quantum number of local orbital exceeds lcut.')
          nindx(l,itype) = nindx(l,itype) + 1
        enddo
      enddo
      maxindx = maxval(nindx)
      write(6,'(/A)')                    'Basic parameters:'
      write(6,'(A,I3)')                  '  Number of spins            =',nspin
      write(6,'(A,I3)')                  '            centers          =',ncent
      write(6,'(A,I3)')                  '            types            =',ntype
      write(6,'(A,'//chr(ntype)//'I3)')  '            equivalent atoms =',neq
      if(latpar<0) then
        write(6,'(A,F10.5,A)')           '  Lattice parameter  =   ',1d0,'  (undefined)'
      else
        write(6,'(A,F10.5)')             '  Lattice parameter  =   ',latpar
      endif
      write(6,'(A,3F10.5,/25x,3F10.5,/25x,3F10.5)')
     &                                   '  Primitive vectors  =   ',lat
      write(6,'(A,F13.5)')               '  Unit-cell volume   =',   vol
      write(6,'(A,3F10.5,/25x,3F10.5,/25x,3F10.5)')
     &                                   '  Reciprocal vectors =   ',rlat
      write(6,'(A,F13.5)')               '  Reciprocal volume  =',   rvol
      write(6,'(A,F10.5)')               '  Reciprocal cutoff  =   ',gcut
      if(invs)  write(6,'(A)')           '  System has inversion symmetry.'
      write(6,'(/A)')                    'Unit-cell geometry:'
      write(6,'(A)')                     '  # Ty  El              Coord.'
      if(latpar<0) then
        if(latpar<0) Warn('Lattice parameter is undefined in DFT input; scaling factor set to 1.')
        latpar = 1
      endif
      allocate(atype(ncent))
      icent = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          icent        = icent + 1
          atype(icent) = itype
          write(6,'(2I3,2X,A2,3F10.5)') icent,itype,zlabel(ztype(itype)),cent(:,icent)
        enddo
      enddo

      ! RESTART
      call getkey(inp,'RESTART',line,status=istat)
      if     (istat==0)              then ; restart = restart_def('02200')
      else if(istat==1.or.line=='1') then ; restart = restart_def('22233')
      else if(            line=='2') then ; restart = restart_def('33333')
      else                                ; restart = restart_def(line)
      endif

      ! ITERATE
      call getkey(inp,'ITERATE', lines, status=istat)
      it_elow = -huge(it_elow)
      if(istat==1) then
        if(l_soc) then ; it_mode = 3
        else           ; it_mode = 2
        endif
      else if(istat==2) then
        i = 2
        select case(lines(1))
          case('NR')   ; it_mode = 1
          case('SR')   ; it_mode = 2
          case('FR')   ; it_mode = 3
          case default ; it_mode = 2 ; i = 1 ; if(l_soc) it_mode = 3
        end select
        if(size(lines)>=i) then
          line = lines(i) ; if(line(len_trim(line)-1:)=='eV') line(len_trim(line)-1:) = ' '
          if(isreal(line)) then
            it_elow = ch2r(lines(i),'*ITERATE')
            i       = i + 1
          endif
          if(size(lines)>=i) then
            if     (lines(i)=='STOP')  then ; it_stop = 1
            else if(lines(i)=='BANDS') then ; it_stop = 2 ; i = i + 1
            else                            ; Error('Unknown argument after ITERATE: '//trim(lines(i)))
            endif
            if(size(lines)>i) Error('Too many arguments after ITERATE.')
          endif
        endif
      endif
      if(allocated(lines)) deallocate(lines)
      ! ignore ITERATE if RESTART requested
      if(iand(restart,R_KS)/=0.and.it_stop==0.and.it_mode/=0) then
        inquire(file=file_wavef,exist=ldum)
        if(ldum) then
          it_mode = 0
          Info('Wavefunctions will be read from restart file as requested (RESTART). ITERATE is ignored.')
        endif
      endif
      ! redefine l_soc
      select case(it_mode)
        case(1) ; if(l_soc) Warn('ITERATE NR used on top of a calculation with SOC. SOC will be disabled.')
                  l_soc = .false.
        case(2) ; if(l_soc) Warn('ITERATE SR used on top of a calculation with SOC. SOC will be disabled.')
                  l_soc = .false.
        case(3) ; if(.not.l_soc) Warn('ITERATE FR used on top of a calculation w/o SOC. SOC will be enabled.')
                  l_soc = .true.
                  sqa   = 0
      end select      

c Write SOC information
      if(l_soc) then
        write(6,'(/A)')      'Spin-orbit coupling included.'
        write(6,'(A,3F8.4)') 'Spin quantization axis:',sin(sqa(1))*cos(sqa(2)),sin(sqa(1))*sin(sqa(2)),cos(sqa(1))
        if(any(sqa/=0)) Warn('Spin quantization axis different from z. Use with care.')
      endif

c Check if sigma_xc exists
      inquire(file='qsgw',exist=l_qsgw)
      if(l_qsgw.and.it_mode>0.and.it_stop==0) Error('No QSGW with ITERATE.')

c
c Read symmetry information. Then prepare the following arrays.
c   sym(i)%rrot   : rotation matrix wrt reciprocal lattice vectors [ = transpose ( sym(sym(i)%inv)%rot ) ]
c   sym(i)%inv    : inverse symmetry operation
c   sym(i)%symmor : true if symmetry operation is symmorphic
c   pcent(c,i)    : / rotation of atom center c through operation i leads to
c   tcent(:,c,i)  : \ cent(:,pcent(c,i)) + tcent(:,c,i)  [ tcent is a lattice vector ]
c If the system has no inversion symmetry, combinations with time reversal symmetry are appended.
      call getkey(inp,'TRSOFF',trsoff,default=.false.)
      if(l_soc.and.nspin==2) trsoff = .true.
      call read_symmetry
      invsym = 0
      do i = 1,nsym
        if(all(sym(i)%rot==reshape([-1,0,0,0,-1,0,0,0,-1],[3,3])).and.all(sym(i)%transl==0)) invsym = i
      enddo
      nsym1 = nsym
      nsymt = nsym
      if(invs.and.invsym==0) Error('Did not find index of inversion symmetry operation.')
      if((invsym==0.or.l_soc).and..not.trsoff) then ! ".not.trsoff" only needed if we change the default above
        allocate   ( sym1(nsym) ) ; sym1         = sym
        deallocate ( sym )        ; nsym         = 2*nsym
        allocate   ( sym(nsym) )  ; sym(:nsym/2) = sym1
      else if(invsym/=0) then
        if(.not.invs) then
          Warn('Found inversion symmetry operation. You seem to use an executable compiled without -DINV')
          write(0,'(A)') '          on a system with inversion symmetry. You should better use an executable compiled with -DINV.'
          write(0,'(A)') '          Anyway, time reversal symmetry will not be used.'
        endif
        trsoff = .true.
      endif
      allocate ( pcent(ncent,nsym),tcent(3,ncent,nsym),rarr(3),symtab(nsym,nsym) )
      if(any(sym(1)%rot-reshape([1,0,0,0,1,0,0,0,1],[3,3])/=0))
     &  Error('First symmetry operation is not the identity.')
      ! Write out properties of symmetry operations
      write(6,'(/A)')   'Symmetry information:'
      write(6,'(A,I3)') '  Number of symmetry operations =',nsym
      write(6,'(A)')    '   #  Det           Axis           Angle   Syp  Inv'
      symtab = 0
      do i = 1,nsym1
        write(6,'(2X,I2,''  '''NoA) i
        call rotprop( rdum1,rdum2,rarr,i )
        if(rdum2==0) then
          write(6,'(I3,1X,3A7,A9'NoA) nint(rdum1),'    -- ','    -- ','    -- ','     --  '
        else
          rarr = rarr / maxval(abs(rarr)) ! originally: rarr = rarr / abs(rarr(minloc(abs(rarr),1,abs(rarr)>1d-10)))
          write(6,'(I3,1X,3F7.2,F9.2'NoA) nint(rdum1),rarr,rdum2/pi*180
        endif
        if(any(sym(i)%transl/=0)) then
          write(6,'(3X,''No '''NoA)
          sym(i)%symmor = .false.
        else
          write(6,'(3X,''Yes'''NoA)
          sym(i)%symmor = .true.
        endif
        sym(i)%inv = 0
        do j = 1,nsym1
          do k = 1,nsym1
            kvec = matmul(sym(i)%rot,sym(j)%transl) + sym(i)%transl - sym(k)%transl
            if(all(matmul(sym(i)%rot,sym(j)%rot)==sym(k)%rot).and.sum(abs(kvec-nint(kvec)))<1d-10) then
              if(symtab(i,j)/=0) Error('Error in symmetry multiplication table.')
              symtab(i,j) = k
              if(k==1) then
                write(6,'(I4)') j
                if(sym(i)%inv/=0) Error('inverse operation already defined.')
                sym(i)%inv  = j
                sym(i)%rrot = transpose_int ( sym(j)%rot ) ! temporary fix for ifc
              endif
            endif
          enddo

          if(symtab(i,j)==0) Error('Could not define symmetry multiplication table.')
          if(all(matmul(sym(i)%rot,sym(j)%rot)==reshape([1,0,0,0,1,0,0,0,1],[3,3]))
     &        .and.all(modulo1r(matmul(sym(i)%rot,sym(j)%transl)+sym(i)%transl)<1d-10)) then
            if(j/=sym(i)%inv) Error('error.')
          endif
        enddo

        if(sym(i)%inv==0) Error('inverse operation not found.')
        ! Check symmetry operations
        icent = 0
        do itype = 1,ntype
          do ineq = 1,neq(itype)
            icent          = icent + 1
            rarr           = matmul(sym(i)%rot,cent(:,icent)) + sym(i)%transl
            icent1         = sum(neq(:itype-1))
            pcent(icent,i) = 0
            do ineq1 = 1,neq(itype)
              icent1 = icent1 + 1
              if( all( abs( rarr-cent(:,icent1)-nint(rarr-cent(:,icent1)) ) <1d-10 ) ) then
                if(pcent(icent,i)/=0) Error('Double atom mapping.')
                pcent(icent,i)   = icent1
                tcent(:,icent,i) = nint ( rarr-cent(:,icent1) )
                if(maxval(abs(cent(:,icent1)+tcent(:,icent,i)-rarr))>1d-6)
     &            Error('Could not define atom mapping.')
              endif
            enddo
            if(pcent(icent,i)==0) then
              write(6,'(A,I3,A)') 'getinput: Symmetry operation',i,' is inconsistent with the crystal structure.'
              write(6,'(A,I3,A)') '          Mapping of atom',icent,' failed.'
              Error('Symmetry operation inconsistent with crystal structure.')
            endif
          enddo
        enddo
      enddo
      deallocate (rarr)

      ! Define symcent (symcent(ic): symmetry that rotates to ic from parent center)
      allocate(symcent(ncent))
      symcent = 0
      icent   = 1
      do itype = 1,ntype
        do isym = 1,nsym1
          icent1 = pcent(icent,isym)
          if(icent1<1.or.icent1>ncent) Error('pcent pointer incorrect.')
          if(symcent(icent1)==0) symcent(icent1) = isym
        enddo
        icent = icent + neq(itype)
      enddo
      if(any(symcent==0)) Error('Definition of array symcent failed.')

      ! Check sym()%rrot
      allocate ( rarr2(3,3),rarr3(3,3,2) )
      rarr3(:,:,1) = matmul(transpose( lat), lat) / (2*pi)
      rarr3(:,:,2) = matmul(transpose(rlat),rlat) / (2*pi)
      do j = 1,nsym1
        rarr2 = matmul( rarr3(:,:,1) ,matmul(sym(j)%rot,rarr3(:,:,2)))
        if(sum(abs(sym(j)%rrot-rarr2))>1d-12) then
          Error('Non-integer matrix element in sym()%rrot. Symmetry broken. Check geometry.')
        endif
      enddo
      deallocate ( rarr2,rarr3 )

      ! define combinations with time reversal symmetry
      if(.not.trsoff) then
        write(6,'(2X,I2,A,I2,A)') nsym/2+1,' - ',nsym,'  Time reversal symmetries'
        do isym = 1,nsym/2
          sym(nsym/2+isym)%transl     =  sym(isym)%transl
          sym(nsym/2+isym)%rot        =  sym(isym)%rot
          sym(nsym/2+isym)%rrot       = -sym(isym)%rrot
          sym(nsym/2+isym)%inv        =  sym(isym)%inv + nsym/2
          pcent(  :,nsym/2+isym)      =  pcent(  :,isym)
          tcent(:,:,nsym/2+isym)      =  tcent(:,:,isym)
          sym(nsym/2+isym)%esoc       =  matmul ( reshape ( [ 0,1,-1,0 ] , [2,2] ) , sym(isym)%esoc )
          symtab(nsym/2+1:,nsym/2+1:) =  symtab(:nsym/2,:nsym/2)
          symtab(nsym/2+1:,:nsym/2)   =  symtab(:nsym/2,:nsym/2) + nsym/2
          symtab(:nsym/2,nsym/2+1:)   =  symtab(:nsym/2,:nsym/2) + nsym/2
        enddo
      endif

      CHKMEM0(sym)

c     Define k-points and irreducible Brillouin zone (IBZ)
c       nkpt           : number of k points defined in BZ (unshifted set)
c       nkpt1          : number of k points actually stored in the DFT output files (nkpti+nkpti2)
c       nkpt2          : number of k points (unshifted + shifted set)
c       nkpti          : number of irred. k points (these are the first in the list)
c       nkpti2         : number of irred. k points of shifted set (starting at nkpt+1)
c       kpt(:,nkpt1)   : k vectors (first nkpti vectors define the IBZ)
c       pkpt(i,j,k)    : points to k-point index of ( i/nkpt3(1) , j/nkpt3(2) , k/nkpt3(3) ) (internal coord.)
c       kptp(i)        : parent of ith k point
c       kindx(i)       : points to k index of cmt/cpw that contains data for ith k point { storeibz/lkptadd = FF:1,..,nkpt TF:kptp FT:1,..,nkpt2 TT:kptp(:nkpt),kptp(nkpt+1:)-nkpt+nkpti }
c       symkpt(i)      : symmetry operation that relates parent k point to ith k point
c       kptsym(i,j)    : \ k' = Rk - G  with  k  = kpt(:,i),           R = sym(j)%rrot,
c       gkptsym(:,i,j) : /                    k' = kpt(:,kptsym(i,j)), G = gkptsym(:,i,j)
c       phase(k)       : phase acquired by rotation to k (only for INV)
      if(it_mode/=0.and.it_stop==2) then ; call getkey(inp,'BZ',nkpt3,mini=1,default=[1,1,1]) ! ITERATE BANDS does not need BZ
      else                               ; call getkey(inp,'BZ',nkpt3,mini=1)
      endif
      nkpt = product(nkpt3)
      allocate ( kpt(3,nkpt),kptp(nkpt),symkpt(nkpt),iarr(3) )
      l = 0
      call ibz(kpt,kptp,symkpt,nkpti,[0d0,0d0,0d0],.false.)
      allocate ( pkpt(0:nkpt3(1),0:nkpt3(2),0:nkpt3(3)) )
      pkpt = 0
      do ikpt = 1,nkpt
        iarr                          = nint ( kpt(:,ikpt) * nkpt3 )
        pkpt(iarr(1),iarr(2),iarr(3)) = ikpt
      enddo
      pkpt(nkpt3(1),    :   ,    :   ) = pkpt(0,:,:)
      pkpt(    :   ,nkpt3(2),    :   ) = pkpt(:,0,:)
      pkpt(    :   ,    :   ,nkpt3(3)) = pkpt(:,:,0)
      if(any(pkpt<=0).or.any(pkpt>nkpt)) Error('Definition of pkpt pointer failed.')
      deallocate ( iarr )

      call getkey(inp,'KPT', lines, default=[' '] )
      allocate ( rarr(3),iarr(3) )

      do
      nkpt_label = 0
      do i = 1,size(lines)
        if(lines(i)==' ') cycle
        call parser(parse_out,lines(i),'1=<<2*||2/3*>><(4,5,6)||[4,5,6]>','KPT')
        if(len_trim(parse_out(1))>10) Error('Maximally ten characters allowed for KPT label.')
        if(isinteger(parse_out(1))) Error('Numbers not allowed as KPT labels.')
        rdum = 1
        if(parse_out(3)/=' ') then ; rdum = ch2r(parse_out(3),'KPT') ; if(parse_out(1)=='+') fac_kptadd = rdum ; endif
        if(parse_out(2)/=' ')        rdum = ch2r(parse_out(2),'KPT')/rdum
        do j = 1,3 ; rarr(j) = ch2r(parse_out(j+3),'KPT') ; enddo ; rarr = rdum * rarr
        if(index(lines(i),'(')==0) then
          if(latpar==1) Error('Lattice parameter unset (=1). [...] definition cannot be used.')
          rarr = matmul(rarr,lat) / latpar ! unit is 2pi/a
        endif
        nkpt_label = nkpt_label + 1
        if(allocated(kpt_label)) then
          do k = 1,nkpt_label
            if(parse_out(1)==ckpt_label(k)) Error('Double definition of KPT label '//trim(parse_out(1)))
          enddo
          kpt_label(:,nkpt_label) = rarr ! stored without back-folding into first BZ
          ckpt_label(nkpt_label)  = parse_out(1)
          if(any(abs(modulo1r(rarr)-rarr)>1d-10)) then
            Info('k point is converted from ('//Chn(rarr,',')//') to ('//Chn(modulo1r(rarr),',')//')')
            rarr = modulo1r(rarr)
          endif
          if(parse_out(1)(1:1)=='+') then
            if(lkptadd) Error('Double definition of additional k point (''+'').')
            if(all(rarr==0)) then
              Info('Zero vector not treated as additional k vector.')
            else
              lkptadd = .true.
              kptadd  = rarr
            endif
          endif
        endif
      enddo
      if(allocated(kpt_label)) exit
      allocate(kpt_label(3,nkpt_label),ckpt_label(nkpt_label))
      enddo

      deallocate ( lines )

      ! Define k-point mapping wrt symmetry operations (->kptsym,gkptsym)
      if(.not.lkptadd) then
        nkpti2 = 0
        nkpt1  = nkpti
        nkpt2  = nkpt
      else
        call reallocate ( kpt, 3,2*nkpt )
        call reallocate ( kptp,  2*nkpt )
        call reallocate ( symkpt,2*nkpt )
        call ibz(kpt(:,nkpt+1:),kptp(nkpt+1:),symkpt(nkpt+1:),nkpti2,kptadd,.false.)
        kptp(nkpt+1:) = kptp(nkpt+1:) + nkpt
        nkpt1         = nkpti+nkpti2
        nkpt2         = 2*nkpt
        allocate ( pkpt1(0:nkpt3(1),0:nkpt3(2),0:nkpt3(3)) )
        pkpt1 = 0
        do ikpt = nkpt+1,nkpt2
          iarr                           = nint( modulo1r(kpt(:,ikpt)-kptadd) * nkpt3 )
          pkpt1(iarr(1),iarr(2),iarr(3)) = ikpt
        enddo
        pkpt1(nkpt3(1),    :   ,    :   ) = pkpt1(0,:,:)
        pkpt1(    :   ,nkpt3(2),    :   ) = pkpt1(:,0,:)
        pkpt1(    :   ,    :   ,nkpt3(3)) = pkpt1(:,:,0)
        if(any(pkpt1<=nkpt).or.any(pkpt1>nkpt2)) Error('Definition of pkpt1 pointer failed.')
      endif

      allocate ( kptsym(nkpt2,nsym),gkptsym(3,nkpt2,nsym) )
      kptsym  = 0
      gkptsym = 0
      do i = 1,nkpt2
        do k = 1,nsym
          if(i<=nkpt) then
            rarr = matmul(sym(k)%rrot,kpt(:,i))
            iarr = nint ( modulo1(rarr) * nkpt3 )
            idum = pkpt(iarr(1),iarr(2),iarr(3))
            kvec = kpt(:,idum)
          else
            rarr = matmul(sym(k)%rrot,kpt(:,i))
            kvec = modulo1r(rarr-kptadd) * nkpt3
            iarr = nint ( kvec ) ; if( any(abs( kvec - iarr )>1d-12) ) cycle
            idum = pkpt1(iarr(1),iarr(2),iarr(3))
            kvec = kpt(:,idum)
          endif
          kptsym(i,k)    = idum
          iarr           = nint ( rarr - kvec )
          gkptsym(:,i,k) = iarr
          if(idum==0.or.any(abs(rarr-iarr-kvec)>1d-8)) then
            write(0,'(A,I3)') 'getinput: K-point mapping failed for symmetry operation',k
            write(0,'(A)')    '          (K-point mesh and symmetry consistent?)'
            Error(' ')
          endif
        enddo
      enddo
      deallocate (rarr,iarr)

c      do i = 1,nkpti
c        do j = 2,nsym
c          if(kptsym(i,j)==i) then
c            if(sum(abs( reshape([0,1,1,0],[2,2]) - abs(sym(j)%esoc) ))<1d-10) then
c              write(*,'(A,$)') ' **'
c            else
c              write(*,'(A,$)') '   '
c            endif
c            write(*,'(2I3,4(2F9.3,2X))') i,j,sym(j)%esoc
c          endif
c        enddo
c      enddo
c      stop

      if(any(kpt<0)) Bug('Found negative component or signed zero in kpoint array.')

      ! Calculate phases
      allocate ( phase(nkpt2) )
      phase = 1
# ifdef INV
      do ikpt = 1,nkpt2
        if(kptp(ikpt)/=ikpt) then
          isym        = symkpt(ikpt)
          phase(ikpt) = exp(-img * 2*pi * dot_product(kpt(:,ikpt),sym(isym)%transl) )
        endif
      enddo
# endif

      ! --------------------------------------------------------------------

      if(wrtkpt) then
        Info('Autoset WRTKPT')
      else
        call getkey(inp,'WRTKPT', wrtkpt, default=.false.)
      endif
      ! Write kpts file if WRTKPT is set
      if(wrtkpt) call write_kpt

      call getkey(inp,'KPTPATH',  lines, status=i)
      if(i==1) then
        Error('Argument missing after KPTPATH.')
      else if(i>1) then
        n = len_trim(lines(1))
        if(size(lines)>2) Error('Expected one or two arguments after KPTPATH.')
        if(lines(1)(:1)//lines(1)(n:n)/='()') Error('Expected (...) as first argument after KPTPATH.')
        line      = lines(1)(2:n-1)
        nkpt_path = count ( [ (line(i:i)==',',i=1,len(line)) ] ) + 1
        allocate ( kpt_path(3,nkpt_path) )
        do i = 1,nkpt_path-1
          ind           = index(line,',') ; if(line=='+') Error('extra point ''+'' not allowed in path.')
          kpt_path(:,i) = getkptr(line(:ind-1))
          line          = line(ind+1:)
          if(i/=1) then
            if(all(kpt_path(:,i)==kpt_path(:,i-1))) Error('Successive points are identical in path.')
          endif
        enddo
        kpt_path(:,nkpt_path) = getkptr(line)
        mkpt_path             = 0 ! default
        if(size(lines)==2) mkpt_path = ch2i(lines(2),'KPTPATH')
        if(mkpt_path<0) Error('KPTPATH: Number of steps must not be negative.')
        ! Write qpts file if WRTKPT is set
        if(wrtkpt) call write_qpts(mkpt_path)
        deallocate(lines)
      endif

      endSingle

      Mpi( call Mcast(wrtkpt) )
      if(wrtkpt) call finish

# ifdef MPI
      ! from read_param
      call Mcast(nspin) ; call Mcast(ncent)    ; call Mcast(ntype) ; call Mcastl(neq)    ; call Mcastl(lcut)
      call Mcast(gcut)  ; call Mcastl(ztype)   ; call Mcastl(cent) ; call Mcast(latpar)  ; call Mcast(lat)
      call Mcast(invs)  ; call Mcast(l_soc)    ; call Mcast(sqa)   ; call Mcastl(nlo)    ; call Mcastl(llo)
      ! ... derived
      call Mcast(rlat)  ; call Mcast(vol)      ; call Mcast(rvol)  ; call Mcast(maxlcut) ; call Mcast(maxindx)
      call Mcastl(nindx); call Mcastl(atype)
      ! from read_symmetry
      call Mcast(nsym)
      MnoR( allocate(sym(nsym)) )
      do i = 1,nsym
        call Mcast(sym(i)%rot) ; call Mcast(sym(i)%transl)
      enddo
      ! ... derived
      do i = 1,nsym
        call Mcast(sym(i)%rrot) ; call Mcast(sym(i)%inv) ; call Mcast(sym(i)%symmor) ; call Mcast(sym(i)%esoc)
      enddo
      call Mcast(invsym)  ; call Mcast(nsymt)   ; call Mcast(trsoff) ; call Mcastl(pcent)  ; call Mcastl(tcent)
      call Mcastl(symtab) ; call Mcastl(symcent)
      ! from spex.inp and derived
      call Mcastl(nkpt3)  ; call Mcast(nkpt)    ; call Mcastl(kpt)   ; call Mcastl(kptp)   ; call Mcastl(symkpt)
      call Mcastl(pkpt)   ; call Mcast(kptadd)  ; call Mcast(nkpti)  ; call Mcast(lkptadd)
      call Mcast(nkpti2)  ; call Mcast(nkpt1)   ; call Mcast(nkpt2)  ; call Mcastl(kptsym) ; call Mcastl(gkptsym)
      call Mcastl(phase)
      call Mcast(nkpt_label) ; MnoR( allocate(ckpt_label(nkpt_label)) )
      call Mcastl(kpt_label)
      do i = 1,nkpt_label
        call mpi_bcast(ckpt_label(i),len(ckpt_label(1)),mpi_character,0,Mcomm,Merr)
      enddo
      if(lkptadd) call Mcastl(pkpt1)
# endif

      MnoR( inp = fopen(trim(inpfile),action='read') )

# ifdef INV
      Rif(l_soc)               Error('Spin-orbit coupling; wrong executable, compile without -DINV!')
      Rif(.not.invs)           Error('No inversion symmetry; wrong executable, compile without -DINV!')
# else
      Rif(invs.and..not.l_soc) Error('Inversion symmetry; wrong executable, compile with -DINV!')
# endif

      Rbegin

      ! IBC
      call getkey(inp,'IBC',         oibc, default=0, mini=0, maxi=2)

      ! BLOECHL
      call getkey(inp,'BLOECHL', obloechl, default=.false.)

      ! NBAND
      call getkey(inp,'NBAND',       line, default=' ')
      maxene    = -huge(0d0)
      nband0    = -1
      ene_extra = 0
      if(line/=' ') then
        ind = index(line,'+')
        if(ind/=0) then
          ene_extra  = ch2r(line(ind+1:),'*NBAND') ; if(ene_extra<=0)
     &                                               Error('shift for extra energy must be positive.')
          line(ind:) = ' '
        endif
        if(index(line,'.')/=0) then ; maxene = ch2r(line,'*NBAND')
        else                        ; nband0 = ch2i(line, 'NBAND')
                                      if(nband0<=0) Error('Argument to NBAND must be positive.')
        endif
      endif
      if(ene_extra/=0) Error('Extrapolar approximation not supported anymore.')

      ! CUTZERO
      call getkey(inp,'CUTZERO', cutzero, default=0d0, mine=0d0)

      ! JOB (see jobtype.h) --------------------------------------------------------------------------------------------
      call getkey(inp,'JOB',lines,status=i)
      if(i==0) allocate(lines(0))
      j = 0
      i = 1
      allocate ( job(0) ) ! dummy allocation
      do while(i<=size(lines))
        select case(lines(i))

          case('HF','GW','GT','GWT','SX','COSX','PBE0','KS')
            line = lines(i)
            call job_next(job,j,.false.)
            job(j)%label = line(:6)
            if(line=='HF')   job(j)%type = J_HF
            if(line=='GW')   job(j)%type = J_GW
            if(line=='GWT')  job(j)%type = J_GW
            if(line=='SX')   job(j)%type = J_SX
            if(line=='COSX') job(j)%type = J_COSX
            if(line=='PBE0') job(j)%type = J_PBE0
            if(line=='KS')   job(j)%type = J_KS
            if(line=='GT')   job(j)%type = J_GT
            allocate ( band1(0),kpt1(0),spin1(0) ) ; l = 0
            do
              i = i + 1 ; if(i>size(lines)) then ; if(lines(i-1)==line) Warn('No job parameters after'//job(j)%label)
                                                   exit
                          endif
                          if(any(lines(i)==jobarg)) exit
              if(lines(i)=='FULL') then
                if(job(j)%full)                   Error('Double definition of FULL in this job.')
                if(any(job(j)%type==[J_GT,J_KS])) Error('FULL calculations not possible with this job type.')
                job(j)%full = .true.
                cycle
              endif
              if(nspin==1.or.l_soc) then
                call parser(parse_out,lines(i),'1:(3)','JOB '//trim(line))
                ispin = 1
              else
                call parser(parse_out,lines(i),'1:<2>(3)','JOB '//trim(line))
                if     (parse_out(2)=='u') then ; ispin = 1
                else if(parse_out(2)=='d') then ; ispin = 2
                else if(parse_out(2)==' ') then ; ispin = 0 ! both spins
                else                            ; Error('Spin not correctly specified; use ''u'' or ''d''.')
                endif
              endif
              call str2iarr(iarr,parse_out(3),'JOB '//trim(line)//'  (band definition)')
              call getkptarr(iarr1,parse_out(1),.false.)
              if(ispin==0) then ; idum = l + size(iarr) * size(iarr1) * 2
              else              ; idum = l + size(iarr) * size(iarr1)
              endif
              call reallocate ( band1,idum )
              call reallocate ( kpt1, idum )
              call reallocate ( spin1,idum )
              do s = 1,2 ; if(ispin/=0.and.ispin/=s) cycle
                do m = 1,size(iarr1)
                  do k = 1,size(iarr)
                    l        = l + 1
                    band1(l) = iarr(k)
                    kpt1(l)  = iarr1(m)
                    spin1(l) = s
                  enddo
                enddo
              enddo
              deallocate ( iarr,iarr1 )
            enddo
            allocate ( job(j)%band(l),job(j)%kpt(l),job(j)%spin(l) )
            job(j)%band = band1
            job(j)%kpt  = kpt1
            job(j)%spin = spin1
            deallocate ( band1,kpt1,spin1 )

          case('SUSCEP','SUSCEPR')
            line = lines(i)
            do
              i = i + 1 ; if(i>size(lines)) then ; if(lines(i-1)==line) then ; Error('no job parameters.')
                                                   else                        ; exit ; endif
                          endif
                          if(any(lines(i)==jobarg)) exit
              call job_next(job,j,.false.) ; job(j)%label = line
              if(line=='SUSCEP')  then
                job(j)%type = J_SUS
                if(index(lines(i),'..')/=0) then
                  call parser(parse_out,lines(i),'1:<2>{3..4,5}</8><,OUT=6>','JOB '//trim(line))
                  Info('".." obsolete; use ":" instead.')
                else
                  call parser(parse_out,lines(i),'1:<2>{3<:4,5>}</8><,OUT=6>','JOB '//trim(line))
                endif
              else if(line=='SUSCEPR') then
                job(j)%type = J_SUSR
                if(index(lines(i),'..')/=0) then
                  call parser(parse_out,lines(i),'1:<2>{3..4,5}</8><,KERN=7><,OUT=6>','JOB '//trim(line))
                  Info('".." obsolete; use ":" instead.')
                else
                  call parser(parse_out,lines(i),'1:<2>{3<:4,5>}</8><,KERN=7><,OUT=6>','JOB '//trim(line))
                endif
                if(parse_out(7)/=' ') then
                  job(j)%kern = parse_out(7)(:4)
                  if(job(j)%kern/='BSE') then
# ifndef TDDFT
                    Error('TDDFT kernel specified; wrong executable, compile with -DTDDFT!')
# endif
                  endif
                endif
              endif
              if(parse_out(4)==' ') then
                job(j)%freqfile = parse_out(3) ; l = fopen(job(j)%freqfile) ; call fclose(l) ! check existence of file job(j)%freqfile
              else
                job(j)%freq(1) = ch2r(parse_out(3),'*JOB '//trim(line))
                job(j)%freq(2) = ch2r(parse_out(4),'*JOB '//trim(line))
                job(j)%freq(3) = ch2r(parse_out(5),'*JOB '//trim(line))
                if(job(j)%freq(1)>job(j)%freq(2)) Error('Lower frequency must be smaller than upper frequency.')
                if(job(j)%freq(3)<=0)             Error('Frequency increment must be positive.')
              endif
              call getkptarr(job(j)%kpt,parse_out(1),.false.)
              allocate ( job(j)%spin(1) )
              if     (parse_out(6)=='BIN') then ; job(j)%out = -1
              else if(parse_out(6)/=' ')   then ; job(j)%out = ch2i(parse_out(6),'JOB '//trim(line))
                                           if(job(j)%out<1) Error('Submatrix for output must have finite size.')
              endif
              if(parse_out(2)/=' ') then
                select case(parse_out(2)) ; case('+')    ; job(j)%spin = 0
                                            case('uu')   ; job(j)%spin = 1
                                            case('dd')   ; job(j)%spin = 2
                                            case('ud')   ; job(j)%spin = 3
                                            case('du')   ; job(j)%spin = 4
                                            case default ; Error('Spin label must be +, uu, dd, ud or du.')
                end select
                if(nspin==1.and.job(j)%spin(1)>1) Error('System is not spin-polarized.')
                if(line=='SUSCEPR') then
                  if(job(j)%spin(1)==1.or. job(j)%spin(1)==2) Error('Spin labels uu/dd not allowed for SUSCEPR.')
                  if(job(j)%kern==' ' .and.job(j)%spin(1)/=0) Error('A spin label but no kernel is given.')
                endif
              else
                job(j)%spin = 0
              endif
c              if(job(j)%kern=='BSE') then
c                call job_next(job,j,.true.)
c                job(j-1)%type      = 4
c                job(j-1)%kpt       = 0
c                job(j-1)%freq(1)   = 0
c                job(j-1)%freq(2)   = 0
c                job(j-1)%label     = 'SCREEN'
c                job(j-1)%wrtdielec = .false.
c                job(j-1)%wrtscreen = .false.
c                job(j-1)%wrtsuscep = .false.
c                job(j-1)%kern      = ' '
c                if(.not.restart) Error('Use RESTART option for BSE calculation.')
c              endif
              if(parse_out(8)/=' ') call cutfreqc(j,parse_out(8))
            enddo

          case('DIELEC','SCREEN')
            line = lines(i)
            do
              i = i + 1 ; if(i>size(lines)) then ; if(lines(i-1)==line) then ; Error('no job parameters.')
                                                   else                      ; exit ; endif
                          endif
                          if(any(lines(i)==jobarg)) exit
              call job_next(job,j,.false.) ; job(j)%label = line(:6)
              if(line=='DIELEC') job(j)%type = J_DIEL
              if(line=='SCREEN') job(j)%type = J_SCR
              if(index(lines(i),'..')/=0) then
                call parser(parse_out,lines(i),'1:{2<..3,4>}</6><,OUT=5>','JOB '//line(:6))
                Info('".." obsolete; use ":" instead.')
              else
                call parser(parse_out,lines(i),'1:{2<:3,4>}</6><,OUT=5>','JOB '//line(:6))
              endif
              if(parse_out(3)==' ') then
                job(j)%freqfile = parse_out(2) ; k = fopen(job(j)%freqfile) ; call fclose(k) ! check existence of file job(j)%freqfile
              else
                job(j)%freq(1) = ch2r(parse_out(2),'*JOB '//trim(line))
                job(j)%freq(2) = ch2r(parse_out(3),'*JOB '//trim(line))
                job(j)%freq(3) = ch2r(parse_out(4),'*JOB '//trim(line))
                if(job(j)%freq(1)>job(j)%freq(2).or.job(j)%freq(1)<0.or.job(j)%freq(3)<=0) then
                  write(6,'(A/A)') 'getinput: Frequency range definition wrong:',
     &                             '          a:b,c with 0 >= a >= b and c > 0'
                  Error('Frequency range definition wrong.')
                endif
              endif
              call getkptarr(job(j)%kpt,parse_out(1),.false.)
              allocate ( job(j)%spin(1) )
              job(j)%spin = 0
              if     (parse_out(5)=='BIN') then ; job(j)%out = -1
              else if(parse_out(5)/=' ')   then ; job(j)%out = ch2i(parse_out(5),'JOB '//trim(line))
                                           if(job(j)%out<1) Error('Submatrix for output must have finite size.')
              endif
            enddo
            if(parse_out(6)/=' ') call cutfreqc(j,parse_out(6))

          case('SCREENW')
            line = lines(i)
            do
              i = i + 1 ; if(i>size(lines)) then ; if(lines(i-1)==line) then ; Error('no job parameters.')
                                                   else                      ; exit ; endif
                                                   endif
                          if(any(lines(i)==jobarg)) exit
              call job_next(job,j,.false.)
              job(j)%label = line(:6)
              job(j)%type  = J_SCRW
              if(index(lines(i),'..')/=0) then
                call parser(parse_out,lines(i),'<2>{3..4,5}</7><,SHELL=6>','JOB '//trim(line))
                Info('".." obsolete; use ":" instead.')
              else
                call parser(parse_out,lines(i),'<2>{3<:4,5>}</7><,SHELL=6>','JOB '//trim(line))
              endif
              if(parse_out(4)==' ') then
                if(parse_out(3)=='0') then
                  job(j)%freq = [ 0d0,0d0,1d0 ]
                else
                  job(j)%freqfile = parse_out(3) ; l = fopen(job(j)%freqfile) ; call fclose(l) ! check existence of file job(j)%freqfile
                endif
              else
                job(j)%freq(1) = ch2r(parse_out(3),'*JOB '//trim(line))
                job(j)%freq(2) = ch2r(parse_out(4),'*JOB '//trim(line))
                job(j)%freq(3) = ch2r(parse_out(5),'*JOB '//trim(line))
                if(job(j)%freq(1)>job(j)%freq(2).or.job(j)%freq(1)<0.or.job(j)%freq(3)<=0) then
                  write(6,'(A/A)') 'getinput: Frequency range definition wrong:',
     &                             '          a:b,c with 0 >= a > b and c > 0'
                  Error('Frequency range definition wrong.')
                endif
              endif
              allocate ( job(j)%spin(1) )  ; job(j)%spin = 1
              if(parse_out(6)/=' ') then ; job(j)%out  = ch2i(parse_out(6),'JOB SCREENW')
              else                         ; job(j)%out  = 0
              endif
              if(parse_out(2)/=' ') then
                select case(parse_out(2)) ; case('+')    ; job(j)%spin = 0
                                            case('uu')   ; job(j)%spin = 1
                                            case('dd')   ; job(j)%spin = 2
                                            case('ud')   ; job(j)%spin = 3
                                            case('du')   ; job(j)%spin = 4
                                            case default ; Error('Spin label must be +, uu, dd, ud or du.')
                end select
              else if(nspin==2) then ; Error('System is spin-polarized; SCREENW spin label must be set.')
              endif
              if(nspin==1.and.job(j)%spin(1)>1) Error('System is not spin-polarized.')
            enddo
            if(parse_out(7)/=' ') call cutfreqc(j,parse_out(7))

          case('RPAENE')
            call job_next(job,j,.false.)
            job(j)%type  = J_RPA
            job(j)%label = lines(i)(:6)
            i            = i + 1

          case('HFENE')
            call job_next(job,j,.false.)
            job(j)%type  = J_HFE
            job(j)%label = lines(i)(:6)
            i            = i + 1

          case default
            Error('Unknown job name: '//trim(lines(i)))

        end select
      enddo
      njob = j
      deallocate ( lines )
      ! -----------------------------------------------------------------

      ! WRITE
      call getkey(inp,'WRITE', lines, default=[' '])
      if(lines(1)/=' ') then
        do i = 1,size(lines)
          if     (lines(i)=='SUSCEP') then ; wrtsuscep = .true.
          else if(lines(i)=='DIELEC') then ; wrtdielec = .true.
          else if(lines(i)=='SCREEN') then ; wrtscreen = .true.
          else if(lines(i)=='INFO'  ) then ; wrtinfo   = .true.
          else if(lines(i)=='EXT'   ) then ; wrtext    = .true.
          else ; Error('Arguments to WRITE must be one of SUSCEP, DIELEC, SCREEN, INFO, EXT.')
          endif
        enddo
      endif
      deallocate ( lines )

      ! GAUSS
      call getkey(inp,'GAUSS',       gauss, default=[0d0,0d0], mine=0d0, allow_eV=.true.)

      ! MEM
      call getkey(inp,'MEM',        maxmem, default=10000d0,    mine=0d0) ; maxmem = maxmem*megabyte

      ! TIMING
      call getkey(inp,'TIMING',     otimer, default=.false.)

      ! STOREBZ
      call getkey(inp,'STOREIBZ', storeibz, default=.false.)
      if(storeibz) Info('STOREIBZ is set by default since version 05.06. (STOREBZ gives the old default behavior.)')
      call getkey(inp,'STOREBZ',  storeibz, default=.true.)

      ! NOSYM
      call getkey(inp,'NOSYM',     use_sym, default=.true.)

      ! ALIGNBD
      call getkey(inp,'ALIGNBD',   alignbd, default=.false.)
      if(alignbd.and..not.lkptadd) Error('ALIGNBD only allowed if an additional k point is used.')

      ! BANDINFO
      call getkey(inp,'BANDINFO', bandinfo, default=.false.)
      allocate ( lines(2) )

      ! BANDOMIT
      call getkey(inp,'BANDOMIT',    lines, default=[' ',' '])
      lomit = .false.
      if(lines(1)/=' ') then
        call parser(parse_out,lines(1),'(1)','BANDOMIT')
        call str2iarr(omit,parse_out(1),'BANDOMIT  (band definition)')
        if(size(omit)==0) Warn('Empty band definition after BANDOMIT.')
        if(len_trim(lines(2))/=3) Error('Second argument must be of the form XXX where X=0/1.')
        do i = 1,3
          if     (lines(2)(i:i)=='1') then ; lomit(i) = .true.
          else if(lines(2)(i:i)/='0') then ; Error('Second argument must be of the form XXX where X=0/1.')
          endif
        enddo
      endif
      deallocate ( lines )

      ! GCUT (MBASIS)
      call getkey(inp,'GCUT', gcutm, section='MBASIS', default=gcut*3/4, mini=0d0)

      ! LCUT (MBASIS)
      allocate(lcutm(ntype))
      if(maxindx>0) then
        call getkey(inp,'LCUT', lcutm, section='MBASIS', default=lcut/2, mini=0) ; maxlcutm = maxval(lcutm)
      else
        lcutm = -1
      endif

      ! TOL (MBASIS)
      if(maxindx>0) call getkey(inp,'TOL', rdum, section='MBASIS', default=0.0001d0, mini=0d0)

      ! OPTIMIZE (MBASIS)
      allocate ( lines(2) )
      call getkey(inp,'OPTIMIZE', lines, section='MBASIS',  default=[' ',' '])
      optm    = 0
      optcutm = 0
      if(lines(1)/=' ') then
        if     (lines(1)=='MB') then ; optm = 1
        else if(lines(1)=='PW') then ; optm = 2
        else                         ; Error('Syntax error; use OPTIMIZE MB or OPTIMIZE PW.')
        endif
        if(index(lines(2),'.')/=0) then
          optcutm = ch2r(lines(2),'OPTIMIZE')
          if(optcutm>0) optcutm = -4*pi/optcutm**2 ! negative sign to distinguish from (positive) integer numbers
        else
          optcutm = ch2i(lines(2),'OPTIMIZE')
          if(optcutm<=0) Error('Argument to OPTIMIZE must be real or a positive integer.')
          if(optcutm<10) then
            if(option==101) then
              Warn('OPTIMIZE: Integer value < 10. Stop overridden by option -o=101')
            else
              Error('OPTIMIZE: Integer value < 10. You probably meant to write "'//Chr(nint(optcutm))//'.0"! Override with -o=101')
            endif
          endif
        endif
      endif
      deallocate ( lines )

      ! NOAPW (MBASIS)
      call getkey(inp,'NOAPW',    noapw, section='MBASIS',  default=.false.)

      ! LCUT (WFPROD)
      allocate(lcutp(ntype))
      call getkey(inp,'LCUT',     lcutp, section='WFPROD',  default=lcut, mini=0)
      if(any(lcutp>lcut)) Error('LCUT(WFPROD) cannot be larger than LCUT(MBASIS).')

      ! MTTHR (WFPROD)
      call getkey(inp,'MTTHR',    mtthr, section='WFPROD',  default=0d0,  mini=0d0)

      ! APPROXPW (WFPROD)
      call getkey(inp,'APPROXPW', fullpw, section='WFPROD', default=.true.)

      ! FFT (WFPROD)
      call getkey(inp,'FFT',      lines, section='WFPROD',  default=[' '])
      if     (size(lines)>1)  then ; Error('At most one argument expected.')
      else if(size(lines)==0) then ; gcutf = ( 2*gcut + gcutm ) * 2d0/3
      else
        if     (lines(1)==' ')     then ; gcutf = 0                                                                     ! default v<=1000
                      if(vol>1200) then ; gcutf = ( 2*gcut + gcutm ) * 2d0/3 ; Info('Autoset FFT '//Chr(gcutf)) ; endif ! default v> 1000
        else if(lines(1)=='EXACT') then ; gcutf = 2*gcut + gcutm
        else if(lines(1)=='OFF')   then ; gcutf = 0
        else                            ; gcutf = ch2r(lines(1),'FFT');if(gcutf<=0) Error('Argument must be positive.')
        endif
      endif
      deallocate ( lines )

      ! NOSTORE (COULOMB)
      call getkey(inp,'NOSTORE', coulstore, section='COULOMB', default=.true.)

      ! LEXP (COULOMB)
      call getkey(inp,'LEXP',         lexp, section='COULOMB', default=max(14,2*maxlcutm), mini=2*maxlcutm, writeout=.false.)
# if defined(MPI) || defined(LOAD)
      if(coulstore) then; coulstore = .false. ; Info('Autoset COULOMB->NOSTORE') ; endif
# endif

      ! TETRAF (SUSCEP)
      call getkey(inp,'TETRAF',     tetraf, section='SUSCEP',  default=.false.)
      if(gauss(1)/=0.and.tetraf) Warn('Keyword TETRAF ignored.')
      if(allocated(lines)) deallocate ( lines )

      ! FREQINT (SENERGY)
      call getkey(inp,'FREQINT',      line, section='SENERGY', default='AUTO')
      if     (line=='AUTO')   then ; freqint = 0
      else if(line=='SPLINE') then ; freqint = 1
      else if(line=='PADE')   then ; freqint = 2
      else if(line=='NONLIN') then ; freqint = 3
      else                         ; Error('Wrong argument to FREQINT; either SPLINE or PADE.')
      endif

      ! LOGDIV (SENERGY)
      call getkey(inp,'LOGDIV',    ologdiv, section='SENERGY', default=.false.)

      ! SMOOTH (SENERGY)
      call getkey(inp,'SMOOTH',       iarr, section='SENERGY', default=[1,1], mini=1)
      smooth = 1
      if     (size(iarr)==0) then ; Error('Missing arguments after SMOOTH.')
      else if(size(iarr)>2)  then ; Error('Too many arguments after SMOOTH.')
      else                        ; smooth(:size(iarr)) = iarr
      endif
      deallocate(iarr)

      ! ZERO (SENERGY)
      call getkey(inp,'ZERO',        ozero, section='SENERGY', default=.false. )

      ! VXC (SENERGY)
      call getkey(inp,'VXC',          line, section='SENERGY', default='CALC' )
      if     (line=='READ') then ; ovxc = 0
      else if(line=='CALC') then ; ovxc = 1
      else                       ; Error('Expecting one of ''READ'' and ''CALC'' (VXC).')
      endif
      if((any(job%type==J_PBE0).or.it_mode/=0).and.ovxc==0) then ; ovxc = 1 ; Info('Autoset VXC CALC') ; endif
      if(allocated(lines)) deallocate ( lines )

      ! WSCALE (WANNIER)
      call getkey(inp,'WSCALE', wscale, section='WANNIER', default=1d0, mini=0d0)
      if(wscale==0) then
        do ijob = 1,njob
          if(job(ijob)%type==J_SUSR.and.job(ijob)%kern=='BSE') then
            rdum = job(ijob)%freq(1) / job(ijob)%freq(3)
            if(job(ijob)%freq(1)>0.or.job(ijob)%freq(2)<0.or.abs(rdum-nint(rdum))>1d-10)
     &                              Error('First SUSCEPR job must contain frequency 0 for WSCALE=0.')
            if(job(ijob)%kpt(1)/=1) Error('Use Gamma point (1) in first SUSCEPR job for WSCALE=0.')
            exit
          endif
        enddo
      endif

      if(allocated(kpt_path)) then
        ! write k points defined in KPTPATH
        write(6,'(/A,I3)') 'K points in path:'
        call getkptarr(iarr,'PATH',.true.)
        write(6,'(I5,I6,''  ('',F7.5,'','',F7.5,'','',F7.5,'')'')') (i,iarr(i),kpt(:,iarr(i)),i=1,size(iarr))
        deallocate ( iarr )
      endif

      CHKMEM0(spex.inp)

      Rend

# ifdef MPI
      call Mcastl(lcutm)    ; call Mcastl(lcutp)    ; call Mcastl(omit)     ; call Mcast(nband0)  ; call Mcast(bandinfo)
      call Mcast(oibc)      ; call Mcast(optm)      ; call Mcast(lexp)      ; call Mcast(freqint) ; call Mcastl(smooth)
      call Mcast(njob)      ; call Mcast(it_mode)   ; call Mcast(it_stop)   ; call Mcast(it_elow) ; call Mcast(cutzero)
      call Mcast(gcutm)     ; call Mcast(tolerance) ; call Mcast(optcutm)   ; call Mcast(mtthr)   ; call Mcast(gcutf)
      call Mcastl(gauss)    ; call Mcast(noapw)     ; call Mcast(fullpw)    ; call Mcast(ozero)   ; call Mcast(maxmem)
      call Mcast(coulstore) ; call Mcast(tetraf)    ; call Mcast(alignbd)   ; call Mcast(lomit)   ; call Mcast(ologdiv)
      call Mcast(storeibz)  ; call Mcast(obloechl)  ; call Mcast(l_soc)     ; call Mcast(maxlcutm)
      call Mcast(wrtsuscep) ; call Mcast(wrtdielec) ; call Mcast(wrtscreen) ; call Mcast(wrtinfo) ; call Mcast(wrtext)
      call Mcast(restart)   ; call Mcast(maxene)    ; call Mcast(use_sym)
      MnoR( allocate(job(njob)) )
      if(njob/=0) then
        call Mcast(job(:)%indx)
        call Mcast(job(:)%full)
        call Mcast(job(:)%type)
        call Mcast(job(:)%out)
        do i = 1,njob
          call Mcastl(job(i)%band)
          call Mcastl(job(i)%kpt)
          call Mcastl(job(i)%spin)
          call Mcast(job(i)%kern)
          call Mcast(job(i)%label)
          call Mcast(job(i)%freqfile)
          call Mcast(job(i)%freq)
        enddo
      endif
      if(nnode/=1) then
        maxmem = maxmem / nnode
        RInfo(Chr(nnode)//' virtual nodes defined by MPISPLIT. MEM set to '//Chf(maxmem/megabyte,'F12.2')//' MB per virtual node.')
      endif
# endif      

c
c     Define k index
      allocate ( kindx(nkpt2) )
      if(storeibz) then
        if(lkptadd) then
          kindx(:nkpt)   = kptp(:nkpt)
          kindx(nkpt+1:) = kptp(nkpt+1:) - nkpt + nkpti
        else
          kindx = kptp
        endif
      else
        kindx = [(i,i=1,nkpt2)]
      endif

c
c     Define numbers of spins
      nspin1 = nspin
      nspin2 = nspin
      nspin3 = 1
      if(l_soc) then
        nspin1 = 1 ! dimension of nband, ene, etc. (spin is not a good quantum number anymore)
        nspin2 = 2 ! dimension of cmt, cpw, etc.   (wavefunctions written as two-component spinors)
        nspin3 = 2 ! spin components per state
      endif

c     Define job indices (trivial)
      do i = 1,njob
        job(i)%indx = i
      enddo

c     Define   fac(i) = i!
c             sfac(i) = sqrt(i!)
c           facfac(i) = (2*i+1)!!
      n    = 0
      rdum = 1
      do while(rdum<=huge(1d0)/(2*n+3.0001))
        n    = n + 1
        rdum = rdum * (2*n+1)
      enddo
      allocate ( fac(0:n),sfac(0:n),facfac(-1:n) )
      fac(0)       = 1
      sfac(0)      = 1
      facfac(-1:0) = 1
      do i = 1,n
        fac(i)    = fac(i-1)    * i
        sfac(i)   = sfac(i-1)   * sqrt(i*1d0)
        facfac(i) = facfac(i-1) * (2*i+1)
      enddo

c     Determine Wigner-D matrices
      ll = max(maxlcut,maxlcutm)
      allocate ( dwgn(-ll:ll,-ll:ll,0:ll,nsym) )
      do isym = 1,nsym
        call dwigner(dwgn(:,:,:,isym),matmul(lat,matmul(sym(isym)%rot,transpose(rlat)))/(2*pi),ll)
        ! time reversal symmetry: multply with D(m1,m2) = delta(-m1,m2) * (-1)**m2
        !                         which fulfills conjg(Y_lm) = SUM(m1) Y_lm1 D(m1,m)
        if(isym>nsymt) then
          do l = 0,ll
            do m1 = -l,l
              do m2 = -l,-1
                cdum                = dwgn(m1, m2,l,isym)
                dwgn(m1, m2,l,isym) = dwgn(m1,-m2,l,isym) * (-1)**m2
                dwgn(m1,-m2,l,isym) = cdum                * (-1)**m2
               enddo
            enddo
          enddo
        endif
      enddo

      allocate(grid(ntype))

      beginSingle

c Read radial grids and functions
      call read_radial
      ! Check allocation state of rgrid(:,:) and gridf(:,:)
      ogrid = 0
      if(allocated(rgrid)) then
        if(allocated(gridf)) then ! integration weights provided
          ogrid = 1
          Info('Radial integration weights provided by DFT code.')
          if(any(abs(rgrid(1,:))>1d-12)) Info('Radial grid does not start at zero. Low-r contribution determined from power law.')
          if(oibc/=0) Error('Using non-exponential grids with IBC not implemented.')
          Error('Sorry, external integration weights not implemented yet.')
        else                      ! interpolation required
          ogrid = 2
          Info('No radial integration weights provided. Radial grid has to be interpolated.')
          write(6,'(/A)') 'Grid interpolation (first, increment):'
          do itype = 1,ntype
            if(abs(rgrid( grid(itype)%number,itype)-grid(itype)%radius)>1d-12) Error('Last radial grid point not on MT boundary.')
            if(any(rgrid(:grid(itype)%number,itype)>grid(itype)%radius+1d-12)) Error('There are grid points beyond MT boundary.')
            if(any(rgrid(:grid(itype)%number,itype)<                  -1d-12)) Error('There are negative radial grid points.')
            grid(itype)%first = rgrid(1,itype)
            if(abs(rgrid(1,itype))<1d-12) then
              Info('Radial grid starts at zero. Zero point removed for standard grid definition.')
              grid(itype)%first = rgrid(2,itype)
            endif
# ifdef IPOLTEST
#   warning IPOLTEST (test of radial grid interpolation)
            grid(itype)%first = grid(itype)%first * 1.5
# endif
            grid(itype)%increment = log( grid(itype)%radius/grid(itype)%first ) / (grid(itype)%number-1)
            write(6,'(A,I3,F15.10,F10.5)') '  type =',itype,grid(itype)%first,grid(itype)%increment
          enddo
          allocate(rgrid0(maxgrid,ntype))
          rgrid0 = rgrid
          deallocate(rgrid)
          do ispin = 1,nspin
            do itype = 1,ntype
              do l = 0,lcut(itype) ; m = max(2-l,1)
                do n = 1,nindx(l,itype)
                  rdum = ubas(n,l,itype,ispin) + grid(itype)%radius * dubas(n,l,itype,ispin)
                  call interpolate_radial(bas1(:,n,l,itype,ispin),rgrid0(:,itype),rdum,m,grid(itype))
                  call interpolate_radial(bas2(:,n,l,itype,ispin),rgrid0(:,itype),0d0, m,grid(itype))
                enddo
              enddo
            enddo
          enddo
        endif
      endif
      ! Trivial copy in case of l_soc and nspin=1
      if(l_soc.and.nspin==1) then
        bas1(:,:,:,:,2) = bas1(:,:,:,:,1)
        bas2(:,:,:,:,2) = bas2(:,:,:,:,1)
        ubas(:,:,:,2)   = ubas(:,:,:,1)
        dubas(:,:,:,2)  = dubas(:,:,:,1)
        ebas(:,:,:,2)   = ebas(:,:,:,1)
      endif

c
c Calculate grid points
c (Adjust first mesh point using nested intervals if last radial grid point deviates from MT radius. "Crack a nut with a sledghammer!")
      if(ogrid/=1) then
        allocate ( rgrid(maxgrid,ntype) )
        do itype = 1,ntype
          rdum = def_get_rgrid(itype,grid(itype)%first)
          rdum = abs(rdum-grid(itype)%radius)
          if(rdum>1d-14) then
            if(rdum>1d-8)
     &      Error('Last radial grid point deviates from MT boundary by too much: '//Chf(rdum,'ES8.1'))            
            Warn( 'Last radial grid point deviates from MT boundary by '//Chf(rdum,'ES8.1')//'. Try to adjust...')
            rdum  = rdum / grid(itype)%radius
            rdum1 = grid(itype)%first * (1-2*rdum) ; f1 = def_get_rgrid(itype,rdum1)
            rdum2 = grid(itype)%first * (1+2*rdum) ; f2 = def_get_rgrid(itype,rdum2)
            i     = 0
            do
              rdum = (rdum1+rdum2) / 2             ; f  = def_get_rgrid(itype,rdum)
              if(max(f1,f2)<grid(itype)%radius.or.min(f1,f2)>grid(itype)%radius)
     &                              Error('Nested interval search failed for radial grid.')
              i = i + 1 ; if(i>100) Error('Nested interval search for radial grid did not converge.')                            
              if(abs(f-grid(itype)%radius)<1d-14) then ; exit
              else if(f>grid(itype)%radius)       then ; f2 = f ; rdum2 = rdum
              else                                     ; f1 = f ; rdum1 = rdum
              endif
            enddo
            grid(itype)%first = rdum
          endif
        enddo
      endif

c
c Fast numerical integration with intgrf after this initialization
      if(.not.allocated(gridf)) call intgrf_init

c Read core electron wavefunctions
      allocate(lcutc(ntype))
      call read_core
      call getkey(inp,'CORESOC',ldum, default=.false.)
      if(maxlcutc>maxlcut)        Error('l(core) > l(basis).')
      if(ldum.and..not.lcore_soc) Error('Core states are not SOC splitted (see read_core). Cannot use CORESOC.')
      nelec = sum(neq(:)*ztype(:))
      do itype = 1,ntype
        do l = 0,lcutc(itype)
          if(lcore_soc.and.l/=0) then ; nelec = nelec - nindxc(l,itype) * neq(itype) * (2*l+1)
          else                        ; nelec = nelec - nindxc(l,itype) * neq(itype) * (2*l+1) * 2
          endif
        enddo
      enddo
      if(lcore_soc) then
        if(ldum) then ! average core states over spin
          if(nspin==2) then
            ecore(:,:,:,1)   = ( ecore(:,:,:,1)   + ecore(:,:,:,2)   ) / 2
            core1(:,:,:,:,1) = ( core1(:,:,:,:,1) + core1(:,:,:,:,2) ) / 2
            core2(:,:,:,:,1) = ( core2(:,:,:,:,1) + core2(:,:,:,:,2) ) / 2
            ecore(:,:,:,2)   = ecore(:,:,:,1)
            core1(:,:,:,:,2) = core1(:,:,:,:,1)
            core2(:,:,:,:,2) = core2(:,:,:,:,1)
          endif
        else          ! average core states over SOC
          do itype = 1,ntype
            do l = 1,lcutc(itype)
              do n = 1,nindxc(l,itype),2
                ecore(  n/2+1,l,itype,:) = ( ecore(  n,l,itype,:) * (2*l+2) + ecore(  n+1,l,itype,:) * (2*l) ) / (4*l+2)
                core1(:,n/2+1,l,itype,:) = ( core1(:,n,l,itype,:) * (2*l+2) + core1(:,n+1,l,itype,:) * (2*l) ) / (4*l+2)
                core2(:,n/2+1,l,itype,:) = ( core2(:,n,l,itype,:) * (2*l+2) + core2(:,n+1,l,itype,:) * (2*l) ) / (4*l+2)
              enddo
              nindxc(l,itype) = nindxc(l,itype) / 2
              maxindxc        = max(maxindxc,nindxc(l,itype))
            enddo
          enddo
          call reallocate ( ecore ,         maxindxc,maxlcutc,ntype,nspin )
          call reallocate ( core1 , maxgrid,maxindxc,maxlcutc,ntype,nspin )
          call reallocate ( core2 , maxgrid,maxindxc,maxlcutc,ntype,nspin )
        endif
      endif
      lcore_soc = ldum

      ! CORES
      allocate ( cores(maxindxc,0:maxlcutc,ntype),lines(ntype) )
      call getkey(inp,'CORES',lines,default=[('()',i=1,ntype)])
      cores = .false.
      do itype = 1,ntype
        if(lines(itype)/='()') then
          call parser(parse_out,lines(itype),'(1)','CORES')
 2        ind = index(parse_out(1),',') - 1
          if(ind==-1) ind = len_trim(parse_out(1))
          if(ind==0)  Error('syntax error in CORES definition.')
          line = parse_out(1)(:ind)
          if(line==' ') Error('syntax error in CORES definition.')
          l = -1
          do ll = 0,12
            if(line(2:2)==lchar(ll)) l = ll
          enddo
          if(l==-1) then
            if(line(2:2)=='') then ; Error('missing label for l quantum number (CORES): '//trim(line))
            else                   ; Error('unknown label for l quantum number (CORES): '//line(2:2))
            endif
          endif
          n = ch2i(line(:1),'CORES') - l
          if(n<1)   Error('wrong principal number.')
          if(    lcore_soc .and. all(line(3:)/=['   ',trim(chr(abs(2*l-1)))//'/2',trim(chr(2*l+1))//'/2']) .or.
     &      .not.lcore_soc .and.     line(3:)/= '   ' ) Error('unknown core state: '//trim(line))
          if(lcore_soc.and.l>0) then
            if(line(3:)=='3/2') then ; n = 2*n-1
            else                     ; n = 2*n
            endif
          endif
          if(l>lcutc(itype).or.n>nindxc(l,itype)) then
            write(0,'(A,A,A,I3)') 'getinput: core state ',trim(line),' does not exist in atom type',itype
            Error(' ')
          endif
          cores(n,l,itype) = .true.
          if(lcore_soc.and.l>0.and.line(3:)==' ') cores(n-1,l,itype) = .true.
          if(ind/=len_trim(parse_out(1))) then
            line         = parse_out(1)(ind+2:)
            parse_out(1) = line
            goto 2
          endif
        endif
      enddo
      deallocate(lines)

c Write out information about radial grids and functions
c      if(any(grid(:)%number<7)) then
c        write(6,'(A)') 'getinput: Number of grid points is smaller than 7.'
c        Error('Number of grid points is smaller than 7.')
c      endif
      if(any(grid(:)%number>0)) then
        if(all(mod(grid(:)%number,6)==1)) then
          write(6,'(/A)') 'Info: Only Simpson integration is used.'
        else
          write(6,'(/A)') 'Info: Simpson and Lagrange integration is used.'
        endif
      endif
      write(6,'(/A)')                   'Radial mesh and basis functions:'
      write(6,'(A,'//chr(ntype)//'I5)') '  Number of radial grid points:    ',grid(:)%number
      write(6,'(A,'//chr(ntype)//'I5)') '  Number of radial basis functions:',(sum(nindx(:,itype)),itype=1,ntype)
      write(6,'(A,'//chr(ntype)//'I5)') '  Angular momentum quantum number: ',lcut
      write(6,'(A,I5)')                 '  Maximal index:                   ',maxindx
      rdum = 4*pi/3 * sum(grid%radius**3)
      write(6,'(A,F15.5,A)')            '  Muffin-tin volume:           ',rdum,'  ('//Chr(nint(rdum/vol*100))//'%)'

c Read in potential, lattice harmonics and stars
c Construct and diagonalize Hamiltonian (only for keyword ITERATE)
      maxlmindx = maxval( [ ( sum( [ (nindx(l,itype)*(2*l+1),l=0,lcut(itype)) ] ),itype=1,ntype) ] )
      call read_pot(1)
# ifdef MPI
      ! Reallocate vpw, vpw_xc, vmt, and vmt_xc as shared memory
      allocate(iarr2(3,4))
      iarr2(1,:) = lbound(vpw)
      iarr2(2,:) = ubound(vpw) ; if(any(iarr2(:2,4)/=[1,nspin])) Bug('4th dimension of vpw is not 1:nspin')
      iarr2(3,:) = iarr2(2,:) - iarr2(1,:) + 1
      deallocate(vpw,vpw_xc,vmt,vmt_xc)
      endSingle
      call Mcastl(iarr2)
      call Mcast(maxgrid)
      call Mcast(maxlh)
      Nallocate ( vmt,    (S_ maxgrid,maxlh,ntype,nspin S_) )
      Nallocate ( vmt_xc, (S_ maxgrid,maxlh,ntype,nspin S_) )
      Nallocate ( vpw_,   (S_ iarr2(3,1),iarr2(3,2),iarr2(3,3),nspin S_) ) ! vpw and vpw_xc have lower bounds different from 1
      Nallocate ( vpw_xc_,(S_ iarr2(3,1),iarr2(3,2),iarr2(3,3),nspin S_) ) 
      vpw   (iarr2(1,1):iarr2(2,1),iarr2(1,2):iarr2(2,2),iarr2(1,3):iarr2(2,3),1:nspin) => vpw_    ; win_vpw    = win_vpw_ 
      vpw_xc(iarr2(1,1):iarr2(2,1),iarr2(1,2):iarr2(2,2),iarr2(1,3):iarr2(2,3),1:nspin) => vpw_xc_ ; win_vpw_xc = win_vpw_xc_
      deallocate(iarr2)
      beginSingle
      call read_pot(1)
# else
      CHKMEM(vpw) ; CHKMEM(vpw_xc) ; CHKMEM(vmt) ; CHKMEM(vmt_xc)
# endif      
      ! Interpolate (if necessary)
      if(ogrid==2) then
        do ispin = 1,nspin
          do itype = 1,ntype
            do l = 1,nlh(itype)
              call interpolate_radial(vmt   (:,l,itype,ispin),rgrid0(:,itype),0d0,3,grid(itype))
              call interpolate_radial(vmt_xc(:,l,itype,ispin),rgrid0(:,itype),0d0,3,grid(itype))
            enddo
            do l = 0,lcutc(itype) ; m = max(2-l,1)
              do i = 1,nindxc(l,itype)
                call interpolate_radial(core1(:,i,l,itype,ispin),rgrid0(:,itype),0d0,m,grid(itype))
                call interpolate_radial(core2(:,i,l,itype,ispin),rgrid0(:,itype),0d0,m,grid(itype))
              enddo
            enddo
          enddo
        enddo
      endif
      ! Spherical part
      do ispin = 1,nspin
        do itype = 1,ntype
          n                        = grid(itype)%number
          vmt   (:n,1,itype,ispin) = vmt   (:n,1,itype,ispin) / rgrid(:n,itype) / sqrt(4*pi) ! undo scaling with r and include Y_00
          vmt_xc(:n,1,itype,ispin) = vmt_xc(:n,1,itype,ispin) / rgrid(:n,itype)              ! undo scaling with r
        enddo
      enddo

c
c Write information about core states
      if(maxindxc==0) then
        write(6,'(/A)') 'No core electrons'
      else if(.not.lcore_soc) then
        write(6,'(/A)') 'Core states:'
        do ispin = 1,nspin
          if(nspin==2.and.ispin==1) write(6,'(A)') '  Spin up'
          if(nspin==2.and.ispin==2) write(6,'(A)') '  Spin down'
          do itype = 1,ntype
            if(all(nindxc(:lcutc(itype),itype)==0)) cycle
            write(6,'(A,I3,A)') '  Atom type',itype,'  ('//zlabel(ztype(itype))//')'
            write(6,'(A)')      '  nl      Energy'
            do l = 0,lcutc(itype)
              do n = 1,nindxc(l,itype)
                write(6,'(I3,A,F12.5)') n+l,lchar(l),ecore(n,l,itype,ispin)
              enddo
            enddo
          enddo
        enddo
      else if(lcore_soc.and.nspin==1) then
        write(6,'(/A)') 'Fine structure of core states:'
        do itype = 1,ntype
          if(all(nindxc(:lcutc(itype),itype)==0)) cycle
          write(6,'(A,I3,A)') '  Atom type',itype,'  ('//zlabel(ztype(itype))//')'
          write(6,'(A)')      '  nlj        Energy'
          do l = 0,lcutc(itype)
            do n = 1,nindxc(l,itype),min(l+1,2)
              if(l==0) then ; write(6,'(I3,A,  ''1/2'',F12.5)') n,        lchar(l),      ecore(n,  l,itype,1)
              else          ; write(6,'(I3,A,I1,''/2'',F12.5)') (n+1)/2+l,lchar(l),2*l-1,ecore(n+1,l,itype,1)
                              write(6,'(I3,A,I1,''/2'',F12.5)') (n+1)/2+l,lchar(l),2*l+1,ecore(n,  l,itype,1)
              endif
            enddo
          enddo
        enddo
      else if(lcore_soc.and.nspin==2) then
        write(6,'(/A)') 'Fine structure of core states with Zeeman splitting:'
        allocate(eig(2),evec(2,2),ham(2,2))
        do itype = 1,ntype
          write(6,'(A,I3,A)') '  Atom type',itype,'  ('//zlabel(ztype(itype))//')'
          write(6,'(A'NoA)    '  nl   mj      Energy' ; if(lcutc(itype)>0) write(6,'(15X,A'NoA) '    Core'
          write(6,*)
          do l = 0,lcutc(itype)
            do n = 1,nindxc(l,itype),min(l+1,2)
              rdum  = intgrf( ( core1(:,n,l,itype,1) * core1(:,n,l,itype,1) +
     &                          core2(:,n,l,itype,1) * core2(:,n,l,itype,1) ) * (vmt(:,1,itype,1)-vmt(:,1,itype,2))/2, itype )
              if(l>0) then
                rdum1 = intgrf( ( core1(:,n+1,l,itype,1) * core1(:,n+1,l,itype,1) +
     &                            core2(:,n+1,l,itype,1) * core2(:,n+1,l,itype,1) ) * (vmt(:,1,itype,1)-vmt(:,1,itype,2))/2, itype )
                rdum2 = intgrf( ( core1(:,n,  l,itype,1) * core1(:,n+1,l,itype,1) +
     &                            core2(:,n,  l,itype,1) * core2(:,n+1,l,itype,1) ) * (vmt(:,1,itype,1)-vmt(:,1,itype,2))/2, itype )
              endif
              do m = -(2*l+1),2*l+1,2
                if(abs(m)==2*l+1) then
                  if(l==0) then ; write(6,'(I3,A,I3,''/2'',F12.5)') n,        lchar(l),m,ecore(n,l,itype,1) + m*rdum / (2*l+1)
                  else          ; write(6,'(I3,A,I3,''/2'',F12.5)') (n+1)/2+l,lchar(l),m,ecore(n,l,itype,1) + m*rdum / (2*l+1)
                  endif
                else
                  ham(1,1) = ecore(n,  l,itype,1) + m*rdum  / (2*l+1)
                  ham(2,2) = ecore(n+1,l,itype,1) - m*rdum1 / (2*l+1)
                  ham(1,2) = -sqrt(1d0*(2*l+1-m)*(2*l+1+m))*rdum2 / (2*l+1)
                  ham(2,1) = ham(1,2)
                  call diagonalize(evec,eig,ham)
                  write(6,'(I3,A,I3,''/2'',2F12.5,3X,2F9.4)') (n+1)/2+l,lchar(l),m,eig(2:1:-1),evec(:,2)
                endif
              enddo
            enddo
          enddo
        enddo
        deallocate(eig,evec,ham)
      endif
      write(6,'(A,I5)') '  Number of valence electrons:',nelec

      ! Average Coulomb potential
      rdum = 0
      do ispin = 1,nspin
        do itype = 1,ntype
          rdum = rdum + neq(itype) * 4*pi * ! vmt and vmt_xc defined in terms of Y_lm except for spherical vmt!
     &                  intgrf( rgrid(:,itype)**2 * (vmt(:,1,itype,ispin)-vmt_xc(:,1,itype,ispin)/sqrt(4*pi)) , itype )
        enddo
        do i = lbound(vpw,1),ubound(vpw,1)
        do j = lbound(vpw,2),ubound(vpw,2)
        do k = lbound(vpw,3),ubound(vpw,3)
          rdum = rdum + (vpw(i,j,k,ispin)-vpw_xc(i,j,k,ispin))
        enddo
        enddo
        enddo
      enddo
      write(6,'(/A,F9.5,''  ('',F9.5,'' eV )'')') 'Average Coulomb potential:',rdum/vol,rdum*hartree/vol

      ! k-point list
      write(6,'(/A)')   'Sampling of Brillouin zone:'
      write(6,'(A,I6)') '  Number of k points:',nkpt
      write(6,'(A,I6)') '              in IBZ:',nkpti
      if(lkptadd) write(6,'(A,I6)')  '             in EIBZ:',nkpti2
      write(6,'(A)')    '  k points of IBZ:'
      write(6,'(I5,''  ('',F7.5,'','',F7.5,'','',F7.5,'')'','//
     &            '''  ['',F8.5,'','',F8.5,'','',F8.5,'' ]'')') (i,kpt(:,i),matmul(rlat,kpt(:,i))*latpar/(2*pi),i=1,nkpti)
      if(lkptadd) then
        write(6,'(A)')  '  k points of EIBZ (shifted set):'
        write(6,'(I5,''  ('',F7.5,'','',F7.5,'','',F7.5,'')'','//
     &              '''  ['',F8.5,'','',F8.5,'','',F8.5,'' ]'')') (i,kpt(:,i),matmul(rlat,kpt(:,i))*latpar/(2*pi),
     &                                                                                               i=nkpt+1,nkpt+nkpti2)
      endif

      endSingle

# ifdef MPI
      ! from read_radial
      call Mcast(grid%number) ; call Mcast(maxgrid) ; call Mcast(grid%first) ; call Mcast(grid%increment) ; call Mcast(grid%radius)
      call Mcastl(bas1)       ; call Mcastl(bas2)   ; call Mcastl(ubas)      ; call Mcastl(dubas)
      call Mcastl(ebas)       ; call Mcast(ogrid)
      ldum = allocated(rgrid0); call Mcast(ldum)
      if(ldum) call Mcastl(rgrid0)
      ! ... derived
      call Mcastl(rgrid)      ; call Mcastl(gridf)
      ! from read_core
      call Mcastl(lcutc)      ; call Mcastl(nindxc) ; call Mcastl(ecore)     ; call Mcastl(core1)         ; call Mcastl(core2)
      call Mcast(maxlcutc)    ; call Mcast(maxindxc); call Mcastl(cores)     ; call Mcast(lcore_soc)
      ! ... derived
      call Mcast(nelec)       ; call Mcast(maxlmindx)
      ! from read_pot
      call Mcastl(nlh)        ; call Mcast(maxlh)   ; call Mcastl(llh)       ; call Mcastl(mlh)           ; call Mcastl(nmlh)
      call Mcastl(clh)
      Nfence(vpw) ; Nfence(vpw_xc) ; Nfence(vmt) ; Nfence(vmt_xc)
      Obegin
      call Mcast(vmt,comm=Ocomm) ; call Mcast(vmt_xc,comm=Ocomm) ; call Mcast(vpw,comm=Ocomm) ; call Mcast(vpw_xc,comm=Ocomm)
      Oend
      Nfence(vpw) ; Nfence(vpw_xc) ; Nfence(vmt) ; Nfence(vmt_xc)
# endif
      CHKMEM(core1) ; CHKMEM(core2) ; CHKMEM(bas1) ; CHKMEM(bas2)

      if(it_mode/=0) then
        NoInv( if(invsym/=0.and.l_soc) then ; call iterate_inv_soc(it_mode,it_stop) ; else )
        call iterate(it_mode,it_stop)
        NoInv( endif )
      else
        allocate ( nband(nkpt2,nspin2),ngpt(nkpt2) ) ; nband = 0 ; ngpt = 0
      endif

      beginSingle

c Read G points
      if(it_mode==0) call read_gpt

c Read energies
      if(it_mode==0) then
        call read_ene
        maxeband = minval(nband(:nkpt1,:nspin1))
        if(maxband<nband0) Error('Did not find enough bands. Reduce NBAND to at least '//Chr(maxband)//' or omit NBAND.')
        do ispin = 1,nspin1
          do ikpt = 1,nkpt1
            do n = 1,nband(ikpt,ispin)-1
              if(ene(n+1,ikpt,ispin)-ene(n,ikpt,ispin)<-1d-12) Error('Energies are not in ascending order. Check read_ene!')
            enddo
            if(nband(ikpt,ispin)<maxband) ene(nband(ikpt,ispin)+1:,ikpt,ispin) = 0
          enddo
        enddo
        if(maxene/=-huge(0d0)) then ! NBAND real
          do ispin = 1,nspin1
            do ikpt = 1,nkpt1
              nband(ikpt,ispin) = count(ene(:nband(ikpt,ispin),ikpt,ispin)<=maxene)
            enddo
          enddo
          maxband  = maxval(nband)
          maxeband = min(maxeband,maxband)
        else if(nband0/=-1) then ! NBAND integer
          do ispin = 1,nspin1
            do ikpt = 1,nkpt1
              if(nband(ikpt,ispin)>nband0) then
                n = nband0
                do while(abs(ene(n,ikpt,ispin)-ene(nband0+1,ikpt,ispin))<edeg)
                  n = n - 1 ; if(n==0) Error('Reached band index zero while cutting degenerate subspace.')
                enddo
                nband(ikpt,ispin) = n
              endif
            enddo
          enddo
          maxband  = maxval(nband)
          maxeband = min(maxeband,maxband)
        endif
        call reallocate( ene , maxband,nkpt2,nspin1 )
        do ikpt = nkpti2,1,-1
          ene(:,nkpt+ikpt,:) = ene(:,nkpti+ikpt,:)
        enddo
      endif
      if(l_soc.and.nspin1==1) nband(:,2) = nband(:,1)

c Distribute G points on related k points
      maxgpt = size(gpt1,2)
      if(lkptadd) then
        do ikpt = nkpti2,1,-1
          gpt1(:,:,nkpt+ikpt) = gpt1(:,:,nkpti+ikpt)
          ngpt(nkpt+ikpt)     = ngpt(nkpti+ikpt)
        enddo
      endif
      do ikpt = 1,nkpt2
        if(kptp(ikpt)/=ikpt) then
          ngpt(ikpt) = ngpt(kptp(ikpt))
          do j = 1,ngpt(ikpt)
            gpt1(:,j,ikpt) = matmul(sym(symkpt(ikpt))%rrot,gpt1(:,j,kptp(ikpt))) +
     &                       gkptsym(:,kptp(ikpt),symkpt(ikpt))
          enddo
        endif
      enddo

c Calculate pointers pgpt and pntgpt, and unique G vectors gpt
      allocate ( pgpt(maxgpt,nkpt2),iarr(3) )
      allocate ( pntgpt(minval(gpt1(1,:,:)):maxval(gpt1(1,:,:)),
     &                  minval(gpt1(2,:,:)):maxval(gpt1(2,:,:)),
     &                  minval(gpt1(3,:,:)):maxval(gpt1(3,:,:))) )
      ! define gpt and pgpt
      ldum   = .false.
 3    pntgpt = 0 ! pntgpt temporary array
      j      = 0
      do ikpt = 1,nkpt2
        do i = 1,ngpt(ikpt)
          iarr = gpt1(:,i,ikpt)
          if(pntgpt(iarr(1),iarr(2),iarr(3))==0) then
            j = j + 1
            if(ldum) gpt(:,j) = iarr
            pntgpt(iarr(1),iarr(2),iarr(3)) = j
          endif
          pgpt(i,ikpt) = pntgpt(iarr(1),iarr(2),iarr(3))
        enddo
      enddo
      if(.not.ldum) then
        ldum    = .true.
        ngptall = j
        allocate ( gpt(3,ngptall) )
        goto 3
      endif
      ! define pntgpt
      pntgpt = 0
      do i = 1,ngptall
        iarr = gpt(:,i)
        pntgpt(iarr(1),iarr(2),iarr(3)) = i
      enddo
      deallocate ( iarr )

c Calculate Fourier series of step function (e.g., for overlap matrices ...)
      if(.not.allocated(cstep)) then
        allocate ( iarr(3) )
        iarr = [ ( maxval(gpt(i,:)) - minval(gpt(i,:)) + 1 , i = 1,3 ) ]
        Allocate_ ( cstep,( -iarr(1):iarr(1) , -iarr(2):iarr(2) , -iarr(3):iarr(3) ) )
        do i = -iarr(1),iarr(1)
        do j = -iarr(2),iarr(2)
        do k = -iarr(3),iarr(3)
          cstep(i,j,k) = stepfunction( [i,j,k] )
        enddo
        enddo
        enddo
        deallocate ( iarr )
      endif

      endSingle

      if(allocated(gpt1)) deallocate( gpt1 )

# ifdef MPI
      ! from read_gpt
      call Mcast(ngpt)
      ! from read_ene (nband, ene might later be redefined and broadcasted again)
      call Mcast(nband) ; call Mcast(maxband) ; call Mcast(maxeband) ; call Mcastl(ene)
      ! ... derived
      call Mcastl(pgpt) ; call Mcastl(pntgpt) ; call Mcast(ngptall)  ; call Mcastl(gpt)
      call Mcast(maxgpt)
# endif

c Allocate arrays cpw and cmt
# ifndef LOAD
      if(it_mode==0) then
        call checkmem('before wavef alloc',0d0)
        Rif(storeibz) then ; n = nkpt1
        else               ; n = nkpt2
        endif
        Mpi( call Mcast(n) )
        Nallocate0 ( cpw,(S_ maxgpt,maxband,n,nspin2 S_) )
        Nallocate0 ( cmt,(S_ maxlmindx,ncent,maxband,n,nspin2 S_) )
        call checkmem('after wavef alloc',0d0)
      endif

c Read plane-wave and muffin-tin coefficients of wavefunctions (->cpw,cmt)
      if(it_mode==0) then
        call cpu_time(time1)
        Rwrite(6,'(/A'NoA) 'Read wavefunctions from harddisc... '
        call read_wavef
        Rcall cpu_done(time1)
      endif
      CHKMEM0(wavef done)
# endif

      beginSingle

      n = (nelec+1) / 2
      if(maxeband<n) then
        if(maxeband==0)               Error('Number of bands is zero.')
        if(maxene/=-huge(0d0)) then ; Error('Number of bands too small. Increase maximal energy (NBAND)!')
        else                        ; Error('Number of bands too small. Increase NBAND to at least '//Chr(n))
        endif
      endif

      if(allocated(omit)) then
        if(any(omit<1).or.any(omit>maxband)) Warn('Band index after BANDOMIT out of range.')
      endif

c Add a "zero" function such that sum(|cpw|**2) becomes minimal.
      call getkey(inp,'MINCPW', mincpw, section='WFPROD',  default=0, mini=0)
      if(mincpw>0) then
        Mpi( Error('MINCPW not implemented for MPI version.') )
        if(fullpw)   Warn('MINCPW should not be used without APPROXPW.')
        if(storeibz) Error('MINCPW only implemented together with STOREBZ.')
        Warn('MINCPW not tested lately!')
        write(6,'(/A)') 'Add "zero" functions such that sum(|cpw|**2) becomes minimal ...'
        write(6,'(A)')  '  Smallest eigenvalues'
        do ikpt = 1,nkpt1
          n = ngpt(ikpt) ; if(mincpw>n) Error('MINCPW value exceeds number of G vectors.')
          allocate ( rarr(n),olappw(n,n),cpwhlp(n,n) )
          if(ikpt<=nkpti) then ; i = ikpt
          else                 ; i = ikpt - nkpti + nkpt
          endif
          call olap_gpt(olappw,n,ikpt,ikpt)              ! calculate IPW overlap
          call diagonalize(cpwhlp,rarr,olappw,-mincpw-1) ! obtain eigenvector of smallest eigenvalue
          if(rarr(1)<0) Warn('Negative eigenvalue of overlap matrix.')
          write(6,'(A,I4,A'NoA) '  K point',i,':'
          write(6,'(ES8.1'NoA)  rarr(:mincpw)
          write(6,*)
          do j = 1,mincpw
            if(abs(rarr(j)-rarr(mincpw+1))>1d-10) then
              write(6,'(15X,A,F12.5'NoA) 'Reduction of sum(|cpw|**2):',
     &          sum(abs(cpw(:n,:maxband,ikpt,:nspin2)**2)) / (maxband*nspin2)
              do ispin = 1,nspin2
                do iband = 1,maxband
                  rdum                     = -dble(dotprod(cpw(:n,iband,ikpt,ispin),cpwhlp(:,j)))
                  cpw(:n,iband,ikpt,ispin) = cpw(:n,iband,ikpt,ispin) + rdum * cpwhlp(:,j)
                enddo
              enddo
              write(6,'(A,F12.5)') '   ->',sum(abs(cpw(:n,:maxband,ikpt,:nspin2)**2)) / (maxband*nspin2)
            else
              write(6,'(15X,A,I3)') 'Number of eigenvectors reduced to',j-1
              exit
            endif
          enddo
          deallocate ( rarr,olappw,cpwhlp )
        enddo
      endif

c Shift or replace KS energies by QP energies from ENERGY file
      call getkey(inp,'ENERGY',line,default=' ')
      if(line/=' ') then
        if(line(:1)=='+'.or.line(:1)=='-') then
          write(6,'(/,A)') 'Energies are shifted by '//trim(line)
          rdum = ch2r(line,'*ENERGY')
          where(ene/=0) ene = ene + rdum
        else
          write(6,'(/,A)') 'Energies are read in from '//trim(line)
          allocate(rarr3(maxband,nkpti,nspin1))
          rarr3 = huge(0d0)
          iunit = fopen(line)
          ldum  = .false.
          read(iunit,'(A)',iostat=ind) line
          if(ind/=0)        Error('ENERGY file is empty.')
          if(line(:1)/='#') Error('ENERGY file does not start with header #...')
          do while(index(line,' Energy')==0)
            read(iunit,'(A)',iostat=ind) line
            if(ind/=0)                      Error('Reached end of ENERGY file while looking for column definition.')
            if(line/=' '.and.line(:1)/='#') Error('Found data in ENERGY file before column definition.')
          enddo
          idum = 0
          if(index(line,' b ')/=0) then ; idum = idum + 1 ; write(6,'(A)') 'File includes block information.' ; endif
          if(index(line,' s ')/=0) then ; idum = idum + 2 ; write(6,'(A)') 'File includes spin index.'        ; endif
          do
            read(iunit,'(A)',iostat=ind) line ; if(ind/=0) exit
            ind = index(line,'#')
            if(ind/=0) then
              if(index(line,'Energy')/=0) then
                ldum = index(line,' b ')/=0
                cycle
              endif
            else
              ind = len(line) + 1
            endif
            if(line(:ind-1)/=' ') then
              ispin = 1
              if     (idum==0) then ; read(line(:ind-1),*) iband,     ikpt,      rdum
              else if(idum==1) then ; read(line(:ind-1),*) iband,ikpt,ikpt,      rdum
              else if(idum==2) then ; read(line(:ind-1),*) iband,     ikpt,ispin,rdum
              else if(idum==3) then ; read(line(:ind-1),*) iband,ikpt,ikpt,ispin,rdum
              endif
              rdum = rdum / hartree
              write(6,'(2X,A,I4,A,I2,A,I4,A,2F12.6)') 'k='//Chr(ikpt)//', s='//Chr(ispin)//', n='//Chr(iband)//': '//
     &          Chr(ene(iband,ikpt,ispin))//' -> '//Chr(rdum)
              if(ispin>nspin1)            Error('ENERGY file: spin index > nspin1.')
              if(ikpt>nkpti)              Error('ENERGY file: k-point index > nkpti.')
              if(iband>nband(ikpt,ispin)) Error('ENERGY file: band index > nband().')
              if(rarr3(iband,ikpt,ispin)/=huge(0d0))
     &                                    Warn('State redefined: k = '//Chr(ikpt)//', s = '//Chr(ispin)//', n = '//Chr(iband))
              rarr3(iband,ikpt,ispin) = rdum
            endif
          enddo
          call fclose(iunit)
          call complete_ene(rarr3,1,maxband,0,1,i,j,k,2,.true.)
          where(rarr3/=huge(0d0)) ene(:,:nkpti,:) = rarr3
          do k = nkpti+1,nkpt
            ene(:,k,:) = ene(:,kptp(k),:)
          enddo
          write(6,'(2X,A)') 'Missing energies; adjusted to KS level: '//Chr(i)//' states, '//Chr(j)//' bands, '//Chr(k)//' kpoints'
          deallocate ( rarr3 )
        endif
      endif

c DELTAEX: Adjust exchange splitting
      call getkey(inp,'DELTAEX',deltaex,default=0d0,allow_eV=.true.)
      if(deltaex/=0d0) then
        if(nspin1==1) Error('DELTAX cannot be used for present system. Use for spin-polarized systems without SOC.')
        write(6,'(/A,F9.5,A)') 'Exchange splitting adjusted by', deltaex*hartree, ' eV.'
        where(ene(:,:,1)/=0) ene(:,:,1) = ene(:,:,1) - deltaex / 2
        where(ene(:,:,2)/=0) ene(:,:,2) = ene(:,:,2) + deltaex / 2
      endif

c Move coefficients, energies, and number of bands from nkpti+1:nkpti+nkpti2 to nkpt+1:nkpt+nkpti2 if lkptadd=.true.
      if(lkptadd) then
        do ikpt = nkpti2,1,-1
# ifndef LOAD
          if(.not.storeibz) then
            cpw(:,:,nkpt+ikpt,:)   = cpw(:,:,nkpti+ikpt,:)
            cmt(:,:,:,nkpt+ikpt,:) = cmt(:,:,:,nkpti+ikpt,:)
          endif
# endif
          ene(:,nkpt+ikpt,:)       = ene(:,nkpti+ikpt,:)
          nband(nkpt+ikpt,:)       = nband(nkpti+ikpt,:)
        enddo
      endif

c Distribute energies and number of bands over equivalent k points
      do ikpt = 1,nkpt2
        if(kptp(ikpt)/=ikpt) then
          ene(:,ikpt,:) = ene(:,kptp(ikpt),:)
          nband(ikpt,:) = nband(kptp(ikpt),:)
        endif
      enddo

c Refine edeg: criterion for degenerate states
      option = 104
      if(use_sym) then
        allocate(rarr(maxband*(nkpti+nkpti2)*nspin1))
        k = 0
        do ispin = 1,nspin1
          do ikpt = 1,nkpt2
            if(kptp(ikpt)==ikpt) then
              do i = 1,nband(ikpt,ispin)-1
                k       = k + 1 ; if(k>size(rarr)) Bug('Index exceeded dimension.')
                rarr(k) = ene(i+1,ikpt,ispin) - ene(i,ikpt,ispin)
              enddo
            endif
          enddo
        enddo
        call rorderf(rarr,k)
        call deg_limit(edeg,rdum1,rarr,k,.false.)
        Rwrite(6,'(/A)') 'Automatic degeneracy limit: '//Chf(edeg, 'ES8.1')
        Rwrite(6,'(A)')  '          degeneracy gap:   '//Chf(rdum1,'ES8.1')
        if( edeg<1d-12 .or. edeg>1d-6 .or. rdum1<1d3 ) then
          if(edeg<1d-12) write(6,'(A)') 'Degeneracy limit too low. Might have fallen into degenerate states.'
          if(edeg>1d-6)  write(6,'(A)') 'Degeneracy limit too high. Bad symmetry or no degeneracies.'
          if(rdum1<1d3)  write(6,'(A)') 'Degeneracy gap between degenerate and non-degenerate states too small. '//
     &                                  'Bad symmetry or no degeneracies.'
          if(rarr(1)>edeg0) then
            write(6,'(A)') 'Minimum energy difference is larger than default degeneracy limit '//Chf(edeg0,'ES8.1')
            write(6,'(A)') 'All states are assumed non-degenerate ...'
            edeg = edeg0
          else if(option==103) then
            write(6,'(A)') 'Option 103 set: Will continue with determined degeneracy limit anyway ...'
            Warn('Automatic degeneracy limit out of reasonable range (see output). Stop overridden with option -o=103.')
          else if(option==104) then
            write(6,'(A)') 'Option 104 set: Will continue with default degeneracy limit '//Chf(edeg0,'ES8.1')
            Warn('Automatic degeneracy limit out of reasonable range (see output). Stop overridden with option -o=104.')
            edeg = edeg0
          else if(option==105) then
            write(6,'(A)') 'Option 105 set: Will continue without degeneracies. Minimum energy distance: '//Chf(rarr(1),'ES8.1')            
            Warn('Automatic degeneracy limit out of reasonable range (see output). Stop overridden with option -o=105.')
            edeg = rarr(1) - 1d-12
          else
            write(6,'(A)')  'Minimum energy distance: '//Chf(rarr(1),'ES8.1')
            write(6,'(/A)') 'Cannot continue.'
            write(6,'(A)')  'Please check if the system''s symmetry is fully taken into account in the input data.' 
            write(6,'(/A)') 'To continue with this input data, use the keyword NOSYM (recommended)'
            write(6,'(A)')  'or override the stop with special option 103, 104, or 105 (e.g., -o=103):'
            write(6,'(A)')  '103 : Continue with determined degeneracy limit'
            write(6,'(A)')  '104 : Continue with default degeneracy limit '//Chf(edeg0,'ES8.1')
            write(6,'(A)')  '105 : Continue with no degeneracies.'
            write(6,'(/A)') 'Please note: This is a new consistency check.'
            write(6,'(A)')  '             Please notify the developers together with the output above.'
            Error('Automatic degeneracy limit out of reasonable range (see output).')
          endif
        endif
        deallocate(rarr)
      endif
      
c Determine degeneracies:
c   deg(i,k,s) = band index of first (last) state in a group of degenerate states if deg(i,k,s) < (>=) i
c such that
c   index of first state = i          if deg(i,k,s) >= i,     deg(i,k,s)      otherwise
c   index of last  state = deg(i,k,s) if deg(i,k,s) >= i, deg(deg(i,k,s),k,s) otherwise
      allocate ( deg(maxband,nkpt2,nspin1) )
      if(use_sym) then
        deg   =  0
        rdum1 = -huge(0d0)
        rdum2 =  huge(0d0)
        do ispin = 1,nspin1
          do ikpt = 1,nkpt2
            i = 1
            do while(i<=nband(ikpt,ispin))
              j = i
              do while(j<nband(ikpt,ispin))
                if(abs(ene(j+1,ikpt,ispin)-ene(j,ikpt,ispin))>edeg) exit
                j = j + 1
              enddo
              !if(i>1) then
              !  rdum = abs(ene(i,ikpt,ispin)-ene(i-1,ikpt,ispin))
              !  if(rdum<rdum2) write(*,'(3F20.10,I3)') rdum,ene(i,ikpt,ispin),ene(i-1,ikpt,ispin),ikpt
              !endif
              !if(j>i) then
              !  rdum = abs(ene(j,ikpt,ispin)-ene(i,ikpt,ispin))
              !  if(rdum>rdum1) write(*,'(80X,3F20.10,I3)') rdum,ene(i,ikpt,ispin),ene(j,ikpt,ispin),ikpt
              !endif
              if(i>1) rdum2       = min(rdum2,abs(ene(i,ikpt,ispin)-ene(i-1,ikpt,ispin)))
              if(j>i) rdum1       = max(rdum1,abs(ene(j,ikpt,ispin)-ene(i,ikpt,ispin)))
              deg(i,ikpt,ispin)   = j
              if(i<j) deg(i+1:j,ikpt,ispin) = i
              rdum                = sum(ene(i:j,ikpt,ispin)) / (j-i+1)
              ene(i:j,ikpt,ispin) = rdum
              i                   = j + 1
            enddo
          enddo
        enddo
        if(rdum1>=0) then
          Rwrite(6,'(/A,ES8.1)') 'Maximal energy difference in deg. subspaces:',rdum1
          Rwrite(6,'( A,ES8.1)') 'Minimal energy difference between subspaces:',rdum2
          if(abs(rdum1-rdum2)<10*edeg) Warn('Small energy distance between degenerate subspaces. Consider using NOSYM!')
        else
          Rwrite(6,'(/A,ES8.1)') 'Minimal energy difference between states:',rdum2
        endif
      else
        do i = 1,maxband
          deg(i,:,:) = i
        enddo
      endif

c Replace ene_extra by extrapolar energy
      if(ene_extra/=0) ene_extra = maxval(ene) + ene_extra

c Check job(:)%band
      do i = 1,njob
        if(.not.allocated(job(i)%band)) cycle
        do j = 1,size(job(i)%band)
          ikpt  = job(i)%kpt(j)
          ispin = job(i)%spin(j)
          nn    = nband(ikpt,ispin)
          if(job(i)%band(j)<1.or.job(i)%band(j)>nn) then
            write(0,'(A)') 'getinput: Band index out of range after keyword JOB '//trim(job(i)%label)//' :'
            if(job(i)%band(j)<1)  write(0,'(10X,I2,A)')    job(i)%band(j),' < 1'
            if(job(i)%band(j)>nn) write(0,'(10X,I3,A,I4)') job(i)%band(j),' >',nn
            Error(' ')
          endif
        enddo
      enddo

      write(6,'(/A)')        'Wavefunctions:'
      write(6,'(A,3(I5,A))') '  Number of bands:',minval(nband),' -',maxband,'  (',maxeband,' )'
      write(6,'(A,3(I5,A))') '        G vectors:',minval(ngpt), ' -',maxgpt, '  (',ngptall,' )'

c
c     Construct wavefunctions at related k points
# ifndef LOAD
      if(.not.storeibz) then
        call cpu_time(time1)
        write(6,'(/A'NoA) 'Generate wavefunctions at equivalent k points... '
        do k = 1,nkpt2
          if(kptp(k)==k) cycle
#   ifndef old_trafo
          do ispin = 1,nspin2
            call waveftrafo2(cmt(:,:,:,k,ispin),cmt(:,:,:,kptp(k),ispin),
     &                       cpw(:,  :,k,ispin),cpw(:,  :,kptp(k),ispin),nband(k,ispin),symkpt(k),kptp(k))
          enddo
#   else
          do ispin = 1,nspin2           
            call waveftrafo2(cmt(:,:,:,k,ispin),cpw(:,:,k,ispin),kptp(k),symkpt(k),ispin)
          enddo
#   endif
#     ifndef INV
          if(l_soc) then
            do i = 1,maxband
              call waveftrafo_soc(cmt(:,:,i,k,:),size(cmt,1)*size(cmt,2),symkpt(k))
              call waveftrafo_soc(cpw(:,i,k,:),  size(cpw,1)            ,symkpt(k))
            enddo
          endif
#     endif
c#   endif
        enddo
        call cpu_time(time2) ; write(6,'(A,F7.2,A)') 'done   ( Timing:',time2-time1,' )'
      endif
      ! clean up
      do ispin = 1,nspin2
        do ikpt = 1,nkpt2
          k = kindx(ikpt)
          if(ngpt(ikpt)       <maxgpt)  cpw(ngpt(ikpt)+1:,         :,k,ispin) = 0
          if(nband(ikpt,ispin)<maxband) cpw(:,  nband(ikpt,ispin)+1:,k,ispin) = 0
          if(nband(ikpt,ispin)<maxband) cmt(:,:,nband(ikpt,ispin)+1:,k,ispin) = 0
        enddo
      enddo

      write(6,'(/A,F8.1," MB")') 'Size of array cpw:',size(cpw)*MBYTES / megabyte
      write(6,'( A,F8.1," MB")') 'Size of array cmt:',size(cmt)*16d0   / megabyte
# endif

c
c     Align band structures (option ALIGNBD)
      if(alignbd) then
        if(storeibz) Error('ALIGNBD only implemented together with STOREBZ.')
        Mpi( Error('ALIGNBD not implemented for MPI version.') )
        call teststop('ALIGNBD not tested lately.')
        inquire ( file='spex.abd',exist=ldum )
        if(ldum) then ; write(6,'(/A)') 'File spex.abd exists and will be read.'
                        iunit = fopen('spex.abd',form='unformatted',status='old')
                        read(iunit) kvec
                        write(6,'(A,3F11.8)') 'K point: ',kvec
        else          ; write(6,'(/A)') 'File spex.abd does not exist and will be created.'
                        iunit = fopen('spex.abd',form='unformatted',status='new')
                        write(iunit) kptadd
        endif
        write(6,'(A'NoA) 'Align bands... '
        if(l_soc) Error('SOC not implemented for ALIGNBD.')
        call cpu_time(time1)
        allocate ( olapmt(maxlmindx,maxlmindx,ncent,nspin),done(maxband) )
        call wfolap_init_mt(olapmt,kptadd)
        do ispin = 1,nspin
          do ikpt = 1,nkpt
            allocate ( olappw(ngpt(ikpt),ngpt(ikpt+nkpt)) )
            call olap_gpt(olappw,ngpt(ikpt),ikpt,ikpt+nkpt)
            done  = .false.
            bandi = 1
            do while(bandi<=nband(ikpt,ispin))
              bandf = deg(bandi,ikpt,ispin)
              if(bandf/=bandi) then
                n = bandf-bandi+1
                allocate ( olap(n,n) )
                if(.not.ldum) then
                  allocate ( iarr(n),iarr1(n),rarr(maxband) ) ; rarr = 0
                  do iband1 = bandi,bandf
                    do iband2 = 1,nband(nkpt+ikpt,ispin)
                      if(done(iband2)) cycle
                      rdum = abs( wfolap( cmt(:,:,iband1,ikpt,ispin),     cpw(:,iband1,ikpt,ispin),
     &                                    cmt(:,:,iband2,ikpt+nkpt,ispin),cpw(:,iband2,ikpt+nkpt,ispin),ikpt,ikpt+nkpt,ispin,
     &                                    olapmt,olappw ) )
                      rarr(iband2) = rarr(iband2) + rdum**2
                    enddo
                  enddo
                  do i = 1,n
                    iarr(i) = maxloc(rarr,1)
                    rarr(iarr(i)) = 0
                  enddo
                  call rorderp(iarr1,iarr*1d0,n) ; iarr = iarr(iarr1)
                  do i = 1,n
                    iband2       = iarr(i)
                    done(iband2) = .true.
                    do j = 1,n
                      iband1 = bandi + j - 1
                      olap(j,i) = wfolap( cmt(:,:,iband1,ikpt,ispin),     cpw(:,iband1,ikpt,ispin),
     &                                    cmt(:,:,iband2,ikpt+nkpt,ispin),cpw(:,iband2,ikpt+nkpt,ispin),ikpt,ikpt+nkpt,ispin,
     &                                    olapmt,olappw )
                    enddo
                  enddo
                  deallocate ( iarr,iarr1,rarr )
                  olap = matmul(olap,invert(sqrtm(matmul(conjg(transpose(olap)),olap)))) ! make orthogonal
                  write(iunit) ispin,ikpt,bandi,bandf
                  write(iunit) olap
                else
                  read(iunit) i,j,k,l
                  if(any([i,j,k,l]/=[ispin,ikpt,bandi,bandf])) then
                    write(6,'(A)')     'getinput: Eigenspace parameters in spex.abd differ from current input.'
                    write(6,'(A,4I5)') '          Current indices:   ',ispin,ikpt,bandi,bandf
                    write(6,'(A,4I5)') '          Read from spex.abd:',i,j,k,l
                    Error('Eigenspace parameters in spex.abd differ from current input.')
                  endif
                  read(iunit) olap
                endif
                cpw(:,bandi:bandf,ikpt,ispin) = matmul(transpose(olap),cpw(:,bandi:bandf,ikpt,ispin)) ! transform cpw
                do i = 1,ncent
                  cmt(:,i,bandi:bandf,ikpt,ispin) = matmul(transpose(olap),cmt(:,i,bandi:bandf,ikpt,ispin)) ! transform cmt
                enddo
                deallocate ( olap )
              endif
              bandi = bandf + 1
            enddo
            deallocate ( olappw )
          enddo
        enddo
        deallocate ( olapmt,done )
        call fclose(iunit)
        call cpu_time(time2) ; write(6,'(A,F7.2,A)') 'done   ( Timing:',time2-time1,' )'
      endif

      endSingle

# ifdef MPI
      call Mcast(nband) ; call Mcast(ene) ; call Mcast(deltaex) ; call Mcastl(deg) ; call Mcast(ene_extra)
      call Mcast(edeg)
#   ifndef LOAD
      Nfence(cpw)
      Nfence(cmt)
      if(Nrank==0) then
        call Mcast(cpw,comm=Ocomm)
        call Mcast(cmt,comm=Ocomm)
      endif
      Nfence(cpw)
      Nfence(cmt)
#   endif
# endif

c     SOC - second variation
      Rcall getkey(inp,'PLUSSOC', ldum, default=.false.) ; Mpi( call Mcast(ldum) )
      if(ldum) call plussoc

      Rcall getkey(inp,'ORBITALS',ldum,section='WANNIER',default=.false.) ; Mpi( call Mcast(ldum) ) ! prepared for parallel Wannier90 library
      if(ldum) then
        call wannier
      else if(any(job%type==J_GT).or.any(job%type==J_SCRW).or.any(job%kern=='BSE')) then
        RError('No Wannier orbitals defined.')
      endif

      beginSingle

c     Check whether system is a metal
      allocate ( wintgr(nkpt2,maxeband,nspin1) )
      call fermi(efermi,egap,metal,wintgr(:nkpt,:,:),.true.)
      write(6,'(A,F12.8,'' Ha'')') 'Maximal energy: ',maxval(ene(:maxeband,:,:))

      if(metal) then ; eip = efermi
      else           ; eip = maxval(ene(:maxeband,:,:),ene(:maxeband,:,:)<=efermi+1d-8)
      endif

      if(lkptadd) then
        if(gauss(1)==0) then ; call tetrahedron_init(wintgr(nkpt+1:,:,:),nkpt,1,maxeband,efermi,0,.true.)
        else                 ; call gauss_init      (wintgr(nkpt+1:,:,:),nkpt,1,maxeband,efermi,0,.true.)
        endif
        write(6,'(A,F20.16)') 'Electron count error of the shifted k-point set:',abs(nelec-sum(wintgr(nkpt+1:,:,:))*2/nspin2)
      endif

c Determine "highest (lowest) occupied (unoccupied) state" -> bando,bandu
      bandu = 0
      do i = 1,maxeband
        if(maxval(wintgr(:,i,:))/=0)                                bando = i
        if(abs(1-minval(wintgr(:nkpt,i,:))*nkpt)>1d-6.and.bandu==0) bandu = i
      enddo

      call reallocate(wintgr,nkpt2,bando,nspin1)

      if(bandinfo) then
        Load( Error('BANDINFO not implemented for Spex compiled with -DLOAD.') )
        call getkey(inp,'BANDINFO',lines,status=i)
        if     (i==0) then ; Bug('BANDINFO not found.')
        else if(i==1) then ; allocate(lines(ntype)) ; lines = '[s,p,d,f,g]'
        endif
        !call band_info_old
        iunit = fopen('spex.binfo',status='unknown')
        call fclose(iunit,status='delete')
        call band_info(lines,title='input data')
        deallocate(lines)
      endif

      call getkey(inp,'DIPOLE', ldum, default=.false.)
      if(ldum) call dipole

      endSingle

c     Determine checksum
      checksum = 0
      rdum     = 0
      rdum2    = 0
      Load( beginSplit(Mrank) )
      Mpi( i = -1 )
      do ispin = 1,nspin1
        do ikpt = 1,nkpt+nkpti2
          if(ikpt>nkpti.and.ikpt<=nkpt) cycle
# ifdef LOAD
          McycleP(i)
          k = 1
          allocate( cmt(maxlmindx,ncent,nband(ikpt,ispin),1,ispin:ispin+nspin3-1) )
          allocate( cpw(maxgpt,         nband(ikpt,ispin),1,ispin:ispin+nspin3-1) )
          call read_wavef2([(iband,iband=1,nband(ikpt,ispin))],nband(ikpt,ispin),ikpt,ispin,cmt,cpw)
# else
          k = kindx(ikpt)
# endif
          do s = ispin,ispin+nspin3-1 ! sum over spins for SOC
            do iband = 1,nband(ikpt,ispin) ; NoLoad( McycleP(i) )
              icent = 0
              do itype = 1,ntype
                n = sum( [((2*l+1)*nindx(l,itype),l=0,lcut(itype))] )
                do ieq = 1,neq(itype)
                  icent    = icent + 1
                  checksum = checksum + sum( real(cmt(:n,icent,iband,k,s))    ) + sum( imag(cmt(:n,icent,iband,k,s))    )
                  rdum2    = rdum2    + sum( real(cmt(:n,icent,iband,k,s))**2 ) + sum( imag(cmt(:n,icent,iband,k,s))**2 )
                  rdum     = rdum + n
                enddo
              enddo
              checksum = checksum + sum( dble(cpw(:ngpt(ikpt),iband,k,s))    ) NoInv( + sum( imag(cpw(:ngpt(ikpt),iband,k,s))    ) )
              rdum2    = rdum2    + sum( dble(cpw(:ngpt(ikpt),iband,k,s))**2 ) NoInv( + sum( imag(cpw(:ngpt(ikpt),iband,k,s))**2 ) )
              rdum     = rdum + ngpt(ikpt)
            enddo
          enddo
          Load( deallocate(cmt,cpw) )
        enddo
      enddo
      Load( endSplit )
      Mpi( call Msum(checksum) ; call Msum(rdum2) ; call Msum(rdum) )
      checksum = checksum / sqrt(rdum)
      rdum2    = rdum2 / rdum
      Rwrite(6,'(/A,2(ES19.12,A))') 'Wavefunction checksum: ',checksum,'   (',rdum2,' )'

# ifdef MPI
      call Mcastl(wintgr) ; call Mcast(efermi)   ; call Mcast(egap)        ; call Mcast(metal)
      call Mcast(eip)     ; call Mcast(bando)    ; call Mcast(bandu)
# endif

c     Initialize FFT
      if(gcutf>0) then
        call determine_imat
        call fft_init(gcutf)
      endif

c     Define job%band for HF energy (job RPAENE)
      do ijob = 1,njob
        if(any(job(ijob)%type==[J_RPA,J_HFE])) then
          allocate ( job(ijob)%band(bando*nkpti*nspin1) )
          allocate ( job(ijob)%kpt (bando*nkpti*nspin1) )
          allocate ( job(ijob)%spin(bando*nkpti*nspin1) )
          i = 0
          do ispin = 1,nspin
            do ikpt = 1,nkpti
              do iband = 1,bando
                i                 = i + 1
                job(ijob)%band(i) = iband
                job(ijob)%kpt(i)  = ikpt
                job(ijob)%spin(i) = ispin
              enddo
            enddo
          enddo
        endif
      enddo

      return ! Irreducible representations not used at the moment.

c
c     Determine irreducible representations
# if 0
      read(*,*) i
      if(i==0) return
c      inquire ( file='spex.irrep',exist=ldum )
      ldum = .false.
      if(ldum) then
        write(6,'(/A)') 'Table of irreducible representation (spex.irrep) is read in.'
        Error('not implemented.')
      else
        write(6,'(/A)') 'Table of irreducible representations (spex.irrep) does not exist.'
        write(6,'(A)')  'It must be set up... '
        if(l_soc) Error('Calculation of irreps not implemented with SOC.')
        call cpu_time(time1)
        iunit = fopen('spex.irrep',status='unknown',form='unformatted')
        allocate ( olapmt(maxlmindx,maxlmindx,ncent,nspin) )
        call wfolap_init_mt(olapmt,[0d0,0d0,0d0])
        do ispin = 1,nspin
          do ikpt = 1,nkpti
            write(6,'(A,I3)') 'k point',ikpt
            ! Look for exact degeneracies
            bandi = 1
 11         bandf = bandi
            do
              if(bandf==nband(ikpt,ispin).or.abs(ene(min(bandf+1,nband(ikpt,ispin)),ikpt,ispin)-ene(bandi,ikpt,ispin))>edeg)
     &          exit
              bandf = bandf + 1
            enddo
            nn = bandf - bandi + 1 ; write(*,*) bandi,bandf
            allocate ( cmthlp(maxlmindx,ncent,nn),cpwhlp(ngpt(ikpt),nn),olap(nn,nn) )
            i = 0
            do isym = 1,nsym
              if(all(abs(matmul(sym(isym)%rrot,kpt(:,ikpt))-kpt(:,ikpt))<1d-10)) i = i + 1
            enddo
            write(iunit) ikpt,bandi,bandf,i

            ! irreducible representations
            allocate ( olappw(ngpt(ikpt),ngpt(ikpt)) )
            call olap_gpt(olappw,ngpt(ikpt),ikpt,ikpt)
            cdum = 0
            do isym = 1,nsym
              kvec = matmul(sym(isym)%rrot,kpt(:,ikpt))-kpt(:,ikpt)
              if(all(kvec-nint(kvec)<1d-10)) then
                do iband = bandi,bandf
                  call waveftrafo(cmthlp(:,:,iband-bandi+1),cpwhlp(:,iband-bandi+1),ikpt,isym,iband,ispin)
                enddo
                do iband2 = 1,nn
                  do iband1 = 1,nn
                    olap(iband1,iband2) = wfolap( cmt(:,:,iband2+bandi-1,ikpt,ispin),cpw(:,iband2+bandi-1,ikpt,ispin),
     &                                            cmthlp(:,:,iband1),cpwhlp(:,iband1),ikpt,ikpt,ispin,
     &                                            cdum,olappw )
                  enddo
                enddo
                olap = olap / exp( img * 2*pi * dot_product(kpt(:,ikpt),sym(isym)%transl) )
                write(iunit) isym,olap
                write(6,'(I3,2F10.5)') isym,sum( [ (olap(i,i),i=1,nn) ] )
                cdum1 = 0
                do iband1 = 1,nn
                  cdum1 = cdum1 + olap(iband1,iband1)
                enddo
                cdum = cdum + abs(cdum1)**2
                olap = matmul(conjg(transpose(olap)),olap)
                do iband = 1,nn
                  olap(iband,iband) = olap(iband,iband) - 1
                enddo
                rdum = sum(abs(olap))
                if(rdum>1d-8) then
                  Warn('Irreducible representation of bands '//Chr(bandi)//' - '//Chr(bandf)//' not orthogonal.')
                  write(0,'(A,F15.10)') '          Error: ',rdum
                endif

              endif
            enddo
            write(*,*) cdum,i ; read(*,*)

            deallocate ( cmthlp,cpwhlp,olap,olappw )
            bandi = bandf + 1
            if(bandi<=nband(ikpt,ispin)) goto 11
          enddo
        enddo
        call fclose(iunit)
        call cpu_time(time2) ; write(6,'(A,F10.5)') 'Timing: ',time2-time1
      endif
      Error(' ')
      deallocate ( olapmt )
# endif

      contains

c     ------------------

      function restart_def(string)
      implicit none
      integer                  :: restart_def
      character(*), intent(in) :: string
      character                :: ch
      integer                  :: i,n
      n           = len_trim(string) ; if(n<4.or.n>5) Error('RESTART argument: 1 or 2 or XXXX or XXXXX (X=0,1,2,3)')
      restart_def = 0
      do i = 0,n-1
        ch = string(n-i:n-i)
        select case(ch)
          case('1')    ; restart_def = restart_def + 2**(2*i)
          case('2')    ; restart_def = restart_def +            2**(2*i+1)
          case('3')    ; restart_def = restart_def + 2**(2*i) + 2**(2*i+1)
          case default ; if(ch/='0') Error('RESTART argument: expected 1 or 2 or XXXX or XXXXX with X=0,1,2,3')
        end select
      enddo
      end function restart_def

c     ------------------

c     Initialize new job with defaults or with the values of the previous job (if copy=.true.).
c     (Increases job counter j and allocation of job(:) by one.)
      subroutine job_next(job,j,copy)
      implicit none
      integer,       intent(inout) :: j
      logical,       intent(in)    :: copy
      type(jobtype), allocatable   :: job(:)
      type(jobtype), allocatable   :: job1(:)
      allocate ( job1(j) )
      j                  = j + 1
      job1               = job ; deallocate(job) ; allocate(job(j))
      job(:j-1)          = job1
      if(copy) then
        call copy_job(j-1,job(j))
      else
        job(j)%indx      = 0
        job(j)%type      = 0
        job(j)%out       = 1
        job(j)%freq      = 0
        job(j)%label     = ' '
        job(j)%freqfile  = ' '
        job(j)%full      = .false.
        job(j)%kern      = ' '
      endif
      deallocate(job1)
      end subroutine job_next

c     ------------------

      subroutine cutfreqc(j,parsed)
      implicit none
      integer,      intent(inout) :: j
      character(*), intent(in)    :: parsed
      real_dp                     :: f1,f2,f3
      integer                     :: n,nfrq,i
      if(job(j)%freqfile/=' ') Error('Subdivision of frequency file not supported.')
      n    = ch2i(parsed,'JOB '//trim(job(j)%label))
      f1   = job(j)%freq(1)
      f2   = job(j)%freq(2)
      f3   = job(j)%freq(3)
      nfrq = nint( (f2-f1) / f3) + 1      
      if(n>nfrq) then
        n = nfrq
        Warn('More subdivisions than frequencies. Number of subdivisions reduced.')
      endif
      do i = 1,n
        if(i>1) call job_next(job,j,.true.)
        job(j)%freq(1) = f1 + f3 * ( (i-1) * nfrq / n     )
        job(j)%freq(2) = f1 + f3 * (    i  * nfrq / n - 1 )
      enddo
      end subroutine cutfreqc

c     ------------------

c     Returns k point defined in string which can be a number (e.g. "1") or a k-point label (e.g. "X").
      function getkptr(string)
      implicit none
      real_dp                  :: getkptr(3)
      character(*), intent(in) :: string
      logical                  :: noninteger
      real_dp                  :: kvec(3)
      integer                  :: k
      if(noninteger(string)) then
        if(len_trim(string)>10) Error('KPT label has more than ten characters.')
        kvec = huge(0d0)
        do k = 1,nkpt_label
          if(string==ckpt_label(k)) kvec = kpt_label(:,k)
        enddo
        if(kvec(1)==huge(0d0)) then
          if(string=='+') then ; Error('no additional k point (+) defined.')
          else                 ; Error('KPT label unknown: '//trim(string))
          endif
        endif
        getkptr = kvec
      else
        read(string,*) k
        if(k<=0.or.k>nkpt) Error('k point index out of range.')
        getkptr = kpt(:,k)
      endif
      end function getkptr

c     ------------------

c     Returns index of the k point defined in string which can be a number (e.g. "1") or a k-point label (e.g. "X").
      function getkpt(string)
      implicit none
      integer                  :: getkpt
      character(*), intent(in) :: string
      logical                  :: noninteger
      real_dp                  :: kvec(3)
      integer                  :: ivec(3),k
      if(noninteger(string)) then
        if(len_trim(string)>10) Error('KPT label has more than ten characters.')
        kvec = huge(0d0)
        do k = 1,nkpt_label
          if(string==ckpt_label(k)) kvec = kpt_label(:,k)
        enddo
        if(kvec(1)==huge(0d0)) then
          if(string=='+') then ; Error('no additional k point (+) defined.')
          else                 ; Error('KPT label unknown: '//trim(string))
          endif
        else if(string=='+') then
          if(lkptadd) then
            getkpt = nkpt + 1
          else if(all(abs(modulo1r(kvec))<1d-12)) then
            getkpt = 1
          else
            Bug('lkptadd not set but kvec set and kvec /=0.')
          endif
        else
          kvec   = modulo1r(kvec)
          ivec   = nint ( kvec*nkpt3 )
          getkpt = pkpt(ivec(1),ivec(2),ivec(3))
          if(any(abs(kvec-kpt(:,getkpt))>1d-12)) Error('Requested k point not found in k-point set: ('//Chfn(kvec,',','F7.5')//')')
        endif
      else
        read(string,*) getkpt
        if(getkpt<=0.or.getkpt>nkpt) Error('k point index out of range: '//chr(getkpt))
      endif
      end function getkpt

c     ------------------

      subroutine getkptarr(kptarr,string,lwrite)
      implicit none
      integer,      allocatable :: kptarr(:)
      character(*), intent(in)  :: string
      logical,      intent(in)  :: lwrite
      integer,      allocatable :: kline(:,:),kpath(:,:)
      integer                   :: iunit_i,iunit_r
      integer                   :: i,n,ivec(3)
      logical                   :: noninteger
      if(noninteger(string)) then
        if(string=='IBZ') then
          allocate ( kptarr(nkpti) )
          kptarr = [(i,i=1,nkpti)]
          return
        else if(string=='BZ') then
          allocate ( kptarr(nkpt) )
          kptarr = [(i,i=1,nkpt)]
          return
        else if(string=='PATH') then
          if(.not.allocated(kpt_path)) Error('KPTPATH not defined.')
          if(lwrite.and.wrtext) then
            iunit_i = fopen('qpts.i',status='unknown')
            iunit_r = fopen('qpts.r',status='unknown')
          endif
          allocate ( kpath(3,1) ) ! dummy allocation
          do i = 1,nkpt_path-1
            call getpath(kline,kpt_path(:,i),kpt_path(:,i+1)) ! kline: line from i to i+1
            if(size(kline,2)>0) then
              call reallocate(kpath,3,size(kpath,2)+size(kline,2)-1)
              kpath(:,size(kpath,2)-size(kline,2)+1:) = kline ! append kline to kpath
            else
              Warn('Empty line in k-point path.')
            endif
            deallocate ( kline )
          enddo
          n = size(kpath,2)
          allocate ( kptarr(n) )
          if(lwrite.and.wrtext) then
            write(iunit_i,'(I4)') n
            write(iunit_r,'(I4)') n
          endif
          do i = 1,n
            if(lwrite.and.wrtext) write(iunit_r,'(3F16.8)') matmul(rlat,kpath(:,i)*1d0/nkpt3)*latpar/(2*pi)
            ivec      = modulo ( kpath(:,i) , nkpt3 )
            kptarr(i) = pkpt(ivec(1),ivec(2),ivec(3))
            if(lwrite.and.wrtext) write(iunit_i,'(2X,F12.10,3X,F12.10,3X,F12.10)') kpt(:,kptarr(i))
          enddo
          deallocate ( kpath )
          if(lwrite.and.wrtext) then
            call fclose(iunit_i)
            call fclose(iunit_r)
          endif
          return
        endif
        if(len_trim(string)/=1) Error('KPT label has more than one character.')
        allocate ( kptarr(1) )
        kptarr = getkpt(string)
        if(kptarr(1)==0) then
          if(string=='+') then ; Error('no additional k point (+) defined.')
          else                 ; Error('KPT label unknown.')
          endif
        endif
      else
        allocate ( kptarr(1) )
        read(string,*) kptarr
        if(kptarr(1)<=0.or.kptarr(1)>nkpt) Error('k point index out of range.')
      endif
      end subroutine getkptarr

c     ------------------

      subroutine getpath(kline,kvec1,kvec2)
      implicit none
      integer, allocatable :: kline(:,:)
      real_dp, intent(in)  :: kvec1(3),kvec2(3)
      real_dp, allocatable :: norm(:)
      integer, allocatable :: pointer(:)
      integer              :: kpt1(3),kpt2(3),i1,i2,i3,n,i
      real_dp              :: kvec(3),dkpt0(3),dkpt1(3),dist,rdum1,rdum2
      logical              :: ldef
      do i = 1,3
        kpt1(i) = floor   ( min(kvec1(i),kvec2(i)) * nkpt3(i) ) ! define parallelepiped to choose
        kpt2(i) = ceiling ( max(kvec1(i),kvec2(i)) * nkpt3(i) ) ! the k points from
      enddo
      if(any(abs(nint(kvec1*nkpt3)-kvec1*nkpt3)>1d-12)) then
        Info('Line end point not in k-point set.')
      endif
      if(any(abs(nint(kvec2*nkpt3)-kvec2*nkpt3)>1d-12)) then
        Info('Line end point not in k-point set.')
      endif
      dkpt0 = matmul ( rlat , kvec2 - kvec1 )
      dist  = sqrt(sum(dkpt0**2))
      ldef  = .false.
 1    n     = 0
      do i1 = kpt1(1),kpt2(1)
        do i2 = kpt1(2),kpt2(2)
          do i3 = kpt1(3),kpt2(3)
            kvec  = 1d0 * [ i1,i2,i3 ] / nkpt3
            dkpt1 = matmul ( rlat , kvec - kvec1 )
            rdum1 = sqrt(sum(dkpt1**2))
            rdum2 = dist - rdum1                            ! not beyond kvec2?
            rdum1 = dot_product(dkpt0,dkpt1) - dist * rdum1 ! on one line?
            if(abs(rdum1)<1d-10.and.rdum2+1d-10>0) then
              n = n + 1
              if(ldef) then
                kline(:,n) = [ i1,i2,i3 ]
              endif
            endif
          enddo
        enddo
      enddo
      if(.not.ldef) then
        ldef = .true.
        allocate ( kline(3,n) )
        goto 1
      else
        allocate ( norm(n),pointer(n) )
        do i = 1,n
          kvec    = 1d0 * kline(:,i) / nkpt3
          norm(i) = sum(matmul(rlat,kvec-kvec1)**2)
        enddo
        call rorderp(pointer,norm,n)
        kline = kline(:,pointer)
        deallocate ( norm,pointer )
      endif
      end subroutine getpath

c     ------------------

c
c Write q points to a file for band structure plot (additional q points, can be read by Fleur)
c   kpt_path(:,1:nkpt_path)
c The increment is the shortest k vector divided by inc (default INC)
# define INC 20
      subroutine write_qpts(inc)
      implicit none
      integer,       intent(in)  :: inc
      integer                    :: iunit,ikpt,nqpt,n,i
      integer                    :: kgv
      real_dp,       allocatable :: qpt(:,:)
      real_dp                    :: scale,qmin,qlen
      character(10), allocatable :: label(:)
      logical                    :: def
      iunit = fopen('qpts',status='unknown')
      scale = kgv(nkpt3,3)
      def   = .false.
      qmin  = rvol**(1d0/3)
      if(inc==0) then ; qmin = qmin / INC
      else            ; qmin = qmin / inc
      endif
      do
        nqpt = 1
        do ikpt = 1,nkpt_path-1
          qlen = sqrt(sum(matmul(rlat,kpt_path(:,ikpt+1)-kpt_path(:,ikpt))**2))
          n    = max(1,nint(qlen/qmin))
          if(def) then
            do i = 0,n ; if(i==0.and.ikpt>1) cycle
              qpt(:,nqpt+i) = ( kpt_path(:,ikpt+1) * i + kpt_path(:,ikpt) * (n-i) ) / n
              write(iunit,'(4F10.5)') qpt(:,nqpt+i)*scale,0d0
            enddo
          endif
          nqpt = nqpt + n
        enddo
        if(def) exit
        def = .true.
        write(iunit,'(I5,F20.10)') nqpt,scale
        allocate(qpt(3,nqpt))
      enddo
      call fclose(iunit)
      allocate(label(nqpt))
      label = ' '
      do i = 1,nqpt
        do j = 1,nkpt_label
          if(all(abs(qpt(:,i)-kpt_label(:,j))<1d-10)) label(i) = ckpt_label(j)
        enddo
      enddo
      call write_qpt(qpt,nqpt,label)
      deallocate(qpt,label)
      write(0,'(A)') 'File qpts written.'
      end subroutine write_qpts

c     ------------------

      function def_get_rgrid(itype,first)
      implicit none
      real_dp             :: def_get_rgrid
      real_dp, intent(in) :: first
      integer, intent(in) :: itype
      integer             :: n
      rgrid(1,itype) = first      
      do n = 2,grid(itype)%number
        rgrid(n,itype) = rgrid(n-1,itype) * exp(grid(itype)%increment)
      enddo
      def_get_rgrid = rgrid(grid(itype)%number,itype)
      end function def_get_rgrid

c     ------------------

      end

c     ------------------

c     Parses string according to syntax definition and stores the result in output,
c     e.g. string = '12:4-7,1'  with  syntax='1:3-2,4' gives output=['12','7','4','1']
c     Optional arguments can be defined by <...>,
c     e.g. syntax = '1<,2>:4-3,5'    gives output=['12',' ','7','4','1'] and
c     Binary 'or' constructs are defined by ...||...
c     e.g. syntax = '<1:||2;>:4-3,5' gives output=['12',' ','7','4','1'].
c     Nesting <...<...>...> is also allowed.
c     Take care in the definition of the syntax! Most mistakes are not detected and lead
c     to unpredictable results, e.g. '<1:2><3;4>', '1:<2<3>>;4', '<1:2>||<3:4>'.

      ! This is the parser front end.
      ! Creates from syntax a hirarchical list of low-level syntaxes without '<>' and '||' which are then interpreted together with string by the back end.
      subroutine parser(output,string,syntax,key)
      use util
      implicit none
      character(*)              :: output(*)
      character(*), intent(in)  :: string,syntax,key
      character,    allocatable :: syntax1(:)*(len(syntax))
      character                 :: syn*(len(syntax)+1),syn1*(len(syntax)+1)
      logical                   :: lor
      integer                   :: i,j,m,n,level,iadd,isyn,nsyn,ind1,ind2,ind3,ind4,error
      integer                   :: ch2i
      logical                   :: isinteger
      ! Determine index numbers in syntax and allocate output accordingly
      i = 0
      n = 0
      do while (i<len_trim(syntax))
        i = i + 1
        if(isinteger(syntax(i:i))) then
          m = ch2i(syntax(i:i),' ')
          if(i/=len_trim(syntax)) then
            if(isinteger(syntax(i+1:i+1))) then ; m = ch2i(syntax(i:i+1),' ') ; i = i + 1 ; endif
          endif
          if(m<=0) Error('Index must be positive.')
          n = max(n,m)
        endif
      enddo
      if(n==0) Error('Nothing to read.')
      allocate ( syntax1(1) )
      ! Create a list of low-level syntaxes ( -> string1(:) )
      syntax1(1) = syntax
      nsyn       = 1
      isyn       = 1
      do while ( isyn<=nsyn )
        do
          syn  = syntax1(isyn)
          ind1 = index(syn,'<')  ; if(ind1==0) ind1 = len(syn)
          ind2 = index(syn,'||') ; if(ind2==0) ind2 = len(syn)
          if(ind1<ind2) then
            ind2  = ind1 ! will contain position of '>'
            level = 1    ! current level of nesting
            do while ( level/=0 )
              ind2 = ind2 + 1 ; if(ind2>len_trim(syn)) Error('Nesting error.')
              if(syn(ind2:ind2)=='<') level = level + 1
              if(syn(ind2:ind2)=='>') level = level - 1
            enddo
            syn1 = syn(ind1+1:ind2-1) ; ind3 = ind1
            iadd = 0
            lor  = .false.
            do while(syn1/=' ')
              ind4  = 1
              level = 0
              do while ( ind4<=len_trim(syn1) )
                if(syn1(ind4:ind4)=='<')        level = level + 1
                if(syn1(ind4:ind4)=='>') then ; level = level - 1 ; if(level<0) Error('Nesting error.') ; endif
                if(syn1(ind4:ind4+1)=='||'.and.level==0) then ; lor = .true. ; exit ; endif
                ind4 = ind4 + 1
              enddo
              nsyn          = nsyn + 1 ; call reallocate(syntax1,nsyn) ; do j = nsyn,isyn+1,-1 ; syntax1(j) = syntax1(j-1) ; enddo
              syntax1(isyn) = syn(:ind1-1)//syn1(:ind4-1)//syn(ind2+1:)
              iadd          = iadd + 1
              syn1          = syn1(ind4+2:) ; ind3 = ind3 + 1
            enddo
            if(lor) then
              nsyn = nsyn - 1 ; do j = isyn+iadd,nsyn ; syntax1(j) = syntax1(j+1) ; enddo ; call reallocate(syntax1,nsyn)
            else
              syntax1(isyn+iadd) = syn(:ind1-1)//syn(ind2+1:)
            endif
          else if(ind2<ind1) then
            nsyn            = nsyn + 1 ; call reallocate(syntax1,nsyn) ; do j = nsyn,isyn+1,-1 ; syntax1(j) = syntax1(j-1) ; enddo
            syntax1(isyn)   = syn(:ind2-1)
            syntax1(isyn+1) = syn(ind2+2:)
          else
            exit
          endif
        enddo
        isyn = isyn + 1
      enddo
      ! Let parser2 interprete the low-level syntaxes and define output
      do i = 1,nsyn
        call parser2(output,string,syntax1(i),error)
        if(error==0)      then ; exit
        else if(error==2) then ; call parse_error(string,syntax,key) ; endif
      enddo
      if(error/=0) call parse_error(string,syntax,key)
      end

      ! This is the parser back end.
      ! Interpretes the low-level syntax strings.
      subroutine parser2(output,string,syntax,error)
      use util
      implicit none
      integer,      intent(out) :: error
      character(*), intent(out) :: output(*)
      character(*), intent(in)  :: string,syntax
      character,    allocatable :: sep(:)*(len(string))
      character                 :: c,c1
      integer,      allocatable :: indx(:)
      integer                   :: i,j,n,ind0,ind1
      integer                   :: ch2i
      logical                   :: noninteger
      j = 1 ; allocate ( sep(len_trim(syntax)),indx(len_trim(syntax)) ) ; sep(1) = ' ' ; indx(1) = 0
      i = 0
      do while (i<len_trim(syntax))
        i = i + 1
        c = syntax(i:i)
        if(noninteger(c)) then
          sep(j) = trim(sep(j))//c
        else
          c1 = ' '
          if(i/=len_trim(syntax)) then
            c1 = syntax(i+1:i+1) ; if(noninteger(c1)) then ; c1 = ' '
                                   else                    ; i  = i + 1 ; endif
          endif
          j               = j + 1
          indx(j)         = ch2i(c//c1,' ') ; if(indx(j)==0) Error('Zero output index. Check syntax'//char(33)) ! char(33) instead of ! to make simplify work
          sep(j)          = ' '
          output(indx(j)) = ' '
        endif
      enddo
      error  = 0
      n      = j
      ind0   = len_trim(sep(1)) + 1
      if(string(:ind0-1)/=sep(1)) error = 1
      do i = 2,n
        if(sep(i)/=' ') then
          ind1 = index(string(ind0:),trim(sep(i))) ; if(ind1==0) error = 1
          ind1 = ind1 + ind0 - 1
        else
          if(i==n) then ; output(indx(i)) = string(ind0:) ; ind0 = len_trim(string) + 1
          else            ; error = 1 ; endif
          exit
        endif
        output(indx(i)) = string(ind0:ind1-1) ; if(output(indx(i))==' ') error = 1
        ind0            = ind1 + len_trim(sep(i))
      enddo
      if(ind0/=len_trim(string)+1) error = 1
      if(error/=0) then ; do i = 2,n ; output(indx(i)) = ' ' ; enddo ; endif
      deallocate ( sep,indx )
      end

      ! Reports error and stops
      subroutine parse_error(string,syntax,key)
      implicit none
      character(*), intent(in)  :: string,syntax,key
      character                 :: c
      integer                   :: i
      logical                   :: isinteger
      write(0,'(A'NoA) 'parser: Could not parse '//trim(string) ; if(key/=' ') write(0,'(A'NoA)' after keyword '//trim(key)
      write(0,'(A/A'NoA) '.','        Expected following syntax: '
      i = 0
      do while ( i<len_trim(syntax) )
        i = i + 1
        c = syntax(i:i)
        if(isinteger(c)) then
          if(i/=len(syntax)) then ; if(isinteger(syntax(i+1:i+1))) i = i + 1 ; endif
            write(0,'(A'NoA) 'XXX'
        else
          write(0,'(A'NoA) c
        endif
      enddo
      write(0,'(A)')
      Error('Syntax error in input file.')
      end

c     ------------------

c     Calculates determinant, rotation axis, and rotation angle of symmetry operation isym.
c     Also defines SU2 rotation matrix (only used if l_soc=.true.).
      subroutine rotprop(det,angle,axis,isym)
      use global,  only: pi,lat,rlat,sym,l_soc,sqa,nspin
      use util,    only: chr
      use wrapper, only: identity
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: isym
      real_dp, intent(out) :: det,angle,axis(3)
      real_dp              :: rot(3,3),rot1(3,3)
      real_dp              :: determinant
      real_dp              :: cosa,sina,sqaxis(3)
      integer              :: i,j,k
      rot = matmul(lat,matmul(sym(isym)%rot,transpose(rlat)))/(2*pi)
      ! Calculate determinant
      det = determinant(rot)
      if(abs(1-abs(det))>1d-10) Error('determinant neither 1 nor -1')
      if( any ( abs( matmul(transpose(rot),rot)-identity(3) )>1d-10 ) ) then
        write(0,'(/A)')'rotprop: Defective symmetry transformation matrix.'
        write(0,'(A)') '         Not a proper rotation matrix:'
        write(0,'(9X,3F20.15)') rot
        Error('Defective symmetry transformation matrix.')
      endif
      rot1 = rot / det ! rot1 is proper rotation

      cosa = ( rot1(1,1)+rot1(2,2)+rot1(3,3) - 1 ) / 2
      if(abs(1-cosa)<1d-10) then
        angle   = 0
      else
        if(abs(1+cosa)<1d-10) cosa = -1
        i       = maxloc ( [ rot1(1,1),rot1(2,2),rot1(3,3) ] ,1)
        j       = mod(i  ,3) + 1
        k       = mod(i+1,3) + 1
        axis(i) = sqrt ( (rot1(i,i)-cosa) / (1-cosa) )
        sina    = (rot1(k,j)-rot1(j,k))/(2*axis(i))
        if(sina<0) then ! rotation angle defined between 0 and 180 degrees
          sina    = -sina
          axis(i) = -axis(i)
        endif
        axis(j) = (rot1(i,j)+rot1(j,i)) / (2*(1-cosa)*axis(i))
        axis(k) = (rot1(i,k)+rot1(k,i)) / (2*(1-cosa)*axis(i))
        angle   = acos(cosa)
        if(abs(1-sum(axis**2))             >1d-12) Error('Axis vector not normalized.')
        if(sum(abs(matmul(rot1,axis)-axis))>1d-12) Error('Determination of rotation axis failed.')
      endif

      if(nspin==2) then
        sqaxis = [ sin(sqa(1)) * cos(sqa(2)) , sin(sqa(1)) * sin(sqa(2)) , cos(sqa(1)) ]
        if(l_soc.and.any(abs(matmul(rot1,sqaxis)-sqaxis)>1d-12))
     &    Error('Symmetry operation '//trim(chr(isym))//' rotates spin quantization axis.')
      endif
      call rot_su2(sym(isym)%esoc,rot)

      end

c     ------------------

c     Returns SU(2) rotation matrix.
      subroutine rot_su2(su2,rot)
      use global, only: img,sqa
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: su2(2,2)
      real_dp,    intent(in)  :: rot(3,3)
      real_dp                 :: rot1(3,3),axis(3)
      real_dp                 :: det,angle,cosa,sina
      real_dp                 :: determinant
      integer                 :: i,j,k
      det = determinant(rot)
      if(abs(1-abs(det))>1d-10) Error('determinant neither 1 nor -1')
      if( any ( abs( matmul(transpose(rot),rot)-identity(3) )>1d-10 ) ) then
        write(0,'(/A)')'rot_su2: Defective symmetry transformation matrix.'
        write(0,'(A)') '         Not a proper rotation matrix:'
        write(0,'(9X,3F20.15)') rot
        Error('Defective symmetry transformation matrix.')
      endif
      if(any(sqa/=0)) then
        call rotation_sph(rot1,sqa(1),sqa(2))
        rot1 = matmul ( transpose(rot1) , matmul( rot / det , rot1) ) ! rot1 is proper rotation (+trafo to align with sqaxis)
      else
        rot1 = rot / det                                              ! rot1 is proper rotation
      endif

      cosa = ( rot1(1,1)+rot1(2,2)+rot1(3,3) - 1 ) / 2
      if(abs(1-cosa)<1d-10) then
        su2     = reshape ( [ 1,0,0,1 ] , [ 2,2 ] )
      else
        if(abs(1+cosa)<1d-10) cosa = -1
        i       = maxloc ( [ rot1(1,1),rot1(2,2),rot1(3,3) ] ,1)
        j       = mod(i  ,3) + 1
        k       = mod(i+1,3) + 1
        axis(i) = sqrt ( (rot1(i,i)-cosa) / (1-cosa) )
        sina    = (rot1(k,j)-rot1(j,k))/(2*axis(i))
        if(sina<0) then ! rotation angle defined between 0 and 180°
          sina    = -sina
          axis(i) = -axis(i)
        endif
        axis(j) = (rot1(i,j)+rot1(j,i)) / (2*(1-cosa)*axis(i))
        axis(k) = (rot1(i,k)+rot1(k,i)) / (2*(1-cosa)*axis(i))
        angle   = acos(cosa)
        if(abs(1-sum(axis**2))             >1d-12) Error('Axis vector not normalized.')
        if(sum(abs(matmul(rot1,axis)-axis))>1d-12) Error('Determination of rotation axis failed.')
        ! rotated spin eigenvectors in terms of original spin frame (as rows)
        cosa = cos(angle/2)
        sina = sin(angle/2)
        su2  = reshape ( (/ cosa - img*axis(3)*sina , (-img*axis(1)-axis(2))*sina ,
     &                  (-img*axis(1)+axis(2))*sina , cosa + img*axis(3)*sina /) , [ 2,2 ] )
      endif
      end

c     ------------------

      subroutine ibz(kpt,kptp,symkpt,nkpti,kptadd,noncoll)
      use global, only: sym,nsym,nkpt,nkpt3,modulo1r,lat,rlat,pi,trsoff,sqa
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(in)  :: kptadd(3)
      real_dp, intent(out) :: kpt(3,nkpt)
      logical, intent(in)  :: noncoll
      integer, intent(out) :: kptp(nkpt),symkpt(nkpt),nkpti
      integer              :: k1,k2,k3,k(3),isym,i1,i2
      logical              :: def(0:nkpt3(1)-1,0:nkpt3(2)-1,0:nkpt3(3)-1)
      real_dp              :: rot(3,3),sqaxis(3)
      real_dp              :: determinant
      if(noncoll) sqaxis = [ sin(sqa(1)) * cos(sqa(2)) , sin(sqa(1)) * sin(sqa(2)) , cos(sqa(1)) ]
      def = .false.
      i1  = 0
      i2  = nkpt+1
      do k1 = 0,nkpt3(1)-1
        do k2 = 0,nkpt3(2)-1
          do k3 = 0,nkpt3(3)-1
            if(.not.def(k1,k2,k3)) then
              i1            = i1 + 1
              kpt(:,i1)     = modulo1r ( dble([k1,k2,k3])/nkpt3 + kptadd )
              kptp(i1)      = i1
              symkpt(i1)    = 1
              def(k1,k2,k3) = .true.
              do isym = 1,nsym
                if(any(modulo1r(matmul(sym(isym)%rrot,kptadd)-kptadd)>1d-12)) cycle ! skip operations that do not leave kptadd invariant
                if(noncoll) then
                  if(.not.trsoff.and.isym>nsym/2) cycle
                  rot = matmul(lat,matmul(sym(isym)%rot,transpose(rlat)))/(2*pi)
                  rot = rot / determinant(rot)
                  if(any(abs(matmul(rot,sqaxis)-sqaxis)>1d-12)) cycle ! skip operations that do not leave the spin quantization axis invariant
                endif
                k = modulo ( matmul(sym(isym)%rrot,[k1,k2,k3]) , nkpt3 )
                if(.not.def(k(1),k(2),k(3))) then
                  i2                  = i2 - 1
                  kpt(:,i2)           = modulo1r ( 1d0*k/nkpt3 + kptadd )
                  kptp(i2)            = i1
                  symkpt(i2)          = isym
                  def(k(1),k(2),k(3)) = .true.
                endif
              enddo
            endif
          enddo
        enddo
      enddo
      nkpti = i1
      if(i2<nkpt) then
        kpt(:,i2:)  = kpt(:, [(i1,i1=nkpt,i2,-1)])
        kptp(i2:)   = kptp(  [(i1,i1=nkpt,i2,-1)])
        symkpt(i2:) = symkpt([(i1,i1=nkpt,i2,-1)])
      endif
      end

c     ------------------

      subroutine interpolate_radial(y,x,dy,m,grd)
      use global, only: gridtype
      use, intrinsic :: iso_fortran_env
      implicit none
      type(gridtype), intent(in)    :: grd
      integer,        intent(in)    :: m
      real_dp,        intent(in)    :: x(*),dy
      real_dp,        intent(inout) :: y(*)
      real_dp,        allocatable   :: a(:),b(:),c(:)
      real_dp                       :: xx
      integer                       :: i,n
      n = grd%number
      allocate(a(n-1),b(n-1),c(n))
      if(x(1)<1d-12) then ; call cubicspline(a,b,y,x,n,0d0,m,dy,1)
      else                ; call cubicspline(a,b,y,x,n,0d0,2,dy,1)
      endif
      xx = grd%first
      do i = 1,n
        call interpolatef(a,b,c(i),xx,y,x,n)
        xx = xx * exp(grd%increment)
      enddo
      y(:n) = c
      deallocate(a,b,c)
      end

c     ------------------

      subroutine fermi(efermi,egap,metal,wintgr,lwrite)
      use global, only: ene,nkpt,nspin1,nspin2,nelec,maxeband,l_soc,gauss,hartree
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: wintgr(nkpt,maxeband,nspin1),efermi,egap
      logical, intent(out) :: metal
      logical, intent(in)  :: lwrite
      real_dp              :: emin,emax,rdum,rdum1,rdum2
      integer              :: m1,m2,j

c     Check whether system is a metal
      wintgr = 0
      if(nspin1==1) then
        if(l_soc) then
          emin = minval(ene(nelec + 1,:nkpt,:))
          emax = maxval(ene(nelec    ,:nkpt,:))
        else
          if(mod(nelec,2)==0) then ! nelec even: nelec/2 electrons in each spin channel
            emin = minval(ene(      nelec/2 + 1,:nkpt,:))
            emax = maxval(ene(max(1,nelec/2)   ,:nkpt,:))
          else                     ! nelec odd: (nelec-1)/2 + 1/2 electrons in each spin channel
            emin = minval(ene(      nelec/2 + 1,:nkpt,:))
            emax = maxval(ene(      nelec/2 + 1,:nkpt,:))
          endif
        endif
      else
        emin = 0
        emax = 1
        do m1 = 0,nelec                                       ! m1 = Number of spin-up   electrons
          m2   = nelec - m1 ; if(max(m1,m2)+1>maxeband) cycle ! m2 = Number of spin-down electrons
          emin = min ( minval(ene(m1+1,:nkpt,1)) , minval(ene(m2+1,:nkpt,2) ) )
          emax = -huge(emax)
          if(m1>0) emax =       maxval(ene(m1,:nkpt,1))
          if(m2>0) emax = max ( maxval(ene(m2,:nkpt,2)) , emax )
          if(emin>emax) exit
        enddo
        if(emin<=emax) then ! safe start values
          emin = minval ( ene(1,       :nkpt,:) )
          emax = maxval ( ene(maxeband,:nkpt,:) )
        endif
      endif
      if(emin>emax) then
        metal  = .false.
        efermi = ( emin + emax ) / 2
        egap   = emin - emax
        if(lwrite) write(6,'(/A/A)') 'The system is not a metal.','The Fermi energy is placed in the center of the band gap.'
        if     (gauss(1)==0) then ; call tetrahedron_init(wintgr,nkpt,1,maxeband,efermi,0,.false.)
        else if(gauss(1)> 0) then ; call gauss_init      (wintgr,nkpt,1,maxeband,efermi,0,.false.)
        else                      ; Error('wrong argument to GAUSS.')
        endif
        if(any(modulo(wintgr*nkpt+1d-12,1d0)>2d-12)) then ! checks whether any of the weights are neither 0 nor 1/nkpt
          if(gauss(1)==0) then
            Bug('Found tetrahedron weights that are neither 0 nor 1/nkpt.')
          else
            Warn('Fermi edge Gauss-broadened. (Increase first GAUSS argument to avoid.)')
          endif
        endif
        if(lwrite) write(6,'(A,F12.8,'' Ha   ('',F12.8,'' eV )'')') 'Energy gap:     ',egap,egap*hartree
        if(lwrite) write(6,'(A,F12.8,'' Ha   ('',F12.8,'' eV )'')') 'Fermi energy:   ',efermi,efermi*hartree
      else
        metal  = .true. ; if(lwrite) write(6,'(/A)') 'The system is a metal.'
        egap   = 0
c       Determine Fermi energy
        ! First interval guess
        if(l_soc) then
          emin = minval(ene(max(1,nelec)   ,:nkpt,:))
          emax = maxval(ene(      nelec + 1,:nkpt,:))
        else
          emin = minval(ene(max(1,nelec/2)   ,:nkpt,:))
          emax = maxval(ene(      nelec/2 + 1,:nkpt,:))
        endif
        rdum1 = 0
        rdum2 = 0
        ! Start iterative procedure (Intervallschachtelung)
        j = 0
        do
          j = j + 1
          if(rdum1==0.or.rdum2==0) then ; efermi = ( emin + emax ) / 2
          else                          ; efermi = emin + (emax-emin) / (rdum2-rdum1) * (nelec-rdum1)
          endif
          efermi = ( emin + emax ) / 2
          if(gauss(1)==0) then     ; call tetrahedron_init(wintgr(:nkpt,:,:),nkpt,1,maxeband,efermi,0,.false.)
          else if(gauss(1)>0) then ; call gauss_init      (wintgr(:nkpt,:,:),nkpt,1,maxeband,efermi,0,.false.)
          else                     ; Error('wrong argument to GAUSS.')
          endif
          rdum = sum(wintgr) * 2/nspin2
          if(abs(rdum-nelec)<1d-14*sqrt(dble(nkpt)*nelec*nspin1)) exit ! the sqrt() takes account of rounding errors
          if(j==100) then
            write(0,'(A)') 'fermienergy: Fermi energy not converged after 100 iterations.'
            if(gauss(1)>0) write(0,'(13X,A)') 'Reduce GAUSS width!'
            Error('Fermi energy not converged after 100 iterations.')
          endif
          if(rdum>nelec) then ; emax = efermi ; rdum2 = rdum
          else                ; emin = efermi ; rdum1 = rdum
          endif
        enddo
        if(all(modulo(wintgr*nkpt+1d-12,1d0)<2d-12)) then ! checks whether any of the weights are neither 0 nor 1/nkpt
          Bug('All weights are either 0 or 1/nkpt.')
        endif
        if(lwrite) write(6,'(A,F12.8,'' Ha   ('',F12.8,'' eV, Iterations:'',I3,'' )'')') 'Fermi energy:   ',efermi,efermi*hartree,j
      endif

      if(lwrite.and.nspin1==2) write(6,'(A,F12.8,'' mu_B'')') 'Magnetic moment:',sum(wintgr(:,:,1))-sum(wintgr(:,:,2))

      end

c     ------------------

c     Complete undefined energies in ene1(n1:n2,:,:) by the ones in ene(n1:n2,:,:).
c     "Undefined" energies have the value huge(0d0).
c
c     If all energies ene1(n) with n<m (n>m) are undefined, the new energies are ene1(n) = ene(n) + ene1(m) - ene(m).
c     If all energies ene1(n) with m<n<p [and ene(m)/=ene(p)] are undefined, the new energies are
c     ene1(n) = ene(n) + [ene1(m)-ene(m)] * [ene(p)-ene(n] + [ene1(p)-ene(p)] * [ene(n)-ene(m)] / [ene(p)-ene(m)]
c
c     spin = 0 (both spins), 1 (spin up), 2 (spin down)
c     dim  = 1 (2) : ene1 is real (complex) array (imaginary part is not modified, only ordered if order>=1)
c
c     If order>=1, the energies are size-ordered.
c     If order =2, the arrays cmt/cpw are reordered accordingly.
c
c     Returns, in addition, ns (number of undefined states), nb (number of undefined bands), and
c     and nk (number of undefined kpoints).
c
      subroutine complete_ene(ene1,n1,n2,spin,dim,ns,nb,nk,order,lwrite)
      use global, only: nkpti,nspin1,nband,maxeband,ene,cmt,cpw,ncent,maxgpt,maxlmindx,edeg
      use util,   only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: n1,n2,spin,dim,order
      logical, intent(in)    :: lwrite
      integer, intent(out)   :: ns,nk,nb
      real_dp, intent(inout) :: ene1(dim,n1:n2,nkpti,*)
      real_dp                :: rdum,rdum1
      character(2)           :: sgn,sgn1
      integer, allocatable   :: pnt(:)
      integer                :: ispin,ispin1,ikpt,icent
      integer                :: i,j,n
      ispin1 = 1
      ns     = 0
      nk     = 0
      do ispin = 1,nspin1 ; if(spin==0) then ; ispin1 = ispin ; else if(ispin/=spin) then ; cycle ; endif
        do ikpt = 1,nkpti
          n  = min(n2,nband(ikpt,ispin))
          ns = ns + count(ene1(1,:n,ikpt,ispin1)==huge(0d0))
          if(all(ene1(1,:n,ikpt,ispin1)==huge(0d0))) nk = nk + 1
        enddo
      enddo
      if(lwrite.and.ns>0) write(6,'(A)') Chr(ns)//' energies are undefined. They are set to KS energies:'
      nb = 0
      do ispin = 1,nspin1 ; if(spin==0) then ; ispin1 = ispin ; else if(ispin/=spin) then ; cycle ; endif
        n = min(n2,maxeband)
        do i = n1,n
          if(all(ene1(1,i,:,ispin1)==huge(0d0))) nb = nb + 1
        enddo
      enddo
      do ispin = 1,nspin1 ; if(spin==0) then ; ispin1 = ispin ; else if(ispin/=spin) then ; cycle ; endif
        do ikpt = 1,nkpti
          n = min(n2,nband(ikpt,ispin))
          i = n1 - 1
          do while(i<n)
            j = i + 1 ; do while(ene1(1,j,ikpt,ispin1)==huge(0d0)) ; j = j + 1 ; if(j>n) exit ; enddo
            if(i==0.and.j>n) then
              ene1(1,:n,ikpt,ispin1)    = ene(n1:n,ikpt,ispin)
              if(lwrite) write(6,'(2X,A)') 'k='//Chr(ikpt)//', s='//Chr(ispin)//':  KS energies'
            else if(i==n1-1.and.j<=n.and.j>1) then
              rdum                      = ene1(1,j,ikpt,ispin1) - ene(j,ikpt,ispin)
              sgn                       = '+ ' ; if(rdum<0) sgn = '- '
              ene1(1,:j-1,ikpt,ispin1)  = ene(n1:j-1,ikpt,ispin) + rdum
              if(lwrite) write(6,'(2X,A)') 'k='//Chr(ikpt)//', s='//Chr(ispin)//', n<='//Chr(j-1)//
     &          ': KS energies '//sgn//Chr(abs(rdum))
            else if(i>n1-1.and.j>n) then
              rdum                      = ene1(1,i,ikpt,ispin1) - ene(i,ikpt,ispin)
              sgn                       = '+ ' ; if(rdum<0) sgn = '- '
              ene1(1,i+1:n,ikpt,ispin1) = ene(i+1:n,ikpt,ispin) + rdum
              if(lwrite) write(6,'(2X,A)') 'k='//Chr(ikpt)//', s='//Chr(ispin)//', n>='//Chr(i+1)//
     &          ': KS energies '//sgn//Chr(abs(rdum))
            else if(j-i>1) then
              if(abs(ene(i,ikpt,ispin)-ene(j,ikpt,ispin))<edeg) then
                rdum                        = ene1(1,i,ikpt,ispin1) - ene(i,ikpt,ispin)
                sgn                         = '+ ' ; if(rdum<0) sgn = '- '
                ene1(1,i+1:j-1,ikpt,ispin1) = ene(i+1:j-1,ikpt,ispin) + rdum
                if(lwrite) write(6,'(2X,A)') 'k='//Chr(ikpt)//', s='//Chr(ispin)//', n='//Chr(i+1)//'->'//Chr(j-1)//
     &            ': KS energies '//sgn//Chr(abs(rdum))
              else
                rdum                        = ene1(1,i,ikpt,ispin1) - ene(i,ikpt,ispin)
                rdum1                       = ene1(1,j,ikpt,ispin1) - ene(j,ikpt,ispin)
                sgn                         = '+ ' ; if(rdum<0) sgn  = '- '
                sgn1                        = '+ ' ; if(rdum<0) sgn1 = '- '
                ene1(1,i+1:j-1,ikpt,ispin1) = ene(i+1:j-1,ikpt,ispin) + (
     &            rdum  * ( ene(j,ikpt,ispin) - ene(i+1:j-1,ikpt,ispin) +
     &            rdum1 * ( ene(i+1:j-1,ikpt,ispin) - ene(i,ikpt,ispin) ) ) ) / ( ene(j,ikpt,ispin) - ene(i,ikpt,ispin) )
                if(lwrite) write(6,'(2X,A)') 'k='//Chr(ikpt)//', s='//Chr(ispin)//', n='//Chr(i+1)//'->'//Chr(j-1)//
     &            ': KS energies '//sgn//Chr(abs(rdum))//' -> '//sgn1//Chr(abs(rdum1))
              endif
            endif
            if(j>=n) exit
            i = j + 1 ; do while(ene1(1,i,ikpt,ispin1)/=huge(0d0)) ; i = i + 1 ; if(i>n) exit ; enddo ; i = i - 1
          enddo
        enddo
      enddo
      if(order==1.or.order==2) then
        do ispin = 1,nspin1 ; if(spin==0) then ; ispin1 = ispin ; else if(ispin/=spin) then ; cycle ; endif
          do ikpt = 1,nkpti
            n = min(n2,nband(ikpt,ispin))
            if(any(ene1(1,:n,ikpt,ispin1)==huge(0d0))) Bug('Could not complete energy array.')
            allocate(pnt(n1:n))
            call rorderp(pnt,ene1(1,:n,ikpt,ispin1),n-n1+1)
            pnt = pnt + n1 - 1
            if(any(pnt/=[(i,i=n1,n)]).and.maxval(abs(ene1(1,:,ikpt,ispin1)-ene1(1,pnt,ikpt,ispin)))>edeg) then
              if(lwrite) write(6,'(A)') 'States have to be re-ordered: k='//Chr(ikpt)//', s='//Chr(ispin)
              ene1(:,:,ikpt,ispin1) = ene1(:,pnt,ikpt,ispin1)
              if(order==2) then
                do i = 1,maxgpt
                  cpw(i,n1:n,ikpt,ispin) = cpw(i,pnt,ikpt,ispin)
                enddo
                do icent = 1,ncent
                  do i = 1,maxlmindx
                    cmt(i,icent,n1:n,ikpt,ispin) = cmt(i,icent,pnt,ikpt,ispin)
                  enddo
                enddo
              endif
            endif
            deallocate(pnt)
          enddo
        enddo
      else if(order/=0) then
        Bug('order parameter out of range.')
      endif
      end
