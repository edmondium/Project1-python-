c Wrapper routines for reading from and writing to HDF5 files.
c      
c hdf_fopen  : Opens HDF5 file.
c hdf_fclose : Closes HDF5 file.
c hdf_rdwr   : Reads or writes data set.
c hdf_rdwr_a : Reads or writes data as attribute.
c hdf_dim    : Queries dimensions of data set.
c
c See below for more details.

# if defined(MPI) && defined(HDF5ser)
#   warning MPI but HDF5 library only serial; this might lead to extra I/O overhead and is not well tested yet.
#   undef MPI
#   define MPI_
# endif

# include "cppmacro.h"

# ifdef MaxChunkBytes
#   define Chunk MaxChunk / (storage_size(buf)/8)      
# else
#   define Chunk MaxChunk
# endif

      module Hwrapper

      use, intrinsic :: iso_fortran_env

      integer :: Hpos ! array dimension that offset refers to (rightmost is zeroth)

      interface hdf_rdwr
      module procedure hdf_rdwr_i1,hdf_rdwr_i2,hdf_rdwr_i3,
     &                 hdf_rdwr_r1,hdf_rdwr_r2,hdf_rdwr_r3,hdf_rdwr_r4,hdf_rdwr_r5,
     &                 hdf_rdwr_c1,hdf_rdwr_c2,hdf_rdwr_c3,hdf_rdwr_c4
      end interface

      interface hdf_rdwr_a
      module procedure hdf_rdwr_a_i, hdf_rdwr_a_r, hdf_rdwr_a_str,
     &                 hdf_rdwr_a_i1,hdf_rdwr_a_r1,hdf_rdwr_a_c1
      end interface

      private :: array_divide

# define ERRSTOP(arg) if(Herr<0) Error('Fatal HDF5 error: '//trim(arg))
# define LARGEATT(name,sub) if(storage_size(buf)/8*size(buf)>64*1024) then ; if(mode>0) Info('Large attribute '//name//' (>64kB) written as dataset.') ; call sub ; return ; endif

      contains

c     --------

c     Opens a HDF5 file "Hname" and returns file handle in Hfile.
c       mode = 0 : readonly           (file opened for reading only)
c              1 : readwrite          (file opened for reading and writing; create file if it does not exist)
c              2 : truncate&readwrite (file opened for reading and writing; always create new file)

      subroutine hdf_fopen(Hfile,Hname,mode)
      use hdf5
# ifdef MPI
      use Mwrapper
      use global, only: Mrank,Mcomm
# endif
      implicit none
      character(*),     intent(in)  :: Hname
      integer,          intent(in)  :: mode
      integer(HID_T),   intent(out) :: Hfile
      integer                       :: Herr,Hmode
      logical                       :: exist,Hpar
      Mpi( integer(HID_T)           :: Hplist )
      Mpi( include 'mpif.h' )
      if(mode==0) then ; Hmode = H5F_ACC_RDONLY_F
      else             ; Hmode = H5F_ACC_RDWR_F
      endif
      ifR inquire(file=Hname,exist=exist)
      Rif(mode==0.and..not.exist) Error('File '//trim(Hname)//' not found.')
# ifdef MPI
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call Mcast(exist)
        call h5pcreate_f(H5P_FILE_ACCESS_F,Hplist,Herr)          ; ERRSTOP(Hname)
        call h5pset_fapl_mpio_f(Hplist,Mcomm,mpi_info_null,Herr) ; ERRSTOP(Hname)
        if(mode==2.or..not.exist) then
          call h5fcreate_f(Hname,H5F_ACC_TRUNC_F,Hfile,Herr,
     &                                        access_prp=Hplist) ; ERRSTOP(Hname)
        else
          call h5fopen_f(Hname,Hmode,Hfile,Herr,Hplist)          ; ERRSTOP(Hname)
        endif
        call h5pclose_f(Hplist,Herr)                             ; ERRSTOP(Hname)
        return
      endif
# endif
      if(mode==2.or..not.exist) then
        call h5fcreate_f(Hname,H5F_ACC_TRUNC_F,Hfile,Herr)       ; ERRSTOP(Hname)
      else
        call h5fopen_f(Hname,Hmode,Hfile,Herr)                   ; ERRSTOP(Hname)
      endif
      end subroutine hdf_fopen

c     --------

c     Closes HDF5 file (with handle Hfile)

      subroutine hdf_fclose(Hfile)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hfile
      integer                    :: Herr
      call h5fclose_f(Hfile,Herr) ; ERRSTOP('file close')
      end subroutine hdf_fclose

c     --------

c     Queries dimensions of data set "Hname" at Hloc

      subroutine hdf_dim(Hloc,Hname,Hdim)
      use hdf5
      implicit none
      character(*),     intent(in)  :: Hname
      integer(HID_T),   intent(in)  :: Hloc
      integer(HID_T)                :: Hset,Hspc
      integer(HSIZE_T), intent(out) :: Hdim(:)
      integer(HSIZE_T)              :: Hmaxdim(size(Hdim))
      integer                       :: Herr
      call h5dopen_f(Hloc,Hname,Hset,Herr)                     ; ERRSTOP(Hname)
      call h5dget_space_f(Hset,Hspc,Herr)                      ; ERRSTOP(Hname)
      call h5sget_simple_extent_dims_f(Hspc,Hdim,Hmaxdim,Herr) ; ERRSTOP(Hname)
      if(Herr/=size(Hdim)) Error('Wrong number of dimensions.')
      call h5sclose_f(Hspc,Herr)                               ; ERRSTOP(Hname)
      call h5dclose_f(Hset,Herr)                               ; ERRSTOP(Hname)
      end subroutine hdf_dim

c     --------

c     Read/Write data attribute "Hname" from/to Hloc
c       mode = 0 : read
c              1 : write
c              2 : write/check (write if attibute does not exist, otherwise check if it is identical to buf)

      recursive subroutine hdf_rdwr_a_i1(Hloc0,Hname0,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      logical                    :: ldum
      integer                    :: buf(:),buf1(size(buf))
      integer(HID_T)             :: Hspc,Hatt,Hloc
      integer(HSIZE_T)           :: Hdim(1)
      integer                    :: Herr
      LARGEATT(Hname0,hdf_rdwr_i1(Hloc0,Hname0,mode,buf))
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      Hdim(1) = size(buf)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)                           ; ERRSTOP(Hname)
        call h5aread_f(Hatt,H5T_NATIVE_INTEGER,buf,Hdim,Herr)          ; ERRSTOP(Hname)
      else
        if(mode==2) then
          call h5aexists_f(Hloc,Hname,ldum,Herr)                       ; ERRSTOP(Hname)
          if(ldum) then
            call hdf_rdwr_a(Hloc,Hname,0,buf1)
            if(any(buf/=buf1)) Error('Attribute(s) '//trim(Hname)//' incorrect.')
            return
          endif
        endif
        call h5screate_simple_f(1,Hdim,Hspc,Herr)                      ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,H5T_NATIVE_INTEGER,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,H5T_NATIVE_INTEGER,buf,Hdim,Herr)         ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                     ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                                       ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)                       ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_i1

      recursive subroutine hdf_rdwr_a_r1(Hloc0,Hname0,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      real_dp                    :: buf(:),buf1(size(buf))
      logical                    :: ldum
      integer(HID_T)             :: Hspc,Hatt,Hloc
      integer(HSIZE_T)           :: Hdim(1)
      integer                    :: Herr
      LARGEATT(Hname0,hdf_rdwr_r1(Hloc0,Hname0,mode,buf))
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      Hdim(1) = size(buf)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)                          ; ERRSTOP(Hname)
        call h5aread_f(Hatt,H5T_NATIVE_DOUBLE,buf,Hdim,Herr)          ; ERRSTOP(Hname)
      else
        if(mode==2) then
          call h5aexists_f(Hloc,Hname,ldum,Herr)                      ; ERRSTOP(Hname)
          if(ldum) then
            call hdf_rdwr_a(Hloc,Hname,0,buf1)
            if(any(abs(buf-buf1)>1d-10)) Error('Attribute(s) '//trim(Hname)//' incorrect.')
            return
          endif
        endif
        call h5screate_simple_f(1,Hdim,Hspc,Herr)                     ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,H5T_NATIVE_DOUBLE,buf,Hdim,Herr)         ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                    ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                                      ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)                      ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_r1

      recursive subroutine hdf_rdwr_a_c1(Hloc0,Hname0,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      complex_dp                 :: buf(:),buf1(size(buf))
      real_dp                    :: rbuf(2,size(buf))
      logical                    :: ldum
      integer(HID_T)             :: Hspc,Hatt,Hloc
      integer(HSIZE_T)           :: Hdim(2)
      integer                    :: Herr
      LARGEATT(Hname0,hdf_rdwr_c1(Hloc0,Hname0,mode,buf))
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      Hdim(1) = 2
      Hdim(2) = size(buf)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)                          ; ERRSTOP(Hname)
        call h5aread_f(Hatt,H5T_NATIVE_DOUBLE,rbuf,Hdim,Herr)         ; ERRSTOP(Hname)
        buf = rbuf(1,:) + (0d0,1d0) * rbuf(2,:)
      else
        if(mode==2) then
          call h5aexists_f(Hloc,Hname,ldum,Herr)                      ; ERRSTOP(Hname)
          if(ldum) then
            call hdf_rdwr_a(Hloc,Hname,0,buf1)
            if(any(abs(buf-buf1)>1d-10)) Error('Attribute(s) '//trim(Hname)//' incorrect.')
            return
          endif
        endif
        rbuf(1,:) = real(buf)
        rbuf(2,:) = imag(buf)
        call h5screate_simple_f(2,Hdim,Hspc,Herr)                     ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,H5T_NATIVE_DOUBLE,rbuf,Hdim,Herr)        ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                    ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                                      ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)                      ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_c1

      subroutine hdf_rdwr_a_i(Hloc,Hname,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname
      integer                    :: buf
      integer                    :: buf1(1)
      if(mode/=0) buf1(1) = buf
      call hdf_rdwr_a_i1(Hloc,Hname,mode,buf1)
      if(mode==0) buf = buf1(1)
      end subroutine hdf_rdwr_a_i

      subroutine hdf_rdwr_a_r(Hloc,Hname,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname
      real_dp                    :: buf
      real_dp                    :: buf1(1)
      if(mode/=0) buf1(1) = buf
      call hdf_rdwr_a_r1(Hloc,Hname,mode,buf1)
      if(mode==0) buf = buf1(1)
      end subroutine hdf_rdwr_a_r

      subroutine hdf_rdwr_a_str(Hloc0,Hname0,mode,str)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      character(*)               :: str
      integer(HID_T)             :: Hspc,Hatt,Htyp,Hloc
      integer(HSIZE_T)           :: Hdim(1)
      integer(SIZE_T)            :: Hlen
      integer                    :: Herr
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      if(mode==2) Error('mode==2 not implemented for this routine.')
      Hlen    = len(str)
      Hdim(1) = 1
      call h5tcopy_f(H5T_NATIVE_CHARACTER,Htyp,Herr)     ; ERRSTOP(Hname)
      call h5tset_size_f(Htyp,Hlen,Herr)                 ; ERRSTOP(Hname)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)             ; ERRSTOP(Hname)
        call h5aread_f(Hatt,Htyp,str,Hdim,Herr)          ; ERRSTOP(Hname)
      else
        call h5screate_simple_f(1,Hdim,Hspc,Herr)        ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,Htyp,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,Htyp,str,Hdim,Herr)         ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                       ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                         ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)         ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_str

c     --------

c     Read/Write data set "Hname" from/to Hloc
c
c     mode = 0 : read
c            1 : write
c            2 : write/check (write data set if it does not exist, otherwise just return)
c
c     off = offset for dimension Hpos (see Hpos above)
c     str = stride for dimension Hpos      

# define call_hdf_rdwr(arg1,arg2,arg3,arg4) if(present(off).and.present(str)) then ; call hdf_rdwr(arg1,arg2,arg3,arg4,off,str) ; else if(present(off)) then ; call hdf_rdwr(arg1,arg2,arg3,arg4,off) ; else if(present(str)) then ; call hdf_rdwr(arg1,arg2,arg3,arg4,str) ; else ; call hdf_rdwr(arg1,arg2,arg3,arg4) ; endif

      subroutine hdf_rdwr_i1(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      integer                       :: buf(:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(1),Hoff(1),Hstr(1)
      integer                       :: Herr
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                 ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                      ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        Hoff(1) = off
        Hstr(1) = 1 ; if(present(str)) Hstr(1) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                             ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                              ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)            ; ERRSTOP(Hname)
        call h5screate_simple_f(1,Hdim,Hmem,Herr)                                        ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) )  ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                       ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                       ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                           ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(1,Hdim,Hspc,Herr)                                      ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_INTEGER,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                     ; ERRSTOP(Hname)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                              ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                         ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_i1

c     --------

      recursive subroutine hdf_rdwr_i2(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      integer                       :: buf(:,:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(2),Hoff(2),Hstr(2)
      integer                       :: Herr
      integer                       :: lb(2),ub(2),d,i
      if(size(buf,kind=int_dp)>Chunk) then
        Info('Dataset '//trim(Hname)//' divided as its size exceeds an internal HDF5 limit.')
        if(any([off,str]/=[0,1])) Error('off=0 or str=1. Dataset division not implemented for this case.')
        lb = lbound(buf)
        ub = ubound(buf)
        call array_divide(d,i,lb,ub)
        ub(d) = i
        call_hdf_rdwr(Hloc,trim(Hname)//'_0',mode,buf(lb(1):ub(1),lb(2):ub(2)))
        ub(d) = ubound(buf,d)
        lb(d) = i+1
        call_hdf_rdwr(Hloc,trim(Hname)//'_1',mode,buf(lb(1):ub(1),lb(2):ub(2)))
        return
      endif
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf,1)
      Hdim(2) = size(buf,2)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                 ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                      ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        if(Hpos<0.or.Hpos>1) Bug('Wrong position.')
        Hoff = 0 ;                  Hoff(2-Hpos) = off
        Hstr = 1 ; if(present(str)) Hstr(2-Hpos) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                             ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                              ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)            ; ERRSTOP(Hname)
        call h5screate_simple_f(2,Hdim,Hmem,Herr)                                        ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) )  ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                       ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                       ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                           ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(2,Hdim,Hspc,Herr)                                      ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_INTEGER,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                     ; ERRSTOP(Hname)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                              ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                         ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_i2

c     --------

      recursive subroutine hdf_rdwr_i3(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      integer                       :: buf(:,:,:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(3),Hoff(3),Hstr(3)
      integer                       :: Herr
      integer                       :: lb(3),ub(3),d,i
      if(size(buf,kind=int_dp)>Chunk) then
        Info('Dataset '//trim(Hname)//' divided as its size exceeds an internal HDF5 limit.')
        if(any([off,str]/=[0,1])) Error('off=0 or str=1. Dataset division not implemented for this case.')
        lb = lbound(buf)
        ub = ubound(buf)
        call array_divide(d,i,lb,ub)
        ub(d) = i
        call_hdf_rdwr(Hloc,trim(Hname)//'_0',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        ub(d) = ubound(buf,d)
        lb(d) = i+1
        call_hdf_rdwr(Hloc,trim(Hname)//'_1',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        return
      endif
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf,1)
      Hdim(2) = size(buf,2)
      Hdim(3) = size(buf,3)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                 ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                      ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        if(Hpos<0.or.Hpos>1) Bug('Wrong position.')
        Hoff = 0 ;                  Hoff(3-Hpos) = off
        Hstr = 1 ; if(present(str)) Hstr(3-Hpos) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                             ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                              ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)            ; ERRSTOP(Hname)
        call h5screate_simple_f(3,Hdim,Hmem,Herr)                                        ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) )  ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                       ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                       ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                           ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(3,Hdim,Hspc,Herr)                                      ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_INTEGER,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_INTEGER,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                     ; ERRSTOP(Hname)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                              ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                         ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_i3

c     --------

      subroutine hdf_rdwr_r1(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      real_dp                       :: buf(:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(1),Hoff(1),Hstr(1)
      integer                       :: Herr
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                     ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        Hoff(1) = off
        Hstr(1) = 1 ; if(present(str)) Hstr(1) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                            ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                             ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)           ; ERRSTOP(Hname)
        call h5screate_simple_f(1,Hdim,Hmem,Herr)                                       ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) )  ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                      ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                      ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                          ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(1,Hdim,Hspc,Herr)                                     ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                    ; ERRSTOP(Hname)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                             ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                        ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_r1

c     --------

      recursive subroutine hdf_rdwr_r2(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      real_dp                       :: buf(:,:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(2),Hoff(2),Hstr(2)
      integer                       :: Herr
      integer                       :: lb(2),ub(2),d,i
      if(size(buf,kind=int_dp)>Chunk) then
        Info('Dataset '//trim(Hname)//' divided as its size exceeds an internal HDF5 limit.')
        if(any([off,str]/=[0,1])) Error('off=0 or str=1. Dataset division not implemented for this case.')
        lb = lbound(buf)
        ub = ubound(buf)
        call array_divide(d,i,lb,ub)
        ub(d) = i
        call_hdf_rdwr(Hloc,trim(Hname)//'_0',mode,buf(lb(1):ub(1),lb(2):ub(2)))
        ub(d) = ubound(buf,d)
        lb(d) = i+1
        call_hdf_rdwr(Hloc,trim(Hname)//'_1',mode,buf(lb(1):ub(1),lb(2):ub(2)))
        return
      endif
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf,1)
      Hdim(2) = size(buf,2)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                     ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        if(Hpos<0.or.Hpos>1) Bug('Wrong position.')
        Hoff = 0 ;                  Hoff(2-Hpos) = off
        Hstr = 1 ; if(present(str)) Hstr(2-Hpos) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                            ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                             ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)           ; ERRSTOP(Hname)
        call h5screate_simple_f(2,Hdim,Hmem,Herr)                                       ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f (Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                      ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                      ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                          ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(2,Hdim,Hspc,Herr)                                     ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                    ; ERRSTOP(Hname)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                             ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                        ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_r2

c     --------

      recursive subroutine hdf_rdwr_r3(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      real_dp                       :: buf(:,:,:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(3),Hoff(3),Hstr(3)
      integer                       :: Herr
      integer                       :: lb(3),ub(3),d,i
      if(size(buf,kind=int_dp)>Chunk) then
        Info('Dataset '//trim(Hname)//' divided as its size exceeds an internal HDF5 limit.')
        if(any([off,str]/=[0,1])) Error('off=0 or str=1. Dataset division not implemented for this case.')
        lb = lbound(buf)
        ub = ubound(buf)
        call array_divide(d,i,lb,ub)
        ub(d) = i
        call_hdf_rdwr(Hloc,trim(Hname)//'_0',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        ub(d) = ubound(buf,d)
        lb(d) = i+1
        call_hdf_rdwr(Hloc,trim(Hname)//'_1',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        return
      endif
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf,1)
      Hdim(2) = size(buf,2)
      Hdim(3) = size(buf,3)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                     ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        if(Hpos<0.or.Hpos>2) Bug('Wrong position.')
        Hoff = 0 ;                  Hoff(3-Hpos) = off
        Hstr = 1 ; if(present(str)) Hstr(3-Hpos) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                            ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                             ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)           ; ERRSTOP(Hname)
        call h5screate_simple_f(3,Hdim,Hmem,Herr)                                       ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f (Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                      ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                      ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                          ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(3,Hdim,Hspc,Herr)                                     ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                             ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                        ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_r3

c     --------

      recursive subroutine hdf_rdwr_r4(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      real_dp                       :: buf(:,:,:,:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(4),Hoff(4),Hstr(4)
      integer                       :: Herr
      integer                       :: lb(4),ub(4),d,i
      if(size(buf,kind=int_dp)>Chunk) then
        Info('Dataset '//trim(Hname)//' divided as its size exceeds an internal HDF5 limit.')
        if(any([off,str]/=[0,1])) Error('off=0 or str=1. Dataset division not implemented for this case.')
        lb = lbound(buf)
        ub = ubound(buf)
        call array_divide(d,i,lb,ub)
        ub(d) = i
        call_hdf_rdwr(Hloc,trim(Hname)//'_0',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)))
        ub(d) = ubound(buf,d)
        lb(d) = i+1
        call_hdf_rdwr(Hloc,trim(Hname)//'_1',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)))
        return
      endif
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf,1)
      Hdim(2) = size(buf,2)
      Hdim(3) = size(buf,3)
      Hdim(4) = size(buf,4)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                     ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        if(Hpos<0.or.Hpos>3) Bug('Wrong position.')
        Hoff = 0 ;                  Hoff(4-Hpos) = off
        Hstr = 1 ; if(present(str)) Hstr(4-Hpos) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                            ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                             ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)           ; ERRSTOP(Hname)
        call h5screate_simple_f(4,Hdim,Hmem,Herr)                                       ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f (Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                      ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                      ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                          ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(4,Hdim,Hspc,Herr)                                     ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                    ; ERRSTOP(Hname)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                             ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                        ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_r4

c     --------

      recursive subroutine hdf_rdwr_r5(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      real_dp                       :: buf(:,:,:,:,:)
      logical                       :: ldum,Hpar
      integer, optional, intent(in) :: off,str
      integer(HID_T)                :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)              :: Hdim(5),Hoff(5),Hstr(5)
      integer                       :: Herr
      integer                       :: lb(5),ub(5),d,i
      if(size(buf,kind=int_dp)>Chunk) then
        Info('Dataset '//trim(Hname)//' divided as its size exceeds an internal HDF5 limit.')
        if(any([off,str]/=[0,1])) Error('off=0 or str=1. Dataset division not implemented for this case.')
        lb = lbound(buf)
        ub = ubound(buf)
        call array_divide(d,i,lb,ub)
        ub(d) = i
        call_hdf_rdwr(Hloc,trim(Hname)//'_0',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)))
        ub(d) = ubound(buf,d)
        lb(d) = i+1
        call_hdf_rdwr(Hloc,trim(Hname)//'_1',mode,buf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)))
        return
      endif
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim(1) = size(buf,1)
      Hdim(2) = size(buf,2)
      Hdim(3) = size(buf,3)
      Hdim(4) = size(buf,4)
      Hdim(5) = size(buf,5)
# ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                ; ERRSTOP(Hname)
        call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                     ; ERRSTOP(Hname)
      endif
# endif
      if(present(off)) then
        if(Hpos<0.or.Hpos>4) Bug('Wrong position.')
        Hoff = 0 ;                  Hoff(5-Hpos) = off
        Hstr = 1 ; if(present(str)) Hstr(5-Hpos) = str
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                            ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                             ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)           ; ERRSTOP(Hname)
        call h5screate_simple_f(5,Hdim,Hmem,Herr)                                       ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f (Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) ) ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                      ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                      ; ERRSTOP(Hname)
      else
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                          ; ERRSTOP(Hname)
          call h5dread_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )   ; ERRSTOP(Hname)
        else
          call h5screate_simple_f(5,Hdim,Hspc,Herr)                                     ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hset,Herr)                 ; ERRSTOP(Hname)
          call h5dwrite_f(Hset,H5T_NATIVE_DOUBLE,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                    ; ERRSTOP(Hname)
        endif
      endif
# ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                             ; ERRSTOP(Hname)
# endif
      call h5dclose_f(Hset,Herr)                                                        ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_r5

c     --------

      subroutine hdf_rdwr_c1(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      complex_dp                    :: buf(:)
      integer, optional, intent(in) :: off,str
      if(present(off)) then
        if(present(str)) then
          call real_hdf_rdwr_c1(Hloc,Hname,mode,buf,size(buf),off,str)
        else
          call real_hdf_rdwr_c1(Hloc,Hname,mode,buf,size(buf),off,-1)
        endif
      else
        call real_hdf_rdwr_c1(Hloc,Hname,mode,buf,size(buf),-1,-1)
      endif
      end subroutine hdf_rdwr_c1

c     --------

      subroutine hdf_rdwr_c2(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      complex_dp                    :: buf(:,:)
      integer, optional, intent(in) :: off,str
      if(present(off)) then
        if(present(str)) then
          call real_hdf_rdwr_c2(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),off,str)
        else
          call real_hdf_rdwr_c2(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),off,-1)
        endif
      else
        call real_hdf_rdwr_c2(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),-1,-1)
      endif
      end subroutine hdf_rdwr_c2

c     --------

      subroutine hdf_rdwr_c3(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      complex_dp                    :: buf(:,:,:)
      integer, optional, intent(in) :: off,str
      if(present(off)) then
        if(present(str)) then
          call real_hdf_rdwr_c3(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),size(buf,3),off,str)
        else
          call real_hdf_rdwr_c3(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),size(buf,3),off,-1)
        endif
      else
        call real_hdf_rdwr_c3(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),size(buf,3),-1,-1)
      endif
      end subroutine hdf_rdwr_c3

c     --------

      subroutine hdf_rdwr_c4(Hloc,Hname,mode,buf,off,str)
      use hdf5
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      complex_dp                    :: buf(:,:,:,:)
      integer, optional, intent(in) :: off,str
      if(present(off)) then
        if(present(str)) then
          call real_hdf_rdwr_c4(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),size(buf,3),size(buf,4),off,str)
        else
          call real_hdf_rdwr_c4(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),size(buf,3),size(buf,4),off,-1)
        endif
      else
        call real_hdf_rdwr_c4(Hloc,Hname,mode,buf,size(buf,1),size(buf,2),size(buf,3),size(buf,4),-1,-1)
      endif
      end subroutine hdf_rdwr_c4

c     --------

      subroutine hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      use hdf5
      implicit none
      integer(HID_T), intent(out) :: Hloc
      integer(HID_T), intent(in)  :: Hloc0
      character(*),   intent(out) :: Hname
      character(*),   intent(in)  :: Hname0
      integer                     :: Herr,ind
      ind = index(Hname0,'/',back=.true.)
      if(ind/=0) then
        if(Hname0(:1)=='/') Bug('Leading "/".')
        call h5gopen_f(Hloc0,Hname0(:ind-1),Hloc,Herr) ; ERRSTOP(Hname)
        Hname = Hname0(ind+1:)
      else
        Hname = Hname0
        Hloc  = Hloc0
      endif
      end subroutine hdf_gopen

c     --------

      subroutine array_divide(d,i,lb,ub)
      implicit none
      integer, intent(out) :: d,i
      integer, intent(in)  :: lb(:),ub(:)
      if(size(lb)/=size(ub)) Bug('Different sizes of lb and ub.')
      do d = size(lb),1,-1
        if(ub(d)-lb(d)>0) then
          i = ( ub(d) - lb(d) ) / 2 + lb(d)
          return
        endif
      enddo
      Bug('Subroutine array_divide failed. Array cannot be divided.')
      end subroutine array_divide

c     --------      

      end module Hwrapper

c     --------

c     Real calls for complex datatype

c     --------

      subroutine real_hdf_rdwr_c1(Hloc,Hname,mode,buf,b1,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b1
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b1)
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

c     --------

      subroutine real_hdf_rdwr_c2(Hloc,Hname,mode,buf,b1,b2,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b1,b2
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b1,b2)
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

c     --------

      subroutine real_hdf_rdwr_c3(Hloc,Hname,mode,buf,b1,b2,b3,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b1,b2,b3
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b1,b2,b3)
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

c     --------

      subroutine real_hdf_rdwr_c4(Hloc,Hname,mode,buf,b1,b2,b3,b4,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b1,b2,b3,b4
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b1,b2,b3,b4)
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

c     --------

# undef ERRSTOP

# if defined(MPI_) && defined(HDF5ser)
#   undef MPI_
#   define MPI
# endif
