c     Included from key.F .
c
# if   ARGTYPE >= 5
      subroutine getkey ARGTYPE (iunit,key,arg,section,default,status,writeout)
# elif ARGTYPE >= 3
      subroutine getkey ARGTYPE (iunit,key,arg,section,default,status,writeout,mini,maxi,mine,maxe,allow_eV)
# else
      subroutine getkey ARGTYPE (iunit,key,arg,section,default,status,writeout,mini,maxi,mine,maxe)
# endif

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,      intent(in)              :: iunit
      character(*), intent(in)              :: key
      character(*), intent(in), optional    :: section
      logical,      intent(in), optional    :: writeout
      integer,      intent(out),optional    :: status
# if   ARGTYPE == 1
# define CTYPE i
      integer,      intent(out)             :: arg
      integer,      intent(in), optional    :: default
      integer,      intent(in), optional    :: mini,mine
      integer,      intent(in), optional    :: maxi,maxe
# elif ARGTYPE == 2
# define CTYPE i
      integer,                  allocatable :: arg(:) ! Ifort (V12) does not like intent(out) together with allocatable.
      integer,      intent(in), optional    :: default(:)
      integer,      intent(in), optional    :: mini,mine
      integer,      intent(in), optional    :: maxi,maxe
      integer                               :: i
# elif ARGTYPE == 3
# define CTYPE r
      real_dp,      intent(out)             :: arg
      real_dp,      intent(in), optional    :: default
      real_dp,      intent(in), optional    :: mini,mine
      real_dp,      intent(in), optional    :: maxi,maxe
      logical,      intent(in), optional    :: allow_eV
# elif ARGTYPE == 4
# define CTYPE r
      real_dp,                  allocatable :: arg(:)
      real_dp,      intent(in), optional    :: default(:)
      real_dp,      intent(in), optional    :: mini,mine
      real_dp,      intent(in), optional    :: maxi,maxe
      logical,      intent(in), optional    :: allow_eV
      integer                               :: i
# elif ARGTYPE == 5
      logical,      intent(out)             :: arg
      logical,      intent(in), optional    :: default
# elif ARGTYPE == 6
      character(*), intent(out)             :: arg
      character(*), intent(in), optional    :: default
# elif ARGTYPE == 7
      character(*),             allocatable :: arg(:)
      character(*), intent(in), optional    :: default(:)
      integer                               :: i
# endif
      integer                               :: ios,ind,ind1,narg
      integer,                  parameter   :: maxlen = 65536
      character(40)                         :: keystr
      character(maxlen)                     :: line
      logical                               :: searchkey,writeout1
# if ARGTYPE != 5
      logical                               :: eV = .false.
# endif
      
# if ARGTYPE == 3 || ARGTYPE == 4
      if(present(allow_eV)) eV = allow_eV
# endif

      writeout1 = .false.
      if(present(writeout)) writeout1 = writeout

      if(.not.present(section)) then
        if(write_section/=' ') then
          write(6,'(A)') 'END'
        endif
      endif

      if(present(section)) then
        searchkey = .false. ; keystr = trim(key)//' (section '//trim(section)//')'
      else
        searchkey = .true.  ; keystr = trim(key)
      endif

# if ARGTYPE == 2 || ARGTYPE == 4 || ARGTYPE == 7
      narg = 0
      if(allocated(arg)) then
        narg = size(arg)
        if(narg==0) deallocate ( arg )
      endif
      if(present(default).and.allocated(arg)) then
        if(size(default)/=narg) Error('Number of arguments and default values do not agree.')
      endif
# endif

# if ARGTYPE == 5
      if(present(default)) then
        arg = default
      else
        arg = .false.
      endif
# endif

      rewind(iunit)

      do
        call getline(iunit,line,ind1,ios,.false.) ! ind1 unused
        if(ios/=0.or.line=='EXIT') exit
        line = adjustl(line) ! Remove leading spaces
        if(line(:8)=='SECTION '.and.key(:min(8,len(key)))/='SECTION ') then
          if(searchkey) searchkey = .false.
          if(present(section)) then
            if(adjustl(line(8:))==section) searchkey = .true.
            if(writeout1.and.write_section/=section) then
              if(write_section/=' ') write(6,'(A)') 'END'
              write(6,'(A)') 'SECTION '//trim(section)
              write_section = trim(section)//'      '
            endif
          endif
          cycle
        endif
        if(line=='END') then
          if(searchkey) searchkey = .false.
          if(.not.present(section)) searchkey = .true.
          cycle
        endif
        if(.not.searchkey) cycle
        ind = len_trim(key) + 1
        if(line(:ind)==key) then
          ind1 = index(line(ind:),'#') + ind - 1 ! truncate comments
          if(ind1==ind-1) ind1 = len(line)       !
# if ARGTYPE != 5
          if(present(status).and.line(ind:ind1)==' ') then
            status = 1
            if(writeout1) then
              if(present(section)) write(6,'(A'NoA) '  '
              write(6,'(A)') key
            endif
            goto 2
          endif
# endif          
# if ARGTYPE == 5
          arg = .not.arg
# else
          call getval ARGTYPE (line(ind:ind1),arg,narg,eV,keystr)
# endif
# if ARGTYPE < 5
          if(present(mini)) call check   ARGTYPE (arg,mini,-1,keystr)
          if(present(maxi)) call check   ARGTYPE (arg,maxi, 1,keystr)
          if(present(mine)) call check_e ARGTYPE (arg,mine,-1,keystr)
          if(present(maxe)) call check_e ARGTYPE (arg,maxe, 1,keystr)
# endif
          goto 3 ! write out and leave
        endif
      enddo

      if(present(status)) status = 0

 2    continue

      if(.not.present(default)) then
        if(present(status)) return ! arg remains undefined on return (status=0)
# if ARGTYPE == 5
        Bug('No default given for logical key '//trim(keystr)//'.')
# endif
        Error('Keyword '//trim(keystr)//' not found.')
      endif

# if ARGTYPE == 2 || ARGTYPE == 4 || ARGTYPE == 7
      if(.not.allocated(arg)) allocate(arg(size(default)))
# endif

      arg = default
      return

 3    if(present(status)) status = 2
      if(writeout1) then
# if   ARGTYPE == 5
        if(present(default)) then
          if(arg.neqv.default) then
            if(present(section)) write(6,'(A'NoA) '  '
            write(6,'(A)') key
          endif
        else
          if(arg) then
            if(present(section)) write(6,'(A'NoA) '  '
            write(6,'(A)') key
          endif
        endif
# else
        if(present(section)) write(6,'(A'NoA) '  '
        write(6,'(A'NoA) key//'       '(:7-len_trim(key))
# endif
# if   ARGTYPE == 1 || ARGTYPE == 3
        write(6,'(''  '''NoA) ; call write_ CTYPE (arg)
# elif ARGTYPE == 2 || ARGTYPE == 4
        do i = 1,narg
          write(6,'(''  '''NoA) ; call write_ CTYPE (arg(i))
        enddo
# elif ARGTYPE == 6
        write(6,'(3X,A'NoA) trim(arg)
# elif ARGTYPE == 7
        write(6,'('' '''NoA)
        do i = 1,narg
          write(6,'(2X,A'NoA) trim(arg(i))
        enddo
# endif
# if    ARGTYPE != 5
        write(6,*)
# endif
      endif

      end subroutine getkey ARGTYPE

c -------------------------
      
# if   ARGTYPE != 5
      subroutine getval ARGTYPE (line,arg,narg,eV,keystr)
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      character(*), intent(in)               :: line,keystr
      logical,      intent(in)               :: eV
      integer                                :: narg
      integer                                :: iword,i,j,ios
      character(80)                          :: word(2000),type
# if   ARGTYPE == 1
      integer,      intent(out)              :: arg
      logical                                :: isinteger
      type = ' integer'
# elif ARGTYPE == 2
      integer,                   allocatable :: arg(:)
      logical                                :: isinteger
      type = ' integer'
# elif ARGTYPE == 3
      real_dp,      intent(out)              :: arg
      type = ' real'
# elif ARGTYPE == 4
      real_dp,                   allocatable :: arg(:)
      type = ' real'
# elif ARGTYPE == 6
      character(*), intent(out)              :: arg
      type = ' character'
# elif ARGTYPE == 7
      character(*),              allocatable :: arg(:)
      type = ' character'
# endif
      word  = ' '
      iword = 1
      do i = 1,len(line)
        if(line(i:i)/=' ') then
          j                = len_trim(word(iword)) + 1 ; if(j>len(word(1))) Error('Increase length of word(:).')
          word(iword)(j:j) = line(i:i)
        else if(word(iword)/=' ') then
          iword = iword + 1 ; if(iword>size(word)) Error('Increase dimension of word(:).')
        endif
      enddo
      if(word(iword)==' ') iword = iword - 1
      narg = iword
# if   ARGTYPE == 2 || ARGTYPE == 4 || ARGTYPE == 7
      if(allocated(arg)) then
        if(iword/=size(arg,1)) then
          if(iword>size(arg,1)) then ; word(1) = 'only'//' '//Chr(size(arg,1))
          else                       ; word(1) = Chr(size(arg,1))
          endif
          Error('Read '//Chr(iword)//' arguments after keyword '//trim(keystr)//' but needed '//trim(word(1)))
        endif
      else
        allocate(arg(iword))
      endif
      do i = 1,size(arg,1)
#   if ARGTYPE == 2
        if(.not.isinteger(word(i))) Error('Argument '//trim(word(i))//' is not an integer.')
#   endif
#   if ARGTYPE < 5
        j = max(1, len_trim(word(i)) - 1 )
        if(eV.and.word(i)(j:)=='eV') then
          word(i)(j:) = ' '
          read(word(i),*,iostat=ios) arg(i)
          arg(i) = arg(i) / hartree
        else
          read(word(i),*,iostat=ios) arg(i)
        endif
        if(ios/=0) Error('Could not read'//trim(type)//' argument from '//trim(word(i)))
#   else
        arg(i) = word(i)
#   endif
      enddo
# else
      if(iword>1) Error('Read '//Chr(iword)//' arguments after keyword '//trim(keystr)//' but needed only one')
      if(iword<1) Error('Could not read'//trim(type)//' argument after keyword '//trim(keystr))
#   if ARGTYPE == 1
      if(.not.isinteger(word(1))) Error('Argument '//trim(word(1))//' after keyword '//trim(keystr)//' is not an integer.')
#   endif
#   if ARGTYPE < 5
      j = max(1, len_trim(word(1)) - 1 )
      if(eV.and.word(1)(j:)=='eV') then
        word(1)(j:) = ' '
        read(word(1),*,iostat=ios) arg
        arg = arg / hartree
      else
        read(word(1),*,iostat=ios) arg
      endif
      if(ios/=0) Error('Could not read'//trim(type)//' argument from '//trim(word(1))//' after keyword '//trim(keystr))
#   else
      arg = word(1)
#   endif
# endif

      end subroutine getval ARGTYPE
# endif

# if   ARGTYPE < 5
      subroutine check ARGTYPE (arg,minmax,mode,keystr)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: mode
      logical             :: ldum
      character(*)        :: keystr
      integer             :: i
# if   ARGTYPE == 1
      integer             :: arg,minmax,arg1
      ldum = arg*mode>minmax*mode      ; arg1 = arg
# elif ARGTYPE == 2
      integer             :: arg(:),minmax,arg1
      ldum = any(arg*mode>minmax*mode) ; arg1 = mode*maxval(mode*arg)
# elif ARGTYPE == 3
      real_dp             :: arg,minmax,arg1
      ldum = arg*mode>minmax*mode      ; arg1 = arg
# elif ARGTYPE == 4
      real_dp             :: arg(:),minmax,arg1
      ldum = any(arg*mode>minmax*mode) ; arg1 = mode*maxval(mode*arg)
# endif
      i = (mode + 1) / 2 + 1
      if(ldum) then
        write(0,'(A'NoA) 'Value out of range after keyword '//trim(keystr)
        write(0,'(A'NoA) ': ' ; call write_ CTYPE (arg1) ; write(0,'(1X,A'NoA) '<>'(i:i)//' ' ; call write_ CTYPE (minmax)
        write(0,'(A)')   ' .'
        Error('Value out of range after keyword '//trim(keystr))
      endif
      end subroutine check ARGTYPE
# endif

# if   ARGTYPE < 5
      subroutine check_e ARGTYPE (arg,minmax,mode,keystr)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: mode
      logical             :: ldum
      character(*)        :: keystr
      integer             :: i
# if   ARGTYPE == 1
      integer             :: arg,minmax,arg1
      ldum = arg*mode>=minmax*mode      ; arg1 = arg
# elif ARGTYPE == 2
      integer             :: arg(:),minmax,arg1
      ldum = any(arg*mode>=minmax*mode) ; arg1 = mode*maxval(mode*arg)
# elif ARGTYPE == 3
      real_dp             :: arg,minmax,arg1
      ldum = arg*mode>=minmax*mode      ; arg1 = arg
# elif ARGTYPE == 4
      real_dp             :: arg(:),minmax,arg1
      ldum = any(arg*mode>=minmax*mode) ; arg1 = mode*maxval(mode*arg)
# endif
      i = mode + 2
      if(ldum) then
        write(0,'(A'NoA) 'Value out of range after keyword '//trim(keystr)
        write(0,'(A'NoA) ': ' ; call write_ CTYPE (arg1) ; write(0,'(1X,A'NoA) '<=>='(i:i+1)//' ' ; call write_ CTYPE (minmax)
        write(0,'(A)')   ' .'
        Error('Value out of range after keyword '//trim(keystr))
      endif
      end subroutine check_e ARGTYPE
# undef CTYPE
# endif
