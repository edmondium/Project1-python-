      subroutine reallocate ARGTYPE
# if   ARGTYPE == 1 || ARGTYPE ==  6 || ARGTYPE == 11 || ARGTYPE == 16
     &  (array,upper1B)
# elif ARGTYPE == 2 || ARGTYPE ==  7 || ARGTYPE == 12 || ARGTYPE == 17
     &  (array,upper1B,upper2B)
# elif ARGTYPE == 3 || ARGTYPE ==  8 || ARGTYPE == 13 || ARGTYPE == 18
     &  (array,upper1B,upper2B,upper3B)
# elif ARGTYPE == 4 || ARGTYPE ==  9 || ARGTYPE == 14 || ARGTYPE == 19
     &  (array,upper1B,upper2B,upper3B,upper4B)
# elif ARGTYPE == 5 || ARGTYPE == 10 || ARGTYPE == 15 || ARGTYPE == 20
     &  (array,upper1B,upper2B,upper3B,upper4B,upper5B)
# endif

      use, intrinsic :: iso_fortran_env
      implicit none

      integer ::
# if   ARGTYPE == 1 || ARGTYPE ==  6 || ARGTYPE == 11 || ARGTYPE == 16
     &  upper1A,upper1B,lower1
# elif ARGTYPE == 2 || ARGTYPE ==  7 || ARGTYPE == 12 || ARGTYPE == 17
     &  upper1A,upper2A,upper1B,upper2B,lower1,lower2
# elif ARGTYPE == 3 || ARGTYPE ==  8 || ARGTYPE == 13 || ARGTYPE == 18
     &  upper1A,upper2A,upper3A,upper1B,upper2B,upper3B,lower1,lower2,lower3
# elif ARGTYPE == 4 || ARGTYPE ==  9 || ARGTYPE == 14 || ARGTYPE == 19
     &  upper1A,upper2A,upper3A,upper4A,upper1B,upper2B,upper3B,upper4B,lower1,lower2,lower3,lower4
# elif ARGTYPE == 5 || ARGTYPE == 10 || ARGTYPE == 15 || ARGTYPE == 20
     &  upper1A,upper2A,upper3A,upper4A,upper5A,upper1B,upper2B,upper3B,upper4B,upper5B,lower1,lower2,lower3,lower4,lower5
# endif

# define ADD
# if   ARGTYPE <= 5
      integer,      allocatable ::
# elif ARGTYPE <= 10
      real_dp,      allocatable ::
# elif ARGTYPE <= 15
      complex_dp,   allocatable ::
# else
# undef  ADD
# define ADD *(len(array))
      character(*), allocatable ::
# endif

# if   ARGTYPE == 1 || ARGTYPE == 6 || ARGTYPE == 11 || ARGTYPE == 16
     &                        array(:),help(:)ADD
      upper1A = ubound(array,1)
      lower1  = lbound(array,1)
      allocate (help(lower1:upper1A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B))
# if   ARGTYPE <= 15
      array   = 0
# else
      array   = ' '
# endif
      upper1A = min(upper1A,upper1B)
      array(lower1:upper1A) =
     &  help(lower1:upper1A)
# elif ARGTYPE == 2 || ARGTYPE == 7 || ARGTYPE == 12 || ARGTYPE == 17
     &                        array(:,:),help(:,:)ADD
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      allocate (help(lower1:upper1A,lower2:upper2A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B))
# if   ARGTYPE <= 15
      array   = 0
# else
      array   = ' '
# endif
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      array(lower1:upper1A,lower2:upper2A) =
     &  help(lower1:upper1A,lower2:upper2A)
# elif ARGTYPE == 3 || ARGTYPE == 8 || ARGTYPE == 13 || ARGTYPE == 18
     &                        array(:,:,:),help(:,:,:)ADD
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      upper3A = ubound(array,3)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      lower3  = lbound(array,3)
      allocate (help(lower1:upper1A,lower2:upper2A,lower3:upper3A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B,lower3:upper3B))
# if   ARGTYPE <= 15
      array   = 0
# else
      array   = ' '
# endif
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      upper3A = min(upper3A,upper3B)
      array(lower1:upper1A,lower2:upper2A,lower3:upper3A) =
     &  help(lower1:upper1A,lower2:upper2A,lower3:upper3A)
# elif ARGTYPE == 4 || ARGTYPE == 9 || ARGTYPE == 14 || ARGTYPE == 19
     &                        array(:,:,:,:),help(:,:,:,:)ADD
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      upper3A = ubound(array,3)
      upper4A = ubound(array,4)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      lower3  = lbound(array,3)
      lower4  = lbound(array,4)
      allocate (help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B,lower3:upper3B,lower4:upper4B))
# if   ARGTYPE <= 15
      array   = 0
# else
      array   = ' '
# endif
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      upper3A = min(upper3A,upper3B)
      upper4A = min(upper4A,upper4B)
      array(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A) =
     &  help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A)
# elif ARGTYPE == 5 || ARGTYPE == 10 || ARGTYPE == 15 || ARGTYPE == 20
     &                        array(:,:,:,:,:),help(:,:,:,:,:)ADD
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      upper3A = ubound(array,3)
      upper4A = ubound(array,4)
      upper5A = ubound(array,5)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      lower3  = lbound(array,3)
      lower4  = lbound(array,4)
      lower5  = lbound(array,5)
      allocate (help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A,lower5:upper5A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B,lower3:upper3B,lower4:upper4B,lower5:upper5B))
# if   ARGTYPE <= 15
      array   = 0
# else
      array   = ' '
# endif
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      upper3A = min(upper3A,upper3B)
      upper4A = min(upper4A,upper4B)
      upper5A = min(upper5A,upper5B)
      array(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A,lower5:upper5B) =
     &  help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A,lower5:upper5B)
# endif
# undef ADD

      end subroutine reallocate ARGTYPE

