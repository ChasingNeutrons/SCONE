module genericProcedures

  use numPrecision
  use endfConstants

  implicit none

  interface removeDuplicates
    module procedure removeDuplicates_Char
  end interface removeDuplicates

  interface linFind
    module procedure linFind_Char
  end interface

  interface findDuplicates
    module procedure findDuplicates_Char
  end interface

  interface binarySearch
    module procedure binaryFloorIdxClosed_Real
  end interface

  interface endfInterpolate
    module procedure RealReal_endf_interpolate
  end interface

  interface interpolate
    module procedure RealReal_linlin_elemental_interpolate
  end interface

  integer(shortInt), parameter :: valueOutsideArray = -1,&
                                 tooManyIter       = -2

  contains

  pure function binaryFloorIdxClosed_Real(array,value) result(idx)
    !! Performes binary search of an real sorted array and returns index of the largest element
    !! smaller-or-equal to the requested value. For the value equalt to the largest element
    !! array(size(array)) it returns size(array)-1. For the value equal to the smallest element
    !! it returns 1. It returns -ve index in case of an error. Specific value is defined as a
    !! paramether. Following errors can happen
    !!   valueOutsideArray -> larger or smaller then array bounds
    !!   tooManyIter       -> algorithm did not convarged in required number of iterations
    real(defReal),dimension(:),intent(in) :: array
    real(defReal),intent(in)              :: value
    integer(shortInt)                     :: idx
    integer(shortInt)                     :: bottom, top, i

    ! Find Top and Bottom Index Array
    bottom = 1
    top = size(array)

    ! Check if the element is in array bounds
    if ( value < array(bottom) .or. value >array(top)) then
      idx = valueOutsideArray
      return
    end if

    do i = 1,70
      !Calculate mid point
      idx = (top + bottom)/2

      ! Termination condition
      if (bottom == idx) return

      ! Binary Step
      if (array(idx) <= value ) then
        bottom = idx
      else
        top = idx
      end if
    end do

    ! Failed to end in 70 steps
    idx = tooManyIter

  end function binaryFloorIdxClosed_Real


  pure function linearFloorIdxClosed_Real(array,value) result (idx)
    !! Performes linear search of an real sorted array and returns index of the largest
    !! element smaller-or-equal to the requested value. For the value equal to the largest element
    !! array(size(array)) it returns size(array)-1. For the value equal to the smallest element
    !! it returns 1. It returns -ve index in case of an error. Specific value is defined as a
    !! paramether. Following errors can happen
    !!   valueOutsideArray -> larger or smaller then array bounds
    real(defReal),dimension(:),intent(in) :: array
    real(defReal),intent(in)              :: value
    integer(shortInt)                     :: idx
    integer(shortInt)                     :: i

    if (value > array(size(array)) .or. value < array(1)) then
      idx = valueOutsideArray
      return
    end if

    do idx = size(array)-1,1,-1
      if( array(idx) <= value) return
    end do

  end function linearFloorIdxClosed_Real


  function linearFloorIdxClosed_shortInt(array,value) result(idx)
    !! Performes linear search of an integer sorted array and returns index of the largest element,
    !! which is smaller-or-equal to the requested value. Returns errors for emelents smaller and larger
    !! than the bounds of the array. For the value equal to the smallest element it returns 1 and
    !! for the value equal to the largest element it returns an error.
    integer(shortInt),dimension(:),intent(in) :: Array
    integer(shortInt),intent(in)              :: Value
    integer(shortInt)                         :: idx
    character(100),parameter                  :: Here='linearFloorIdxClosed_shortInt (genericProcedures.f90)'

    ! Check if the value is above the bounds of an array
    if ( Value >= array(size(array))) call fatalError(Here,'Value is above upper bound of the array')

    do idx=size(array),1,-1
      if ( array(idx) <= value ) return
    end do

    call fatalError(Here,'Value is below lower bound of the array')

  end function linearFloorIdxClosed_shortInt

  pure function linearCeilingIdxOpen_shortInt(array,value) result(idx)
    !! Performes linear search of an integer sorted array and returns index of the smallest element,
    !! which is greater-or-equal to the requested value. Returns errors for elements larger than
    !! the upper bound of the array. Returns 1 for values smaller or equal to the lower bound of the
    !! array. Following errors can happen:
    !!   valueOutsideArray -> larger then the upper bound of array
    integer(shortInt),dimension(:),intent(in) :: Array
    integer(shortInt),intent(in)              :: Value
    integer(shortInt)                         :: idx
    character(100),parameter                  :: Here='linearCeilingIdxOpen_shortInt (genericProcedures.f90)'

    do idx=1,size(array)
      if ( array(idx) >= value ) return
    end do

    ! Value is larger than the upper bound of the array
    idx = valueOutsideArray

  end function linearCeilingIdxOpen_shortInt

  subroutine searchError(idx,Here)
    !! Subroutine that checks whether there was an error during search and returns approperiate
    !! message.
    integer(shortInt),intent(in)  :: idx
    character(*),intent(in)       :: Here

    if (idx < 0) then

      select case (idx)
        case (valueOutsideArray)
          call fatalError(Here,'The requested value was outide the array bounds')
        case (tooManyIter)
          call fatalError(Here,'Search did not terminate in hardcoded number of iterations')
        case default
          call fatalError(Here,'Search returned unknown error flag (negative index)')
      end select

    end if

  end subroutine

  subroutine fatalError(Where,Why)
    character(*), intent(in)    :: Why, Where
    character(100)              :: Line, locWhy, locWhere
    character(20)               :: format
    integer(shortInt)           :: i

    Line = repeat('*',100)
    format = '(A100)'
    locWhere = adjustR(where)
    locWhy = adjustR(why)

    print format, Line
    print format, 'Fatal Error has occured in:'
    print format, locWhere
    print *
    print format, 'For the following reason:'
    print format, locWhy
    print *
    print format, Line
    stop
  end subroutine fatalError

  subroutine openToRead(unitNum,File)
    integer(kind=shortInt), intent(in)    :: unitNum
    character(len=*), intent(in)          :: File
    integer(kind=shortInt)                :: errorStat
    character(len=99)                     :: errorMsg

        open ( unit   = unitNum,   &
               file   = File,      &
               status = "old",     &
               action = "read",    &
               iostat = errorStat, &
               iomsg  = errorMsg)

        !errorMsg=adjustR(errorMsg)

        if (errorStat > 0) call fatalError('openToRead subroutine (genericProcedures.f90)', &
                                           errorMsg )
  end subroutine openToRead

  function removeDuplicates_Char(charArray) result(out)
    !! Function that removes duplicates from input character array. It returns array of equal or
    !! smaller size. Unfortunatly Fortran requires output character to have specified length. Length
    !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
    !! length greater then 100
    character(len=*),dimension(:),intent(in)       :: charArray
    character(len=100),dimension(:),allocatable    :: out
    logical(kind=defBool),dimension(:),allocatable :: unique
    integer(kind=shortInt)                         :: i,j

    if (len(charArray) > len(out)) call fatalError('removeDuplicates_Char (genericProcedures.f90)',&
                                                   'Maximum length of input character is 100 ')
    if (size(charArray) == 1) then
      out = charArray
    else
      allocate(unique(size(charArray)))
      unique = .true.
      ! For every element search if it matches any previous element. Change uniqe to false if it is
      ! repeted.
        do i = 1,size(charArray)
          search: &
          do j = 1, i-1
            if( trim(charArray(i)) == trim(charArray(j)) ) then
              unique(i) = .false.
              exit search
            end if
          end do search
        end do
        ! Select elements from charArray for which unique == .true.
        out = pack(charArray, unique)
    end if
  end function

  function findDuplicates_Char(charArray) result(out)
    !! Function that finds duplicates in array of characters. Returns array that contains repeted
    !! element. Unfortunatly Fortran requires output character to have specified length. Length
    !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
    !! length greater then 100
    character(len=*),dimension(:),intent(in)       :: charArray
    character(len=100),dimension(:),allocatable    :: out
    logical(kind=defBool),dimension(:),allocatable :: unique
    integer(kind=shortInt)                         :: i,j

    if (len(charArray) > len(out)) call fatalError('removeDuplicates_Char (genericProcedures.f90)',&
                                                   'Maximum length of input character is 100 ')
    if (size(charArray) == 1) then
      out = charArray
    else
      allocate(unique(size(charArray)))
      unique = .true.
      ! For every element search if it matches any previous element. Change uniqe to false if it is
      ! repeted.
        do i = 1,size(charArray)
          search: &
          do j = 1, i-1
            if( trim(charArray(i)) == trim(charArray(j)) ) then
              unique(i) = .false.
              exit search
            end if
          end do search
        end do
        ! Select elements from charArray for which unique == .true.
        out = pack(charArray, .not.unique)
        out = removeDuplicates(out)
    end if
  end function

  function linFind_Char(charArray,target) result(index)
    !! Searches linearly for the occurance of target in charArray. Returns index of -1 if target
    !! is not found. The index assumes that array begins at 1 (i.e. charArray(1:N)). If array begins
    !! with diffrent index (i.e. A(-5:N))the returned value needs to be approperiatly translated.
    character(len=*),dimension(:),intent(in) :: charArray
    character(len=*),intent(in)              :: target
    integer(kind=shortInt)                   :: index

    do index=1,size(charArray)
      if( trim(charArray(index)) == trim(target) ) return
    end do
    index = -1
  end function

  function arrayConcat(charArray) result(out)
    !! Concatenate strings from an array into a single long character. Trims elements of char Array
    !! and ads on blank between them for separation.
    character(*),dimension(:),intent(in)       :: charArray
    character(:),allocatable                   :: out
    integer(shortInt)                          :: trimLen , i

    ! Find total trim length of elements of charArray
    trimLen=0
    do i=1,size(charArray)
      trimLen = trimLen + len(trim(charArray(i)))
    end do

    allocate(character(trimLen+size(charArray)):: out)
    out = ''
    do i=1,size(charArray)
      out = out // trim(charArray(i)) // ' '
    end do
  end function arrayConcat

  elemental function RealReal_linlin_elemental_interpolate(xMin,xMax,yMin,yMax,x) result(y)
    real(defReal), intent(in) :: xMin, xMax, yMin, yMax, x
    real(defReal)             :: y
    real(defReal)             :: interFactor

    interFactor = (x-xMin)/(xMax-xMin)
    y = yMax * interFactor + (1-interFactor)*yMin
  end function RealReal_linlin_elemental_interpolate

  function RealReal_endf_interpolate(xMin,xMax,yMin,yMax,x,endfNum) result(y)
    real(defReal), intent(in)     :: xMin, xMax, yMin, yMax, x
    integer(shortInt), intent(in) :: endfNum
    real(defReal)                 :: y
    character(100),parameter      :: Here='RealReal_endf_interpolate (genericProcedures.f90)'

    select case (endfNum) ! Naming Convention for ENDF interpolation (inY-inX) i.e. log-lin => logarithmic in y; linear in x
      case (histogramInterpolation)
        y = yMin
      case (linLinInterpolation)
        y = interpolate(xMin,xMax,yMin,yMax,x)
      case (linLogInterpolation)
        y = interpolate(log(xMin),log(xMax),yMin,yMax,log(x))
      case (logLinInterpolation)
        y = interpolate(xMin,xMax,log(yMin),log(yMax),x)
        y = exp(y)
      case (loglogInterpolation)
        y = interpolate(log(xMin),log(xMax),log(yMin),log(yMax),log(x))
        y = exp(y)
      case (chargedParticleInterpolation)
        ! Not implemented
        call fatalError(Here, 'ENDF interpolation law for charged Particles is not implemented')
      case default
        call fatalError(Here, 'Unknown ENDF interpolation number')
    end select

  end function RealReal_endf_interpolate



end module genericProcedures
