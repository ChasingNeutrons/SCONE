module fractal_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  !!
  !! Fractal surface
  !!
  !!
  !! Sample Dictionary Input:
  !!   aab { type fractal; id 92; origin (0.0 0.0 9.0); 
  !!         c (1.0 -3.0 ); it 20; tol 1E-6; limit 3; }
  !!
  !! Boundary Conditions:
  !!   Does not support boundary conditions.
  !!
  !! Private Members:
  !!   origin -> poosition of the middle of the box
  !!   halfwidth -> Halfwidths (half-length) of the box in each direction (must be > 0.0)
  !!   BC -> Boundary conditions - not supported
  !!   limit -> if exceeded, iterations terminate
  !!   tolerance -> if absolute value is less, iterations terminate
  !!   maxIt -> maximum number of iterations to termination
  !!   c -> complex parameter, determining the julia set
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: fractal
    private
    real(defReal), dimension(3)     :: origin     = ZERO
    real(defReal)                   :: limit      = 2.0_defReal
    real(defReal)                   :: tolerance  = 1.0E-3
    complex(defReal)                :: c          = ZERO
    real(defReal)                   :: scale      = ONE
    integer(shortInt)               :: maxIt      = 20
    integer(shortInt), dimension(6) :: BC = VACUUM_BC

  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
  end type fractal


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(fractal), intent(in) :: self
    character(:), allocatable  :: str

    str = 'fractal'

  end function myType

  !!
  !! Initialise box from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!
  subroutine init(self, dict)
    class(fractal), intent(inout)     :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)                        :: id, N
    real(defReal), dimension(:), allocatable :: temp
    character(100), parameter :: Here = 'init (fractal_class.f90)'

    ! Load id
    call dict % get(id,'id')
    if (id <= 0) call fatalError(Here, 'ID must be <=0. Is: '//numToChar(id))
    call self % setID(id)

    ! Load origin
    call dict % get(temp,'origin')
    N = size(temp)
    if (N /= 3) call fatalError(Here,'origin must have size 3. Has: '//numToChar(N))
    self % origin = temp

    ! Load complex parameter
    call dict % get(temp,'c')
    N = size(temp)
    if (N /= 2) call fatalError(Here,'c must have size 2. Has: '//numToChar(N))
    self % c % re = temp(1)
    self % c % im = temp(2)

    ! Load maximum number of iterations
    call dict % getOrDefault(self % maxIt, 'it', 20)
    if (self % maxIt < 1) call fatalError(Here,'it must have a value greater than 0')
    
    ! Load tolerance
    call dict % getOrDefault(self % tolerance, 'tol', 1.0E-3_defReal)
    if (self % tolerance < ZERO) call fatalError(Here,'tolerance must have a value greater than 0')
    
    ! Load overpass value - terminate iterations if this value is exceeded
    call dict % getOrDefault(self % limit, 'limit', 2.0_defReal)
    if (self % limit < ZERO) call fatalError(Here,'limit must have a value greater than zero')

    ! Load scale value to make the fractal bigger without changing its shape
    call dict % getOrDefault(self % scale, 'scale',ONE)
    if (self % scale < ZERO) call fatalError(Here,'scale must have a value greater than zero')

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(fractal), intent(in)  :: self
    real(defReal), dimension(6) :: aabb

    aabb(1:3) = self % origin - INF
    aabb(4:6) = self % origin + INF

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(f)
    class(fractal), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: rb
    real(defReal)                           :: f
    real(defReal)                           :: x, y, xtemp, r2
    integer(shortInt)                       :: i

    ! Move to origin-frame and evaluate
    rb = (r - self % origin) * self % scale

    x = rb(1)
    y = rb(2)
    r2 = self % limit * self % limit

    ! Recursively evaluate the surface equation
    i = 0
    do while(i < self % maxIt .and. x*x + y*y < r2)

      xtemp = x*x - y*y
      y = 2 * x * y + self % c % im
      x = xtemp + self % c % re

      i = i + 1

    end do

    if (i == self % maxIt) then
      f = -ONE
    else
      f = ONE
    end if
    
  end function evaluate

  !!
  !! Return distance to the surface
  !! Not defined for the fractal surface which should only be used with delta tracking.
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(fractal), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d

    d = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !! Not defined for the fractal surface which should only be used with delta tracking.
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(halfspace)
    class(fractal), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace

    halfspace = .false.

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(fractal), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin = ZERO
    self % tolerance = ZERO
    self % limit = 2.0_defReal
    self % c = ZERO
    self % maxIt = 20

  end subroutine kill

end module fractal_class
