module cellUniverse_class

  use numPrecision
  use universalVariables, only : UNDEF_MAT, NUDGE
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill
  implicit none
  private

  !!
  !! Local helper class to group cell data
  !!
  !! Public Members:
  !!   idx -> cellIdx of the cell in cellShelf
  !!   ptr -> Pointer to the cell
  !!
  type, private :: localCell
    integer(shortInt)    :: idx = 0
    class(cell), pointer :: ptr => null()
  end type localCell

  !!
  !! Representation of a universe via cells
  !!
  !! Each local cell in the universe corespondes to a cell given by an ID.
  !! An extra local cell is always defined inside the cellUniverse with UNDEF_MAT
  !! (undefined material) filling. If position is not in any user-defined cell, it is in this
  !! extra cell. Extra cell exists to enable plotting of geometry without fatalErrors.
  !!
  !! After the elements of a universe have been visited 5000 times, the search order of the cells
  !! is changed to be more optimal according to the frequency of visits. 
  !!
  !! Sample Input Dictionary:
  !!   uni { type cellUniverse;
  !!         id 7;
  !!         # origin (2.0 0.0 0.0);    #
  !!         # rotation (23.0 0.0 0.0); #
  !!         cells ( 1 3 4);         }
  !!
  !! Note:
  !!   - Local IDs are assigned in order as in definition. In the example above local id would map
  !!   to following cell  [localID: cellID] {1: 1, 2: 3, 3: 4, 4: UNDEF }
  !!   - Cell overlaps are forbidden, but there is no check to find overlaps.
  !!
  !! Public Members:
  !!   cells -> Structure that stores cellIdx and pointers to the cells
  !!
  !! Private Members:
  !!   visits   -> Number of times a given cell in the universe has been visited. Used for reordering.
  !!   doUpdate -> Whether to update the cell ordering
  !!   permute  -> Array of permutation indices to change the cell search order
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: cellUniverse
    type(localCell), dimension(:), allocatable            :: cells
    integer(longInt), private, dimension(:), allocatable  :: visits
    logical(defBool), private                             :: doUpdate = .true.
    integer(shortInt), private, dimension(:), allocatable :: permute
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
  end type cellUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(cellUniverse), intent(inout)                        :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), dimension(:), allocatable  :: cellTemp
    real(defReal), dimension(:), allocatable      :: temp
    integer(shortInt)                             :: id, N, i
    character(100), parameter :: Here = 'init (cellUniverse_class.f90)'

    ! Load basic data
    call dict % get(id, 'id')
    if (id <= 0) call fatalError(Here, 'Universe ID must be +ve. Is: '//numToChar(id))
    call self % setId(id)

    ! Load origin
    if (dict % isPresent('origin')) then
      call dict % get(temp, 'origin')

      if (size(temp) /= 3) then
        call fatalError(Here, 'Origin must have size 3. Has: '//numToChar(size(temp)))
      end if
      call self % setTransform(origin=temp)

    end if

    ! Load rotation
    if (dict % isPresent('rotation')) then
      call dict % get(temp, 'rotation')

      if (size(temp) /= 3) then
        call fatalError(Here, '3 rotation angles must be given. Has only: '//numToChar(size(temp)))
      end if
      call self % setTransform(rotation=temp)
    end if

    ! Load Cells by ID
    call dict % get(cellTemp, 'cells')
    N = size(cellTemp)

    allocate(self % cells(N))
    self % cells % idx = cellTemp

    ! Create array of times a cell has been visited
    allocate(self % visits(N))
    self % visits = 0
    self % doUpdate = .true.

    ! Create a permutation array - initially boring, stores new cell indices after permutation
    allocate(self % permute(N))
    do i = 1, N
      self % permute(i) = i
    end do

    ! Load pointers and convert cell ID to IDX
    do i = 1, N
      ! Convert cell ID to IDX
      self % cells(i) % idx = cells % getIdx(self % cells(i) % idx)
      self % cells(i) % ptr => cells % getPtr(self % cells(i) % idx)
    end do

    ! Create fill array
    allocate(fill(N+1))

    do i = 1, N
      fill(i) = cells % getFill(self % cells(i) % idx)
    end do
    fill(N+1) = UNDEF_MAT

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(cellUniverse), intent(inout)         :: self
    integer(shortInt), intent(out)             :: localID
    integer(shortInt), intent(out)             :: cellIdx
    real(defReal), dimension(3), intent(in)    :: r
    real(defReal), dimension(3), intent(in)    :: u
    integer(shortInt)                          :: idx

    ! Reorder cells after maximum number of samples
    if (self % doUpdate .and. sum(self % visits) >= 5000) then
      
      !$omp critical
      self % doUpdate = .false.
      do idx = 1, size(self % cells)
        localID = maxloc(self % visits, 1)
        if (self % visits(localID) == -1) exit
        self % permute(idx) = localID
        self % visits(localID) = -1
      end do      
      !$omp end critical
    
    end if

    ! Search all cells
    do idx = 1, size(self % cells)
      localID = self % permute(idx)
      if (self % cells(localID) % ptr % inside(r, u)) then
        cellIdx = self % cells(localID) % idx
        if (self % doUpdate) then
          !$omp atomic
          self % visits(idx) = self % visits(idx) + 1
        end if
        return

      end if
    end do
    
    ! Previous code for a naive search without reordering
    !do localID = 1, size(self % cells)
    !  if (self % cells(localID) % ptr % inside(r, u)) then
    !    cellIdx = self % cells(localID) % idx
    !    return
    !
    !  end if
    !end do

    ! If not found return undefined cell
    ! Already set to localID (== size(self % cells) + 1) by the do loop
    cellIdx = 0

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if in UNDEFINED cell
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(cellUniverse), intent(inout) :: self
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    type(coord), intent(in)            :: coords
    integer(shortInt)                  :: localID
    character(100), parameter :: Here = 'distance (cellUniverse_class.f90)'

    localID = coords % localID

    if (localID > size(self % cells)) then
      call fatalError(Here, 'Particle is in undefined local cell. Local ID: '//numToChar(localID))
    end if

    ! Calculate distance
    call self % cells(localID) % ptr % distance(d, surfIdx, coords % r, coords % dir)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Note: Introduces extra movment to the particle to push it over boundary
  !!   for more efficent search. Distance is NUGDE.
  !!
  subroutine cross(self, coords, surfIdx)
    class(cellUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    integer(shortInt), intent(in)      :: surfIdx

    ! NUDGE position slightly forward to escape surface tolerance
    ! and avoid calculating normal and extra dot-products
    coords % r = coords % r + coords % dir * NUDGE

    ! Find cell
    ! TODO: Some cell neighbout list
    call self % findCell(coords % localID, &
                         coords % cellIdx, &
                         coords % r,       &
                         coords % dir)

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(cellUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cellUniverse), intent(inout) :: self

    ! SUperclass
    call kill_super(self)

    ! Local
    if(allocated(self % cells)) deallocate(self % cells)

  end subroutine kill


end module cellUniverse_class
