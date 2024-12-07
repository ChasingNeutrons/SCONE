module fissionMatrixClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill
  use tallyResult_class,          only : tallyResult

  ! Nuclear Data
  use nuclearDatabase_inter,      only : nuclearDatabase
  use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Response
  use macroResponse_class,        only : macroResponse

  implicit none
  private

  !!
  !! 1-D fission matrix
  !! Generates the fission matrix by considering the birth position of neutrons
  !! in the dungeon at the end of cycle
  !!
  !! Contains only a single map for discretisation
  !!
  !! Notes:
  !!    -> FM is stored in column-major order [prodBin, startBin]
  !!
  !! Private Members:
  !!   map      -> Map to divide phase-space into bins
  !!   startWgt -> Starting Weigths in each bin
  !!   N        -> Number of Bins
  !!
  !! Interface:
  !!   tallyClerk Interface
  !!
  !! Sample dictionary input:
  !!
  !!  clerkName {
  !!      type simpleFMClerk;
  !!      map { <TallyMapDef> }
  !!  }
  !!
  type, public, extends(tallyClerk) :: fissionMatrixClerk
    private
    !! Map defining the discretisation
    class(tallyMap), allocatable             :: map
    real(defReal),dimension(:),allocatable   :: startWgt
    real(defReal),dimension(:),allocatable   :: endWgt
    real(defReal),dimension(:),allocatable   :: eigVec
    integer(shortInt)                        :: N = 0 !! Number of bins
    real(defReal), dimension(:,:), allocatable :: matrix

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportCycleStart
    procedure  :: reportCycleEnd

    ! Overwrite default run-time result procedure
    procedure  :: getResult

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill

    ! Solve for the FM eigenvector
    procedure  :: solve 

  end type fissionMatrixClerk

  !!
  !! Fission matrix result class
  !!    dim1 -> target bin
  !!    dim2 -> orgin bin
  !!
  type,public, extends( tallyResult) :: FMresult
    integer(shortInt)                       :: N  = 0 ! Size of FM
    real(defReal), dimension(:),allocatable :: eigVec ! FM eigenvector
  end type FMResult

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(fissionMatrixClerk), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    character(nameLen), intent(in)      :: name

    ! Assign name
    call self % setName(name)

    ! Read map
    call new_tallyMap(self % map, dict % getDictPtr('map'))

    ! Read size of the map
    self % N = self % map % bins(0)

    ! Allocate space for starting weights
    allocate(self % startWgt(self % N))
    
    ! Allocate end weights for scaling particles
    allocate(self % endWgt(self % N))
    
    ! Allocate fundamental eigenvector for scaling particles
    allocate(self % eigVec(self % N))
    self % eigVec = ONE

    ! Allocate space for the matrix
    allocate(self % matrix(self % N, self % N))

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(fissionMatrixClerk),intent(in)            :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleStart_Code, cycleEnd_Code]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(fissionMatrixClerk), intent(in) :: self
    integer(shortInt)                     :: S

    S = self % N * self % N

  end function getSize

  !!
  !! Process start of the cycle
  !! Calculate starting weights in each bin
  !! Note: done using current positions of particles
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleStart(self, start, mem)
    class(fissionMatrixClerk), intent(inout) :: self
    class(particleDungeon), intent(in)  :: start
    type(scoreMemory), intent(inout)    :: mem
    integer(shortInt)                   :: idx, i

    self % startWgt = ZERO
    self % endWgt = ZERO
    self % matrix = ZERO

    ! Loop through a population and calculate starting weight in each bin
    do i = 1,start % popSize()

      associate (state => start % get(i))
        idx = self % map % map(state)
        if (idx > 0) self % startWgt(idx) = self % startWgt(idx) + state % wgt
      end associate

    end do

  end subroutine reportCycleStart

  !!
  !! Process cycle end
  !! Note: done using current positions of particles and birth positions
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(fissionMatrixClerk), intent(inout) :: self
    class(particleDungeon), intent(in)  :: end
    type(scoreMemory), intent(inout)    :: mem
    integer(shortInt)                   :: p, j
    integer(shortInt), save             :: i, idxOut, idxIn
    real(defReal), save                 :: normFactor
    type(particleState), save           :: birthState
    !$omp threadprivate(normFactor, birthState, i, idxOut, idxIn)

    ! Cycle through bank, identifying particle's birth and death position
    !$omp parallel do
    do p = 1, end % popSize()

      associate(state => end % get(p))
        idxOut = self % map % map(state)
        birthState % r = state % rBirth
        idxIn = self % map % map(birthState)

        ! Score end weight
        if (idxOut > 0) then
          !$omp atomic
          self % endWgt(idxOut) = self % endWgt(idxOut) + state % wgt
        end if

        ! Score element of the matrix
        if((idxOut > 0) .and. (idxIn > 0)) then
          !$omp atomic
          self % matrix(idxOut, idxIn) = self % matrix(idxOut,IdxIn) + state % wgt
        end if
      end associate

    end do
    !$omp end parallel do

    ! Construct fission matrix proper
    !$omp parallel do
    do j = 1, self % N
        
      ! Calculate normalisation factor
      normFactor = self % startWgt(j)
      if (normFactor /= ZERO) normFactor = ONE / normFactor
      
      do i = 1, self % N
        self % matrix(i,j) = self % matrix(i,j) * normFactor
      end do
    end do
    !$omp end parallel do

    ! Obtain the fission matrix eigenvector
    call self % solve()

    ! Modify the eigenvector to scale the dungeon weights
    !$omp parallel do
    do j = 1, self % N
      normFactor = self % endWgt(j)
      if (normFactor /= ZERO) normFactor = ONE / normFactor
      self % eigVec(j) = self % eigVec(j) * normFactor
    end do
    !$omp end parallel do

    self % eigVec = self % eigVec / sum(self % eigVec)
    print *,'Here it is'
    print *, self % eigVec


  end subroutine reportCycleEnd
  
  !!
  !! Solve the fission matrix eigenvalue problem
  !! by power iteration
  !!
  subroutine solve(self)
    class(fissionMatrixClerk), intent(inout) :: self
    real(defReal), dimension(:), allocatable :: b
    real(defReal)                            :: tol, err
    integer(shortInt)                        :: it, i
    integer(shortInt), save                  :: j
    !$omp threadprivate(j)

    tol = 1.0E-3
    err = ONE
    it = 0
    allocate(b(self % N))

    do it = 1, 1000 

      b = self % eigVec 
      self % eigVec = ZERO

      ! Matrix-vector multiply
      !$omp parallel do
      do i = 1, self % N
        do j = 1, self % N
          self % eigVec(i) = self % eigVec(i) + self % matrix(i,j) * b(j)
        end do
      end do
      !$omp end parallel do

      ! Normalise appropriately
      self % eigVec = self % eigVec / norm2(self % eigVec)
      
      err = norm2(self % eigVec - b) / norm2(self % eigVec)
      if (err < tol) exit

    end do

    if (it >= 1000) print *,'FM iterations did not finish'

    self % eigVec = self % eigVec * sum(self % startWgt)

  end subroutine solve

  !!
  !! Return result from the clerk for interaction with Physics Package
  !! Returns FMresult defined in this module
  !! If res is already allocated to a FM of fitting size it reuses already allocated space
  !! This should improve performance when updating estimate of FM each cycle
  !!
  !! See tallyClerk_inter for details
  !!
  pure subroutine getResult(self, res, mem)
    class(fissionMatrixClerk), intent(in)          :: self
    class(tallyResult),allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                  :: mem
    integer(shortInt)                              :: i

    ! Allocate result to FMresult
    ! Do not deallocate if already allocated to FMresult
    ! Its not to nice -> clean up
    if (allocated(res)) then
      select type(res)
        class is (FMresult)
          ! Do nothing

        class default ! Reallocate
          deallocate(res)
          allocate( FMresult :: res)
     end select

    else
      allocate( FMresult :: res)

    end if

    ! Load data inti the FM
    select type(res)
      class is(FMresult)
        ! Check size and reallocate space if needed
        ! This is horrible. Hove no time to polish. Blame me (MAK)
        if (allocated(res % eigVec)) then
          if (any(shape(res % eigVec) /= self % N)) then
            deallocate(res % eigVec)
            allocate(res % eigVec(self % N))
          end if
        else
          allocate(res % eigVec(self % N))
        end if

        ! Set size of the FM
        res % N = self % N

        ! Load entries
        do i = 1,self % N
          res % eigVec(i) = self % eigVec(i)
        end do

    end select

  end subroutine getResult

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(fissionMatrixClerk), intent(in) :: self
    type(scoreMemory), intent(in)    :: mem

    print *, 'fissionMatrixClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(fissionMatrixClerk), intent(in) :: self
    class(outputFile), intent(inout) :: outFile
    type(scoreMemory), intent(in)    :: mem
    integer(shortInt)                :: i, j
    integer(longInt)                 :: addr
    real(defReal)                    :: val, std
    character(nameLen)               :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print map information
    call self % map % print(outFile)

    ! Print fission matrix
    name = 'matrix'
    call outFile % startArray(name, [self % N, self % N])

    do i = 1,self % N
      do j = 1,self % N
        val = self % matrix(i,j)
        call outFile % addResult(val, ZERO)
      end do
    end do
    call outFile % endArray()

    name = 'eig'
    call outFile % startArray(name, [self % N, 1])

    do i = 1,self % N
      val = self % eigVec(i)
      call outFile % addResult(val, ZERO)
    end do
    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  !! See tallyClerk_inter for details
  !!
  elemental subroutine kill(self)
    class(fissionMatrixClerk), intent(inout) :: self

    ! Call superclass
    call kill_super(self)

    if (allocated(self % map)) deallocate(self % map)
    if (allocated(self % startWgt)) deallocate(self % startWgt)
    if (allocated(self % matrix)) deallocate(self % matrix)

    self % N = 0

  end subroutine kill

end module fissionMatrixClerk_class
