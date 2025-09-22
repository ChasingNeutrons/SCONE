module arraysRR_class

  use numPrecision
  use universalVariables
  use constantsRR
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Data
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use dataRR_class,                   only : dataRR

  ! Geometry
  use coord_class,                    only : coordList
  use geometryStd_class,              only : geometryStd

  ! Visualisation
  use visualiser_class,               only : visualiser
  
  ! Tallying
  use tallyMap_inter,                 only : tallyMap
  use particle_class,                 only : particle, particleState
  use tallyAdmin_class,               only : tallyAdmin
 
  ! Random ray info
  use cellRR_class,                   only : cellRR 

  ! For locks
  use omp_lib

  implicit none
  private
 
  
  !!
  !! Object to store all necessary data in random ray
  !! By default, will have the array of cells which contain flux and geometry info.
  !! Also contains nuclear data and random ray settings
  !!
  !! Private Members
  !!   nG             -> Number of energy groups, kept for convenience.
  !!   nCells         -> Number of unique cells in the geometry, kept for convenience.
  !!   lengthPerIt    -> RR active length per iteration, kept for convenience
  !!   XSData         -> Pointer to nuclear data, for convenience.
  !!   geom           -> Pointer to geometry, for convenience.
  !!   rho            -> Stabilisation factor: 0 is no stabilisation, 1 is aggressive stabilisation
  !!   simulationType -> Identifies which simulation to perform: flat/linear, isotropic/anisotropic
  !!   set2D          -> Stabilises LS in 2D problems if true (zeros Z moment)
  !!
  type, public :: arraysRR
    private
    ! Components
    class(cellRR), dimension(:), allocatable, public :: cells
    class(geometryStd), pointer              :: geom => null()
    type(dataRR)                             :: XSData      
    integer(shortInt)                        :: nG = 0
    integer(shortInt)                        :: nCells = 0
    real(defReal)                            :: lengthPerIt = ZERO
    real(defFlt)                             :: rho = 0.0_defFlt
    integer(shortInt)                        :: simulationType = 0
    real(defReal)                            :: totalVolume = ONE
    integer(shortInt)                        :: volPolicy   = hybrid !simAverage
    integer(shortInt)                        :: missPolicy  = hybrid !srcPolicy
    logical(defBool)                         :: set2D       = .false. ! Stabilises LS in 2D problems
    logical(defBool)                         :: hasFixedSource = .false.
    
    ! Other data
    real(defReal)     :: averageHit = ZERO
    integer(shortInt) :: iterations = 0

  contains
    
    ! Public procedures
    procedure :: init
    procedure :: initAdjoint
    procedure :: kill

    ! Access procedures
    procedure :: getDataPointer
    procedure :: getGeomPointer
    procedure :: getNG
    procedure :: getCellHitRate
    procedure :: getSimulationType
    procedure :: getFluxAtAPoint
    procedure :: countFound
    
    ! Change individual elements of the type
    ! Predominantly for use in the transport sweep
    procedure :: getAverageHitRate
    procedure :: wipeCellHits
    procedure :: setActiveLength

    ! Basic RR procedures
    procedure :: resetFluxes
    procedure :: normaliseFluxAndVolume
    procedure :: updateSource
    procedure :: accumulateFluxScores
    procedure :: finaliseFluxScores
    procedure :: calculateKeff
    procedure :: zeroPrevFlux

    ! Tally results for use with MC tally machinery
    procedure :: tallyResults

    ! Output procedures
    procedure :: outputToVTK
    procedure :: outputPointFluxes
    
    ! Private procedures
    procedure, private :: initialiseFixedSource
    
    procedure, private :: normaliseFluxAndVolumeFlatIso
    procedure, private :: normaliseFluxAndVolumeLinearIso
    procedure, private :: normaliseFluxAndVolumeLIFA
    procedure, private :: normaliseFluxAndVolumeFlatAni
    
    procedure, private :: sourceUpdateKernelFlatIso
    procedure, private :: sourceUpdateKernelLinearIso
    procedure, private :: sourceUpdateKernelLIFA
    procedure, private :: sourceUpdateKernelFlatAni
    
    procedure, private :: calculateKeffKernel

    procedure, private :: countHit

  end type arraysRR

contains

  !!
  !! Initialise the arrays object
  !!
  !! The object is fed sizes and requirements by the physics package.
  !! This will allocate the necessary arrays
  !!
  subroutine init(self, db, geom, lengthPerIt, rho, lin, doKinetics, loud, &
                  dictFS, volPolicy, missPolicy, set2D)
    class(arraysRR), intent(inout)                      :: self
    class(baseMgNeutronDatabase), pointer, intent(in)   :: db
    class(geometryStd), pointer, intent(in)             :: geom
    real(defReal), intent(in)                           :: lengthPerIt
    real(defReal), intent(in)                           :: rho
    logical(defBool), intent(in)                        :: lin
    logical(defBool), intent(in)                        :: doKinetics
    logical(defBool), intent(in)                        :: loud
    class(dictionary), pointer, intent(inout), optional :: dictFS
    integer(shortInt), intent(in), optional             :: volPolicy, missPolicy
    logical(defBool), intent(in), optional              :: set2D
    integer(shortInt)                                   :: ani, i
    logical(defBool)                                    :: do2D
    real(defReal), dimension(6)                         :: bb
    character(100), parameter :: Here = 'init (arraysRR_class.f90)'

    call self % XSData % init(db, doKinetics, loud)
    self % nG          = self % XSdata % getNG()
    self % geom        => geom
    self % nCells      = self % geom % numberOfCells()
    
    self % lengthPerIt = lengthPerIt
    self % rho = real(rho, defFlt)

    if (present(volPolicy)) then
      self % volPolicy = volPolicy
    else
      self % volPolicy = simAverage
    end if
    if (present(missPolicy)) then
      self % missPolicy = missPolicy
    else
      self % missPolicy = srcPolicy
    end if

    if (present(set2D)) then
      do2D = set2D
    else
      do2D = .false.
    end if

    ! Assume bounding box of the geometry is filled (and a box)
    ! Can this be relaxed in future?
    bb = self % geom % bounds()
    self % totalVolume = (bb(4) - bb(1)) * (bb(5) - bb(2)) * (bb(6) - bb(3))

    ! Set simulation type
    ! TODO: read ani from nuclear data
    ani = 0
    if (.not. lin .and. ani == 0) then
      self % simulationType = flatIso
    elseif (lin .and. ani == 0) then
      self % simulationType = linearIso
    elseif (.not. lin .and. ani > 0) then
      self % simulationType = flatAni
    else
      self % simulationType = linearAni
    end if

    ! Allocate and initialise cells
    allocate(self % cells(self % nCells))
    do i = 1, self % nCells
      call self % cells(i) % init(self % nG, lin, ani, do2D)
    end do

    ! Initialise the fixed source if present
    if (present(dictFS)) then
      call self % initialiseFixedSource(dictFS)
    end if

  end subroutine init
  
  !!
  !! Overwrite the active length
  !!
  subroutine setActiveLength(self, lengthPerIt)
    class(arraysRR), intent(inout) :: self
    real(defReal), intent(in)      :: lengthPerIt

    self % lengthPerIt = lengthPerIt

  end subroutine setActiveLength
  
  !!
  !! Initialise the adjoint source and update nuclear data.
  !! For now, assumes the adjoint is for global variance reduction.
  !!
  subroutine initAdjoint(self)
    class(arraysRR), intent(inout)         :: self
    integer(shortInt)                      :: cIdx

    call self % xsData % setAdjointXS()

    !$omp parallel do
    do cIdx = 1, self % nCells
      call self % cells(cIdx) % initAdjoint()
    end do
    !$omp end parallel do

  end subroutine initAdjoint

  !!
  !! Initialises fixed sources to be used in the simulation.
  !! Takes a dictionary containing names of materials in the geometry and
  !! source strengths in each energy group and places these in cells containing
  !! these materials.
  !!
  subroutine initialiseFixedSource(self, dict)
    class(arraysRR), intent(inout)               :: self
    class(dictionary), intent(inout)             :: dict
    character(nameLen),dimension(:), allocatable :: names
    real(defReal), dimension(:), allocatable     :: sourceStrength
    integer(shortInt)                            :: i, nSource, cIdx
    integer(shortInt), save                      :: matIdx
    logical(defBool)                             :: found
    character(nameLen)                           :: sourceName
    character(nameLen), save                     :: localName
    character(100), parameter :: Here = 'initialiseFixedSource (arraysRR_class.f90)'
    !$omp threadprivate(matIdx, localName)

    call dict % keys(names)

    nSource = size(names)

    ! Cycle through entries of the dictionary
    do i = 1, nSource

      sourceName = names(i)
      call dict % get(sourceStrength, sourceName)

      ! Ensure correct number of energy groups
      if (size(sourceStrength) /= self % nG) call fatalError(Here,'Source '//sourceName//&
              ' has '//numToChar(size(sourceStrength))//' groups rather than '//numToChar(self % nG))

      ! Make sure that the source corresponds to a material present in the geometry
      found = .false.
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells

        matIdx    = self % geom % geom % graph % getMatFromUID(cIdx)
        localName = self % XSData % getName(matIdx)

        if (localName == sourceName) then

          found = .true.
          call self % cells(cIdx) % initialiseFixedSource(sourceStrength)

        end if

      end do
      !$omp end parallel do

      if (.not. found) call fatalError(Here,'The source '//trim(sourceName)//' does not correspond to '//&
              'any material found in the geometry.')

    end do
    self % hasFixedSource = .true.

  end subroutine initialiseFixedSource

  !!
  !! Return a pointer to the nuclear data object
  !!
  function getDataPointer(self) result(dataPtr)
    class(arraysRR), intent(in), target :: self
    class(dataRR), pointer              :: dataPtr

    dataPtr => self % XSData

  end function getDataPointer
  
  !!
  !! Return a pointer to the geometry object
  !!
  function getGeomPointer(self) result(geomPtr)
    class(arraysRR), intent(in), target :: self
    class(geometryStd), pointer         :: geomPtr

    geomPtr => self % geom

  end function getGeomPointer
  
  !!
  !! Return the simulation type
  !!
  function getSimulationType(self) result(simType)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: simType

    simType = self % simulationType

  end function getSimulationType

  !!
  !! Return the cell hit rate for the given iteration
  !! Also accumulate to average hit rate.
  !!
  !! Only averages after 20 iterations to account for
  !! requiring several iterations to determine which cells 
  !! are present in the geometry.
  !!
  function getCellHitRate(self, it) result(hitRate)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)              :: totalHit, realCells
    real(defReal)                  :: hitRate

    ! Reset averages after iteration 20
    if (it == 21) then
      self % averageHit = ZERO
      self % iterations = 0
    end if

    totalHit = self % countHit()

    if (it > 20) then
      realCells = self % countFound()
    else
      realCells = self % nCells
    end if
    hitRate = real(totalHit,defReal) / realCells

    self % averageHit = self % averageHit + hitRate
    self % iterations = self % iterations + 1

  end function getCellHitRate

  !!
  !! Return the simulation average cell hit rate
  !!
  function getAverageHitRate(self) result(hitRate)
    class(arraysRR), intent(in)   :: self
    real(defReal)                 :: hitRate

    hitRate = self % averageHit / self % iterations

  end function getAverageHitRate

  !! 
  !! Wipe cell hits
  !!
  subroutine wipeCellHits(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: i

    !$omp parallel do
    do i = 1, self % nCells
      call self % cells(i) % unhitCell()
    end do
    !$omp end parallel do

  end subroutine wipeCellHits
  
  !!
  !! Return number of cells found
  !!
  function countFound(self) result(n)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: n
    integer(shortInt)           :: i

    n = 0
    !$omp parallel do reduction(+:n)
    do i = 1, self % nCells
      n = n + merge(1, 0, self % cells(i) % cellTotalHit > 0)
    end do
    !$omp end parallel do

  end function countFound
  
  !!
  !! Return number of cells hit this iteration
  !!
  function countHit(self) result(n)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: n
    integer(shortInt)           :: i

    n = 0
    !$omp parallel do reduction(+:n)
    do i = 1, self % nCells
      n = n + merge(1, 0, self % cells(i) % wasHit())
    end do
    !$omp end parallel do

  end function countHit

  !!
  !! Return number of energy groups used
  !!
  elemental function getNG(self) result(nG)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: nG

    nG = self % nG

  end function getNG
  
  !!
  !! Calls appropriate normalise flux and volume subroutines
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    character(100), parameter :: Here = 'normaliseFluxAndVolume (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso)
        call self % normaliseFluxAndVolumeFlatIso(it)
      case(linearIso)
        call self % normaliseFluxAndVolumeLinearIso(it)
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select
    
  end subroutine normaliseFluxAndVolume

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolumeFlatIso(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    real(defReal)                             :: norm, normIt
    real(defReal), save                       :: vol, volAve, volNaive, D
    real(defFlt), save                        :: sigGG, tot
    real(defFlt), dimension(:), pointer, contiguous, save :: total
    integer(shortInt), save                   :: g, matIdx
    integer(shortInt)                         :: cIdx
    logical(defBool), save                    :: hit, isSrc, smallCell
    character(100), parameter :: Here = 'normaliseFluxAndVolumeFlatIso (arraysRR_class.f90)'
    !$omp threadprivate(total, vol, g, matIdx, sigGG, D, hit, isSrc, volAve, volNaive, tot, smallCell)

    norm = ONE / self % lengthPerIt
    normIt = ONE / (self % lengthPerIt * it)
    
    !$omp parallel do 
    cellLoop: do cIdx = 1, self % nCells
      associate(cell => self % cells(cIdx))
      if (.not. cell % wasFound()) cycle cellLoop
      matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
      
      hit = cell % wasHit()
      isSrc = cell % hasFixedSource()

      ! Is the cell hit frequently?
      smallCell = (real(cell % cellTotalHit / it, defReal) < 1.5)
      
      ! Compute various volume types
      call cell % scaleVolume(norm, normIt, volAve, volNaive)

      ! Decide volume to use
      select case(self % volPolicy)
        case(simAverage)
          vol = volAve
        case(naive)
          vol = volNaive
        case(hybrid)
          if (isSrc) then
            vol = volNaive
          else
            vol = volAve
          end if
        case default
          call fatalError(Here,'Unsupported volume handling requested')
      end select

      if (smallCell) vol = volNaive
      
      ! Reset cycle-wise estimator
      cell % volumeTracks = ZERO
            
      call self % XSData % getTotalPointer(matIdx, total)

      groupLoop: do g = 1, self % nG

        tot = total(g)

        ! Route for non-void materials
        if (matIdx <= self % XSData % getNMat() .and. tot > 0) then
          if (hit) then
     
            ! Can hit a cell but with a tiny volume, such that 
            ! things break a bit - would rather remove this arbitrary
            ! check in future
            if (vol < volume_tolerance) then
              cell % scalarFlux(g) = ZERO
              cycle groupLoop
            end if

            cell % scalarFlux(g) = cell % scalarFlux(g) * norm / (vol * tot)

            ! Presumes non-zero total XS
            sigGG = self % XSData % getScatterXS(matIdx, g, g)
            if ((sigGG < 0) .and. (total(g) > 0)) then
              D = -real(self % rho * sigGG / tot, defReal)
            else
              D = ZERO
            end if

            cell % scalarFlux(g) =  (cell % scalarFlux(g) + &
                    cell % source(g) / tot + D * cell % prevFlux(g) ) / (1 + D)

          else
            
            ! Decide flux treatment to use on missing a cell
            select case(self % missPolicy)
              case(srcPolicy)
                cell % scalarFlux(g) = cell % source(g) / tot
              case(prevPolicy)
                cell % scalarFlux(g) = cell % prevFlux(g)
              case(hybrid)
                if (isSrc) then
                  cell % scalarFlux(g) = cell % prevFlux(g)
                else
                  cell % scalarFlux(g) = cell % source(g) / tot
                end if
              case default
                call fatalError(Here,'Unsupported miss handling requested')
            end select

          end if
        
        ! Alternatively, handle unidentified/void regions
        else
          
          if (vol < volume_tolerance) then
            cell % scalarFlux(g) = ZERO
            cycle groupLoop
          end if

          if (hit) then
            cell % scalarFlux(g) = cell % scalarFlux(g) * norm / vol
          else
            cell % scalarFlux(g) = cell % prevFlux(g)
          end if
        end if

      end do groupLoop
      end associate

    end do cellLoop
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeFlatIso
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for linear isotropic sources
  !!
  subroutine normaliseFluxAndVolumeLinearIso(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    real(defReal)                             :: norm, normVol
    real(defReal), save                       :: vol, volNaive, volAve, norm_V
    real(defFlt), save                        :: D, sigGG
    real(defFlt), dimension(:), pointer, contiguous, save :: total
    integer(shortInt)                         :: cIdx
    integer(shortInt), save                   :: g, matIdx
    logical(defBool), save                    :: hit, isSrc, smallCell
    character(100), parameter :: Here = 'normaliseFluxAndVolumeLinearIso (arraysRR_class.f90)'
    !$omp threadprivate(total, vol, g, matIdx, norm_V, D, sigGG, hit, isSrc, volNaive, volAve, smallCell)

    norm = ONE / self % lengthPerIt
    normVol = ONE / (self % lengthPerIt * it)


    !$omp parallel do schedule(static)
    cellLoop: do cIdx = 1, self % nCells
      associate(cell => self % cells(cIdx))
      if (.not. cell % wasFound()) cycle cellLoop
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx)

      hit = cell % wasHit()
      isSrc = cell % hasFixedSource()
      ! Is the cell hit frequently?
      smallCell = (real(cell % cellTotalHit / it, defReal) < 1.5)

      ! Compute various volume types
      call cell % scaleVolume(norm, normVol, volAve, volNaive)
      
      ! Decide volume to use
      select case(self % volPolicy)
        case(naive)
          vol = volNaive
          norm_V = ONE / cell % volumeTracks
        case(simAverage)
          vol = volAve
          norm_V = it / cell % allVolumeTracks
        case(hybrid)
          if (isSrc) then
            vol = volNaive
            norm_V = ONE / cell % volumeTracks
          else
            vol = volAve
            norm_V = it / cell % allVolumeTracks
          end if
        case default
          call fatalError(Here,'Unsupported volume handling requested')
      end select

      if (smallCell) then
        vol = volNaive
        norm_V = ONE / cell % volumeTracks
      end if
      
      ! Reset cycle-wise estimator
      cell % volumeTracks = ZERO
       
      ! Update centroid and moment matrix information provided volume has been visited
      call cell % scaleCentroidAndMomMat()

      call self % XSData % getTotalPointer(matIdx, total)

      groupLoop: do g = 1, self % nG

        if (matIdx <= self % XSData % getNMat() .and. total(g) > 0) then
          if (hit) then
        
            ! Can hit a cell but with a tiny volume, such that 
            ! things break a bit - would rather remove this arbitrary
            ! check in future
            if (vol < volume_tolerance) then
              cell % scalarFlux(g) = ZERO
              cell % scalarX(g) = ZERO
              cell % scalarY(g) = ZERO
              cell % scalarZ(g) = ZERO
              cycle groupLoop
            end if
      
            cell % scalarFlux(g) = cell % scalarFlux(g) * norm_V / total(g)
            cell % scalarX(g) = cell % scalarX(g) * norm_V / total(g)
            cell % scalarY(g) = cell % scalarY(g) * norm_V / total(g)
            cell % scalarZ(g) = cell % scalarZ(g) * norm_V / total(g)
          
            ! Apply the standard MoC post-sweep treatment and
            ! stabilisation for negative XSs
            ! Presumes non-zero total XS
            sigGG = self % XSData % getScatterXS(matIdx, g, g)
            if (sigGG < 0) then
              D = -self % rho * sigGG / total(g)
            else
              D = 0.0_defFlt
            end if
           
            cell % scalarFlux(g) =  (cell % scalarFlux(g) + cell % source(g)/total(g) &
                    + D * cell % prevFlux(g)) / (1 + D)

          else
            ! Decide flux treatment to use
            associate(mat => cell % momMat)
            select case(self % missPolicy)
              ! Note: this is policy to use the source, not policy for hitting a fixed source
              case(srcPolicy)
                cell % scalarFlux(g) = cell % source(g) / total(g)
                ! OPENMC SETS MOMENTS TO ZERO
                !self % scalarX(idx) = 0.0_defFlt
                !self % scalarY(idx) = 0.0_defFlt
                !self % scalarZ(idx) = 0.0_defFlt
                ! Need to multiply source gradients by moment matrix
                cell % scalarX(g) = real(mat(xx) * cell % sourceX(g) + &
                        mat(xy) * cell % sourceY(g) + mat(xz) * cell % sourceZ(g),defFlt)/ total(g)
                cell % scalarY(g) = real(mat(xy) * cell % sourceX(g) + &
                        mat(yy) * cell % sourceY(g) + mat(yz) * cell % sourceZ(g),defFlt)/ total(g)
                cell % scalarZ(g) = real(mat(xz) * cell % sourceX(g) + &
                        mat(yz) * cell % sourceY(g) + mat(zz) * cell % sourceZ(g),defFlt)/ total(g)
              case(prevPolicy)
                cell % scalarFlux(g) = cell % prevFlux(g)
                cell % scalarX(g) = cell % prevX(g)
                cell % scalarY(g) = cell % prevY(g)
                cell % scalarZ(g) = cell % prevZ(g)
              case(hybrid)
                if (isSrc) then
                  cell % scalarFlux(g) = cell % prevFlux(g)
                  cell % scalarX(g) = cell % prevX(g)
                  cell % scalarY(g) = cell % prevY(g)
                  cell % scalarZ(g) = cell % prevZ(g)
                else
                  cell % scalarFlux(g) = cell % source(g) / total(g)
                  ! OPENMC SETS MOMENTS TO ZERO
                  !self % scalarX(idx) = 0.0_defFlt
                  !self % scalarY(idx) = 0.0_defFlt
                  !self % scalarZ(idx) = 0.0_defFlt
                  ! Need to multiply source gradients by moment matrix
                  cell % scalarX(g) = real(mat(xx) * cell % sourceX(g) + &
                        mat(xy) * cell % sourceY(g) + mat(xz) * cell % sourceZ(g),defFlt)/ total(g)
                  cell % scalarY(g) = real(mat(xy) * cell % sourceX(g) + &
                        mat(yy) * cell % sourceY(g) + mat(yz) * cell % sourceZ(g),defFlt)/ total(g)
                  cell % scalarZ(g) = real(mat(xz) * cell % sourceX(g) + &
                        mat(yz) * cell % sourceY(g) + mat(zz) * cell % sourceZ(g),defFlt)/ total(g)
                end if
              case default
                call fatalError(Here,'Unsupported miss handling requested')
            end select
            end associate
          end if

        else

          ! Apply void treatment
          if (vol < volume_tolerance) then
            cell % scalarFlux(g) = ZERO
            cell % scalarX(g) = ZERO
            cell % scalarY(g) = ZERO
            cell % scalarZ(g) = ZERO
            cycle groupLoop
          end if

          if (hit) then
            cell % scalarFlux(g) = cell % scalarFlux(g) * norm / vol
            cell % scalarX(g) = ZERO
            cell % scalarY(g) = ZERO
            cell % scalarZ(g) = ZERO
          else
            cell % scalarFlux(g) = cell % prevFlux(g)
            cell % scalarX(g) = ZERO
            cell % scalarY(g) = ZERO
            cell % scalarZ(g) = ZERO
          end if

        end if
        
        ! For stability while still accumulating geometric info
        if (it < 10) then
          cell % scalarX(g) = ZERO
          cell % scalarY(g) = ZERO
          cell % scalarZ(g) = ZERO
        end if

      end do groupLoop
      end associate

    end do cellLoop
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeLinearIso
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for flat anisotropic sources
  !!
  subroutine normaliseFluxAndVolumeFlatAni(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it

  end subroutine normaliseFluxAndVolumeFlatAni
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for Linear sources with flat
  !! anisotropic sources
  !!
  subroutine normaliseFluxAndVolumeLIFA(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it

  end subroutine normaliseFluxAndVolumeLIFA
  
  !!
  !! Update all sources given a prevFlux.
  !! Uses ONE_KEFF to scale the fission source.
  !! Uses it to check whether cells which have never been hit
  !! can be neglected.
  !! This nesting allows using combined OMP + SIMD
  !!
  subroutine updateSource(self, ONE_KEFF, it)
    class(arraysRR), intent(inout) :: self
    real(defReal), intent(in)      :: ONE_KEFF
    integer(shortInt), intent(in)  :: it
    real(defFlt)                   :: ONE_K
    integer(shortInt)              :: cIdx
    character(100), parameter      :: Here = 'updateSource (arraysRR_class.f90)'

    ONE_K = real(ONE_KEFF, defFlt)

    select case(self % simulationType)
      case(flatIso)
        !$omp parallel do 
        do cIdx = 1, self % nCells
          call self % sourceUpdateKernelFlatIso(cIdx, ONE_K, it)
        end do
        !$omp end parallel do
      case(linearIso)
        !$omp parallel do 
        do cIdx = 1, self % nCells
          call self % sourceUpdateKernelLinearIso(cIdx, ONE_K, it)
        end do
        !$omp end parallel do
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine updateSource

  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernelFlatIso(self, cIdx, ONE_KEFF, it)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF
    integer(shortInt), intent(in)            :: it
    logical(defBool)                         :: notFound
    real(defFlt)                             :: scatter, fission
    real(defFlt), dimension(self % nG)       :: fluxFlt
    real(defFlt), dimension(:), pointer, contiguous :: nuFission, chi, scatterXS, scatterVec
    integer(shortInt)                        :: matIdx, g, gIn, sIdx1, sIdx2

    associate(cell => self % cells(cIdx))
    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    ! Also check whether cell has even been visited
    notFound = (it > 1 .and. .not. cell % wasFound())
    if (matIdx > self % XSData % getNMat() .or. notFound) then
      do g = 1, self % nG
        cell % source(g) = 0.0_defFlt
        if (allocated(cell % fixedSource)) then
          cell % source(g) = cell % fixedSource(g)
        end if
      end do
      return
    end if

    ! Obtain XSs
    call self % XSData % getProdPointers(matIdx, nuFission, scatterXS, chi)

    fluxFlt = real(cell % prevFlux, defFlt)

    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission)
    do gIn = 1, self % nG
      fission = fission + fluxFlt(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF

    do g = 1, self % nG

      sIdx1 = self % nG * (g - 1) + 1
      sIdx2 = self % nG * g
      scatterVec => scatterXS(sIdx1:sIdx2)

      ! Calculate scattering source
      scatter = 0.0_defFlt
      !$omp simd reduction(+:scatter)
      do gIn = 1, self % nG
        scatter = scatter + fluxFlt(gIn) * scatterVec(gIn)
      end do

      ! Output index
      cell % source(g) = chi(g) * fission + scatter
      if (allocated(cell % fixedSource)) then
        cell % source(g) = cell % source(g) + cell % fixedSource(g)
      end if

    end do
    end associate

  end subroutine sourceUpdateKernelFlatIso
  
  !!
  !! Kernel to update sources given a cell index for linear sources
  !! with isotropic scattering
  !!
  subroutine sourceUpdateKernelLinearIso(self, cIdx, ONE_KEFF, it)
    class(arraysRR), target, intent(inout)  :: self
    integer(shortInt), intent(in)           :: cIdx
    real(defFlt), intent(in)                :: ONE_KEFF
    integer(shortInt), intent(in)           :: it
    logical(defBool)                        :: notFound
    real(defFlt)                            :: scatter, xScatter, yScatter, zScatter, &
                                               fission, xFission, yFission, zFission, &
                                               xSource, ySource, zSource
    real(defFlt), dimension(matSize)        :: invM
    real(defFlt), dimension(self % nG)      :: fluxFlt, xFlt, yFlt, zFlt
    real(defFlt), dimension(:), pointer, contiguous :: nuFission, chi, scatterXS, scatterVec
    integer(shortInt)                       :: matIdx, g, gIn, sIdx1, sIdx2

    associate(cell => self % cells(cIdx))
    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx)
    
    ! Guard against void cells
    ! Also check whether cell has even been visited
    notFound = (it > 1 .and. .not. cell % wasFound())
    if (matIdx > self % XSData % getNMat() .or. notFound) then
      do g = 1, self % nG
        cell % source(g) = 0.0_defFlt
        if (allocated(cell % fixedSource)) then
          cell % source(g) = cell % fixedSource(g)
        end if
        cell % sourceX(g) = 0.0_defFlt
        cell % sourceY(g) = 0.0_defFlt
        cell % sourceZ(g) = 0.0_defFlt
      end do
      return
    end if
    
    ! Invert moment matrix
    invM = cell % invertMatrix()
     
    ! Obtain XSs
    call self % XSData % getProdPointers(matIdx, nuFission, scatterXS, chi)

    fluxFlt  = real(cell % prevFlux, defFlt)
    xFlt = real(cell % prevX, defFlt)
    yFlt = real(cell % prevY, defFlt)
    zFlt = real(cell % prevZ, defFlt)
    
    ! Calculate fission source
    fission = 0.0_defFlt
    xFission = 0.0_defFlt
    yFission = 0.0_defFlt
    zFission = 0.0_defFlt

    !$omp simd reduction(+:fission, xFission, yFission, zFission)
    do gIn = 1, self % nG
      fission = fission + fluxFlt(gIn) * nuFission(gIn)
      xFission = xFission + xFlt(gIn) * nuFission(gIn)
      yFission = yFission + yFlt(gIn) * nuFission(gIn)
      zFission = zFission + zFlt(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF
    xFission = xFission * ONE_KEFF
    yFission = yFission * ONE_KEFF
    zFission = zFission * ONE_KEFF

    do g = 1, self % nG

      sIdx1 = self % nG * (g - 1) + 1
      sIdx2 = self % nG * g
      scatterVec => scatterXS(sIdx1:sIdx2)

      ! Calculate scattering source
      scatter = 0.0_defFlt
      xScatter = 0.0_defFlt
      yScatter = 0.0_defFlt
      zScatter = 0.0_defFlt
      !$omp simd reduction(+:scatter, xScatter, yScatter, zScatter)
      do gIn = 1, self % nG
        scatter = scatter + fluxFlt(gIn) * scatterVec(gIn)
        xScatter = xScatter + xFlt(gIn) * scatterVec(gIn)
        yScatter = yScatter + yFlt(gIn) * scatterVec(gIn)
        zScatter = zScatter + zFlt(gIn) * scatterVec(gIn)
      end do

      cell % source(g) = chi(g) * fission + scatter
      if (allocated(cell % fixedSource)) then
        cell % source(g) = cell % source(g) + cell % fixedSource(g)
      end if
      
      xSource = chi(g) * xFission + xScatter
      ySource = chi(g) * yFission + yScatter
      zSource = chi(g) * zFission + zScatter
      
      if (allocated(cell % fixedX)) then
        xSource = xSource + cell % fixedX(g)
        ySource = ySource + cell % fixedY(g)
        zSource = zSource + cell % fixedZ(g)
      end if
      
      ! Calculate source gradients by inverting the moment matrix
      cell % sourceX(g) = invM(xx) * xSource + &
              invM(xy) * ySource + invM(xz) * zSource
      cell % sourceY(g) = invM(xy) * xSource + &
              invM(yy) * ySource + invM(yz) * zSource
      cell % sourceZ(g) = invM(xz) * xSource + &
           invM(yz) * ySource + invM(zz) * zSource
      
    end do
    end associate

  end subroutine sourceUpdateKernelLinearIso
  
  !!
  !! Kernel to update sources given a cell index for flat sources
  !! with anisotropic scattering
  !!
  subroutine sourceUpdateKernelFlatAni(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF

  end subroutine sourceUpdateKernelFlatAni
  
  !!
  !! Kernel to update sources given a cell index for linear sources
  !! with flat anisotropic scattering
  !!
  subroutine sourceUpdateKernelLIFA(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF

  end subroutine sourceUpdateKernelLIFA
  
  !! 
  !! Calculate keff
  !! Wraps the main kernel call to allow for OMP + SIMD (thanks Fortran)
  !!
  function calculateKeff(self, k0) result(k1)
    class(arraysRR), intent(in)           :: self
    real(defReal), intent(in)             :: k0
    real(defReal)                         :: k1
    integer(shortInt)                     :: cIdx
    real(defReal)                         :: fissTotal, prevFissTotal
    real(defReal), save                   :: fissLocal, prevFissLocal
    character(100), parameter             :: Here = 'calculateKeff (arraysRR_class.f90)'
    !$omp threadprivate (fissLocal, prevFissLocal)

    fissTotal     = ZERO
    prevFissTotal = ZERO
    !$omp parallel do reduction(+:fissTotal, prevFissTotal)
    do cIdx = 1, self % nCells
      call self % calculateKeffKernel(cIdx, fissLocal, prevFissLocal)
      fissTotal     = fissTotal + fissLocal
      prevFissTotal = prevFissTotal + prevFissLocal
    end do 
    !$omp end parallel do

    k1 = k0 * fissTotal / prevFissTotal 
    if ((k1 <= 0) .or. (k1 > 5)) call fatalError(Here, 'Unphysical keff: '//numToChar(k1))
    if (k1 /= k1) call fatalError(Here, 'NaN keff')

  end function calculateKeff
  
  !!
  !! Calculate keff for a single cell
  !!
  subroutine calculateKeffKernel(self, cIdx, fissionRate, prevFissionRate)
    class(arraysRR), target, intent(in)  :: self
    integer(shortInt), intent (in)       :: cIdx
    real(defReal), intent(out)           :: fissionRate, prevFissionRate
    real(defReal)                        :: vol
    integer(shortInt)                    :: g, matIdx
    real(defFlt), dimension(:), pointer, contiguous :: nuSigmaF
    real(defReal), dimension(:), pointer, contiguous :: flux, prevFlux

    fissionRate     = ZERO
    prevFissionRate = ZERO

    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Check whether to continue in this cell
    if (matIdx > self % XSData % getNMat()) return
    if (.not. self % XSData % isFissile(matIdx)) return
      
    call self % XSData % getNuFissPointer(matIdx, nuSigmaF)
    
    if (.not. self % cells(cIdx) % wasFound()) return
    vol = self % cells(cIdx) % volume

    flux => self % cells(cIdx) % scalarFlux
    prevFlux => self % cells(cIdx) % prevFlux

    !$omp simd reduction(+: fissionRate, prevFissionRate)
    do g = 1, self % nG
      fissionRate     = fissionRate     + real(flux(g) * nuSigmaF(g), defReal)
      prevFissionRate = prevFissionRate + real(prevFlux(g) * nuSigmaF(g), defReal)
    end do

    fissionRate     = fissionRate * vol
    prevFissionRate = prevFissionRate * vol

  end subroutine calculateKeffKernel
  
  !!
  !! Zero the previous-step flux
  !!
  subroutine zeroPrevFlux(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: i

    !$omp parallel do
    do i = 1, self % nCells
      call self % cells(i) % zeroPrevFlux()
    end do
    !$omp end parallel do

  end subroutine zeroPrevFlux
  
  !!
  !! Reset fluxes
  !!
  subroutine resetFluxes(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: i
    character(100), parameter      :: Here = 'resetFluxes (arraysRR_class.f90)'

    if (.not. allocated(self % cells)) call fatalError(Here, 'Cells not allocated')
    
    !$omp parallel do
    do i = 1, self % nCells
      if (.not. self % cells(i) % wasFound()) cycle
      call self % cells(i) % resetFluxes()
    end do
    !$omp end parallel do

  end subroutine resetFluxes
  
  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScores(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: i
    character(100), parameter      :: Here = 'accumulateFluxScores (arraysRR_class.f90)'

    if (.not. allocated(self % cells)) call fatalError(Here, 'Cells not allocated')
    
    !$omp parallel do
    do i = 1, self % nCells
      if (.not. self % cells(i) % wasFound()) cycle
      call self % cells(i) % accumulateFluxScores()
    end do
    !$omp end parallel do

  end subroutine accumulateFluxScores
  
  !!
  !! Finalise results
  !!
  subroutine finaliseFluxScores(self, it)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)              :: i
    character(100), parameter      :: Here = 'finaliseFluxScores (arraysRR_class.f90)'

    if (.not. allocated(self % cells)) call fatalError(Here, 'Cells not allocated')
    
    !$omp parallel do
    do i = 1, self % nCells
      if (.not. self % cells(i) % wasFound()) cycle
      call self % cells(i) % finaliseFluxScores(it)
    end do
    !$omp end parallel do

  end subroutine finaliseFluxScores

  !!
  !! Tallies flux-related results, using an input tallyAdmin.
  !! Allows for use of MC tally machinery with RR, although limited
  !! to relatively simple estimators.
  !!
  !! Loops over all phase space points, contributing each to tallies.
  !!
  !! Assumes the cycle will be ended afterwards.
  !!
  subroutine tallyResults(self, tally)
    class(arraysRR), intent(in)                :: self
    type(tallyAdmin), pointer, intent(inout)   :: tally
    type(particle), save                       :: p
    real(defReal), save                        :: vol
    real(defReal), dimension(3), save          :: pos
    integer(shortInt), save                    :: g, matIdx
    integer(shortInt)                          :: i
    !$omp threadprivate(p, vol, pos, g, matIdx)

    !$omp parallel
    call p % build([-INFINITY, -INFINITY, -INFINITY], [ONE, ZERO, ZERO], 1, ZERO, ZERO)
    !$omp end parallel

    !$omp parallel do
    do i = 1, self % nCells

      if (.not. self % cells(i) % wasFound()) cycle
      vol = self % cells(i) % getVolume()
      pos = self % cells(i) % getCellPos()
      call p % teleport(pos)

      matIdx = self % geom % geom % graph % getMatFromUID(i) 
      p % coords % matIdx = matIdx 

      do g = 1, self % nG

        ! The weight should be flux * V, which, for a single score in a tally
        ! will produce SigmaX * flux * V or volume-integrated reaction rate.
        ! For now only call at the end and use finalised flux scores
        !p % w = self % cells(i) % scalarFlux(g) * vol
        p % w = self % cells(i) % getFluxScore(g) * vol
        p % G = g

        call tally % reportInColl(p, .false.)

      end do

    end do
    !$omp end parallel do

  end subroutine tallyResults
  
  !!
  !! Send all arrays of interest to VTK output
  !!
  subroutine outputToVTK(self, viz)
    class(arraysRR), intent(in)               :: self
    class(visualiser), intent(inout)          :: viz
    real(defReal), dimension(:), allocatable  :: resVec
    character(nameLen)                        :: name
    integer(shortInt)                         :: cIdx, g

    allocate(resVec(self % nCells))

    ! Output all fluxes (assuming finalisation of scores happened)
    do g = 1, self % nG
      name = 'flux_g'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = self % cells(cIdx) % getFluxScore(g)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output all flux uncertainties
    do g = 1, self % nG
      name = 'std_g'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = self % cells(cIdx) % getFluxSD(g) /self % cells(cIdx) % getFluxScore(g)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output final iteration sources
    do g = 1, self % nG
      name = 'source_'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = real(self % cells(cIdx) % getSource(g),defReal)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do
    
    if (self % hasFixedSource) then
      do g = 1, self % nG
        name = 'fixedSource_'//numToChar(g)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          resVec(cIdx) = real(self % cells(cIdx) % getFixedSource(g),defReal)
        end do
        !$omp end parallel do
        call viz % addVTKData(resVec,name)
      end do
    end if

    ! Output final volume estimates
    name = 'volume'
    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      resVec(cIdx) = self % cells(cIdx) % volume * self % totalVolume
    end do
    !$omp end parallel do
    call viz % addVTKData(resVec,name)

    ! Output material IDs
    name = 'material'
    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      resVec(cIdx) = self % geom % geom % graph % getMatFromUID(cIdx) 
    end do
    !$omp end parallel do
    call viz % addVTKData(resVec,name)

    call viz % finaliseVTK()

  end subroutine outputToVTK

  !!
  !! Output fluxes at given points
  !!
  subroutine outputPointFluxes(self, out, points, names)
    class(arraysRR), intent(in)                  :: self
    class(outputFile), intent(inout)             :: out
    real(defReal), dimension(:,:), intent(in)    :: points
    character(nameLen), dimension(:), intent(in) :: names
    integer(shortInt)                            :: i
    character(nameLen)                           :: name
    real(defReal), dimension(self % nG)          :: flux, fluxSD
    character(100), parameter                    :: Here = 'outputPointFluxes (arraysRR_class.f90)'

    name = 'pointFlux'
    call out % startBlock(name)

    ! Ensure points and names have the correction dimensions
    if (size(points,1) /= 3) call fatalError(Here, 'Points are not 3D.')
    if (size(points,2) /= size(names)) call fatalError(Here, &
            'Different numbers of sample points to sample names.')

    do i = 1, size(names)
    
      call out % startArray(names(i), [self % nG])
      call self % getFluxAtAPoint(points(:, i), flux, fluxSD)
      call out % addResult(flux, fluxSD)
      call out % endArray()

    end do
    call out % endBlock()

  end subroutine outputPointFluxes

  !!
  !! Returns the flux vector at a point in space
  !!
  subroutine getFluxAtAPoint(self, r, flux, fluxSD)
    class(arraysRR), intent(in)                      :: self
    real(defReal), dimension(3), intent(in)          :: r
    real(defReal), dimension(self % nG), intent(out) :: flux
    real(defReal), dimension(self % nG), intent(out) :: fluxSD
    integer(shortInt)                                :: g, matIdx, cIdx, i
    real(defReal), dimension(3)                      :: mom, momSD, centroid, fluxGrad
    real(defFlt), dimension(matSize)                 :: invM

    ! Identify cell at the given point
    call self % geom % whatIsAt(matIdx, cIdx, r)

    if (cIdx > 0) then
      do g = 1, self % nG
      
        flux(g) = self % cells(cIdx) % getFluxScore(g)
        fluxSD(g) = self % cells(cIdx) % getFluxSD(g)

        ! Include linear moments if available
        if ((self % simulationType == linearIso) .or. &
                (self % simulationType == linearAni)) then

          fluxSD(g) = fluxSD(g) * fluxSD(g)

          mom = self % cells(cIdx) % getFluxMoments(g)
          momSD = self % cells(cIdx) % getFluxMomentSDs(g)
          centroid = self % cells(cIdx) % getCentroid()
          invM = self % cells(cIdx) % invertMatrix()
      
          ! Get flux gradients
          fluxGrad(x) = real(invM(xx) * mom(x) + invM(xy) * mom(y) + invM(xz) * mom(z), defReal)
          fluxGrad(y) = real(invM(xy) * mom(x) + invM(yy) * mom(y) + invM(yz) * mom(z), defReal)
          fluxGrad(z) = real(invM(xz) * mom(x) + invM(yz) * mom(y) + invM(zz) * mom(z), defReal)

          ! Note this will not correctly estimate uncertainty as moments are covariant
          do i = 1, 3
            flux(g) = flux(g) + fluxGrad(i) * (r(i) - centroid(i))
            ! Not sure exactly how to propagate uncertainties - need to do some maths
            fluxSD(g) = fluxSD(g) + momSD(i)**2 * (r(i) - centroid(i))**2
          end do

          if (fluxSD(g) > ZERO) fluxSD(g) = sqrt(fluxSD(g))

        end if

      end do

    else
      print *,'WARNING: No cell found at position '//numToChar(r)
      flux = -ONE
      fluxSD = -ONE
    end if

  end subroutine getFluxAtAPoint
  
  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(arraysRR), intent(inout) :: self

    ! Clean standard contents
    if(allocated(self % cells)) deallocate(self % cells)
    self % geom   => null()
    call self % XSData % kill()
    self % nG     = 0
    self % nCells = 0
    self % lengthPerIt = ZERO
    self % rho         = 0.0_defFlt
    self % simulationType = 0
    self % totalVolume = ONE
    self % averageHit = ZERO
    self % iterations = 0
    self % volPolicy = hybrid
    self % missPolicy = hybrid
    self % set2D = .false.

  end subroutine kill

end module arraysRR_class
