module cellRR_class

  use numPrecision
  use universalVariables
  use constantsRR
  use genericProcedures,              only : fatalError, numToChar
  use dictionary_class,               only : dictionary

  ! For locks
  use omp_lib

  implicit none
  private
 
  !!
  !! Object that holds all random ray information about a cell.
  !! Holds geometric and flux information.
  !! Optionally can also hold flux moments.
  !!
  !! Private Members
  !!   nG             -> Number of energy groups, kept for convenience.
  !!   ani            -> Order of anisotropic flux moments to be stored
  !!   set2D          -> Neglects z-moment for LS stability in 2D calculations
  !!
  !!   scalarFlux     -> Array of scalar flux values of length [nG]
  !!   prevFlux       -> Array of previous scalar flux values of length [nG]
  !!
  !!   source      -> Array of sources [nG]
  !!   fixedSource -> Array of fixed sources [nG]
  !!
  !!   volumeTracks  -> Sum of track lengths for computing volumes
  !!   volume        -> Dimensionless cell volume
  !!
  !!   cellHit      -> Whether the cell was visited this iteration
  !!   cellTotalHit -> Total number of hits for a cell over all iterations
  !!   cellFound    -> Logical of whether the cell was found
  !!   cellPos      -> Array of cell coordinates [3]
  !!
  !!   scalarX   -> X-spatial moments of scalar flux [nG]
  !!   scalarY   -> Y-spatial moments of scalar flux [nG]
  !!   scalarZ   -> Z-spatial moments of scalar flux [nG]
  !!   prevX     -> Previous x-spatial moments of scalar flux [nG]
  !!   prevY     -> Previous y-spatial moments of scalar flux [nG]
  !!   prevZ     -> Previous z-spatial moments of scalar flux [nG]
  !!   sourceX   -> Source x-spatial gradients of scalar flux [nG]
  !!   sourceY   -> Source y-spatial gradients of scalar flux [nG]
  !!   sourceZ   -> Source z-spatial gradients of scalar flux [nG]
  !!   momMat    -> Symmetric spatial moment matrix [matSize]
  !!   momTracks -> Weighted tracks used to computer spatial moment matrices [matSize]
  !!   centroid       -> Cell centroid [3]
  !!   centroidTracks -> Weighted tracks used to computer centroids [3]
  !!
  !!   locks -> OpenMP lock for the cell
  !!
  type, public :: cellRR
    
    ! Data is ordered  in terms of frequency of access during the sweep
    ! Geometry info
    logical(defBool)  :: cellHit = .false.
    integer(longInt)  :: cellTotalHit = 0
    logical(defBool)  :: cellFound = .false.
    
    ! Flux arrays
    real(defFlt), dimension(:), allocatable :: source
    
    ! OMP lock
    integer(kind=omp_lock_kind) :: lock
    real(defReal)               :: volumeTracks = ZERO
    
    real(defReal), dimension(:), allocatable :: scalarFlux
    
    ! Linear source arrays
    real(defFlt), dimension(:), allocatable   :: sourceX
    real(defFlt), dimension(:), allocatable   :: sourceY
    real(defFlt), dimension(:), allocatable   :: sourceZ
    real(defReal), dimension(:), allocatable  :: scalarX
    real(defReal), dimension(:), allocatable  :: scalarY
    real(defReal), dimension(:), allocatable  :: scalarZ
    
    real(defReal), dimension(matSize)          :: momMat = ZERO
    real(defReal), dimension(matSize)          :: momTracks = ZERO
    real(defReal), dimension(3)                :: centroid = ZERO
    real(defReal), dimension(3)                :: centroidTracks = ZERO

    real(defReal), dimension(:), allocatable   :: prevFlux
    real(defReal), dimension(:,:), allocatable :: fluxScores
    
    ! Source arrays
    real(defFlt), dimension(:), allocatable    :: fixedSource

    real(defReal), dimension(:), allocatable   :: prevX
    real(defReal), dimension(:), allocatable   :: prevY
    real(defReal), dimension(:), allocatable   :: prevZ
    real(defFlt), dimension(:), allocatable    :: fixedX
    real(defFlt), dimension(:), allocatable    :: fixedY
    real(defFlt), dimension(:), allocatable    :: fixedZ
    
    real(defReal), dimension(:,:), allocatable :: xScores
    real(defReal), dimension(:,:), allocatable :: yScores
    real(defReal), dimension(:,:), allocatable :: zScores
    
    real(defReal)     :: volume = ZERO
    real(defReal)     :: allVolumeTracks = ZERO
    real(defReal), dimension(3)   :: cellPos = -INFINITY
    
    ! Components
    integer(shortInt) :: nG    = 0
    integer(shortInt) :: ani   = 0
    logical(defBool)  :: set2D = .false.
    
  contains
    
    ! Public procedures
    procedure :: init
    procedure :: kill
    procedure :: initialiseFixedSource
    procedure :: initAdjoint

    ! Access procedures
    procedure :: getSource
    procedure :: getFixedSource
    procedure :: getPrevFlux
    procedure :: getFluxScore
    procedure :: getFluxSD
    procedure :: getVolume
    procedure :: getCellPos
    procedure :: wasHit
    procedure :: wasFound
    procedure :: hasFixedSource
    
    procedure :: getCentroid
    procedure :: getMomentMatrix
    procedure :: getFluxMoments
    procedure :: getFluxMomentSDs

    procedure :: getSourcePointer
    procedure :: getFluxPointer
    procedure :: getSourceMomentPointers
    procedure :: getFluxMomentPointers

    ! Predominantly for use in the transport sweep
    procedure :: incrementVolume
    procedure :: incrementCentroid
    procedure :: incrementMoments
    procedure :: hitCell
    procedure :: unhitCell
    procedure :: newFound
    procedure :: setLock
    procedure :: unsetLock

    ! Pre- and post-sweep procedures
    procedure :: scaleVolume
    procedure :: scaleCentroidAndMomMat
    procedure :: resetFluxes
    procedure :: accumulateFluxScores
    procedure :: finaliseFluxScores
    procedure :: zeroPrevFlux

    ! Private procedures
    procedure :: invertMatrix

  end type cellRR

contains

  !!
  !! Initialise the arrays object
  !!
  !! The object is fed sizes and requirements by the physics package.
  !! This will allocate the necessary arrays
  !!
  subroutine init(self, nG, lin, ani, set2D)
    class(cellRR), intent(inout)  :: self
    integer(shortInt), intent(in) :: nG       
    logical(defBool), intent(in)  :: lin
    integer(shortInt), intent(in) :: ani
    logical(defBool), intent(in)  :: set2D
    character(100), parameter :: Here = 'init (cellRR_class.f90)'

    self % nG = nG
    self % set2D = set2D

    ! Allocate and initialise arrays
    allocate(self % source(self % nG))
    allocate(self % scalarFlux(self % nG))
    allocate(self % prevFlux(self % nG))
    allocate(self % fluxScores(2, self % nG))
    
    self % source        = 0.0_defFlt
    self % scalarFlux    = ZERO
    self % prevFlux      = ONE
    self % fluxScores    = ZERO
    self % volumeTracks  = ZERO
    self % allVolumeTracks = ZERO
    self % volume        = ZERO
    self % cellHit       = .false.
    self % cellTotalHit  = 0
    self % cellPos       = -INFINITY

    ! Allocate linear components, if present
    if (lin) then
    
      allocate(self % sourceX(self % nG))
      allocate(self % sourceY(self % nG))
      allocate(self % sourceZ(self % nG))
      allocate(self % scalarX(self % nG))
      allocate(self % scalarY(self % nG))
      allocate(self % scalarZ(self % nG))
      allocate(self % prevX(self % nG))
      allocate(self % prevY(self % nG))
      allocate(self % prevZ(self % nG))
      allocate(self % xScores(2, self % nG))
      allocate(self % yScores(2, self % nG))
      allocate(self % zScores(2, self % nG))
      
      self % scalarX        = ZERO
      self % scalarY        = ZERO
      self % scalarZ        = ZERO
      self % prevX          = ZERO
      self % prevY          = ZERO
      self % prevZ          = ZERO
      self % sourceX        = 0.0_defFlt
      self % sourceY        = 0.0_defFlt
      self % sourceZ        = 0.0_defFlt
      self % momMat         = ZERO
      self % momTracks      = ZERO
      self % centroid       = ZERO
      self % centroidTracks = ZERO
      self % xScores        = ZERO
      self % yScores        = ZERO
      self % zScores        = ZERO

    end if

    ! TODO: allocate anisotropic components, if present
    if (ani > 0) then
      call fatalError(Here, 'Anisotropic scattering not yet supported')
    end if
    
    ! Initialise OMP locks
#ifdef _OPENMP
    call OMP_init_lock(self % lock)
#endif

  end subroutine init
  
  !!
  !! Initialise the adjoint source and update nuclear data.
  !! For now, assumes the adjoint is for global variance reduction.
  !!
  subroutine initAdjoint(self)
    class(cellRR), intent(inout)     :: self
    logical(defBool)                 :: doLinear
    integer(shortInt)                :: g
    real(defFlt), dimension(matSize) :: invM
    real(defFlt), save               :: xMom, yMom, zMom

    doLinear = .false.
    invM = 0.0_defFlt

    if (.not. allocated(self % fixedSource)) then
      allocate(self % fixedSource(self % nG))
    end if
    self % fixedSource = 0.0_defFlt

    if (allocated(self % scalarX)) then
      doLinear = .true.
      if (.not. allocated(self % fixedX)) then
        allocate(self % fixedX(self % nG))
        allocate(self % fixedY(self % nG))
        allocate(self % fixedZ(self % nG))
      end if
      self % fixedX = 0.0_defFlt
      self % fixedY = 0.0_defFlt
      self % fixedZ = 0.0_defFlt
    end if 

    ! Create fixed source from the flux scores
    ! Presently assumes the response of interest is global flux
    if (.not. self % wasFound()) return
    if (doLinear) invM = self % invertMatrix()

    do g = 1, self % nG

      ! Check for inordinately small flux values.
      ! Note, these can have arbitrarily low magnitude.
      ! Maybe should be something more robust.
      if (self % fluxScores(1, g) == ZERO) return
      self % fixedSource(g) = real(ONE / self % fluxScores(1, g), defFlt)

      ! Linear source treatment relies on performing a Taylor expansion
      ! of q' = 1/phi = 1/(phi_0 + <gradPhi , (r - r0)>) 
      ! = 1/phi_0 - 1/phi^2_0 * <gradPhi , (r - r0)>
      if (doLinear .and. (self % fixedSource(g) > 0)) then
        self % fixedX(g) = invM(xx) * xMom + invM(xy) * yMom + invM(xz) * zMom 
        self % fixedY(g) = invM(xy) * xMom + invM(yy) * yMom + invM(yz) * zMom 
        self % fixedZ(g) = invM(xz) * xMom + invM(yz) * yMom + invM(zz) * zMom 

        self % fixedX(g) = -self % fixedX(g) * self % fixedSource(g) ** 2
        self % fixedY(g) = -self % fixedY(g) * self % fixedSource(g) ** 2
        self % fixedZ(g) = -self % fixedZ(g) * self % fixedSource(g) ** 2
          
      end if

    end do
    
    ! Reinitialise arrays to be used during transport
    self % scalarFlux      = ZERO
    self % prevFlux        = ZERO
    self % fluxScores      = ZERO
    self % source          = 0.0_defFlt
    
    ! Ideally we would have a way of reusing the volume estimators
    ! No compact ideas at the moment, so these will simply be reinitialised
    self % volumeTracks    = ZERO
    self % allVolumeTracks = ZERO
    self % volume          = ZERO

    if (doLinear) then
      self % scalarX        = ZERO
      self % scalarY        = ZERO
      self % scalarZ        = ZERO
      self % prevX          = ZERO
      self % prevY          = ZERO
      self % prevZ          = ZERO
      self % sourceX        = 0.0_defFlt
      self % sourceY        = 0.0_defFlt
      self % sourceZ        = 0.0_defFlt
      self % momMat         = ZERO
      self % momTracks      = ZERO
      self % centroid       = ZERO
      self % centroidTracks = ZERO
      self % xScores        = ZERO
      self % yScores        = ZERO
      self % zScores        = ZERO
    end if

  end subroutine initAdjoint

  !!
  !! Initialises fixed source.
  !!
  subroutine initialiseFixedSource(self, strength)
    class(cellRR), intent(inout)               :: self
    real(defReal), dimension(:)                  :: strength
    character(100), parameter :: Here = 'initialiseFixedSource (cellRR_class.f90)'

    allocate(self % fixedSource(self % nG))

    if (size(strength) /= self % nG) call fatalError(Here,'Source strength is not size nG')

    self % fixedSource = real(strength, defFlt)

  end subroutine initialiseFixedSource

  !!
  !! Return source value given group
  !!
  elemental function getSource(self, g) result(src)
    class(cellRR), intent(in)     :: self
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: src

    src = self % source(g)

  end function getSource
  
  !!
  !! Return fixed source value given group
  !!
  elemental function getFixedSource(self, g) result(src)
    class(cellRR), intent(in)     :: self
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: src

    if (allocated(self % fixedSource)) then
      src = self % fixedSource(g)
    else
      src = 0.0_defFlt
    end if

  end function getFixedSource

  !!
  !! Return pointer to source
  !!
  function getSourcePointer(self) result(src)
    class(cellRR), intent(in), target   :: self
    real(defFlt), dimension(:), pointer, contiguous :: src

    src => self % source

  end function getSourcePointer
  
  !!
  !! Return pointer to flux
  !!
  function getFluxPointer(self) result(flx)
    class(cellRR), intent(in), target    :: self
    real(defReal), dimension(:), pointer, contiguous :: flx

    flx => self % scalarFlux

  end function getFluxPointer
  
  !!
  !! Return pointers to source moments
  !!
  subroutine getSourceMomentPointers(self, src, srcX, srcY, srcZ)
    class(cellRR), intent(in), target :: self
    real(defFlt), dimension(:), pointer, contiguous, intent(out) :: src, srcX, srcY, srcZ

    src => self % source
    srcX => self % sourceX
    srcY => self % sourceY
    srcZ => self % sourceZ

  end subroutine getSourceMomentPointers
  
  !!
  !! Return pointers to flux moments
  !!
  subroutine getFluxMomentPointers(self, flx, flxX, flxY, flxZ)
    class(cellRR), intent(in), target :: self
    real(defReal), dimension(:), pointer, contiguous, intent(out) :: flx, flxX, flxY, flxZ

    flx  => self % scalarFlux
    flxX => self % scalarX
    flxY => self % scalarY
    flxZ => self % scalarZ

  end subroutine getFluxMomentPointers
  
  !!
  !! Return previous flux value given group
  !!
  elemental function getPrevFlux(self, g) result(flux)
    class(cellRR), intent(in)     :: self
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % prevFlux(g)

  end function getPrevFlux
  
  !!
  !! Return final flux value given group
  !!
  elemental function getFluxScore(self, g) result(flux)
    class(cellRR), intent(in)     :: self
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % fluxScores(1, g)

  end function getFluxScore
  
  !!
  !! Return final flux standard deviation given group.
  !! Will return square of flux scores if called before finaliseFluxScores
  !!
  elemental function getFluxSD(self, g) result(flux)
    class(cellRR), intent(in)     :: self
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % fluxScores(2, g)

  end function getFluxSD
  
  !!
  !! Return final flux moment values given group
  !!
  pure function getFluxMoments(self, g) result(flux)
    class(cellRR), intent(in)     :: self
    integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: flux

    flux(1) = self % xScores(1, g)
    flux(2) = self % yScores(1, g)
    flux(3) = self % zScores(1, g)

  end function getFluxMoments
  
  !!
  !! Return final flux moment standard deviations given group
  !! Will return square of moment scores if called before finaliseFluxScores
  !!
  pure function getFluxMomentSDs(self, g) result(fluxSD)
    class(cellRR), intent(in)     :: self
    integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: fluxSD

    fluxSD(1) = self % xScores(2, g)
    fluxSD(2) = self % yScores(2, g)
    fluxSD(3) = self % zScores(2, g)

  end function getFluxMomentSDs
  
  !!
  !! Return volume
  !!
  pure function getVolume(self) result(vol)
    class(cellRR), intent(in)   :: self
    real(defReal)               :: vol

    vol = self % volume

  end function getVolume

  !!
  !! Return cell position
  !!
  pure function getCellPos(self) result(pos)
    class(cellRR), intent(in)      :: self
    real(defReal), dimension(nDim) :: pos

    pos = self % cellPos(1:nDim)

  end function getCellPos
  
  !!
  !! Return cell centroid
  !!
  pure function getCentroid(self) result(cent)
    class(cellRR), intent(in)      :: self
    real(defReal), dimension(nDim) :: cent

    cent = self % centroid

  end function getCentroid
  
  !!
  !! Return moment matrix
  !!
  pure function getMomentMatrix(self) result(mat)
    class(cellRR), intent(in)         :: self
    real(defReal), dimension(matSize) :: mat

    mat = self % momMat

  end function getMomentMatrix
  
  !!
  !! Increment the local volume estimate.
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementVolume(self, length)
    class(cellRR), intent(inout) :: self
    real(defReal), intent(in)    :: length     
    
    self % volumeTracks = self % volumeTracks + length
  
  end subroutine incrementVolume
  
  !!
  !! Increment the local centroid estimate.
  !! rL is the tracklength-weighted centroid
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementCentroid(self, rL)
    class(cellRR), intent(inout)               :: self
    real(defReal), dimension(nDim), intent(in) :: rL

    self % centroidTracks = self % centroidTracks + rL

  end subroutine incrementCentroid
  
  !!
  !! Increment the local moment matrix estimate.
  !! mat is the tracklength-weighted matrix
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementMoments(self, mat)
    class(cellRR), intent(inout)                  :: self
    real(defReal), dimension(matSize), intent(in) :: mat

    self % momTracks = self % momTracks + mat
  
  end subroutine incrementMoments
 
  !!
  !! Check if a cell has been hit
  !!
  elemental function wasHit(self) result (hit)
    class(cellRR), intent(in) :: self
    logical(defBool)          :: hit
    
    hit = self % cellHit
  
  end function wasHit
  
  !!
  !! Hit a cell.
  !! Should only be called in a lock.
  !!
  subroutine hitCell(self)
    class(cellRR), intent(inout) :: self
    
    self % cellHit = .true.
    self % cellTotalHit = self % cellTotalHit + 1
  
  end subroutine hitCell

  !! 
  !! Undo cell hit
  !!
  subroutine unhitCell(self)
    class(cellRR), intent(inout) :: self

    self % cellHit = .false.

  end subroutine unhitCell
  
  !!
  !! Has a cell ever been found?
  !!
  elemental function wasFound(self) result(found)
    class(cellRR), intent(in) :: self
    logical(defBool)          :: found
    
    found = self % cellFound
  
  end function wasFound

  !!
  !! Note that a new cell has been found, saving its position
  !!
  subroutine newFound(self, r)
    class(cellRR), intent(inout)            :: self
    real(defReal), dimension(3), intent(in) :: r

    ! Remove critical if this is to go in the lock
    !$omp critical 
    self % cellPos(:) = r
    self % cellFound = .true.
    !$omp end critical

  end subroutine newFound

  !!
  !! Check if a cell has an inhomogeneous source
  !!
  pure function hasFixedSource(self) result (hasSrc)
    class(cellRR), intent(in)   :: self
    logical(defBool)            :: hasSrc
    
    if (allocated(self % fixedSource)) then
      ! Take an absolute value in case of (possibly desirable?) negative sources
      hasSrc = any(abs(self % fixedSource) > 0.0_defFlt)
    else
      hasSrc = .false.
    end if
  
  end function hasFixedSource

  !!
  !! Set the OMP lock
  !!
  subroutine setLock(self) 
    class(cellRR), intent(inout) :: self

#ifdef _OPENMP
    call OMP_set_lock(self % lock)
#endif

  end subroutine setLock
  
  !!
  !! Unset the OMP lock
  !!
  subroutine unsetLock(self) 
    class(cellRR), intent(inout) :: self

#ifdef _OPENMP
    call OMP_unset_lock(self % lock)
#endif

  end subroutine unsetLock

  !!
  !! Computes cell volume given total track lengths.
  !! Norm is the inverse of the active track length on one iteration.
  !! NormIt is the inverse of the active track length over all iterations.
  !!
  subroutine scaleVolume(self, norm, normIt, volAve, volNaive)
    class(cellRR), intent(inout) :: self
    real(defReal), intent(in)    :: norm
    real(defReal), intent(in)    :: normIt
    real(defReal), intent(out)   :: volAve
    real(defReal), intent(out)   :: volNaive

    ! Actual integral volume
    self % allVolumeTracks = self % allVolumeTracks + self % volumeTracks
    self % volume = self % allVolumeTracks * normIt
      
    ! Simulation-average volume
    volAve = self % volume
      
    ! Cycle-wise/naive volume
    volNaive = self % volumeTracks * norm

  end subroutine scaleVolume
  
  !!
  !! Produce centroid and moment matrix by scaling the weighted estimators
  !! by the total track length in the cell.
  !!
  subroutine scaleCentroidAndMomMat(self)
    class(cellRR), intent(inout) :: self
    real(defReal)                :: invVol

    if (self % allVolumeTracks > ZERO) then
      
      invVol = ONE / self % allVolumeTracks
        
      ! Update centroids
      self % centroid(x) = self % centroidTracks(x) * invVol
      self % centroid(y) = self % centroidTracks(y) * invVol
      self % centroid(z) = self % centroidTracks(z) * invVol
      
      ! Update spatial moments
      self % momMat(xx) = self % momTracks(xx) * invVol
      self % momMat(xy) = self % momTracks(xy) * invVol
      self % momMat(xz) = self % momTracks(xz) * invVol
      self % momMat(yy) = self % momTracks(yy) * invVol
      self % momMat(yz) = self % momTracks(yz) * invVol
      self % momMat(zz) = self % momTracks(zz) * invVol

    end if

  end subroutine scaleCentroidAndMomMat

  !!
  !! Inverts the spatial moment matrix for use in linear source calculations.
  !!
  pure function invertMatrix(self) result(invM)
    class(cellRR), intent(in)         :: self
    real(defFlt), dimension(matSize)  :: invM
    real(defReal), dimension(matSize) :: momVec
    integer(shortInt)                 :: condX, condY, condZ, inversionTest
    real(defReal)                     :: det
    real(defFlt)                      :: one_det  

    momVec = self % momMat

    ! Pre-invert the moment matrix
    ! Need to check for poor conditioning by evaluating the
    ! diagonal elements of the matrix
    condX = 0
    condY = 0
    condZ = 0

    ! Trying out simpler matrix test
    if (momVec(xx) > condition_tolerance) condX = 1
    if (momVec(yy) > condition_tolerance) condY = 1
    if (momVec(zz) > condition_tolerance) condZ = 1
    ! Significantly stabilises 2D linear source problems. Z moments can vary wildly.
    if (self % set2D) condZ = 0

    ! Map conditions to test variable
    inversionTest = condX * 4 + condY * 2 + condZ
    invM = 0.0_defFlt

    select case(inversionTest)
    case(invertXYZ)
      det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
            - momVec(yy) * momVec(xz) * momVec(xz) &
            - momVec(zz) * momVec(xy) * momVec(xy) &
            + 2 * momVec(xy) * momVec(xz) * momVec(yz)
      invM(xx) = real(momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz),defFlt)
      invM(xy) = real(momVec(xz) * momVec(yz) - momVec(xy) * momVec(zz),defFlt)
      invM(xz) = real(momVec(xy) * momVec(yz) - momVec(yy) * momVec(xz),defFlt)
      invM(yy) = real(momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz),defFlt)
      invM(yz) = real(momVec(xy) * momVec(xz) - momVec(xx) * momVec(yz),defFlt)
      invM(zz) = real(momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy),defFlt)

    case(invertYZ)
      det = momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)
      invM(yy) = real(momVec(zz),defFlt)
      invM(yz) = real(-momVec(yz),defFlt)
      invM(zz) = real(momVec(yy),defFlt)

    case(invertXY)
      det = momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)
      invM(xx) = real(momVec(yy),defFlt)
      invM(xy) = real(-momVec(xy),defFlt)
      invM(yy) = real(momVec(xx),defFlt)

    case(invertXZ)
      det = momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)
      invM(xx) = real(momVec(zz),defFlt)
      invM(xz) = real(-momVec(xz),defFlt)
      invM(zz) = real(momVec(xx),defFlt)

    case(invertX)
      det = momVec(xx)
      invM(xx) = 1.0_defFlt

    case(invertY)
      det = momVec(yy)
      invM(yy) = 1.0_defFlt

    case(invertZ)
      det = momVec(zz)
      invM(zz) = 1.0_defFlt

    case default
      det = ONE
    end select

    one_det = real(ONE/det, defFlt)
    invM = invM * one_det
    
    ! Check for zero determinant
    if (abs(det) < det_tolerance) invM = 0.0_defFlt

  end function invertMatrix

  !!
  !! Zero the previous-step flux
  !!
  subroutine zeroPrevFlux(self)
    class(cellRR), intent(inout) :: self
    character(100), parameter    :: Here = 'zeroPrevFlux (cellRR_class.f90)'

    if (allocated(self % prevFlux)) then
      self % prevFlux = ZERO
    else
      call fatalError(Here,'prevFlux has not been initialised')
    end if

  end subroutine zeroPrevFlux
  
  !!
  !! Reset fluxes
  !!
  subroutine resetFluxes(self)
    class(cellRR), intent(inout) :: self
    integer(shortInt)            :: idx
    character(100), parameter    :: Here = 'resetFluxes (cellRR_class.f90)'

    !$omp simd
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = ZERO
    end do
    
    if (allocated(self % scalarX)) then
      !$omp simd
      do idx = 1, size(self % scalarFlux)
        self % prevX(idx) = self % scalarX(idx)
        self % scalarX(idx) = ZERO
        self % prevY(idx) = self % scalarY(idx)
        self % scalarY(idx) = ZERO
        self % prevZ(idx) = self % scalarZ(idx)
        self % scalarZ(idx) = ZERO
      end do
    end if

  end subroutine resetFluxes
  
  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScores(self)
    class(cellRR), intent(inout) :: self
    real(defReal)                :: flux
    integer(shortInt)            :: idx

    !$omp simd
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarFlux(idx)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) + flux
      self % fluxScores(2, idx) = self % fluxScores(2, idx) + flux * flux
    end do

    if (allocated(self % scalarX)) then
      !$omp simd
      do idx = 1, size(self % scalarFlux)
        flux = self % scalarX(idx)
        self % xScores(1, idx) = self % xScores(1, idx) + flux
        self % xScores(2, idx) = self % xScores(2, idx) + flux * flux
      end do
    
      !$omp simd
      do idx = 1, size(self % scalarFlux)
        flux = self % scalarY(idx)
        self % yScores(1, idx) = self % yScores(1, idx) + flux
        self % yScores(2, idx) = self % yScores(2, idx) + flux * flux
      end do
    
      !$omp simd
      do idx = 1, size(self % scalarFlux)
        flux = self % scalarZ(idx)
        self % zScores(1, idx) = self % zScores(1, idx) + flux
        self % zScores(2, idx) = self % zScores(2, idx) + flux * flux
      end do

    end if

  end subroutine accumulateFluxScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScores(self,it)
    class(cellRR), intent(inout)  :: self
    integer(shortInt), intent(in) :: it
    integer(shortInt)             :: idx
    real(defReal)                 :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !$omp simd
    do idx = 1, size(self % scalarFlux)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) * N1
      self % fluxScores(2, idx) = self % fluxScores(2, idx) * N1
      self % fluxScores(2, idx) = Nm1 * (self % fluxScores(2, idx) - &
            self % fluxScores(1, idx) * self % fluxScores(1, idx)) 
      if (self % fluxScores(2, idx) <= ZERO) then
        self % fluxScores(2, idx) = ZERO
      else
        self % fluxScores(2, idx) = sqrt(self % fluxScores(2, idx))
      end if
    end do

    if (allocated(self % scalarX)) then
      do idx = 1, size(self % scalarFlux)
        self % xScores(1, idx) = self % xScores(1, idx) * N1
        self % xScores(2, idx) = self % xScores(2, idx) * N1
        self % xScores(2, idx) = Nm1 * (self % xScores(2, idx) - &
              self % xScores(1, idx) * self % xScores(1, idx)) 
        if (self % xScores(2, idx) <= ZERO) then
          self % xScores(2, idx) = ZERO
        else
          self % xScores(2, idx) = sqrt(self % xScores(2, idx))
        end if
      
        self % yScores(1, idx) = self % yScores(1, idx) * N1
        self % yScores(2, idx) = self % yScores(2, idx) * N1
        self % yScores(2, idx) = Nm1 * (self % yScores(2, idx) - &
              self % yScores(1, idx) * self % yScores(1, idx)) 
        if (self % yScores(2, idx) <= ZERO) then
          self % yScores(2, idx) = ZERO
        else
          self % yScores(2, idx) = sqrt(self % yScores(2, idx))
        end if
      
        self % zScores(1, idx) = self % zScores(1, idx) * N1
        self % zScores(2, idx) = self % zScores(2, idx) * N1
        self % zScores(2, idx) = Nm1 * (self % zScores(2, idx) - &
              self % zScores(1, idx) * self % zScores(1, idx)) 
        if (self % zScores(2, idx) <= ZERO) then
          self % zScores(2, idx) = ZERO
        else
          self % zScores(2, idx) = sqrt(self % zScores(2, idx))
        end if
      end do
    end if

  end subroutine finaliseFluxScores
  
  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(cellRR), intent(inout) :: self

    ! Clean standard contents
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % fixedSource)) deallocate(self % fixedSource)
    
    ! Clean LS contents
    if(allocated(self % scalarX)) deallocate(self % scalarX)
    if(allocated(self % scalarX)) deallocate(self % scalarY)
    if(allocated(self % scalarX)) deallocate(self % scalarZ)
    if(allocated(self % prevX)) deallocate(self % prevX)
    if(allocated(self % prevY)) deallocate(self % prevY)
    if(allocated(self % prevZ)) deallocate(self % prevZ)
    if(allocated(self % sourceX)) deallocate(self % sourceX)
    if(allocated(self % sourceY)) deallocate(self % sourceY)
    if(allocated(self % sourceZ)) deallocate(self % sourceZ)
    if(allocated(self % xScores)) deallocate(self % xScores)
    if(allocated(self % yScores)) deallocate(self % yScores)
    if(allocated(self % zScores)) deallocate(self % zScores)

#ifdef _OPENMP
    call OMP_destroy_lock(self % lock)
#endif
    
    self % nG  = 0
    self % ani = 0
    self % volume = ZERO
    self % volumeTracks = ZERO
    self % allVolumeTracks = ZERO
    self % centroid = ZERO
    self % centroidTracks = ZERO
    self % momMat = ZERO
    self % momTracks = ZERO
    self % cellPos = -INFINITY
    self % cellHit = .false.
    self % cellTotalHit = 0
    self % cellFound = .false.

  end subroutine kill

end module cellRR_class
