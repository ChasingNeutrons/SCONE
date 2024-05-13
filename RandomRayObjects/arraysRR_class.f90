module arraysRR_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryStd_class,              only : geometryStd
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat            => nMat, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init         => init, &
                                             ndReg_getMatNames  => getMatNames, &
                                             ndReg_activate     => activate, &
                                             ndReg_kill         => kill, &
                                             ndReg_getNeutronMG => getNeutronMG
  use materialHandle_inter,           only : materialHandle
  use mgNeutronDatabase_inter,        only : mgNeutronDatabase
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class,    only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast
  
  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Object to store all arrays in random ray
  !! By default, will have the current flux, previous flux, accumulated flux,
  !! source arrays, and geometric arrays.
  !!
  !! Can be extended to having corresponding flux and source moments, as well
  !! as the fixed source.
  !!
  !! Private Members
  !!   nG          -> Number of energy groups, kept for convenience.
  !!   nCells      -> Number of unique cells in the geometry, kept for convenience.
  !!   lengthPerIt -> RR active length per iteration, kept for convenience
  !!   mgData      -> Pointer to nuclear data, for convenience.
  !!   geom        -> Pointer to geometry, for convenience.
  !!   cache       -> Logical check whether to use distance caching
  !!   rho         -> Stabilisation factor: 0 is no stabilisation, 1 is aggressive stabilisation
  !!   linear      -> Are linear flux moments stored?
  !!   aniso       -> Order of ansitropic flux moments to be stored
  !!   linAni      -> Are linear and anisotropic flux moments stored?
  !!   flux        -> Array of scalar flux values of length = nG * nCells
  !!   prevFlux    -> Array of previous scalar flux values of length = nG * nCells
  !!   fluxScore   -> Array of scalar flux values and squared values to be reported 
  !!                  in results, dimension =  [nG * nCells, 2]
  !!
  !!   source      -> Array of sources
  !!   fixedSource -> Array of fixed sources
  !!
  type, public :: arraysRR
    private
    ! Components
    class(baseMgNeutronDatabase), pointer :: mgData      => null()
    class(geometryStd), pointer           :: geom        => null()
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0
    real(defReal)                         :: lengthPerIt = ZERO
    real(defReal)                         :: rho         = ZERO
    logical(defBool)                      :: linear      = .false.
    logical(defBool)                      :: linAni      = .false.
    integer(shortInt)                     :: aniso       = 0
    
    ! Flux arrays
    real(defFlt), dimension(:), allocatable    :: scalarFlux
    real(defFlt), dimension(:), allocatable    :: prevFlux
    real(defReal), dimension(:,:), allocatable :: fluxScores
    
    ! Source arrays
    real(defFlt), dimension(:), allocatable    :: source
    real(defFlt), dimension(:), allocatable    :: fixedSource

    ! Geometry arrays
    real(defReal), dimension(:), allocatable     :: volTracks
    real(defReal), dimension(:), allocatable     :: volume
    integer(shortInt), dimension(:), allocatable :: cellHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos

  contains
    
    procedure :: init
    procedure :: initialiseFixedSource
    procedure :: getSourcePointer
    procedure :: getFluxPointer
    procedure :: incrementVolume
    procedure :: hasHit
    procedure :: hitCell
    procedure :: cellHitRate
    procedure :: wipeCellHits
    procedure :: found
    procedure :: newFound
    procedure :: getNG
    procedure :: setLock
    procedure :: unsetLock
    procedure :: resetFluxes
    procedure :: resetFluxesFlatIso
    procedure :: resetFluxesLinearIso
    procedure :: resetFluxesLinearAni
    procedure :: resetFluxesFlatAni
    procedure :: normaliseFluxAndVolume
    procedure :: normaliseFluxAndVolumeFlatIso
    procedure :: normaliseFluxAndVolumeLinearIso
    procedure :: normaliseFluxAndVolumeLinearAni
    procedure :: normaliseFluxAndVolumeFlatAni
    procedure :: updateSource
    procedure :: sourceUpdateKernelFlatIso
    procedure :: sourceUpdateKernelLinearIso
    procedure :: sourceUpdateKernelLinearAni
    procedure :: sourceUpdateKernelFlatAni
    procedure :: accumulateFluxScores
    procedure :: finaliseFluxAndKeffScores
    procedure :: calculateKeff
    procedure :: calculateKeffKernel
    procedure :: kill

    ! Private procedures

  end type arraysRR

contains

  !!
  !! Initialise the arrays object
  !!
  !! The object is fed sizes and requirements by the physics package.
  !! This will allocate the necessary arrays
  !!
  subroutine init(self, db, geom, nG, nCells, rho, lin, ani, dictFS)
    class(arraysRR), intent(inout) :: self
    class(mgNeutronDatabase), pointer, intent(in)       :: db
    class(geometryStd), pointer, intent(in)             :: geom
    integer(shortInt), intent(in)                       :: nG
    integer(shortInt), intent(in)                       :: nCells
    real(defReal), intent(in)                           :: rho
    logical(defBool), intent(in)                        :: lin
    integer(shortInt), intent(in)                       :: ani
    class(dictionary), pointer, intent(inout), optional :: dictFS
    character(100), parameter :: Here = 'init (arraysRR_class.f90)'

    self % nG          = nG
    self % nCells      = nCells
    self % nEl         = nG * nCells
    self % lengthPerIt = lengthPerIt

    self % rho = rho
    self % lin = lin
    if (ani >= 0) self % ani = ani
    if (self % lin .and. self % ani > 0) self % linAni = .true.

    ! Allocate and initialise arrays
    allocate(self % scalarFlux(self % nEl))
    allocate(self % prevFlux(self % nEl))
    allocate(self % fluxScores(self % nEl, 2))
    allocate(self % source(self % nEl))
    allocate(self % volTracks(self % nCells))
    allocate(self % volume(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))
    
    self % scalarFlux = 0.0_defFlt
    self % prevFlux   = 1.0_defFlt
    self % fluxScores = ZERO
    self % source     = 0.0_defFlt
    self % volTracks  = ZERO
    self % volume     = ZERO
    self % cellHit    = 0
    self % cellFound  = .false.
    self % cellPos    = -INFINITY

    ! Initialise the fixed source if present
    if (present(dictFS)) then
      allocate(self % fixedSource(self % nEl))
      call self % initialiseFixedSource(dictFS)
    end if

    ! TODO: allocate linear and anisotropic components, if present
    
    ! Initialise OMP locks
    allocate(self % locks(self % nCells))
    do i = 1, self % nCells
      call omp_init_lock(self % locks(i))
    end do

  end subroutine init

  !!
  !! Return a pointer to the flux vector for a given cell
  !!
  subroutine getFluxPointer(self, cIdx, fluxVec)
    class(arraysRR), intent(in)         :: self
    integer(shortInt), intent(in)       :: cIdx
    real(defFlt), dimension(:), pointer :: fluxVec
    integer(shortInt)                   :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    fluxVec => self % scalarFlux(baseIdx1:baseIdx2)

  end subroutine getFluxPointer
  
  !!
  !! Increment the local volume estimate in cell cIdx
  !!
  subroutine incrementVolume(self, cIdx, length)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    real(defReal), intent(in)     :: length     
    
    self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
  
  end subroutine incrementVolume
  
  !!
  !! Check if a cell has been hit
  !!
  function hasHit(self, cIdx) result (hit)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt)             :: hit
    
    hit = self % cellHit(cIdx)
  
  end subroutine hasHit
  
  !!
  !! Hit a cell 
  !!
  subroutine hitCell(self, cIdx)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    
    self % cellHit(cIdx) = 1
  
  end subroutine hitCell

  !!
  !! Return the cell hit rate for the given iteration
  !!
  function cellHitRate(self) result(hitRate)
    class(arraysRR), intent(in)  :: self
    integer(shortInt)            :: cIdx, totalHit
    real(defReal)                :: hitRate

    totalHit = 0
    !$omp parallel do reduction(+:totalHit)
    do cIdx = 1, self % nCells
      totalHit = totalHit + self % cellHit(cIdx)
    end do
    !$omp end parallel do
    hitRate = real(hitRate / self % nCells, defReal)

  end function cellHitRate

  !! 
  !! Wipe cell hits
  !!
  subroutine wipeCellHits(self)
    class(arraysRR), intent(in)  :: self
    integer(shortInt)            :: cIdx

    !$omp parallel do
    do cIdx = 1, self % nCells
      self % cellHit(cIdx) = 0
    end do
    !$omp end parallel do

  end subroutine wipeCellHits
  
  !!
  !! Has a cell ever been found?
  !!
  function found(self, cIdx) result(wasFound)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: wasFound
    
    wasFound = self % cellFound(cIdx)
  
  end function found

  !!
  !! Note that a new cell has been found
  !!
  subroutine newFound(self, cIdx, r)
    class(arraysRR), intent(inout)          :: self
    integer(shortInt), intent(in)           :: cIdx
    real(defReal), dimension(3), intent(in) :: r

    !$omp critical 
    self % cellFound(cIdx) = .true.
    self % cellPos(cIdx,:) = r
    !$omp end critical

  end subroutine newFound

  !!
  !! Return number of energy groups used
  !!
  function getNG(self) result(nG)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: nG

    nG = self % nG

  end function getNG

  !!
  !! Return a pointer to the source vector for a given cell
  !!
  subroutine getSourcePointer(self, cIdx, sourceVec)
    class(arraysRR), intent(in)         :: self
    integer(shortInt), intent(in)       :: cIdx
    real(defFlt), dimension(:), pointer :: sourceVec
    integer(shortInt)                   :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    sourceVec => self % source(baseIdx1:baseIdx2)

  end subroutine getSourcePointer

  !!
  !! Set the OMP lock in a given cell
  !!
  subroutine setLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

    call OMP_set_lock(self % locks(cIdx))

  end subroutine setLock
  
  !!
  !! Unset the OMP lock in a given cell
  !!
  subroutine unsetLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

    call OMP_unset_lock(self % locks(cIdx))

  end subroutine unsetLock


  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolumeFlatIso(self, XSData, geom, it)
    class(arraysRR), intent(inout)            :: self
    class(dataRR), pointer, intent(in)        :: XSData
    class(geometrdyStd), pointer, intent(in)  :: geom
    integer(shortInt), intent(in)             :: it
    real(defFlt)                              :: norm
    real(defReal)                             :: normVol
    real(defFlt), save                        :: total, vol, sigGG, D
    real(defFlt), dimension(:), pointer, save :: total
    integer(shortInt), save                   :: g, matIdx, idx
    integer(shortInt)                         :: cIdx
    !$omp threadprivate(total, vol, idx, g, matIdx, sigGG, D)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  geom % geom % graph % getMatFromUID(cIdx) 
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = real(self % volume(cIdx), defFlt)

      call XSData % getTotalPointer(matIdx, total)

      do g = 1, self % nG

        idx   = self % nG * (cIdx - 1) + g

        if (vol > volume_tolerance) then
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm / ( total(g) * vol)
        end if
        
        ! Apply stabilisation for negative XSs
        if (matIdx < VOID_MAT) then
          sigGG = XSData % getSigmaS(matIdx, g, g)

          ! Presumes non-zero total XS
          if ((sigGG < 0) .and. (total > 0)) then
            D = -real(self % rho, defFlt) * sigGG / total(g)
          else
            D = 0.0_defFlt
          end if
        else
          D = 0.0_defFlt
        end if

        self % scalarFlux(idx) =  (self % scalarFlux(idx) + self % source(idx) + D * self % prevFlux(idx) ) / (1 + D)

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeFlatIso
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernelFlatIso(self, cIdx, ONE_KEFF, XSData, geom)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    class(dataRR), pointer, intent(in)       :: XSData
    class(geometrdyStd), pointer, intent(in) :: geom
    real(defReal), intent(in)                :: ONE_KEFF
    real(defFlt)                             :: scatter, fission, ONE_K
    real(defFlt), dimension(:), pointer      :: nuFission, total, chi, scatterXS, &
                                                scatterVec, fluxVec 
    integer(shortInt)                        :: matIdx, g, gIn, baseIdx, idx, sIdx1, sIdx2

    ONE_K = real(ONE_KEFF,defFlt)
    ! Identify material
    matIdx  =  geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx >= VOID_MAT) then
      baseIdx = self % nG * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
      end do
      return
    end if

    ! Obtain XSs
    matIdx = (matIdx - 1) * self % nG
    call XSData % getPointers(matIdx, total, nuFission, scatterXS, chi)

    baseIdx = self % nG * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx + 1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_K

    do g = 1, self % nG

      sIdx1 = self % nG * (g - 1) + 1
      sIdx2 = self % nG * g
      scatterVec => scatterXS(sIdx1:sIdx2)

      ! Calculate scattering source
      scatter = 0.0_defFlt
      
      ! Sum contributions from all energies
      !$omp simd reduction(+:scatter)
      do gIn = 1, self % nG
        scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      self % source(idx) = chi(g) * fission + scatter
      if (allocated(self % fixedSource)) then
        self % source(idx) = self % source(idx) + self % fixedSource(idx)
      end if
      self % source(idx) = self % source(idx) / total(g)

    end do

  end subroutine sourceUpdateKernelFlatIso

  
  !! 
  !! Calculate keff
  !! Wraps the main kernel call to allow for OMP + SIMD (thanks Fortran)
  !!
  function calculateKeff(self, k0, XSData, geom) result(k1)
    class(arraysRR), intent(inout)           :: self
    real(defReal), intent(in)                :: k0
    class(dataRR), pointer, intent(in)       :: XSData
    class(geometrdyStd), pointer, intent(in) :: geom
    real(defReal)                            :: k1
    integer(shortInt)                        :: cIdx
    real(defReal), save                      :: fissLocal, prevFissLocal
    real(defReal)                            :: fissTotal, prevFissTotal
    !$omp threadprivate (fissLocal, prevFissLocal)

    !$omp parallel do schedule(dynamic) reduction(+: fissTotal, prevFissTotal
    do cIdx = 1, self % nCells
      call self % calculateKeffKernel(cIdx, XSData, geom, fissLocal, prevFissLocal)
      fissTotal     = fissTotal + fissLocal
      prevFissTotal = prevFissTotal + prevFissLocal
    end do 
    !$omp parallel end do

    k1 = k0 * fissTotal / prevFissTotal 

  end function calculateKeff
  
  !!
  !! Calculate keff for a single cell
  !!
  subroutine calculateKeffKernel(self, cIdx, XSData, geom, fissionRate, prevFissionRate)
    class(arraysRR), intent(inout)              :: self
    integer(shortInt), intent (in)              :: cIdx
    class(dataRR), pointer, intent(in)          :: XSData
    class(geometrdyStd), pointer, intent(in)    :: geom
    real(defReal), intent(out)                  :: fissionRate, prevFissionRate
    real(defReal)                               :: vol
    integer(shortInt)                           :: g, matIdx
    real(defFlt), dimension(:), pointer         :: nuFission, flux, prevFlux

    fissionRate     = ZERO
    prevFissionRate = ZERO

    ! Identify material
    matIdx =  geom % geom % graph % getMatFromUID(cIdx) 
      
    ! Check whether to continue in this cell
    if (matIdx >= VOID_MAT) return
    if (.not. XSData % isFissile(matIdx)) return
    vol = self % volume(cIdx)
    if (vol <= volume_tolerance) return

    call XSData % getNuFissPointer(matIdx, nuSigmaF)
    flux => self % scalarFlux((self % nG * (cIdx - 1) + 1):(self % nG * cIdx))
    prevFlux => self % prevFlux((self % nG * (cIdx - 1) + 1):(self % nG * cIdx))

    !$omp simd reduction (+: fissionRate, prevFissionRate)
    do g = 1, self % nG
      fissionRate     = fissionRate     + real(flux(g) * nuSigmaF(g), defReal)
      prevFissionRate = prevFissionRate + real(prevFlux(g) * nuSigmaF(g), defReal)
    end do

    fissionRate     = fissionRate * vol
    prevFissionRate = prevFissionRate * vol

  end subroutine calculateKeffKernel
  
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !!
  subroutine resetFluxesFlatIso(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do schedule(static)
    do idx = 1, self % nEl
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesFlatIso

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScores(self)
    class(arraysRR), intent(inout) :: self
    real(defReal), save            :: flux
    integer(shortInt)              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, self % nEl
      flux = real(self % scalarFlux(idx),defReal)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux * flux
    end do
    !$omp end parallel do

  end subroutine accumulateFluxScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScores(self,it)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)              :: idx
    real(defReal)                  :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !$omp parallel do schedule(static)
    do idx = 1, self % nEl
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) * N1
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) * N1
      self % fluxScores(idx,2) = Nm1 *(self % fluxScores(idx,2) - &
            self % fluxScores(idx,1) * self % fluxScores(idx,1)) 
      if (self % fluxScores(idx,2) <= ZERO) then
        self % fluxScores(idx,2) = ZERO
      else
        self % fluxScores(idx,2) = sqrt(self % fluxScores(idx,2))
      end if
    end do
    !$omp end parallel do

  end subroutine finaliseFluxScores
  
  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean contents
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % fixedSource)) deallocate(self % fixedSource)
    if(allocated(self % volTracks)) deallocate(self % volTracks)
    if(allocated(self % volumes)) deallocate(self % volumes)
    if(allocated(self % cellHit)) deallocate(self % cellHit)
    if(allocated(self % cellFound)) deallocate(self % cellFound)
    if(allocated(self % cellPos)) deallocate(self % cellPos)
    
    if(allocated(self % locks)) then
      do i = 1, self % nCells
        call OMP_destroy_lock(self % locks(i))
      end do
      deallocate(self % locks)
    end if
    
    self % geom   => null()
    self % mgData => null()
    self % nG     = 0
    self % nCells = 0
    self % lengthPerIt = ZERO
    self % rho         = ZERO

  end subroutine kill

end module arraysRR_class
