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
    procedure :: getFluxPointer
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
    class(fluxRR), intent(inout) :: self
    class(mgNeutronDatabase), pointer, intent(inout)    :: db
    class(geometryStd), pointer, intent(inout)          :: geom
    integer(shortInt), intent(in)                       :: nG
    integer(shortInt), intent(in)                       :: nCells
    real(defReal), intent(in)                           :: rho
    logical(defBool), intent(in)                        :: lin
    integer(shortInt), intent(in)                       :: ani
    class(dictionary), pointer, intent(inout), optional :: dictFS
    character(100), parameter :: Here = 'init (fluxRR_class.f90)'

    self % mgData => db
    self % geom   => geom

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

    end if

    ! TODO: allocate linear and anisotropic components, if present
    
    ! Initialise OMP locks
    allocate(self % locks(self % nCells))
    do i = 1, self % nCells
      call omp_init_lock(self % locks(i))
    end do

  end subroutine init

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolumeFlatIso(self, it)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    real(defFlt)                   :: norm
    real(defReal)                  :: normVol
    real(defFlt), save             :: total, vol, sigGG, D
    integer(shortInt), save        :: g, matIdx, idx
    integer(shortInt)              :: cIdx
    !$omp threadprivate(total, vol, idx, g, matIdx, sigGG, D)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = real(self % volume(cIdx), defFlt)

      do g = 1, self % nG

        total = self % XSdata % getTotalXS(matIdx, g)
        idx   = self % nG * (cIdx - 1) + g

        if (vol > volume_tolerance) then
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm / ( total * vol)
        end if
        
        ! Apply stabilisation for negative XSs
        if (matIdx < VOID_MAT) then
          sigGG = self % XSdata % getSigmaS(matIdx, g, g)

          ! Presumes non-zero total XS
          if ((sigGG < 0) .and. (total > 0)) then
            D = -real(self % rho, defFlt) * sigGG / total
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
  subroutine sourceUpdateKernelFlatIso(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout) :: self
    integer(shortInt), intent(in)          :: cIdx
    real(defFlt), intent(in)               :: ONE_KEFF
    real(defFlt)                           :: scatter, fission
    real(defFlt), dimension(:), pointer    :: nuFission, total, chi, scatterXS 
    integer(shortInt)                      :: matIdx, g, gIn, baseIdx, idx
    real(defFlt), pointer, dimension(:)    :: fluxVec, scatterVec

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx >= VOID_MAT) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
      end do
      return
    end if

    ! Obtain XSs
    matIdx = (matIdx - 1) * self % nG
    total => self % sigmaT(matIdx + (1):(self % nG))
    scatterXS => self % sigmaS(matIdx * self % nG + (1):(self % nG*self % nG))
    nuFission => self % nuSigmaF(matIdx + (1):(self % nG))
    chi => self % chi(matIdx + (1):(self % nG))

    baseIdx = self % ng * (cIdx - 1)
    fluxVec => self % prevFlux(baseIdx+(1):(self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF

    do g = 1, self % nG

      scatterVec => scatterXS(self % nG * (g - 1) + (1):self % nG)

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
  !!
  function calculateKeff(self, k0) result(k1)
    class(arraysRR), intent(inout)              :: self
    real(defReal), intent(in)                   :: k0
    real(defReal)                               :: k1
    real(defFlt)                                :: fissionRate, prevFissionRate
    real(defFlt), save                          :: fissLocal, prevFissLocal, vol
    integer(shortInt), save                     :: matIdx, g, idx, mIdx
    integer(shortInt)                           :: cIdx
    class(baseMgNeutronMaterial), pointer, save :: mat
    class(materialHandle), pointer, save        :: matPtr
    !$omp threadprivate(mat, matPtr, fissLocal, prevFissLocal, matIdx, g, idx, mIdx, vol)

    fissionRate     = 0.0_defFlt
    prevFissionRate = 0.0_defFlt
    !$omp parallel do schedule(static) reduction(+: fissionRate, prevFissionRate)
    do cIdx = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      if (matIdx >= VOID_MAT) cycle

      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)
      if (.not. mat % isFissile()) cycle

      vol = real(self % volume(cIdx), defFlt)

      if (vol <= volume_tolerance) cycle

      fissLocal = 0.0_defFlt
      prevFissLocal = 0.0_defFlt
      mIdx = (matIdx - 1) * self % nG
      do g = 1, self % nG
        
        ! Source index
        idx = self % nG * (cIdx - 1) + g
        fissLocal     = fissLocal     + self % scalarFlux(idx) * self % nuSigmaF(mIdx + g)
        prevFissLocal = prevFissLocal + self % prevFlux(idx) * self % nuSigmaF(mIdx + g)

      end do

      fissionRate     = fissionRate     + fissLocal * vol
      prevFissionRate = prevFissionRate + prevFissLocal * vol

    end do
    !$omp end parallel do

    ! Update k
    k1 = k0 * fissionRate / prevFissionRate

  end subroutine calculateKeff
  
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
