module fluxRR_class

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
  !! Object to store all flux-related information in random ray
  !! By default, will have the current flux, previous flux, and accumulated flux.
  !! Can be extended to having corresponding flux moments.
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
  type, public :: fluxRR
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
    
    ! Results space
    real(defFlt), dimension(:), allocatable    :: scalarFlux
    real(defFlt), dimension(:), allocatable    :: prevFlux
    real(defReal), dimension(:,:), allocatable :: fluxScores

  contains
    
    procedure :: init
    procedure :: resetFluxes
    procedure :: normaliseFluxAndVolume
    procedure :: accumulateFluxScores
    procedure :: finaliseFluxAndKeffScores
    procedure :: calculateKeff
    procedure :: getFluxPointer
    procedure :: kill

    ! Private procedures

  end type fluxRR

contains

  !!
  !! Initialise the flux object
  !!
  !! The object is fed sizes and requirements by the physics package.
  !! This will allocate the necessary arrays
  !!
  subroutine init(self, db, geom, nG, nCells, rho, lin, ani)
    class(fluxRR), intent(inout) :: self
    class(mgNeutronDatabase), pointer, intent(inout) :: db
    class(geometryStd), pointer, intent(inout)       :: geom
    integer(shortInt), intent(in)                    :: nG
    integer(shortInt), intent(in)                    :: nCells
    real(defReal), intent(in)                        :: rho
    logical(defBool), intent(in)                     :: lin
    integer(shortInt), intent(in)                    :: ani
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


    ! Allocate and initialise results space
    allocate(self % scalarFlux(self % nEl))
    allocate(self % prevFlux(self % nEl))
    allocate(self % fluxScores(self % nEl, 2))
    
    self % scalarFlux = 0.0_defFlt
    self % prevFlux   = 1.0_defFlt
    self % fluxScores = ZERO

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
  subroutine normaliseFluxAndVolume(self, it)
    class(fluxRR), intent(inout)   :: self
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
      call self % geomObjects % normaliseTracks(cIdx, normVol)
      vol = real(self % geomObjects % getVol(cIdx))

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

        self % scalarFlux(idx) =  (self % scalarflux(idx) + self % source(idx) + D * self % prevFlux(idx) ) / (1 + D)

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume
  
  !!
  !! Calculate keff
  !!
  function calculateKeff(self, k0) result(k1)
    class(fluxRR), intent(inout)                :: self
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

      vol = real(self % geomObjects % getVolume(cIdx), defFlt)

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
  subroutine resetFluxes(self)
    class(fluxRR), intent(inout) :: self
    integer(shortInt)            :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxes

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScores(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    real(defReal), save                           :: flux
    integer(shortInt)                             :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
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
    class(fluxRR), intent(inout)  :: self
    integer(shortInt), intent(in) :: it
    integer(shortInt)             :: idx
    real(defReal)                 :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
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
    class(fluxRR), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean Nuclear Data, Geometry and visualisation
    call ndreg_kill()
    call self % viz % kill()

    ! Clean contents
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    
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

end module fluxRR_class
