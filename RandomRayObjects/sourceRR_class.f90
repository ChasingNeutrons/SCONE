module sourceRR_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use dictionary_class,               only : dictionary
  use rng_class,                      only : RNG

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
  
  implicit none
  private

  !!
  !! Source class for random ray calculations.
  !! Will later generalise to having spatial and angular moments.
  !!
  !! Private Members
  !!   geom        -> Pointer to the geometry.
  !!   mgData      -> MG database. Calculation obviously cannot be run in CE.
  !!   nG          -> Number of energy groups, kept for convenience.
  !!   nCells      -> Number of unique cells in the geometry, kept for convenience.
  !!
  !!   source      -> Array of neutron source values of length = nG * nCells
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public :: sourceRR
    private
    ! Components
    class(geometryStd), pointer           :: geom
    class(baseMgNeutronDatabase), pointer :: mgData      => null()
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0

    ! Results space
    real(defFlt), dimension(:), allocatable :: source
    real(defFlt), dimension(:), allocatable :: fixedSource

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: initFixedSource
    procedure :: sourceUpdateKernel
    procedure :: kill

  end type sourceRR

contains

  !!
  !! Initialise source.
  !! If provided the appropriate dictionary, will set a fixed source.
  !!
  subroutine init(self, dictFS)
    class(sourceRR), intent(inout) :: self
    class(dictionary), intent(inout), optional :: dictFS
    class(mgNeutronDatabase),pointer              :: db
    class(geometry), pointer                      :: geom
    class(baseMgNeutronMaterial), pointer         :: mat
    class(materialHandle), pointer                :: matPtr
    character(100), parameter :: Here = 'init (sourceRR_class.f90)'

    ! Ensure that nuclear data is multi-group
    db => ndReg_getNeutronMG()
    if (.not. associated(db)) call fatalError(Here,&
            'No MG nuclear database was constructed')

    ! Ensure nuclear data is baseMgNeutronDatabase
    select type(db)
      type is (baseMgNeutronDatabase)
        self % mgData => db
      class default
        call fatalError(Here,'Unrecognised MG database type')
    end select

    ! Store number of energy groups for convenience
    self % nG = self % mgData % nGroups()

    ! Store number of cells in geometry for convenience
    self % nCells = self % geom % numberOfCells()
    allocate(self % source(self % nCells * self % nG))
    

  end subroutine init

  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF)
    class(sourceRR), target, intent(inout) :: self
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
      self % source(idx) = self % source(idx) / total(g)

    end do

  end subroutine sourceUpdateKernel

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(sourceRR), intent(inout) :: self

    ! Clean contents
    self % geom    => null()
    self % nG        = 0
    self % nCells    = 0
    if(allocated(self % source)) deallocate(self % source)

  end subroutine kill

end module sourceRR_class
