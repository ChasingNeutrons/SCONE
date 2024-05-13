module dataRR_class

  use numPrecision
  use universalVariables
  use rng_class,                      only : RNG

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
  !! Nuclear data in a random ray-friendly format.
  !!
  !! Stores data and provides access in a manner which is more performant
  !! than is done for MC MG data at present.
  !!
  !! TODO: Add kinetic data and higher-order scattering matrices
  !!
  type, public :: dataRR
    private
    ! Components
    class(baseMgNeutronDatabase), pointer :: mgData => null()
    integer(shortInt)                     :: nG     = 0
    integer(shortInt)                     :: nG2    = 0
    integer(shortInt)                     :: nMat   = 0

    ! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable     :: sigmaT
    real(defFlt), dimension(:), allocatable     :: nuSigmaF
    real(defFlt), dimension(:), allocatable     :: sigmaS
    real(defFlt), dimension(:), allocatable     :: chi
    logical(defBool), dimension(:), allocatable :: fissile

    ! Optional kinetic parameters
    logical(defBool)                        :: doKinetics  = .false.
    integer(shortInt)                       :: nP = 0
    real(defFlt), dimension(:), allocatable :: chiD
    real(defFlt), dimension(:), allocatable :: chiP
    real(defFlt), dimension(:), allocatable :: beta
    real(defFlt), dimension(:), allocatable :: invSpeed

    ! Optional higher-order scattering matrices up to P3
    real(defFlt), dimension(:), allocatable :: sigmaS1
    real(defFlt), dimension(:), allocatable :: sigmaS2
    real(defFlt), dimension(:), allocatable :: sigmaS3

  contains
    
    procedure :: init
    procedure :: kill

    ! Access procedures
    procedure :: getPointers
    !procedure :: getAllPointers
    procedure :: getTotalPointer
    procedure :: getNuFissPointer
    procedure :: getTotalXS
    procedure :: getScatterXS
    procedure :: getNumPrec
    procedure :: isFissile

    ! Private procedures
    procedure, private :: getIdxs
    procedure, private :: getScatterIdxs
    !procedure, private :: getKineticIdxs


  end type dataRR

contains

  !!
  !! Initialise necessary nuclear data.
  !! Can optionally include kinetic parameters.
  !!
  subroutine init(self, db, doKinetics, aniOrder)
    class(dataRR), intent(inout)                    :: self
    class(mgNeutronDatabase),pointer, intent(inout) :: db
    logical(defBool), intent(in), optional          :: doKinetics
    integer(shortInt), intent(in), optional         :: aniOrder
    integer(shortInt)                               :: g, g1, m
    class(RNG)                                      :: rand
    class(baseMgNeutronMaterial), pointer           :: mat
    class(materialHandle), pointer                  :: matPtr
    character(100), parameter :: Here = 'init (dataRR_class.f90)'

    ! Ensure nuclear data is baseMgNeutronDatabase
    select type(db)
      type is (baseMgNeutronDatabase)
        self % mgData => db
      class default
        call fatalError(Here,'Unrecognised MG database type')
    end select

    if (present(doKinetics)) then
      self % doKinetics = doKinetics
    else
      self % doKinetics = .false.
    end if

    ! Store number of energy groups for convenience
    self % nG = self % mgData % nGroups()
    self % nG2 = self % nG * self % nG

    ! Initialise local nuclear data
    ! TODO: clean nuclear database afterwards! It is no longer used
    !       and takes up memory.
    self % nMat = mm_nMat()
    allocate(self % sigmaT(self % nMat * self % nG))
    allocate(self % nuSigmaF(self % nMat * self % nG))
    allocate(self % chi(self % nMat * self % nG))
    allocate(self % sigmaS(self % nMat * self % nG * self % nG))
    allocate(self % fissile(self % nMat))

    ! Create a dummy RNG to satisfy the mgDatabase access interface
    call rand % init(1_longInt)

    do m = 1, self % nMat
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      fiss = .false.
      do g = 1, self % nG
        self % sigmaT(self % nG * (m - 1) + g) = real(mat % getTotalXS(g, rand),defFlt)
        self % nuSigmaF(self % nG * (m - 1) + g) = real(mat % getNuFissionXS(g, rand),defFlt)
        if (self % nuSigmaF(self % nG * (m - 1) + g) > 0) fiss = .true.
        self % chi(self % nG * (m - 1) + g) = real(mat % getChi(g, rand),defFlt)
        ! Include scattering multiplicity
        do g1 = 1, self % nG
          self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)  = &
                  real(mat % getScatterXS(g1, g, rand) * mat % scatter % prod(g1, g) , defFlt)
        end do
      end do
      self % fissile(m) = fiss
    end do

    ! Initialise data necessary for kinetic/noise calculations
    if (self % doKinetics) then


    end if

    ! Initialise higher-order scattering matrices
    if (aniOrder > 0) then


    end if

  end subroutine init

  !!
  !! Calculate the lower and upper indices for accessing the XS array
  !! (excluding scattering and kinetic data)
  !!
  subroutine getIdxs(self, matIdx, idx1, idx2)
    class(dataRR), intent(in)      :: self
    integer(shortInt), intent(in)  :: matIdx
    integer(shortInt), intent(out) :: idx1, idx2

    idx1 = (matIdx - 1) * self % nG + 1
    idx2 = matIdx  * self % nG 

  end subroutine getIdxs

  !!
  !! Calculate the lower and upper indices for accessing the scattering XS array
  !!
  subroutine getScatterIdxs(self, matIdx, idx1, idx2)
    class(dataRR), intent(in)      :: self
    integer(shortInt), intent(in)  :: matIdx
    integer(shortInt), intent(out) :: idx1, idx2

    idx1 = (matIdx - 1) * self % nG2 + 1
    idx2 = matIdx * self % nG2 

  end subroutine getScatterIdxs

  !!
  !! Return if a material is fissile
  !!
  elemental function isFissile(self, matIdx) result(isFiss)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx
    logical(defBool)              :: isFiss

    isFiss = self % fissile(matIdx)

  end function isFissile

  !!
  !! Return the number of precursors
  !!
  elemental function getNumPrec(self) result(nP)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nP

    nP = self % nP

  end function getNumPrec

  !!
  !! Return pointers to all commonly used XS data sets
  !! This is done for a given material, across all energies
  !!
  subroutine getPointers(self, matIdx, sigT, nuSigF, sigS, chi)
    class(dataRR), intent(inout)                       :: self
    integer(shortInt), intent(in)                      :: matIdx
    real(defFlt), dimension(:), pointer, intent(inout) :: sigT, nuSigF, sigS, chi

    call self % getIdxs(matIdx, idx1, idx2)
    call self % getScatterIdxs(matIdx, idx1s, idx2s)
    sigT   => self % sigmaT(idx1:idx2)
    nuSigF => self % nuSigmaF(idx1:idx2)
    chi    => self % chi(idx1:idx2)
    sigS   => self % sigmaS(idx1s:idx2s)

  end subroutine getPointers
  
  !!
  !! Return pointers to only the total XS
  !! This is done for a given material, across all energies
  !!
  subroutine getTotalPointer(self, matIdx, sigT)
    class(dataRR), intent(inout)                       :: self
    integer(shortInt), intent(in)                      :: matIdx
    real(defFlt), dimension(:), pointer, intent(inout) :: sigT

    call self % getIdxs(matIdx, idx1, idx2)
    sigT   => self % sigmaT(idx1:idx2)

  end subroutine getTotalPointer
  
  !!
  !! Return pointers to only the nuFission XS
  !! This is done for a given material, across all energies
  !!
  subroutine getNuFissPointer(self, matIdx, nuFiss)
    class(dataRR), intent(inout)                       :: self
    integer(shortInt), intent(in)                      :: matIdx
    real(defFlt), dimension(:), pointer, intent(inout) :: nuFiss

    call self % getIdxs(matIdx, idx1, idx2)
    sigT   => self % nuFission(idx1:idx2)

  end subroutine getNuFissPointer


  !!
  !! Return total XS in a given material and group
  !!
  elemental function getTotalXS(self, matIdx, g) result(sigT)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, g
    real(defFlt)                  :: sigT

    sigT = self % sigmaT((matIdx - 1) * self % nG + g)

  end function getTotalXS
  
  !!
  !! Return scatter XS in a given material, ingoing group, and outgoing group
  !!
  elemental function getScatterXS(self, matIdx, gIn, gOut) result(sigS)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, gIn, gOut
    real(defFlt)                  :: sigS

    sigS = self % sigmaS((matIdx - 1) * self % nG2 + g)

  end function getScatterXS


  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(dataRR), intent(inout) :: self

    ! Clean contents
    self % mgData     => null()
    self % nG         = 0
    self % nG2        = 0
    self % nMat       = 0
    self % nP         = 0
    self % doKinetics = .false.
    if(allocated(self % sigmaT)) deallocate(self % sigmaT)
    if(allocated(self % sigmaS)) deallocate(self % sigmaS)
    if(allocated(self % nusigmaF)) deallocate(self % nuSigmaF)
    if(allocated(self % chi)) deallocate(self % chi)
    if(allocated(self % fissile)) deallocate(self % fissile)
    if(allocated(self % chiD)) deallocate(self % chiD)
    if(allocated(self % chiP)) deallocate(self % chiP)
    if(allocated(self % beta)) deallocate(self % beta)
    if(allocated(self % invSpeed)) deallocate(self % invSpeed)
    if(allocated(self % sigmaS1)) deallocate(self % sigmaS1)
    if(allocated(self % sigmaS2)) deallocate(self % sigmaS2)
    if(allocated(self % sigmaS3)) deallocate(self % sigmaS3)

  end subroutine kill

end module dataRR_class
