module rayHandling_func

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use randomRayMaths_func,            only : exponential
  use dictionary_class,               only : dictionary
  use rng_class,                      only : RNG

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryStd_class,              only : geometryStd

  ! Random ray - or a standard particle
  ! Also particleState for easier output
  use particle_class,                      only : ray => particle, particleState

  implicit none
  private

  !!
  !! Set of functions and subroutines to handle everything to do with rays
  !! in random ray
  !!
  public :: moveRay
  public :: initialiseRay
  public :: transportSweepFlatIso
  public :: transportSweepLinearIso
  public :: transportSweepLinearAni
  public :: transportSweepFlatAni

contains

  subroutine initialiseRay(geom, lb, ub, r, arrays)
    class(geometryStd), pointer, intent(in)                :: geom
    real(defReal), dimension(3), intent(in)                :: lb
    real(defReal), dimension(3), intent(in)                :: ub
    type(ray), intent(inout)                               :: r
    class(arraysRR), pointer, intent(inout)                :: arrays
    real(defReal)                                          :: mu, phi
    real(defReal), dimension(3)                            :: u, rand3, x
    integer(shortInt)                                      :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (rayHandling_func.f90)'

    i = 0
    mu = TWO * r % pRNG % get() - ONE
    phi = TWO_PI * r % pRNG % get()
    u = rotateVector([ONE, ZERO, ZERO], mu, phi)

    rejection : do
      rand3(1) = r % pRNG % get()
      rand3(2) = r % pRNG % get()
      rand3(3) = r % pRNG % get()
      x = lb + (ub - lb) * rand3

      ! Exit if point is inside the geometry
      call geom % whatIsAt(matIdx, cIdx, x, u)
      if (matIdx /= OUTSIDE_MAT) exit rejection

      i = i + 1
      if (i > 5000) then
        call fatalError(Here, 'Infinite loop when searching for ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(x, u, 1, ONE)
    call self % geom % placeCoord(r % coords)

    if (.not. arrays % found(cIdx)) call arrays % newFound(cIdx, x)

  end subroutine initialiseRay
  
  !!
  !! Move ray across a cell, into the next  
  !! Use distance caching or standard ray tracing
  !! Distance caching seems a little bit more unstable
  !! due to FP error accumulation, but is faster.
  !! This can be fixed by resetting the cache after X number
  !! of distance calculations.
  !!
  subroutine moveRay(r, doCache, ints, geom, length, event, cache, hitVacuum)
    type(ray), intent(inout)                 :: r
    logical(defBool), intent(in)             :: doCache
    integer(longInt), intent(inout)          :: ints
    class(geometryStd), pointer, intent(in)  :: geom 
    real(defReal), intent(inout)             :: length
    integer(shortInt), intent(out)           :: event  
    class(distCache), pointer, intent(inout) :: cache
    logical(defBool), intent(out)            :: hitVacuum

    if (doCache) then
      if (mod(ints,20_longInt) == 0)  cache % lvl = 0
      call geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)
    else
      call geom % moveRay_noCache(r % coords, length, event, hitVacuum)
    end if
    ints = ints + 1

  end subroutine moveRay
  
  !!
  !! Set maximum flight distance and ensure ray is active
  !!
  subroutine checkRayLength(totalLength, dead, termination, activeRay, length)
    real(defReal), intent(in)    :: totalLength
    real(defReal), intent(in)    :: dead
    real(defReal), intent(in)    :: termination
    real(defBool), intent(inout) :: activeRay
    real(defReal), intent(inout) :: length
      
    if (totalLength >= dead) then
      length = termination - totalLength 
      activeRay = .true.
    else
      length = dead - totalLength
    end if

  end subroutine checkRayLength

  !!
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweepFlatIso(r, ints, nG, doCache, dead, termination, geom, XSData, arrays)
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: nG
    logical(defBool), intent(in)                          :: doCache
    real(defReal), intent(in)                             :: dead
    real(defReal), intent(in)                             :: termination
    class(geometryStd), pointer, intent(in)               :: geom
    class(dataRR), pointer, intent(in)                    :: XSData
    class(arraysRR), pointer, intent(inout)               :: arrays
    integer(shortInt)                                     :: matIdx, g, cIdx, event, matIdx0
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt
    real(defFlt), dimension(nG)                           :: attenuate, delta, fluxVec
    real(defFlt), pointer, dimension(:)                   :: scalarVec, sourceVec, totVec
    real(defReal), dimension(3)                           :: r0, mu0
    
    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    do g = 1, nG
      fluxVec(g) = arrays % getSourceValue(cIdx,g)
    end do

    ints = 0
    matIdx0 = 0
    totalLength = ZERO
    activeRay = .false.
    do while (totalLength < self % termination)

      ! Get material and cell the ray is moving through
      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        
        ! Cache total cross section
        call XSData % getTotalPointer(matIdx, totVec)
      end if

      ! Set maximum flight distance and ensure ray is active
      call checkRayLength(totalLength, dead, termination, activeRay, length)

      ! Move ray
      call moveRay(r, doCache, ints, geom, length, event, cache, hitVacuum)
      totalLength = totalLength + length
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. arrays % found(cIdx)) call arrays % newFound(cIdx, r % rGlobal + length * HALF * r % dirGlobal())

      lenFlt = real(length,defFlt)
      call arrays % getSourcePointer(cIdx, sourceVec)

      !$omp simd
      do g = 1, nG
        attenuate(g) = exponential(totVec(g) * lenFlt)
        delta(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g)
        fluxVec(g) = fluxVec(g) - delta(g)
      end do

      ! Accumulate to scalar flux
      if (activeRay) then
      
        call arrays % setLock(cIdx)
          call arrays % getFluxPointer(cIdx, scalarVec)
          !$omp simd
          do g = 1, nG
            scalarVec(g) = scalarVec(g) + delta(g) 
          end do
          call arrays % incrementVolume(cIdx, length)
        call arrays % unsetLock(cIdx)

        if (arrays % hasHit(cIdx) == 0) call arrays % hitCell(cIdx)
      
      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, nG
          fluxVec(g) = 0.0_defFlt
        end do
      end if

    end do

  end subroutine transportSweepFlatIso

end module rayHandling_func
