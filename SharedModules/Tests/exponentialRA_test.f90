module exponentialRA_test

  use numPrecision
  use exponentialRA_func, only : exponential
  use pfUnit_mod

  implicit none

contains

  !!
  !! Test exponential rational approximation on a few values
  !! ExponentialRA evaluated (1 - exp(-x))
  !! This is compared against the analytic equivalent
  !!
@Test
  subroutine testExponentialRA()
    real(defFlt)                :: x, res, resRA
    real(defFlt), parameter     :: tol = 1E-5

    x = 0.5
    res   = ONE - exp(-x)
    resRA = exponential(x) 

    @assertEqual(res, resRA, tol)

    x = 0.2_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

    @assertEqual(res, resRA, tol)

    x = 0.03_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

    @assertEqual(res, resRA, tol)

    x = 3.0_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

    @assertEqual(res, resRA, tol)


    x = 0.0001_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

    @assertEqual(res, resRA, tol)
  end subroutine testExponentialRA


end module exponentialRA_test
