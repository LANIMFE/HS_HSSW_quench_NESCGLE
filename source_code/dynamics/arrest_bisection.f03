Module arrest_bisection
  Use module structure_function_selector
  Use scgle
  Implicit None
  
  contains


  Recursive Subroutine hsay_gsvw_scgle_arrest_bisection(Tu,Tl,phi,swlam,lamk,error,Ta,gam)
    Implicit none
    Real * 8, intent(in) :: Tu, Tl, phi, swlam
    Real * 8, intent(in) ::
    Real * 8, intent(inout) :: error, Ta, gam
    Real * 8, parameter :: tol = 1.d-7, gam_max = 1.d20
  End Subroutine
End Module
