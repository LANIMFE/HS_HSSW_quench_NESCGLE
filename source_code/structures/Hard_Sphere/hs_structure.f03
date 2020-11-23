Module hs_structure
  Use math_constants
  Implicit None
  Contains

  !Function that calculates the direct correlation function of a Hard Sphere
  !System using Percus-Yevick Closure relationship
  ! c(k,ϕ)
  !Depends only on the wave-vector "k" for "k>=0", and density "ϕ" for "ϕ>0"
  Function c_hs_py(phi,k) result(c)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: c
    Real * 8 :: k_func,ckdum,dum1,dum2,dum3,dum4,dumsin,dumcos
    Real * 8 :: k2,k3,k4,k6
    dum1 = - ( ( 1.d0 + 2.d0 * phi ) ** 2 )
    dum2 = 6.d0 * phi * ( ( 1.d0 + ( phi / 2.d0 ) ) ** 2 )
    dum3 = - 0.5d0 * phi * ( ( 1.d0 + ( 2.d0 * phi ) ) ** 2 )
    dum4 = ( 1.d0 - phi )
    dum4 = dum4 ** 4
    dum1 = dum1 / dum4
    dum2 = dum2 / dum4
    dum3 = dum3 / dum4
    k2 = k * k
    If (k > 0.075d0 ) Then
      dumsin = sin( k )
      dumcos = cos( k )
      k3 = k2 * k
      k4 = k3 * k
      k6 = k4 * k2
      c = ( dum1 * ( dumsin - k * dumcos ) / ( k3 ) ) + &
      & ( dum2 * ( ( ( 2.d0 * k ) * dumsin ) + ( ( - ( k2 ) + 2.d0 ) * dumcos ) - 2.d0 ) / &
      & ( k4 ) ) + ( dum3 * ( ( ( 4.d0 * ( k3 ) - 24.d0 * k ) * dumsin )&
      & + ( ( - ( k4 ) + 12.d0 * ( k2 ) - 24.d0 ) * dumcos ) + 24.d0 ) / ( k6 ))
    Else
      c = dum1 * ( (1.d0 / 3.d0) - ( k2 / 30.d0 ) ) + dum2 * ( 0.25d0 -  ( k2 / 36.d0 ) ) &
      & + dum3 * ( ( 1.d0 / 6.d0 ) - ( k2 / 48.d0 ) )
    End If
    c = c * 4.d0 * pi
    return
  End Function

  !Function that calculates the inverse structure factor of a Hard Sphere system
  !Using Percus-Yevick Closure relationship
  ! S(k,ϕ)
  !Depends only on the wave-vector "k" for "k>=0", and density "ϕ" for "ϕ>0"
  Function is_hs_py(phi,k) result(is)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: is
    is = c_hs_py(phi,k)
    is = 1.d0 - 6.d0 * phi * is / pi
    return
  End Function

  !Function that calculates the structure factor of a Hard Sphere system
  !Using Percus-Yevick Closure relationship
  ! S(k,ϕ)
  !Depends only on the wave-vector "k" for "k>=0", and density "ϕ" for "ϕ>0"
  Function s_hs_py(phi,k) result(s)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: s
    s = 1.d0 / is_hs_py(phi,k)
    return
  End Function


  !Function that calculates the direct correlation function of a Hard Sphere
  !System using Verlet-Weiss correction
  ! c(k,ϕ)
  !Depends only on the wave-vector "k" for "k>=0", and density "ϕ" for "ϕ>0"
  Function c_hs_vw(phi,k) result(c)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: c
    Real * 8 :: k_vw,phi_vw
    !Verlet-Weiss correction for the volume fraction
    phi_vw = phi*(1.d0 - (phi / 16.d0))
    !Verlet-Weiss correction for the wave-vector
    k_vw = k * ( ( phi_vw / phi ) ** ( 1.d0 / 3.d0 ) )
    c = c_hs_py(phi_vw,k_vw)
  End Function

  !Function that calculates the inverse structure factor of a Hard Sphere system
  !Using Verlet-Weiss correction
  ! S(k,ϕ)
  !Depends only on the wave-vector "k" for "k>=0", and density "ϕ" for "ϕ>0"
  Function is_hs_vw(phi,k) result(is)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: is
    Real * 8 :: phi_vw, k_vw
    !Verlet-Weiss correction for the volume fraction
    phi_vw = phi*(1.d0 - (phi / 16.d0))
    !Verlet-Weiss correction for the wave-vector
    k_vw = k * ( ( phi_vw / phi ) ** ( 1.d0 / 3.d0 ) )
    is = is_hs_py(phi_vw,k_vw)
    return
  End Function

  !Function that calculates the structure factor of a Hard Sphere system
  !Using Verlet-Weiss correction
  ! S(k,ϕ)
  !Depends only on the wave-vector "k" for "k>=0", and density "ϕ" for "ϕ>0"
  Function s_hs_vw(phi,k) result(s)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: s
    s = 1.d0 / is_hs_vw(phi,k)
    return
  End Function

  !This subroutine writes down the parameters related with structural calculation for
	!a Hard Disk system
	Subroutine writting_hs_p(sys_cond,approx,i_unit)
		Real * 8, dimension(:), intent(in) :: sys_cond
		Integer :: i_unit
		Character(*), intent(in) :: approx
		Real * 8 :: phi
		Integer :: i1
		!System conditions
		phi = sys_cond(1)
		Write(i_unit,*) "### Hard Sphere ###"
		Write(i_unit,*) "Volume Fraction phi = ", phi
		Write(i_unit,*) "###############################"
		Write(i_unit,*) "### Structure Approximation ###"
		Select case (approx)
		case ("PY00")
			Write(i_unit,*) "				Percus-Yevick Wertheim Solution"
		case ("VW00")
			Write(i_unit,*) "				Verlet Weis correction to PY Wertheim solution"
		case Default
			Write(i_unit,*) "				"//approx
		End Select
    Write(i_unit,*) "###############################"
	End Subroutine

  Function P_HS_CS(T,phi) result(P)
    Implicit None
    Real * 8, intent(in) :: T, phi
    Real * 8 :: P
    P = 6.d0 * phi * T * (1.d0 + phi**2 - phi **3) / ( pi * (1.d0 - phi)**3 )
  End Function

End Module
