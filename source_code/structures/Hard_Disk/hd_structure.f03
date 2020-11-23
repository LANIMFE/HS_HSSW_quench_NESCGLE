Module hd_structure
	Use math_constants
	Use quadratures
	Implicit None
	contains

	!Function that calculates the direct correlation function of a Hard Disk Function
	!using Rosenfeld Fundamental Measure Theory structure
	Function ck_HD_ros(phi,k) result(ck)
		Implicit None
		Real * 8, intent(in) :: phi, k
		Real * 8 :: ck
		Real * 8 :: a0,b0,g0,k2
		Real * 8 ::  phip, phim, J10, J00,J11
		phip = 1.d0 + phi
		phim = 1.d0 - phi
		a0 = 1.d0 + ( (2.d0 * phi - 1.d0) * phip / (phim ** 3) ) + (2.d0 * phi/phim)
		a0 = a0 / phi
		b0 = (phip / (phim ** 2)) -1.d0 - (3.d0 * phi / phim)
		b0 = b0 / phi
		g0 = 1.d0 / phim
		k2 = k * k
		If (k>0.075d0) Then
			J10 = BESSEL_JN(1,k/2.d0)
			J11 = BESSEL_JN(1,k)
			J00 = BESSEL_JN(0,k/2.d0)
			ck = (4.d0 * a0 * J10 * J10/ k2) + (2.d0 * b0 * J00 * J10 / k) + (2.d0 * g0 * J11 / k)
		Else
      ck = a0 * ( (0.25d0 - k2 / 64.d0) ) + b0 * (0.5d0 - 3.d0 * k2 / 64.d0) + g0 * (1.d0 - k2 / 8.d0)
		End If
		ck = - pi * ck
		return
	End Function

	!Function that calculates the inverse of the structure factor of a Hard Disk Function
	!using Rosenfeld Fundamental Measure Theory structure
	Function isk_HD_ros(phi,k) result(isk)
		Implicit None
		Real * 8, intent(in) :: phi, k
		Real * 8 :: isk
		Real * 8 :: rho
		rho = 4.d0 * phi / pi
		isk = 1.d0 - rho * ck_HD_ros(phi,k)
		return
	End Function

	!Function that calculates the structure factor of a Hard Disk Function
	!using Rosenfeld Fundamental Measure Theory structure
	Function sk_HD_ros(phi,k) result(sk)
		Implicit None
		Real * 8, intent(in) :: phi, k
		Real * 8 :: sk
		sk = 1.d0 / isk_HD_ros(phi,k)
		return
	End Function

	!Auxiliary Function for the Baus-Colot direct correlation function
	Function aux_F_BC(q) result(F)
		Implicit None
		Real * 8, intent(in) :: q
		Real * 8 :: F
		If (q>0.075d0) Then
			F = 2.d0 * BESSEL_JN(1,q) / q
		Else
			F = 1.d0 - (q * q / 8.d0)
		End If
	End Function

	!Function that calculates the direct correlation function of a Hard Disk Function
	!using the Modified Baus-Colot structure by Guo and Riebel
	Function ck_HD_BCGR(phi,k) result(ck)
		Implicit None
		Real * 8, intent(in) :: phi, k
		Real * 8 :: ck
		Real * 8 :: c0, a, d, h, Fq, a2, x2, aqx, f2ak2, phi2,phi3,phi4
		Integer, parameter :: n_integrand = 2**10
		Real * 8, dimension(n_integrand) :: h_integrand, xp, xw
		Integer :: i1
		phi2 = phi * phi
		phi3 = phi2 * phi
		phi4 = phi2 * phi2
		d = 0.98946d0
		a = 0.3699d0 * phi4 - 1.2511d0 * phi3 + 2.0199d0 * phi2 - 2.2373d0 * phi + 2.1d0
		c0 = - (1.d0 - d * phi2) / ((1.d0 - 2.d0 * phi + d * phi2) ** 2)
		call Clenshaw_Curtis_quadrature(n_integrand,xp,xw)
		call rescale_OW(n_integrand,xp,xw,-1.d0,1.d0,1.d0/a,1.d0)
		Fq = aux_F_BC(k)
		a2 = a * a
		Do i1 = 1, n_integrand
			x2 = xp(i1) * xp(i1)
			aqx = a * k * xp(i1)
			h_integrand(i1) = sqrt( 1.d0 - x2 ) * ( Fq - a2 * x2 * aux_F_BC(aqx) )
		End Do
		h = 16.d0 * sum(h_integrand * xw) / pi
		f2ak2 = aux_F_BC(a * k / 2.d0)
		f2ak2 = f2ak2 * f2ak2
		ck = (4.d0 * (1.d0 - a2 * phi) * Fq) + (a2 * phi * (f2ak2 + h))
		ck = pi * c0 * ck / 4.d0
		return
	End Function

	!Function that calculates the inverse of the structure factor of a Hard Disk Function
	!using the Modified Baus-Colot structure by Guo and Riebel
	Function isk_HD_BCGR(phi,k) result(isk)
		Implicit None
		Real * 8, intent(in) :: phi, k
		Real * 8 :: isk
		Real * 8 :: rho
		rho = 4.d0 * phi / pi
		isk = 1.d0 - rho * ck_HD_BCGR(phi,k)
		return
	End Function

	!Function that calculates the structure factor of a Hard Disk Function
	!using the Modified Baus-Colot structure by Guo and Riebel
	Function sk_HD_BCGR(phi,k) result(sk)
		Implicit None
		Real * 8, intent(in) :: phi, k
		Real * 8 :: sk
		sk = 1.d0 / isk_HD_BCGR(phi,k)
		return
	End Function

	!This subroutine writes down the parameters related with structural calculation for
	!a Hard Disk system
	Subroutine writting_hd_p(sys_cond,approx,i_unit)
		Real * 8, dimension(:), intent(in) :: sys_cond
		Integer :: i_unit
		Character(*), intent(in) :: approx
		Real * 8 :: phi, T, lambda
		Integer :: i1
		!System conditions
		phi = sys_cond(1)
		Write(i_unit,*) "### Hard Disk ###"
		Write(i_unit,*) "Area Fraction phi = ", phi
		Write(i_unit,*) "###############################"
		Write(i_unit,*) "### Structure Approximation ###"
		Select case (approx)
		case ("ROSF")
			Write(i_unit,*) "				Rosenfeld Fundamental Measure Theory"
		case ("BCGR")
			Write(i_unit,*) "				Modified Baus-Colot structure by Guo and Riebel"
		case Default
			Write(i_unit,*) "				"//approx
		End Select
    Write(i_unit,*) "###############################"
	End Subroutine


End Module
