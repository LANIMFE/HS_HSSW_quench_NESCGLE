Module structure_function_selector
	Use hs_structure
	Use hssw_structure
	Use hsay_structure
	Use hd_structure
	Use quadratures
	Implicit None
	contains

	!This functions serves the purpose of selecting a structural function
	!with a system potential, an approximation and parameters:
	!"F" options:
	!		"c" = Direct correlation function
	!		"is" = Inverse of the Structure Factor
	!		"s" = Structure Factor
	!"potential" options:
	!		"HS" = Hard Sphere
	!								"Approx" options for Hard Sphere:
	!									"PY" = Percus-Yevick
	! "fp" are the floating point parameters
	! "ip" are the integer point parameters
	!	"qp" is the quadrature array
	! "qw" is the weight array
	!needed by the function and are grouped inside a
	!vector, look inside the function for the potential and approximation
	!parameters needed and their organization
	Function	structure(F,potential,Approx,fp,ip,qp,qw) result(FR)
		Implicit None
		character(len=2), intent(in) :: F
		character(len=4), intent(in) :: potential,Approx	!Cases selectors
		real * 8, dimension(:), intent(in) :: fp	!Function real parameters
		real * 8, dimension(:), intent(inout), optional :: qp, qw !Quadrature arrays for points and weights
		integer, dimension(:), intent(inout), optional :: ip !Integer parameters
		real * 8 :: FR	!Function Result
		real * 8 :: phi, k, T, lambda, z, rmin, rmax
		integer :: qnp,selector1
		real * 8, dimension(2) :: CCfp
		integer, dimension(1) :: CCip
		!Potentials:
		select case(potential)
		!Hard Sphere:
		case("HS00")
			!Parameters:
			!volume fraction "phi" -> params(1)
			!wave vector "k" -> params(2)
			phi = fp(1)
			k = fp(2)
			!Approximations:
			select case(Approx)
			!Percus-Yevick
			case("PY00")
				!Functions:
				select case(F)
				case("c0")
					FR = c_hs_py(phi,k)
				case("is")
					FR = is_hs_py(phi,k)
				case("s0")
					FR = s_hs_py(phi,k)
				end select
			case("VW00")
				!Functions:
				select case(F)
				case("c0")
					FR = c_hs_vw(phi,k)
				case("is")
					FR = is_hs_vw(phi,k)
				case("s0")
					FR = s_hs_vw(phi,k)
				end select
			end select
		!Hard Disk:
		case("HD00")
			!Parameters:
			!volume fraction "phi" -> params(1)
			!wave vector "k" -> params(2)
			phi = fp(1)
			k = fp(2)
			!Approximations:
			select case(Approx)
			!Rosenfeld FMT
			case("ROSF")
				!Functions:
				select case(F)
				case("c0")
					FR = ck_HD_ros(phi,k)
				case("is")
					FR = isk_HD_ros(phi,k)
				case("s0")
					FR = sk_HD_ros(phi,k)
				end select
			!Modified Baus-Colot by Guo and Riebel
			case("BCGR")
				!Functions:
				select case(F)
				case("c0")
					FR = ck_HD_BCGR(phi,k)
				case("is")
					FR = isk_HD_BCGR(phi,k)
				case("s0")
					FR = sk_HD_BCGR(phi,k)
				end select
		end select
		!Hard Sphere + Square Well
		case("HSSW")
			!Parameters:
			!volume fraction "phi" -> params(1)
			!temperature "T" -> params(2)
			!square well length "lambda" -> params(3)
			!wave vector "k" -> params(4)
			phi = fp(1)
			T = fp(2)
			lambda = fp(3)
			k = fp(4)
			!Approximations:
			select case(Approx)
			!Sharma-Sharma with Verlet Weiss
			case("SHVW")
				!Functions:
				select case(F)
				case("c0")
					FR = c_hssw_shvw(phi,T,lambda,k)
				case("is")
					FR = is_hssw_shvw(phi,T,lambda,k)
				case("s0")
					FR = s_hssw_shvw(phi,T,lambda,k)
				end select
				!Generalized Sharma-Sharma with Verlet Weiss
			case("GSVW")
				!Functions:
				select case(F)
				case("c0")
					FR = c_hssw_gsvw(phi,T,lambda,k)
				case("is")
					FR = is_hssw_gsvw(phi,T,lambda,k)
				case("s0")
					FR = s_hssw_gsvw(phi,T,lambda,k)
				end select
			case("BOLT")
				!Functions:
				select case(F)
				case("s0")
					FR = s_hssw_bolt(phi,T,lambda,k)
				end select
			case("STOS")
				!Functions:
				select case(F)
				case("s0")
					FR = s_hssw_stos(phi,T,lambda,k)
				end select
			end select
		!Hard Sphere + Attractive Yukawa
		case("HSAY")
			!Parameters:
			!volume fraction "phi" -> params(1)
			!temperature "T" -> params(2)
			!Yukawa potential length parameter "z" -> params(3)
			!wave vector "k" -> params(4)
			phi = fp(1)
			T = fp(2)
			z = fp(3)
			k = fp(4)
			!Approximations:
			select case(Approx)
			!Sharma-Sharma with Verlet Weiss
			case("SHVW")
				!Functions:
				select case(F)
				case("c0")
					FR = c_hsay_shvw(phi,T,z,k)
				case("is")
					FR = is_hsay_shvw(phi,T,z,k)
				case("s0")
					FR = s_hsay_shvw(phi,T,z,k)
				end select
			!Generalized Sharma-Sharma with Verlet Weiss
			case("GSVW")
				!Parameters
				!Number of points -> ip(1)
				!Selector to compute Clenshaw Curtis quadrature -> ip(2) (0 for calculation)
				CCfp(1) = fp(5)
				CCfp(2) = fp(6)
				CCip(1) = ip(1)
				selector1 = ip(2)
				If (selector1 == 0) Then
					!Call quad_select("CLCU",CCip,CCfp,qp,qw)
					ip(2) = 1
				End If

				select case(F)
				case("c0")
					FR = c_hsay_gsvw(phi,T,z,k,qp,qw)
				case("is")
					FR = is_hsay_gsvw(phi,T,z,k,qp,qw)
				case("s0")
					FR = s_hsay_gsvw(phi,T,z,k,qp,qw)
				end select
			end select
		end select
		return
	End Function
End Module
