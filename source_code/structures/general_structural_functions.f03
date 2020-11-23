Module general_structural_functions
	use math_constants
	use quadratures
	use hssw_structure
	use hd_structure
	Implicit None
	contains

	!Function that calculates the radial distribution function g(r) through the
	!structure factor. The function needs:
	! "r" -> value of the radial length
	! "rho" -> particle density of the system
	! "k" -> wave vector array
	! "kw" -> wave vector weights
	! "sk" -> structure factor
	Function radial_dist_3D(r,rho,k,kw,s) result(g)
		Implicit None
		Real * 8, intent(in) :: r, rho !radial length and density
		Real * 8, intent(in), dimension(:) :: k, kw, s !wave vector, weights and structure factor
		Real * 8 :: g
		Real * 8, dimension(size(k)) :: integrand
		Integer :: i1
		do i1 = 1, size(k)
			if ( k(i1) * r >= 0.075d0 ) integrand(i1) = dsin( k(i1) * r ) / r
			if ( k(i1) * r <  0.075d0 ) integrand(i1) = k(i1) - (( r * r * (k(i1) ** 3))/6.d0)
			!integrand(i1) = sin( k(i1) * r ) / r
			integrand(i1) = k(i1) * kw(i1) * (S(i1) - 1.d0) * integrand(i1)
		end do
		g = 1.d0 + ( 1.d0 / ( 2.d0 * rho *  pi * pi ) ) * sum(integrand)
		return
	End function


	!Function that calculates the radial distribution function g(r) through the
	!structure factor. The function needs:
	! "r" -> value of the radial length
	! "rho" -> particle density of the system
	! "k" -> wave vector array
	! "kw" -> wave vector weights
	! "sk" -> structure factor
	Function radial_dist_2D(r,rho,k,kw,s) result(g)
		Implicit None
		Real * 8, intent(in) :: r, rho !radial length and density
		Real * 8, intent(in), dimension(:) :: k, kw, s !wave vector, weights and structure factor
		Real * 8 :: g
		Real * 8, dimension(size(k)) :: integrand
		Integer :: i1
		do i1 = 1, size(k)
			if ( k(i1) * r >= 0.075d0 ) integrand(i1) = dsin( k(i1) * r ) / r
			if ( k(i1) * r <  0.075d0 ) integrand(i1) = k(i1) - (( r * r * (k(i1) ** 3))/6.d0)
			!integrand(i1) = sin( k(i1) * r ) / r
			integrand(i1) = k(i1) * kw(i1) * (S(i1) - 1.d0) * integrand(i1)
		end do
		g = 1.d0 + ( 1.d0 / ( 2.d0 * rho *  pi * pi ) ) * sum(integrand)
		return
	End function

	!Function that calculates the coordination number integral for a given set of points
	!given the value up to the starting point. The function needs:
	! "rho" -> system density
	! "coor_last" -> the value of the coordination number up to the starting value of the integral
	! "r" -> radial distance vector
	! "rw" -> weights of the radial distance
	! "gr" -> radial distribution function evaluated in the "r" distances
	Function coordination_num_3D(rho,coor_last,r,rw,gr) result (coor_n)
		Implicit None
		Real * 8, intent(in) :: rho,coor_last
		Real * 8, intent(in), dimension(:) :: r, rw, gr
		Real * 8 :: coor_n
		Real * 8, dimension(size(r)) :: integrand
		integrand = r * r  * gr * rw
		coor_n = coor_last + 4.d0 * pi * rho * sum(integrand)
		return
	End Function

	!Function that calculates the potential energy of the system
	! "rho" -> system density
	! "r" -> radial distance vector
	! "rw" -> weights of the radial distance
	! "gr" -> radial distribution function evaluated in the "r" distances
	! "ur" -> pair potential array for each distance
	! "energy" -> output value of the energy
	Function pot_energy_3D(rho,r,rw,gr,ur) result (energy)
		Implicit None
		Real * 8, intent(in) :: rho
		Real * 8, intent(in), dimension(:) :: r, rw, gr, ur
		Real * 8 :: energy
		Real * 8, dimension(size(r)) :: integrand
		integrand = r * r  * gr * ur * rw
		energy = 2.d0 * pi * rho * sum(integrand)
		return
	End Function

	!This subroutine calculates the radial distribution function as well as some other
	!properties dependent of this function. This function needs:
	! "rho" -> system density
	! "k" -> Wave-vectors array to calculate the Fourier transform
	! "kw" -> Wave-vectors weights
	! "Sk" -> Structure Factor array
	! "r" -> radial distance array in which we want to compute "g(r)"
	! "dr" -> differential array of radial distance for each distance in the array "r"
	! "dr_quadp" -> quadrature points associated with the integrals over "g(r)"
	! "dr_quadw" -> quadrature weights associated with the integrals over "g(r)"
	!The outputs are:
	! "g_r" -> Radial distribution function "g(r)"
	! "coor_n_r" -> Coordination number depending on the cutting value "r"
	! "d" -> dimension of the system
	Subroutine calc_gr_dependent_functions(rho,k,kw,Sk,r,dr_quadp,dr_quadw,g_r,coor_n_r,d)
		Implicit None
		Real * 8, intent(in), dimension(:) :: k, kw, Sk, r, dr_quadp, dr_quadw
		Real * 8, intent(out), dimension(:) :: g_r,coor_n_r
		Real * 8, intent(in) :: rho
		Integer, intent(in) :: d
		Real * 8, dimension(size(dr_quadp)) :: drqp, drqw, dg_r
		Real * 8 :: dr,coor_n_last
		Integer :: i1,i2
		!open(unit=11, file="gr_test_F.dat", status="replace")
	  coor_n_r = 0.d0
		coor_n_last = 0.d0
		drqw = dr_quadw
		dr = r(2)-r(1)
		if (d == 3) then
			Do i1=1, size(r)
				drqp = dr_quadp + r(i1) - dr
				g_r(i1)=  radial_dist_3D(r(i1),rho,k,kw,sk)
			   do i2=1, size(drqp)
			     dg_r(i2) = radial_dist_3D(drqp(i2),rho,k,kw,sk)
			   end do
			   coor_n_r(i1) = coordination_num_3D(rho,coor_n_last,drqp,drqw,dg_r)
			  ! write(11,*) r, gr,grtail, gr+grtail,ncoor, (ncoor)/(4.d0 * pi *rho* (r ** 3)/3.d0)
			   coor_n_last = coor_n_r(i1)
		  End do
		Else if(d == 2) Then
		End If
	  !close(11)
	End Subroutine

	!This subroutine writes down the information of the system used, as well as the
	!approximation scheme
	Subroutine sys_writting(sys,sysp,approx,i_unit)
		Implicit None
		Real * 8, dimension(:), intent(in) :: sysp
		Character(len=4), intent(in) :: sys,approx
		Integer :: i_unit
		select case(sys)
		case("HSSW")
			call writting_hssw_p(sysp,approx,i_unit)
		case("HD00")
			call writting_hd_p(sysp,approx,i_unit)
		case("HS00")
			call writting_hs_p(sysp,approx,i_unit)
		end select
	End Subroutine

	!This subroutine is for the writting of the structure factor as well as the
	!weights of the quadrature and the information file
	Subroutine s_writting(k,sk,kw,i_unit)
		Implicit None
		Real * 8, dimension(:), intent(in) :: k, sk, kw
		Integer :: i1
		Integer, intent(in) :: i_unit
		!k S(k) kw WRITTING
		Write (i_unit,*) "# 			k 								S(k)								Weights"
		Do i1=1, size(k)
			Write(i_unit,*) k(i1),sk(i1),kw(i1)
		End Do
	End Subroutine

	!This subroutine is for the writting of the radial distribution function
	!as well as other radial distance dependent functions
	Subroutine g_writting(r,gr,coor_n_r,i_unit)
		Implicit None
		Real * 8, dimension(:), intent(in) :: r, gr,coor_n_r
		Integer, Intent (in) :: i_unit
		Integer :: i1
		Write (i_unit,*) "# 			r 								g(r)						Coordination Number"
		Do i1=1, size(r)
			Write(i_unit,*) r(i1),gr(i1),coor_n_r(i1)
		End Do
	End Subroutine

	!This subroutine is for the information with respect the wave-vectors used
	Subroutine k_info_writting(k,quadp,k_quad,i_unit)
		Implicit None
		Real * 8, dimension(:), intent(in) :: k,quadp
		Character(*), intent(in) :: k_quad
		Integer :: i1
		Integer, intent(in) :: i_unit
		!INFO FILE WRITTING
		Write(i_unit,*) "### Wave-Vectors Quadrature ###"
		Call quad_writting(k_quad,i_unit,quadp)
		Write(i_unit,*) "### Wave-Vectors General INFO ###"
		Write(i_unit,*) "number of points     = ", size(k)
		Write(i_unit,*) "k minimum value      = ", k(1)
		Write(i_unit,*) "k maximum value      = ", k(size(k))
	End Subroutine

	!This subroutine is for the information with respect the radial distance used
	Subroutine r_info_writting(r,r_params,r_quad,i_unit)
		Implicit None
		Real * 8, dimension(:), intent(in) :: r, r_params
		Character(*), intent(in) :: r_quad
		Integer, intent(in) :: i_unit
		!INFO FILE WRITTING
		Write(i_unit,*) "### Radial distance Quadrature ###"
		Call quad_writting(r_quad,i_unit,r_params)
		Write(i_unit,*) "### Radial distance General INFO ###"
		Write(i_unit,*) "number of points     = ", size(r)
		Write(i_unit,*) "r minimum value      = ", r(1)
		Write(i_unit,*) "r maximum value      = ", r(size(r))
	End Subroutine

	!This subroutine is for the information with respect the potential
	Subroutine Energy_info_writting(energy,r_params,r_quad,i_unit)
		Implicit None
		Real * 8 :: energy
		Real * 8, dimension(:), intent(in) :: r_params
		Character(*), intent(in) :: r_quad
		Integer, intent(in) :: i_unit
		!INFO FILE WRITTING
		Write(i_unit,*) "### Potential Energy info ###"
		Write(i_unit,*) "Potential Energy     = ", energy
		Call quad_writting(r_quad,i_unit,r_params)
	End Subroutine

	!This subroutine is for the information with respect the radial distance differential used
	Subroutine dr_info_writting(dr,dr_w,dr_params,dr_quad,i_unit)
		Implicit None
		Real * 8, dimension(:), intent(in) :: dr, dr_w, dr_params
		Character(*), intent(in) :: dr_quad
		Integer, intent(in) :: i_unit
		Integer :: i1
		Write(i_unit,*) "### Radial distance Differencial General INFO ###"
		Write(i_unit,*) "number of points     = ", size(dr)
		Write(i_unit,*) "dr minimum value      = ", dr(1)
		Write(i_unit,*) "dr maximum value      = ", dr(size(dr))
		Write(i_unit,*) "###############################"
		Write(i_unit,*) "### Radial distance Differencial Quadrature ###"
		Call quad_writting(dr_quad,i_unit,dr_params)
		Write(i_unit,*) "### Radial distance Differencial Quadrature values ###"
		Write(i_unit,*) "						r												Weights"
		Do i1=1, size(dr)
			Write(i_unit,*) dr(i1), dr_w(i1)
		End Do
		Write(i_unit,*) "###############################"
	End Subroutine

	Subroutine file_name_constructor(sys,sysp,folder,preffix,id,suffix,file_name)
		Implicit None
		Real * 8, dimension(:), intent(in) :: sysp
		Character(*), intent(in) :: sys,folder,preffix,suffix
		Integer, intent(in) :: id
		Real * 8 :: phi, T, lambda, z
		Integer :: i_phi, di_phi, i_T, di_T, i_L, di_L, i_z, di_z
		Integer, parameter :: i_phi_m=1, di_phi_m=6, i_T_m =2, di_T_m = 6, id_m=2
		Integer, parameter :: i_L_m = 2, di_L_m = 2, i_z_m = 2, di_z_m = 2
		Character (len=i_phi_m+di_phi_m+1) :: phi_char
		Character (len=i_T_m+di_T_m+1) :: T_char
		Character (len=i_L_m+di_L_m+1) :: L_char
		Character (len=i_z_m+di_z_m+1) :: z_char
		Character (len=id_m) :: id_char
		Character (len=100) :: sys_ID
		Character (*) :: file_name
		Character (len=50) :: FMT0,FMT1,FMT2,FMTID
		Integer :: i1
		!Global Formats
		Write(FMT0,"(A,I1,A1,I1,A5,I1,A1,I1,A1)") "(I",i_phi_m,".",i_phi_m,",A1,I",di_phi_m,".",di_phi_m,")"
		Write(FMTID,"(A2,I1,A1,I1,A1)") "(I",id_m,".",id_m,")"
		!Global character variables
		phi = sysp(1)
		i_phi = int(phi)
		di_phi = nint((phi - dble(i_phi)) * (10.d0 ** di_phi_m) )
		Write(id_char,trim(FMTID)) id
		Write(phi_char,trim(FMT0)) i_phi,"_",di_phi
		!System selector
		select case(sys)
		case("HSSW")
			!Formats
			Write(FMT1,"(A2,I1,A1,I1,A5,I1,A1,I1,A1)") "(I",i_T_m,".",i_T_m,",A1,I",di_T_m,".",di_T_m,")"
			Write(FMT2,"(A2,I1,A1,I1,A5,I1,A1,I1,A1)") "(I",i_L_m,".",i_L_m,",A1,I",di_L_m,".",di_L_m,")"
			!System conditions
			lambda = sysp(3)
			T = sysp(2)
			!Calculating the associated integers for phi, T and lambda:
			i_T = int(T)
			di_T = nint( (T- dble(i_T)) * (10.d0 ** di_T_m) )
			i_L = int(lambda)
			di_L = nint( (lambda- dble(i_L)) * (10.d0 ** di_L_m) )
			!Writting character variables
			Write(T_char,FMT1) i_T,"_",di_T
			Write(L_char,FMT2) i_L,"_",di_L
			Write(id_char,"(I2.2)") id
			!Writting System ID
			sys_ID = trim("_HSSW_phi_"//phi_char//"_T_"//T_char//"_L_"//L_char//"_id_"//id_char)
		case("HD00")
			!Writting System ID
			sys_ID = trim("_HD_phi_"//phi_char//"_id_"//id_char)
		case("HS00")
			!Writting System ID
			sys_ID = trim("_HS_phi_"//phi_char//"_id_"//id_char)
		case("HSAY")
			!Formats
			Write(FMT1,"(A2,I1,A1,I1,A5,I1,A1,I1,A1)") "(I",i_T_m,".",i_T_m,",A1,I",di_T_m,".",di_T_m,")"
			Write(FMT2,"(A2,I1,A1,I1,A5,I1,A1,I1,A1)") "(I",i_z_m,".",i_z_m,",A1,I",di_z_m,".",di_z_m,")"
			!System conditions
			z = sysp(3)
			T = sysp(2)
			!Calculating the associated integers for phi, T and lambda:
			i_T = int(T)
			di_T = nint( (T- dble(i_T)) * (10.d0 ** di_T_m) )
			i_z = int(z)
			di_z = nint( (z- dble(i_z)) * (10.d0 ** di_z_m) )
			!Writting character variables
			Write(T_char,FMT1) i_T,"_",di_T
			Write(z_char,FMT2) i_z,"_",di_z
			Write(id_char,"(I2.2)") id

			!Writting System ID
			sys_ID = trim("_HSAY_phi_"//phi_char//"_T_"//T_char//"_z_"//z_char//"_id_"//id_char)
		case default
			sys_ID = trim(sys//"_phi_"//phi_char//"_id_"//id_char)
		end select
		file_name = trim(trim(folder)//trim(preffix)//trim(sys_ID)//trim(suffix))
	End Subroutine
End Module
