Module nescgle
Use math_constants
Use omp_lib
Use scgle
Implicit None

contains

!Subroutine that calculates the alpha vector used for the calculation of
!the S(k,u(t)) structure
Subroutine Calc_alphaSku(k_v,SkF_v,D0_val,alpha_v)
	implicit none
	real * 8, dimension(:), intent(in) :: k_v, SkF_v
	real * 8, dimension(:), intent(out) :: alpha_v
	real * 8, intent (in) :: D0_val
	alpha_v = 2.d0 * (k_v ** 2) * D0_val / SkF_v
End Subroutine

!Subroutine that calculates the structure factor of the system at a given
!u material time through an instantaneous isochoric quench
Function sku_quench(Sk0_v,SkF_v,alpha_v,u_val) result(sku)
	implicit none
	real * 8, dimension(:), intent(in) :: Sk0_v, SkF_v, alpha_v
	real * 8, intent(in) :: u_val
	real * 8, dimension(size(Sk0_v)) :: Sku
	integer :: i1
	If (u_val > 0.d0) Then
		sku = SkF_v + ( Sk0_v - SkF_v ) * exp( - alpha_v * u_val )
	Else
		sku= Sk0_v
	End If
	return
End Function

!Subroutine that calculates the value of the material time u in which a system arrest
!calculated up to a max value of u and gamma when an instantaneous isochoric quench
!is applied to the system
Function calc_ua_isoch_inst_quench(Sk0,SkF,alpha,lam,k,kw,rho,d) result(ua_val)
	implicit none
	real * 8, dimension(:), intent(in) :: Sk0, SkF, alpha, lam, k, kw
	real * 8, dimension(size(Sk0)) :: sku
	real * 8, intent(in) :: rho
	real * 8 :: gam_val
	integer, intent(in) :: d
	real * 8 :: ua_val
	logical :: convergence
	real * 8 :: u_low, u_up, gam_test_u, gam_test_l, error, u_test, gam_test
	integer :: i1
	Real * 8, parameter :: gam_max_val = 1.d20, tol = 1.d-8
	Real * 8, parameter :: ua_max_val = 1.d20, ua_min_val = 1.d-20
	u_low = 0.d0
	u_up = 1.d-4
	sku = sku_quench( Sk0,SkF,alpha,u_up)
	gam_test_u = gamma(d,rho,k,kw,sku,lam)
	sku = sku_quench(Sk0,SkF,alpha,u_low)
	gam_test_l = gamma(d,rho,k,kw,sku,lam)
	If ( gam_test_u > gam_max_val .OR. gam_test_u /= gam_test_u ) Then
		convergence = .false.
		do while (convergence .eqv. .false.)
			u_low = u_up
			u_up  = u_up * 2.d0
			sku = sku_quench( Sk0,SkF,alpha,u_up)
			gam_test = gamma(d,rho,k,kw,sku,lam)
			If ( gam_test < gam_max_val .or. u_up > ua_max_val )	convergence = .true.
		end do
	else if (gam_test_l < gam_max_val ) then
		u_up = u_low
		return
	end if
	ua_val  = u_up
	gam_val = gam_test
	error =  abs(u_up - u_low) / u_up
	if ( error > tol ) convergence = .false.
	if ( u_up > ua_max_val ) convergence = .true.
	do while (convergence .eqv. .false.)
		u_test = (u_up + u_low) / 2.d0
		sku = sku_quench( Sk0,SkF,alpha,u_test)
		gam_test = gamma(d,rho,k,kw,sku,lam)
		if ( gam_test < gam_max_val ) then
			u_up = u_test
			gam_val = gam_test
		else
			u_low = u_test
		end if
		error =  abs(u_up - u_low) / u_up
		if ( error < tol ) convergence = .true.
	end do
	ua_val = u_up
	print*, ua_val, gam_val, error
	return
End Function

!Subroutine that search in an isochore the arrest temperature and then

Real * 8 Function gamma_u(Sk0_v,SkF_v,alpha_v,u_val,lam_v,k_v,kw_v,phi_val,d)
	implicit none
	real * 8, dimension(:), intent(in) :: Sk0_v, SkF_v, alpha_v, lam_v, k_v, kw_v
	real * 8, intent(in) :: u_val, phi_val
	integer, intent(in) :: d
	real * 8, dimension(size(Sk0_v)) :: Sku_v
	Real * 8 :: rho
	rho = phi_val * 6.d0 / pi
	!Call Calc_Sk_u(Sk0_v,SkF_v,alpha_v,Sku_v,u_val)
	gamma_u = gamma(d,rho,k_v,kw_v,Sku_v,lam_v)
	return
End Function

!This subroutine is for the writting of the structure factor as well as the
!weights of the quadrature and the information file
Subroutine su_si_sf_writting(k,sku,ski,skf,kw,i_unit)
	Implicit None
	Real * 8, dimension(:), intent(in) :: k, sku, ski, skf, kw
	Integer :: i1
	Integer, intent(in) :: i_unit
	!k S(k) kw WRITTING
	Write (i_unit,*) "# //	k //	S(k;u)	//	Si(k)	//	Sf(k)	//	Weights		//"
	Do i1=1, size(k)
		Write(i_unit,*) k(i1),sku(i1),ski(i1),skf(i1),kw(i1)
	End Do
End Subroutine



!Subroutine diagram_phi_T(phi,Ti,Tfl,system,k,kw,lamk,)
!	Implicit None
	!!! Volume/Area fraction and # density !!!
!	Real * 8 :: phi, rho
	!!! Quench Initial Temperature // Quench Final Temperature !!!
!	Real * 8 :: Ti, Tf
	!!! Structure factors !!!
!	Real * 8, dimension(:) :: sku, ski, skf
	!!! Wave vectors !!!
!	Real * 8, dimension(:) :: k, kw
	!!! Radial distribution functions !!!
!	Real * 8, dimension(:) :: gru, gri, grf
	!!! Radial vectors !!!
!	Real * 8, dimension(:) :: r, rw
	!!! Differential Radial vectors !!!
!	Real * 8, dimension(:) :: dr, drw
	!!! Integer for the wave-vector parameter !!!
!	Integer :: ikp
	!!! System and system approximation !!!
!	Character ( len = 4 ) :: sys, sys_approx
	!!! System parameters !!!
!	Real * 8, dimension(:) :: s_params
	!!! Writting Folder !!!
!	Character ( len = 100 ) :: folder
	!!! Writting Wave vectors !!!
!	Real * 8, dimension(:) :: k_wr, kw_wr
	!!! Writting Structure factors !!!
!	Real * 8, dimension(:) :: sku_wr, ski_wr, skf_wr
	!!! Writting File Names !!!
!	Character(len=100) :: gam_file_name, S_file_name, g_file_name, S_inf_file_name
!	Character(len=100) :: gam_inf_file_name
	!!! Dummy Variables !!!
!	Real * 8 :: ua
!	Integer :: i1,i2,i3
!	Do i1 = 1, 4000
!		phi = i1 * 1.d-4
!		rho = phi * 6.d0 / pi !3D
!		Ti = 4.d0
!		swlambda = 1.5d0
!		s_params(1) = phi
!		s_params(2) = Ti
!		Do i2 = 1, knp
!			s_params(ikp) = kv(i2)
!			SkI(i2) = structure("s0",sys,sys_approx,s_params)
!		End Do
!		Do i2 = 1, knpgr
!			s_params(ikp) = kvgr(i2)
!			skigr(i2) = structure("s0",sys,sys_approx,s_params)
!		End Do
!		i3 = 800
!		ua_val = 0.d0
!		call file_name_constructor(sys,s_params,folder,"gamma",id,".dat",S_file_name)
!		call file_name_constructor(sys,s_params,folder,"info_ua",id,".inf",inf_file_name)
!		open (unit=104, file=trim(S_file_name), status="Replace")
!		open (unit=105, file=trim(inf_file_name), status="Replace")
!		call k_info_writting(kv,k_quadp,k_approx,105)
!		call lambda_info_writting(kc,105)
!		write(104,*) "#phi,   Tf,   gamma,    ua"
!		close(105)
!		Do while (ua_val < 1.d20)
!			T_low = i3 * 1.d-4
!			s_params(2) = T_low
!			Print*, phi, T_low
!			call file_name_constructor(sys,s_params,folder,"s_ua",id,".dat",S_file_name)
!			call file_name_constructor(sys,s_params,folder,"g_ua",id,".dat",g_file_name)
!			call file_name_constructor(sys,s_params,folder,"info_ua",id,".inf",inf_file_name)
!			Print*, S_file_name
!			open (unit=101, file=trim(S_file_name), status="Replace")
!			open (unit=102, file=trim(g_file_name), status="Replace")
!			open (unit=103, file=trim(inf_file_name), status="Replace")
!			Do i2 = 1, knp
!				s_params(ikp) = kv(i2)
!				SkF(i2) = structure("s0",sys,sys_approx,s_params)
!			End Do
!			Print*, "Skf Done"
!			call Calc_alphaSku(kv,SkF,D0,alphaku)
!			Print*, "alphaku Done"
!			ua_val = calc_ua_isoch_inst_quench(SkI,SkF,alphaku,lamkv,kv,kwv,rho,d)
!			skv = sku_quench(SkI,SkF,alphaku,ua_val)
!			if (ua_val < 1.d20) Then
!			 gam_val = gamma(d,rho,kv,kwv,skv,lamkv)
!			 ua_writting = ua_val
!			 Ta = T_low
!			 Print*, "gamma=",gam_val, "ua=",ua_val
!			 Write(104,*) phi, T_low, gam_val, ua_val
!			End If
!			Do i2 = 1, knpgr
!				s_params(ikp) = kvgr(i2)
!				skfgr(i2) = structure("s0",sys,sys_approx,s_params)
!			End Do
!			call Calc_alphaSku(kvgr,skfgr,D0,alphakgr)
!			skvgru = sku_quench(skigr,skfgr,alphakgr,ua_val)
!			call calc_gr_dependent_functions(rho,kvgr,kwvgr,skvgru,rv,drv,drw,grv,coor_n_rv,d)
!			call s_writting(kvgr,skvgru,kwvgr,101)
!			call g_writting(rv,grv,coor_n_rv,102)
!			call k_info_writting(kvgr,kgr_quadp,k_approx,103)
!			call r_info_writting(rv,r_quadp,r_approx,103)
!			call dr_info_writting(drv,drw,dr_quadp,dr_approx,103)
!			flush(101)
!			flush(102)
!			flush(103)
!			Close(101)
!			Close(102)
!			Close(103)
!			i3=i3+1
!		End Do
!		Write(201,*) phi, Ta, gam_val, ua_writting
!		flush(201)
!	End Do
!	Close(201)
!End Subroutine


End Module
