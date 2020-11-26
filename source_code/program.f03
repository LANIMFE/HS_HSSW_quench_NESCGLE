! Main Program !
!It will only call a subroutine to do something
Program Main
  Use math_constants
  Use hs_structure
  Use structure_function_selector
  Use hssw_structure
  Use hsay_structure
  Use hd_structure
  Use scgle
  Use nescgle
  Use general_structural_functions
  Use omp_lib
  Use quadratures
  Implicit None
  Real * 8 :: phi, Tf, swlambda
  integer :: i1
  logical :: writting_op
  Character :: writting_choice
  Write(*,*) "IMPORTANT MESSAGE:"
  Write(*,*) "This program computes the mobility of a Hard Sphere (HS) + Square-Well (SW) system"
  Write(*,*) "when quenching from an infinite temperature (HS limit) to a given temperature. It"
  Write(*,*) "requires as an input the interaction length lambda with range (1:5] of the SW part,"
  Write(*,*) "the HS density as a volume fraction (0:0.582] and the final quenching temperature (0:10]."
  Write(*,*) "Selecting any parameter outside of the expected ranges may result in the program"
  Write(*,*) "behaving unexpectedly. The mobility is saved in the folder '*/data/'.  Additionally,"
  Write(*,*) "the user can select to save all the other kinetic properties given by the theory in the"
  Write(*,*) "subfolder '*/data/kinetics/'."
  Write(*,*) "**********************************************************************************"
  Write(*,*) "Please state the value of lambda (1 : 5]"
  Read(*,*) swlambda
  Write(*,*) "lambda=",swlambda
  Write(*,*) "Please state the value of the volume fraction (0 : 0.582]"
  Read(*,*)  phi
  Write(*,*) "volume fraction=",phi
  Write(*,*) "Please state the value of the final quenching temperature (0 : 10]"
  Read(*,*)  Tf
  Write(*,*) "Final temperature=",Tf
  Write(*,*) "Do you desire to write the additional kinetic properties?"
  Write(*,*) "Please input 'y' for yes and 'n' for no"
  Read(*,*) writting_choice
  If( writting_choice == 'y' ) then
    writting_op=.true.
    Write(*,*) "Additional kinetic properties will be saved in '*/data/kinetics/'"
  else
    writting_op=.false.
    Write(*,*) "No additional kinetic properties will be saved"
  endif
  call HS_HSSW_quench_mobility_dl_convergence_no_gr(phi,Tf,swlambda,writting_op)

Contains

  Subroutine HS_HSSW_quench_mobility_dl_convergence_no_gr(phi,Tf,swlambda,writting_op)
    Implicit None
    real * 8, intent(in) :: phi, Tf,swlambda
    logical, intent(in) :: writting_op
    real * 8, dimension(:), Allocatable :: s_params_i, s_params_f
    real * 8, dimension(2**10) :: qp,qw
    real * 8 :: phi_B_l, phi_B_u,T_c,phi_c,T_sharma
    real * 8 ::  T, z, rho, D0,gam_val,ua_val,ua_writting
    integer, parameter :: knp = 2**10, tnp = 2**9, nst = 10, decim = 32

    integer, parameter :: ktnp = 12
    Real * 8, dimension(7) :: dyn_p
    real * 8, dimension(knp) :: kv, kwv, skv, lamkv,SkI,SkF,alphaku
    real * 8, dimension(knp/2) :: kv1, kwv1, kv2, kwv2
    real * 8, dimension(ktnp) :: kv_t, skv_t,lamkv_t, SkI_t, SkF_t, alphaku_t

    real * 8, dimension(4,knp) :: s_arrays
    real * 8, dimension(3,ktnp) :: s_w_arrays
    real * 8, dimension(3) :: k_quadp, k_quadp1, k_quadp2

    real * 8, dimension(3) :: quad_ch_p
    real * 8, parameter :: k_max = 40.96d0, k_max_1=4.d0, k_min = 0.d0, dt = 1.d-7
    real * 8, parameter :: kgr_max = (2**6) * 40.96d0, kgr_min = 0.d0
    real * 8 :: dk, kc, k, Dl,gam, T_up, T_low,Ts,Ta, Dl_I, Dl_F
    real * 8 :: t_val, dt_linear, db,kgr_min1, kgr_max1
    logical :: convergence
    character (len=100) :: folder,folder_kinetics,folder_mobility
    character (len=4) :: sys, k_approx, dr_approx, r_approx,sys_approx
    character (len=4) :: sys_i, sys_f, sys_approx_i, sys_approx_f
    character (len=4) :: u_char
    character (len=200) :: u_file_name
    integer :: i1,i2,i3,i4,i5,id,ikp,d,u_count,ikp_i,ikp_f
    integer, parameter :: quad_ch_np = 2**10

    real * 8, dimension(quad_ch_np) :: quad_ch_v, quad_ch_w, quad_ch_v2, quad_ch_w2
    real * 8 :: ur_min, ur_max
    integer, dimension(5) :: dyn_units
    integer, dimension(6) :: units
    Real * 8 :: du,u_val,bu,bu_p,rel_dif_b, rel_dif_bF,dr,ncoor_last
    Real * 8, parameter :: db_threshold = 5.0d-3, bF_threshold=1.d-4
    Real * 8, parameter :: dtw_wthreshold = 0.1d0, tw_wstart = 1.d-4
    Real * 8 :: tw_w_1, tw_w, u_w_1,u_w, du_w_1, dtw_w_1
    Real * 8 :: tw_w_p, u_p, dtw_w,Dl_tw,du_w
    Real * 8 :: uw, energy, press, energy_ua, press_ua
    folder= "data/"
    folder_mobility="data/"
    folder_kinetics= "data/kinetics/"
    !swlambda = 1.5d0
    d = 3
    ikp = 4
    !!!!!
    rho = 6.d0 * phi / pi

    k = 6.1d0
    z = 2.d0
    dr = 1.d-2
    d = 3
    ikp = 4
    ikp_i=2
    ikp_f=4

    Allocate (s_params_i(ikp_i))
    Allocate (s_params_f(ikp_f))

    k_quadp(1) = k_min
    k_quadp(2) = k_max
    k_quadp(3) = dble(knp)

    k_quadp1(1) = k_min
    k_quadp1(2) = k_max_1
    k_quadp1(3) = dble(knp/2)

    k_quadp2(1) = k_max_1
    k_quadp2(2) = k_max
    k_quadp2(3) = dble(knp/2)

    quad_ch_p(1) = -1.d0
    quad_ch_p(1) = 1.d0
    quad_ch_p(3) = quad_ch_np

    Write(k_approx,"(A4)") "CLCU"
    Write(dr_approx,"(A4)") "CLCU"
    Write(r_approx,"(A4)") "RECT"
    Write(sys_f,"(A4)") "HSSW"
    Write(sys_approx_f,"(A4)") "GSVW"
    Write(sys_i,"(A4)") "HS00"
    Write(sys_approx_i,"(A4)") "VW00"
    D0 = 1.d0
    dk = 1.d-2
    convergence = .false.
    !NOTA: estoy partiendo en dos intervalos la integracion
    call quad_select(k_approx,k_quadp1,kv1,kwv1)
    call quad_select(k_approx,k_quadp2,kv2,kwv2)
    call Clenshaw_Curtis_quadrature(quad_ch_np,quad_ch_v,quad_ch_w)
    Do i1=1, knp/2
      kv(i1)=kv1(i1)
      kv(i1 + (knp/2)) = kv2(i1)
      kwv(i1)=kwv1(i1)
      kwv(i1 + (knp/2)) = kwv2(i1)
    End Do
    id = 1
    kc = 2.d0 * pi * 1.305d0
    do i1=1, knp
      lamkv(i1) = lambda(kv(i1),kc)
    end do
    !Writting wave-vectors
    kv_t(1) = 2.d-1
    kv_t(2) = 4.d-1
    kv_t(3) = 6.d-1
    kv_t(4) = 8.d-1
    kv_t(5) = 1.d0
    kv_t(6) = 2.d0
    kv_t(7) = 3.d0
    kv_t(8) = 4.d0
    kv_t(9) = 5.d0
    kv_t(10) = 6.18d0
    kv_t(11) = 2 * pi
    kv_t(12) = 7.18
    do i1=1, ktnp
      lamkv_t(i1) = lambda(kv_t(i1),kc)
    end do
    s_arrays(1,:) = kv
    s_arrays(2,:) = kwv
    s_arrays(4,:) = lamkv
    s_w_arrays(1,:) = kv_t
    s_w_arrays(3,:) = lamkv_t
    !!!!!
    !Dynamic units for writting


    units(1)=101
    units(2)=103
    units(3)=104
    units(4)=105
    units(5)=106
    units(6)=107
    dyn_units(1) = units(3)
    dyn_units(2) = units(4)
    dyn_units(3) = units(5)
    dyn_units(4) = units(2)
    dyn_units(5) = units(6)

    s_params_i(1) = phi !Initial phi
    !s_params_i(2) = Ti !Initial Temperature
    !s_params_i(3) = swlambda !Initial length of square well
    !!!!
    s_params_f(1) = phi !Final phi
    s_params_f(2) = Tf !Final Temperature
    s_params_f(3) = swlambda !Final length of square well
    !!!!
    dyn_p(1) = dble (  d  )
    dyn_p(2) = dble ( nst )
    dyn_p(3) = dble ( tnp )
    dyn_p(4) = dble (decim)
    dyn_p(5) = rho
    dyn_p(6) = D0
    dyn_p(7) = dt
    !!!!
    call file_name_constructor(sys_f,s_params_f,folder_mobility,"bu",id,".dat",u_file_name)
    Open (unit = 21, file=trim(u_file_name), status="Replace")
    Write (21,*) "#u // bu // tw "

    !time identifier of kinetic files
    If (writting_op .eqv. .True.) Then
      call file_name_constructor(sys_f,s_params_f,folder_kinetics,"t_u_files",id,".dat",u_file_name)
      Open (unit = 22, file=trim(u_file_name), status="Replace")
      Write (22,*) "#u digit //  tw  // u // bu "
    endif

    Do i1 = 1, knp
      s_params_i(ikp_i) = kv(i1)
      s_params_f(ikp_f) = kv(i1)
      SkI(i1) = structure("s0",sys_i,sys_approx_i,s_params_i)
      SkF(i1) = structure("s0",sys_f,sys_approx_f,s_params_f)
    End Do
    Do i1 = 1, ktnp
      s_params_i(ikp_i) = kv_t(i1)
      s_params_f(ikp_f) = kv_t(i1)
      SkI_t(i1) = structure("s0",sys_i,sys_approx_i,s_params_i)
      SkF_t(i1) = structure("s0",sys_f,sys_approx_f,s_params_f)
    End Do
    call Calc_alphaSku(kv,SkF,D0,alphaku)
    call Calc_alphaSku(kv_t,SkF_t,D0,alphaku_t)
    Print*, "Calculating ua:"
    ua_val = calc_ua_isoch_inst_quench(SkI,SkF,alphaku,lamkv,kv,kwv,rho,d)
    Print*, "ua = ", ua_val
    skv   = sku_quench(SkI,SkF,alphaku,ua_val)
    skv_t = sku_quench(SkI_t,SkF_t,alphaku_t,ua_val)
    !Calculating initial diffussion
    If (writting_op .eqv. .True.) then
      call file_open_HSSW_quench_mobility_1(sys_i,s_params_i,folder_kinetics,id,units,0)
    End If
    s_arrays(3,:)   = SkI
    s_w_arrays(2,:) = SkI_t
    !writting_op = .True.
    Print*, "Calculating Initial state dynamics"
    call Long_t_dynamics_Dl_convergence(dyn_p,s_arrays,Dl_I,writting_op,s_w_arrays,dyn_units)
    If (writting_op .eqv. .True.) Then
      call efh_HSSW_quench_mobility(kv,skI,kwv,units,sys_i,s_params_i,sys_approx_i,Dl_I,gam_val,kc)
    End If
    If (writting_op .eqv. .True.) Then
      call file_open_HSSW_quench_mobility_1(sys_f,s_params_f,folder_kinetics,id,units,999)
    End If
    s_arrays(3,:)   = skv
    s_w_arrays(2,:) = skv_t
    Print*, "Calculating Final state dynamics"
    !call SCGLE_dynamics(dyn_p,s_arrays,Dl_F,writting_op,s_w_arrays,dyn_w_arrays,dyn_units)
    call Long_t_dynamics_Dl_convergence(dyn_p,s_arrays,Dl_F,writting_op,s_w_arrays,dyn_units)
    If (Dl_F < 1.d-6) Dl_F = 1.d-6
    If (writting_op .eqv. .True.) Then
      call efh_HSSW_quench_mobility(kv,skv,kwv,units,sys_f,s_params_f,sys_approx,Dl_F,gam_val,kc)
    End If
    gam_val = gamma(d,rho,kv,kwv,Skv,lamkv)
    Print*, "gamma=", gam_val
    !Calculating g(r;t=0) and g_a(r)
    du = 1.d-7
    u_val = 0.d0
    t_val = 0.d0
    Write (21,*) u_val, Dl_I, t_val

    flush(21)
    bu_p = Dl_I
    convergence = .false.
    u_count = 0
    tw_w_p = tw_wstart!dtw_wthreshold
    do while (convergence .eqv. .false.)
      u_val = u_val + du
      skv   = sku_quench(SkI,SkF,alphaku,u_val)
      skv_t = sku_quench(SkI_t,SkF_t,alphaku_t,u_val)
      s_arrays(3,:)   = Skv
      s_w_arrays(2,:) = Skv_t
      call Long_t_dynamics_Dl_convergence(dyn_p,s_arrays,bu,writting_op,s_w_arrays,dyn_units)
      If (bu < 1.d-6) bu = 1.d-6
      db = abs(bu_p - bu)
      dt_linear = du / bu
      t_val = t_val + dt_linear
      Write (21,*) u_val, bu, t_val
      flush(21)
      !Writting wanted times
      dtw_w = log10(t_val) - log10(tw_w_p)
      print*, dtw_w, t_val, tw_w_p
      If ( t_val > tw_wstart .and. dtw_w > dtw_wthreshold ) Then
        tw_w_1 = t_val - dt_linear
        u_w_1 = u_val - du
        If (tw_w_1 < tw_wstart) Then
          dtw_w_1  = tw_wstart - tw_w_1
          du_w_1   = bu * dtw_w_1
        Else If (dtw_w > dtw_wthreshold) Then
          dtw_w_1  = tw_w_p * ( 10**(dtw_wthreshold) ) - tw_w_1
          du_w_1   = bu * dtw_w_1
        Else
        End If
        u_w = u_w_1 + du_w_1
        skv   = sku_quench(SkI,SkF,alphaku,u_w)
        skv_t = sku_quench(SkI_t,SkF_t,alphaku_t,u_w)
        !File writting
        u_count=u_count+1
        Print*, "Calculating writting log equispaced, count:", u_count
        call file_open_HSSW_quench_mobility_1(sys_f,s_params_f,folder_kinetics,id,units,u_count)
        s_arrays(3,:)   = Skv
        s_w_arrays(2,:) = Skv_t
        call Long_t_dynamics_Dl_convergence(dyn_p,s_arrays,Dl_tw,writting_op,s_w_arrays,dyn_units)
        call efh_HSSW_quench_mobility(kv,skv,kwv,units,sys_f,s_params_f,sys_approx,Dl_tw,gam_val,kc)
        tw_w_p = tw_w_1 + dtw_w_1
        Print*, "For the values: u_w=",u_w,"tw_approx=", tw_w_p
        If (writting_op .eqv. .True.) Then
          Write (22,*) u_count, tw_w_p, u_w, Dl_tw
          flush(22)
        End if
      End If
      !Changing du and out while options
      rel_dif_b  = abs((bu_p / bu) - 1.d0)
      rel_dif_bF = abs((bu / Dl_F) - 1.d0)
      If ( rel_dif_b < db_threshold ) Then
        du = du * 2.d0
      Else
        du = du / 2.d0
      End If
      If (u_val + du >= ua_val ) convergence = .True.
      If (rel_dif_bF <= bF_threshold ) convergence = .True.
      bu_p = bu
    End Do
    Write (21,*) ua_val, Dl_F
    flush(21)
    close(21)
    If (writting_op .eqv. .True.) Then
      close(22)
    endif
  End Subroutine

  Subroutine file_open_HSSW_quench_mobility_1(sys,s_params,folder,id,units,u_count)
    Implicit None
    integer, Intent(in) :: id, u_count
    character (len=4), Intent(in) :: sys
    character (len=100), Intent(in) :: folder
    integer, dimension(:), Intent(in) :: units
    real * 8, dimension(:), Intent(in) :: s_params
    character (len=200) :: S_file_name, inf_file_name, u_file_name, g_file_name
    character (len=200) :: Fc_file_name, Fs_file_name, dyn_file_name, tau_file_name
    character(len=5) :: u_char
    Write(u_char,"(A2,I3.3)") "u_",u_count
    if (u_count==999) Write(u_char,"(A5)") "ua000"
    call file_name_constructor(sys,s_params,folder,u_char//"_s",id,".dat",S_file_name)
    call file_name_constructor(sys,s_params,folder,u_char//"_tau",id,".dat",tau_file_name)
    call file_name_constructor(sys,s_params,folder,u_char//"_Fc",id,".dat",Fc_file_name)
    call file_name_constructor(sys,s_params,folder,u_char//"_Fs",id,".dat",Fs_file_name)
    call file_name_constructor(sys,s_params,folder,u_char//"_dyn",id,".dat",dyn_file_name)
    call file_name_constructor(sys,s_params,folder,u_char//"_info",id,".inf",inf_file_name)
    open (unit=units(1), file=trim(S_file_name), status="Replace")
    open (unit=units(2), file=trim(inf_file_name), status="Replace")
    open (unit=units(3), file=trim(Fc_file_name), status="Replace")
    open (unit=units(4), file=trim(Fs_file_name), status="Replace")
    open (unit=units(5), file=trim(dyn_file_name), status="Replace")
    open (unit=units(6), file=trim(tau_file_name), status="Replace")
  End Subroutine

  Subroutine close_units(units)
    Implicit None
    Integer, Dimension(:) :: units
    Integer :: i1
    Do i1=1, size(units)
      close(units(i1))
    End Do
  End Subroutine

  Subroutine efh_HSSW_quench_mobility(kv,sk,kwv,units,sys,s_params,sys_approx,Dl,gam_val,kc)
    Implicit None
    real * 8, Dimension(:), Intent(in) :: kv, sk, kwv, s_params
    integer, Dimension(:), Intent(in) :: units
    real * 8, Intent(in) :: Dl, gam_val, kc
    character (len=4), Intent(in) :: sys, sys_approx
    call s_writting(kv,sk,kwv,units(1))
    call sys_writting(sys,s_params,sys_approx,units(2))
    call lambda_info_writting(kc,units(2))
    call dynamics_Long_t_info_writting(Dl,gam_val,units(2))
    call close_units(units)
  End Subroutine

End program
