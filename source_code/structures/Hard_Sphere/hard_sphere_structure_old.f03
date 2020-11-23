Module hard_sphere_structure_module
  Implicit None
  Real * 8, Dimension(:), Allocatable :: Int_absci_1, Int_weight_1
  Integer, Parameter :: Int_points=100
  Real * 8, Dimension(:) , Allocatable :: musj
  Contains

  !Function for the calculation of the structure factor for a Hard Sphere System
  Real * 8 Function Sk_HS_py(vw_option,phi_val,k_val)
    Implicit None
    Logical, Intent(in) :: vw_option  !Logical variable for Verlet-Weiss correction if True
    Real * 8, Intent(in) :: phi_val !Volume fraction value
    Real * 8, Intent (in) :: k_val !Value of the wave vector
    Real * 8 :: k_vw,k_func,ckdum,phi_vw,dum1,dum2,dum3,dum4,dumsin,dumcos
    Real * 8 :: k_vw2,k_vw3,k_vw4,k_vw6
    If (vw_option .EQV. .TRUE.) Then
      !Verlet-Weiss correction for the volume fraction
      phi_vw=phi_val*(1.d0 - (phi_val / 16.d0))
      !Verlet-Weiss correction for the wave-vector
      k_vw= k_val * ( ( phi_vw / phi_val ) ** ( 1.d0 / 3.d0 ) )
    Else
      phi_vw=phi_val
      k_vw= k_val
    End If
    dum1= - ( ( 1.d0 + 2.d0 * phi_vw ) ** 2 )
    dum2= 6.d0 * phi_vw * ( ( 1.d0 + ( phi_vw / 2.d0 ) ) ** 2 )
    dum3= - 0.5d0 * phi_vw * ( ( 1.d0 + ( 2.d0 * phi_vw ) ) ** 2 )
    dum4= ( 1.d0 - phi_vw )
    dum4= dum4 ** 4
    dum1= dum1 / dum4
    dum2= dum2 / dum4
    dum3= dum3 / dum4
    k_vw2=k_vw * k_vw
    If (k_vw > 0.075d0 ) Then
      dumsin= sin( k_vw )
      dumcos= cos( k_vw )
      k_vw3=k_vw2 * k_vw
      k_vw4=k_vw3 * k_vw
      k_vw6=k_vw4 * k_vw2
      ckdum= ( dum1 * ( dumsin - k_vw * dumcos ) / ( k_vw3 ) ) + &
      & ( dum2 * ( ( ( 2.d0 * k_vw ) * dumsin ) + ( ( - ( k_vw2 ) + 2.d0 ) * dumcos ) - 2.d0 ) / &
      & ( k_vw4 ) ) + ( dum3 * ( ( ( 4.d0 * ( k_vw3 ) - 24.d0 * k_vw ) * dumsin )&
      & + ( ( - ( k_vw4 ) + 12.d0 * ( k_vw2 ) - 24.d0 ) * dumcos ) + 24.d0 ) / ( k_vw6 ))
    Else
      ckdum = dum1 * ( (1.d0 / 3.d0) - ( k_vw2 / 30.d0 ) ) + dum2 * ( 0.25d0 -  ( k_vw2 / 36.d0 ) ) &
      & + dum3 * ( ( 1.d0 / 6.d0 ) - ( k_vw2 / 48.d0 ) )
    End If
    ckdum =  24.d0 * phi_vw * ckdum
    Sk_HS_py = 1.d0 / ( 1.d0 - ckdum )
    return
    !*********!
  End Function

  !Subroutine for the calculation of the Structure factor of a Hard Sphere System
  !for an array of wave vectors "k_val" of dimension "k_points" that goes into
  !a predefined array "sk_val"
  !vw_option Logical value, True for Verlet-Weiss Correction, False Otherwise
  Subroutine Calc_Sk_HS_py_mono(vw_option,phi_val,k_points,k_val,sk_val)
    Implicit None
    Logical, Intent(in) :: vw_option  !Logical variable for activating Verlet-Weiss correction if True
    Real * 8, Intent(in) :: phi_val !Volume fraction of the system
    Real * 8, Dimension(:), Intent(in) :: k_val !Wave vector array
    Real * 8, Dimension(:), Intent(out) :: sk_val !Structure factor array
    Integer, Intent(in) :: k_points
    Integer :: i1
    Do i1=1, k_points
      sk_val(i1)=Sk_HS_py(vw_option,phi_val,k_val(i1))
    End Do
  End Subroutine

  !Function for the calculation of the structure factor of Hard Disk System using Rosenfeld Approximation
  Real * 8 Function Sk_HD_ros(phi_val,k_val)
    Implicit None
    Real * 8, Intent(in) :: phi_val
    Real * 8, Intent(in) :: k_val
    Integer :: i1
    Real * 8 :: A,B,G,J0,J1,J12,Z1,Z2,X, phi_val2, phi_val3, phi_val4,dum1
    Real * 8 :: k_val2
    phi_val2=phi_val * phi_val
    phi_val3=phi_val2 * phi_val
    phi_val4=phi_val3 * phi_val
    Z1= ( 1.d0 + 0.128d0 * ( phi_val2 ) + 0.027d0 * ( phi_val3 ) + 0.06d0 * &
    & ( phi_val4 ) ) / ( 1.d0 - phi_val ) ** 2
    Z2= ( ( 2.d0 * 0.128d0 * phi_val ) + ( 3.d0 * 0.027d0 * ( phi_val2 ) ) + &
    & ( 4.d0 * 0.06d0 * ( phi_val3 ) ) ) * ( ( 1.d0 - phi_val ) ** 2 )
    Z2= ( Z2 / ( 1.d0 - phi_val ) ** 4 ) + ( 2.d0 * ( Z1 / ( 1.d0 - phi_val ) ) )
    !Z1=1.0d0/(1.d0-eta(1))**2.d0
    !Z2=2.0d0/(1.d0-eta(1))**3.d0
    X= Z1 + phi_val * Z2
    G= Sqrt( Z2 / 2.d0 )
    A= ( 1.d0 + ( ( 2.d0 * phi_val - 1.d0 ) * X ) + ( 2.d0 * phi_val * G ) ) / phi_val
    B= ( ( ( 1.d0 - phi_val ) * X ) - 1.d0 - ( 3.d0 * phi_val * G ) ) / phi_val
    If ( k_val > 0.075d0 ) Then
    	J0= BESSEL_JN ( 0 , k_val/ 2.d0 )
    	J1= BESSEL_JN ( 1 , k_val/ 2.d0 )
    	J12= BESSEL_JN ( 1 , k_val )
      dum1= A * ( ( 2.d0 * J1 / k_val ) ** 2 ) + ( B * 2.d0 * J0 * J1 /k_val ) + G * 2.d0 * J12 /k_val
    Else
      k_val2= k_val * k_val
      dum1= A * ( ( 0.25d0 - k_val2 / 64.d0 ) ) + ( B  * ( 0.5d0 - 3.d0 * k_val2 / 64.d0 )  ) &
      & + G * ( 1.d0 - k_val2 / 8.d0  )
    End If
    dum1=  dum1 * 4.d0 * phi_val
    Sk_HD_ros= 1.d0 / ( 1.d0 + dum1 )
    return
  End Function

  !Subroutine for the calculation of the Structure factor of a Hard Disk System
  !for an array of wave vectors "k_val" of dimension "k_points" that goes into
  !a predefined array "sk_val" using the approximation derived by Rosenfeld
  Subroutine Calc_Sk_hd_ros_mono(phi_val,k_points,k_val,sk_val)
    Implicit None
    Real * 8, Intent(in) :: phi_val !Volume fraction of the system
    Real * 8, Dimension(:), Intent(in) :: k_val !Wave vector array
    Real * 8, Dimension(:), Intent(out) :: sk_val !Structure factor array
    Integer, Intent(in) :: k_points
    Integer :: i1
    Do i1=1, k_points
      sk_val(i1)=Sk_HD_ros(phi_val,k_val(i1))
    End Do
  End Subroutine

  !!Function that calculates the structure factor of a HARD Sphere + Atractive Yukawa
  !System with particle size equal to the unit using Sharma Sharma approximation
  Real * 8 Function  Sk_HSAY(vw_option,phi_val,Temp_val,z_yuk,k_val)
    Implicit None
    Logical, Intent(in) :: vw_option  !Logical variable for Verlet-Weiss correction if True
    Real * 8, Intent (in) :: phi_val !Volume fraction value
    Real * 8, Intent (in) :: k_val !Value of the wave vector
    Real * 8, Intent (in) :: Temp_val !Value of the Temperature
    Real * 8, Intent (in) :: z_yuk !Value of the screening constant
    Real * 8 :: dum1,ck_dum
    ck_dum= ( 1.d0 / Sk_HS_py(vw_option,phi_val,k_val) )
    If ( k_val > 0.075d0 ) Then
      dum1= 24.d0 * phi_val * ( ( k_val * cos(k_val) ) + ( z_yuk * sin(k_val) ) ) &
      & / ( Temp_val * k_val * ( ( k_val ** 2 ) + ( z_yuk ** 2) ) )
    Else
      dum1= 24.d0 * phi_val * ( ( 1.d0 - 0.5d0 * k_val ** 2 ) + ( z_yuk * ( 1.d0 - &
      & k_val ** 2 / 6.d0 ) ) ) / ( Temp_val * ( ( k_val ** 2 ) + ( z_yuk ** 2) ) )
    End If
    ck_dum= ck_dum-dum1
    sk_HSAY= 1.d0 / ck_dum
    return
  End Function

  !Subroutine for the calculation of the Structure factor of a Hard Sphere + Atractive
  !Yukawa System for an array of wave vectors "k_val" of dimension "k_points" that goes into
  !a predefined array "sk_val", with "vw_option" Logical value which for True selects Verlet-Weiss Correction
  Subroutine Calc_Sk_HSAY_sh_sh(vw_option,phi_val,Temp_val,z_yuk,k_points,k_val,sk_val)
    Implicit None
    Logical, Intent(in) :: vw_option  !Logical variable for activating Verlet-Weiss correction if True
    Real * 8, Intent(in) :: phi_val !Volume fraction of the system
    Real * 8, Intent(in) :: Temp_val !Value of the temperature
    Real * 8, Intent(in) :: z_yuk !Value of the Yukawa screening reach
    Real * 8, Dimension(:), Intent(in) :: k_val !Wave vector array
    Real * 8, Dimension(:), Intent(out) :: sk_val !Structure factor array
    Integer, Intent(in) :: k_points
    Integer :: i1
    Do i1=1, k_points
      sk_val(i1)=Sk_HSAY(vw_option,phi_val,Temp_val,z_yuk,k_val(i1))
    End Do
  End Subroutine

  !Function that calculates the structure factor of a HARD Disk + Atractive Yukawa
  !System with particle size equal to the unit using Sharma Sharma approximation
  Real * 8 Function  Sk_HDAY(phi_val,Temp_val,z_yuk,k_val)
    Implicit None
    Real * 8, Intent (in) :: phi_val !Volume fraction value
    Real * 8, Intent (in) :: k_val !Value of the wave vector
    Real * 8, Intent (in) :: Temp_val !Value of the Temperature
    Real * 8, Intent (in) :: z_yuk !Value of the screening constant
    Real * 8 :: dum1,ck_dum,y, J0
    Real * 8, Parameter :: dy=1.d-5
    ck_dum= ( 1.d0 / Sk_HD_ros(phi_val,k_val) )
    !FOURIER TRANSFORM OF beta u(r)
    y= 1.d0 !initial value of integration -> Particle size =1.d0
    dum1= 0.d0
    Do While ( exp(- z_yuk * y ) > dy )
    	J0= BESSEl_JN(0, k_val * y)
    	dum1= dum1 + J0 * exp(- y * z_yuk  )
    	y= y + dy
    End Do
    dum1 = 8.d0 * phi_val * dum1 * ( dy ) * exp( z_yuk ) / Temp_val !!!!!!!!!!!!!!!!!!!!!!!!!!!rho beta u
    ck_dum=ck_dum-dum1
    sk_HDAY= 1.d0 / ck_dum
    return
  End Function

  !Function that calculates the structure factor of a HARD Disk + Atractive Yukawa
  !System with particle size equal to the unit using Sharma Sharma approximation
  Real * 8 Function  Sk_HDAY2(phi_val,Temp_val,z_yuk,k_val)
    Implicit None
    Real * 8, Intent (in) :: phi_val !Volume fraction value
    Real * 8, Intent (in) :: k_val !Value of the wave vector
    Real * 8, Intent (in) :: Temp_val !Value of the Temperature
    Real * 8, Intent (in) :: z_yuk !Value of the screening constant
    Real * 8 :: dum1,ck_dum,y, J0, weight
    Integer :: EOFint,i1
    ck_dum= ( 1.d0 / Sk_HD_ros(phi_val,k_val) )
    !FOURIER TRANSFORM OF beta u(r)
    dum1= 0.d0
    EOFint=0
    i1=0
    Do i1=1,Int_points
      y= Int_absci_1(i1)
      weight=Int_weight_1(i1)
      J0= BESSEl_JN(0, k_val * y)
    	dum1= dum1 + weight *J0 * exp(- y * z_yuk  )
    End Do
    dum1 = 8.d0 * phi_val * ( 1.d0 / ( SQRT( z_yuk ** 2 + k_val ** 2) ) - dum1  ) * exp( z_yuk ) / Temp_val !!!!!!!!!!!!!!!!!!!!!!!!!!!rho beta u
    ck_dum=ck_dum-dum1
    sk_HDAY2= 1.d0 / ck_dum
    return
  End Function

  !Obtention of the abscissas and weights for the integral used in Sk_HDAY2
  Subroutine Sk_HDAY2_integral_constants(Allocation)
    Implicit None
    Integer :: EOFint,i1
    Logical, Intent(in) :: Allocation
    If (Allocation .eqv. .True.) Then
      Allocate(Int_absci_1(Int_points))
      Allocate(Int_weight_1(Int_points))
      Open(unit=33, File="./source_code/operations/Integral_method1_abscisas_weights.dat", Status="Old" )
      Do i1=1,Int_points
        Read (33, *, IOSTAT=EOFint) Int_absci_1(i1), Int_weight_1(i1)
      End Do
      Close(33)
    Else
      Deallocate(Int_absci_1)
      Deallocate(Int_weight_1)
    End If
  End Subroutine

  !Subroutine for the calculation of the Structure factor of a Hard Disk + Atractive
  !Yukawa System for an array of wave vectors "k_val" of dimension "k_points" that goes into
  !a predefined array "sk_val", with "vw_option" Logical value which for True selects Verlet-Weiss Correction
  Subroutine Calc_Sk_HDAY_sh_sh(phi_val,Temp_val,z_yuk,k_points,k_val,sk_val)
    Implicit None
    Real * 8, Intent(in) :: phi_val !Volume fraction of the system
    Real * 8, Intent(in) :: Temp_val !Value of the temperature
    Real * 8, Intent(in) :: z_yuk !Value of the Yukawa screening reach
    Real * 8, Dimension(:), Intent(in) :: k_val !Wave vector array
    Real * 8, Dimension(:), Intent(out) :: sk_val !Structure factor array
    Integer, Intent(in) :: k_points
    Integer :: i1
    Do i1=1, k_points
      !sk_val(i1)=Sk_HDAY(phi_val,Temp_val,z_yuk,k_val(i1))
      sk_val(i1)=Sk_HDAY2(phi_val,Temp_val,z_yuk,k_val(i1))
    End Do
  End Subroutine

  Real * 8 Function spinodal_T_HDAY(phi_val,z_yuk)
    Implicit None
    Real * 8, Intent(in) :: phi_val, z_yuk
    spinodal_T_HDAY= 8.d0 * phi_val * Sk_HD_ros(phi_val,0.d0) / z_yuk
  End

  Real * 8 Function spinodal_T_HSAY(phi_val,z_yuk,vw_option)
    Implicit None
    Real * 8, Intent(in) :: phi_val, z_yuk
    Logical, intent(in) :: vw_option
    Real * 8 :: phi_dum
    if (vw_option .eqv. .true. ) then
      phi_dum = phi_val*(1.d0 - (phi_val / 16.d0))
    else
      phi_dum = phi_val
    end if
    spinodal_T_HSAY = 24.d0 * phi_dum * (z_yuk + 1.d0) * Sk_HS_py(vw_option,phi_val,0.d0) / (z_yuk ** 2)
  end

  Real * 8 Function blip_TLJ(Temp_val,dimension,nu_val)
    Use integration_module
    Implicit None
    Real * 8, Intent(in) :: Temp_val,dimension,nu_val
    Real * 8 :: dx, x
    Real * 8, Dimension(:), Allocatable :: Integrand
    Integer :: nx, i1
    dx=1.d-6
    nx= Int(1.d0/dx)
    Allocate(Integrand(nx))
    Do i1=1, nx
      x=i1*dx
      Integrand(i1) = ( (x) ** (dimension-1.d0) ) * exp(- ( (x ** - (2 * nu_val ) ) &
      &- 2.d0 * (x ** (-nu_val)) + 1.d0) / Temp_val )
    End Do
    blip_TLJ=1.d0 - dimension * DIntegral_simpson_3_8(Integrand,nx,dx)
    blip_TLJ=blip_TLJ**( 1.d0 / dimension )
    Deallocate(Integrand)
  End

  Subroutine Calc_gr_HDAY(grs,rs,kmax_val, phi_val, Temp_val, z_yuk, np, k)
    Implicit None
    Real * 8, Dimension(:), Intent(out) :: grs
    Real * 8, Dimension(:), Intent(out) :: rs
    Real * 8, Intent(in) :: kmax_val
    REal * 8, Intent(in) :: phi_val, Temp_val, z_yuk
    Real * 8, Dimension(:), Intent(out) :: k
    Integer, Intent(in) :: np
    Real * 8 :: K_max, R_max, ri,kj, J0,J1
    Integer :: i1, i2
    Do i1=1, np
      k(i1)= musj(i1) * kmax_val /musj(np)
      rs(i1) = musj(i1) / kmax_val
    End Do

    Do i1=1,np
      grs(i1)=0.d0
      Do i2=1,np-1
        J0=BESSEl_JN(0,(k(i2)*rs(i1)))
        J1=BESSEl_JN(1,(k(i2)*rs(np)))
        J1= J1 * J1
        grs(i1)= grs(i1) + (Sk_HDAY2(phi_val,Temp_val,z_yuk,k(i1)) * J0 / J1 )
      End Do
      grs(i1) = grs(i1) / (4.d0 * Atan(1.d0) * rs(np) * rs(np) )
    End Do
  End Subroutine

  Subroutine Calc_Sk_HDAY_low_rho(np,musj,rs,ks, R_max, phi_val, Temp_val, Z_yuk, grV, SkV)
    Implicit None
    Real * 8, Dimension(:), Intent(out) :: grV, SkV
    Real * 8, Dimension(:) :: rs, ks
    Real * 8, Dimension(:), Intent(in) :: musj
    Real * 8, Intent(in) :: phi_val, Temp_val, Z_yuk
    Integer, Intent(in) :: np
    Real * 8 :: R_max, J0 , J1, kn, delta_rho
    Integer :: i1, i2

    !Calculation of the wave vectors, radial distances and radial distribution
    Do i1=1, np
      ks(i1) = musj(i1) / R_max
      rs(i1) = musj(i1) * R_max / musj(np)
      grV(i1) = low_rho_gr_HDAY(rs(i1), Temp_val, Z_yuk)
    End Do
    kn = ks(np)
    !Calculation of the structure factor for each wave vector
    !S(k) = 1 + (ρ4π/k) ∫ g(r) Sin(kr) dr
    Do i1=1,np
      SkV(i1) = 0.d0
      !J0=rs(1)*grV(1)*BESSEl_JN(0,(ks(i1) * rs(1)))
      Do i2=2,np!-1
        delta_rho= rs(i2) - rs(i2-1)
        J0=BESSEl_JN(0,(ks(i1) * rs(i2)))
        J1=BESSEl_JN(1,(kn * rs(i2)))
        J1= J1 * J1
        SkV(i1) = SkV(i1) + ( (grV(i2)-1.d0) * J0 / J1 )
        !J1=rs(i2)*(grV(i2)-1.d0)*BESSEl_JN(0,(ks(i1) * rs(i2)))
        !SkV(i1) = SkV(i1) + 0.5d0 * delta_rho * (J0+J1)
        !J0=J1
      End Do
      SkV(i1) = 1.d0 + (16.d0 * phi_val * SkV(i1) / (kn * kn ))
      !SkV(i1) = 1.d0 + (8.d0 * phi_val * SkV(i1))
    End Do

  End Subroutine

 Real * 8 Function low_rho_gr_HDAY(X_val,Temp_val,Z_val)
   Implicit None
   Real * 8 :: X_val, Temp_val, Z_val
   If (X_val < 1.d0 ) Then
     low_rho_gr_HDAY =0.d0
   Else
     low_rho_gr_HDAY = exp( exp(- Z_val * (X_val-1) ) / (Temp_val * X_val) )
   End If
 End

 Subroutine Calc_hk_ck(Sk_val,knp_val,rho_val,hk_v,ck_v)
    Implicit None
    Real * 8, Dimension(:), Intent(in) :: Sk_val
    Real * 8, Dimension(:), Intent(out) :: hk_v, ck_v
    Integer, Intent(in) :: knp_val
    Real * 8 :: rho_val
    hk_v = (Sk_val - 1.d0 ) / rho_val
    ck_v = (1.d0 - (1.d0 / Sk_val)  )/rho_val
  End Subroutine

  Subroutine Calc_FT_2D_radial(rnp,Fr,Fk,rv,kv)
    Implicit None
    Real * 8, Dimension(:) :: Fr, Fk, rv, kv
    Integer :: rnp
    Integer :: i1,i2
    Real * 8 :: kn,kn2, int_const, J0, J1
    kn=kv(rnp)
    kn2=kn*kn
    int_const =  4.d0* 4.d0 * Atan(1.d0) / kn2
    Do i1=1,rnp
      Fk(i1) = 0.d0
      Do i2=1,rnp-1
        J0=BESSEl_JN(0,(kv(i1) * rv(i2)))
        J1=BESSEl_JN(1,(kn * rv(i2)))
        J1= J1 * J1
        Fk(i1) = Fk(i1) + ( (Fr(i2)) * J0 / J1 )
      End Do
      Fk(i1) = int_const * Fk(i1)
    End Do

  End Subroutine

  Subroutine Calc_IFT_2D_radial(rnp,Fr,Fk,rv,kv)
    Implicit None
    Real * 8, Dimension(:) :: Fr, Fk, rv, kv
    Integer :: rnp
    Integer :: i1,i2
    Real * 8 :: rn,rn2, int_const, J0, J1
    rn=rv(rnp)
    rn2=rn*rn
    int_const = 1.d0 / (rn2 * 4.d0 * Atan(1.d0))
    Do i1=1,rnp
      Fr(i1) = 0.d0
      Do i2=1,rnp-1
        J0=BESSEl_JN(0,(kv(i2) * rv(i1)))
        J1=BESSEl_JN(1,(rn * kv(i2)))
        J1= J1 * J1
        Fr(i1) = Fr(i1) + ( (Fk(i2)) * J0 / J1 )
      End Do
      Fr(i1) = int_const * Fr(i1)
    End Do

  End Subroutine

  Subroutine Calc_sk_HDAY_Sharma_Sharma_V2(knp_val,rv,kv, phi_val, Temp_val, Z_yuk, grV, SkV)
   Implicit None
   Real * 8, Dimension(:) :: grV, SkV
   Real * 8, Dimension(:) :: rv, kv
   Real * 8, Intent(in) :: phi_val, Temp_val, Z_yuk
   Integer, Intent(in) :: knp_val
   Real * 8, Dimension(:), Allocatable :: ck_v, hk_v, cr_v, hr_v, Sk_aux, Yk_v,Yr_v
   Real * 8 :: rho_val,dum1
   Integer :: i1, i2
   Allocate(ck_v(knp_val))
   Allocate(hk_v(knp_val))
   Allocate(cr_v(knp_val))
   Allocate(hr_v(knp_val))
   Allocate(Sk_aux(knp_val))
   rho_val= phi_val / Atan(1.d0)
   Call Calc_Sk_hd_ros_mono(phi_val, knp_val,kv,SkV)
   Call Calc_hk_ck(Skv,knp_val,rho_val,hk_v,ck_v)
   cr_v=0.d0
   Do i1=1, knp_val
     if( rv(i1) >= 1.d0 ) cr_v(i1) = exp(exp(-Z_yuk * (rv(i1) - 1.d0) ) / (Temp_val*rv(i1))) - 1.d0
   End Do
   Call Calc_FT_2D_radial(knp_val,cr_v,Sk_aux,rv,kv)
   ck_v = ck_v + Sk_aux
   SkV= 1.d0 / (1.d0 - rho_val * ck_v)
   !Open(unit=10, file="ck_PY.dat", Status="Replace")
   !Do i1=1, knp_val
   ! Write(10,*), kv(i1), ck_v(i1), Sk_aux(i1)
   !End Do
   !Close(10)
   Deallocate(ck_v)
   Deallocate(hk_v)
   Deallocate(cr_v)
   Deallocate(hr_v)
   Deallocate(Sk_aux)
 End Subroutine

End Module
