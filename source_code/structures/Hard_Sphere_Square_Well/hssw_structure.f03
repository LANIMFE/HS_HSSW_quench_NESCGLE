Module hssw_structure
  Use hs_structure
  Use math_constants
  Use quadratures
  Implicit None
  Contains

    !Function that calculates the direct correlation function of a
    !Square Well System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    !for the hard core part
    ! c(k,ϕ,λ)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !square well range "λ>1"
    Function c_hssw_shvw_aux(phi,T,lambda,k) result(c_aux)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: lambda !square well range
      Real * 8, Intent (in) :: k !Value of the wave vector
      Real * 8 :: c_aux
      Real * 8 :: k_func,ckdum,chs
      Real * 8 :: k2,lk,l5,l3
      chs = c_hs_vw(phi,k)
      k2 = k * k
      if (k > 0.075d0) then
        lk = lambda * k
        c_aux = ((cos(k) - lambda * cos(lk)) / k) + ((sin(lk) - sin(k)) / k2)
        c_aux = c_aux / k
      else
        l3 = lambda ** 3
        l5 = lambda ** 5
        c_aux = (1.d0/3.d0) * (l3 - 1.d0) - (1.d0/30.d0) * (l5 - 1.d0) * k2
      end if
      c_aux = 4.d0 * pi * c_aux / T
      return
    End Function


  !Function that calculates the direct correlation function of a Hard Sphere +
  !Square Well System using Sharma-Sharma approximation Scheme and Verlet-Weiss
  !for the hard core part
  ! c(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function c_hssw_shvw(phi,T,lambda,k) result(c)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: c
    Real * 8 :: chs
    chs = c_hs_vw(phi,k)
    c = chs + c_hssw_shvw_aux(phi,T,lambda,k)
    return
  End Function

  !Function that calculates the inverse structure factor of a Hard Sphere +
  !Square Well System using Sharma-Sharma approximation Scheme and Verlet-Weiss
  ! 1/S(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function is_hssw_shvw(phi,T,lambda,k) result(is)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: is,phi_vw, chs
    phi_vw = phi*(1.d0 - (phi / 16.d0))
    chs = c_hs_vw(phi,k)
    is = (phi_vw * chs) + (phi * c_hssw_shvw_aux(phi,T,lambda,k))
    is = 1.d0 - 6.d0 * is / pi
    return
  End Function

  !Function that calculates the structure factor of a Hard Sphere +
  !Square Well System using Sharma-Sharma approximation Scheme and Verlet-Weiss
  ! S(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function s_hssw_shvw(phi,T,lambda,k) result(s)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: s
    s = 1.d0 / is_hssw_shvw(phi,T,lambda,k)
    return
  End Function

  !Function that calculates the spinodal temperature of a Hard sphere +
  !Square Well system using Sharma-Sharma approximation, it needs the information of:
  ! "phi" -> volume fraction of the system
  ! " lambda" -> length of the square well potential
  Function Ts_hssw_shvw(phi,lambda) result(T)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: lambda  !length of the square well potential
    Real * 8 :: T
    T = 8.d0 * phi * (lambda ** 3 - 1) * s_hs_vw(phi,0.d0)
  End Function


  !Function that calculates the direct correlation function of a
  !Square Well System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
  !for the hard core part
  ! c(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function c_hssw_gsvw_aux(phi,T,lambda,k) result(c_aux)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: c_aux
    Real * 8 :: k_func,ckdum,chs
    Real * 8 :: k2,lk,l5,l3
    k2 = k * k
    if (k > 0.075d0) then
      lk = lambda * k
      c_aux = ((cos(k) - lambda * cos(lk)) / k) + ((sin(lk) - sin(k))/k2)
      c_aux = c_aux / k
    else
      l3 = lambda ** 3
      l5 = lambda ** 5
      c_aux = (1.d0/3.d0) * (l3 - 1.d0) - (1.d0/30.d0) * (l5 - 1.d0) * k2
    end if
    c_aux = 4.d0 * pi * c_aux * ( exp(1.d0/T) - 1.d0)
    return
  End Function

  !Function that calculates the direct correlation function of a Hard Sphere +
  !Square Well System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
  !for the hard core part
  ! c(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function c_hssw_gsvw(phi,T,lambda,k) result(c)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: c
    Real * 8 :: k_func,ckdum,chs
    Real * 8 :: k2,lk,l5,l3
    chs = c_hs_vw(phi,k)
    c = c_hs_vw(phi,k) + c_hssw_gsvw_aux(phi,T,lambda,k)
    return
  End Function

  !Function that calculates the inverse structure factor of a Hard Sphere +
  !Square Well System using Sharma-Sharma approximation Scheme and Verlet-Weiss
  ! 1/S(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function is_hssw_gsvw(phi,T,lambda,k) result(is)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: is,phi_vw, chs
    phi_vw = phi*(1.d0 - (phi / 16.d0))
    chs = c_hs_vw(phi,k)
    is = (phi_vw * chs) + (phi *  c_hssw_gsvw_aux(phi,T,lambda,k))
    is = 1.d0 - 6.d0  * is / pi
    return
  End Function

  !Function that calculates the structure factor of a Hard Sphere +
  !Square Well System using Sharma-Sharma approximation Scheme and Verlet-Weiss
  ! S(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function s_hssw_gsvw(phi,T,lambda,k) result(s)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: s
    s = 1.d0 / is_hssw_gsvw(phi,T,lambda,k)
    return
  End Function

  !Function that calculates the structure factor of a Hard Sphere +
  !Square Well System using Boltzmann approximation Scheme
  ! S(k,ϕ,λ)
  !wave-vector "k" for "k>=0"
  !density "ϕ" for "ϕ>0"
  !temperature "T" for "T>0"
  !square well range "λ>1"
  Function s_hssw_bolt(phi,T,lambda,k) result(s)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: T !Temperature
    Real * 8, Intent(in) :: lambda !square well range
    Real * 8, Intent (in) :: k !Value of the wave vector
    Real * 8 :: s, x, lk, k3, l3
    k3=k * k * k
    x = exp(1.d0/T) - 1.d0
    lk = lambda * k
    s = 0.d0
    If (k > 0.075d0) Then
      s = sin(lk) - sin(k) + k * cos(k) - lk * cos(lk)
      s = s * x
      s = s + k * cos(k) - sin(k)
      s = s / k3
    Else
      l3 = lambda * lambda * lambda
      s = (x * (l3 - 1.d0) - 1.d0) / 3.d0
    End If
    s = 1.d0 + (24.d0 * phi * s )
    return
  End Function

  !Santos SW S(k)
  Function s_hssw_stos (phi,T,lambda,k) Result(s)
    Implicit None
    Integer :: i1
    Real * 8, Intent(in) :: phi, T, lambda, k
    Real * 8 :: s
    Real * 8 :: A0,A01,L1,L2,S1,S2,S3,sw_tau !F(t) Parameters
    Real * 8 :: alpha1,alpha2,alpha3,beta1,beta2,beta3,beta4,gamma0,gamma1,gamma2 !L2 Parameters
    Real * 8 :: Lam1,Lam2,Lam3,Lam4,Lam5,Lam6,tau1,tau2 !convenient dummies variables
    Complex * 16 :: Ftp,Ftm,Gtp,Gtm !F(t),F(-t),G(t),G(-t)
    Real * 8 :: Ftdum1,Ftdum2,Ftdum3
    real * 8 dum1,dum2,dum3
    Complex * 16 :: im
    Real * 8 dumeta
    dumeta=phi
    im=(0.d0,1.d0)
    !Allocate (Ftp(kpoints))
    !Allocate (Ftm(kpoints))
    !Allocate (Gtp(kpoints))
    !Allocate (Gtm(kpoints))
    Lam1=lambda
    Lam2=Lam1**2
    Lam3=Lam1**3
    Lam4=Lam1**4
    Lam5=Lam1**5
    Lam6=Lam1**6
    !A0
    A0=exp(1.d0/T)-1.d0
    !A0 prima
    A01=A0*((Lam1-1.d0)**2)
    !print*, "A01=",A01
    !square well tau
    sw_tau=1.d0/(12.d0*A0*Lam1*(Lam1-1.d0))
    tau1=1.d0/sw_tau
    !print*, "tau1=",tau1
    tau2=1.d0/(sw_tau**2)
    !print*, "tau2=",tau2
    !alpha 1
    alpha1=2.d0*A01*(2.d0+Lam1)+((1.d0+(4.d0*Lam1)+Lam2-3.d0*A01)*tau1/6.d0)
    !print*, "alpha1", alpha1
    !alpha2
    alpha2=3.d0+A01*(7.d0+2.d0*Lam1-3.d0*Lam2)-((7.d0+Lam1+16.d0*Lam2+A01*&
    &(23.d0+15.d0*Lam1+15.d0*Lam2+7.d0*Lam3))*tau1/12.d0)
    !print*, "alpha2", alpha2
    !alpha3
    alpha3=-2.d0-2.d0*A01*(1.d0+2.d0*Lam1+3.d0*Lam2)+((7.d0+Lam1-2.d0*Lam2+A01*&
    &(7.d0+15.d0*Lam1+21.d0*Lam2+11.d0*Lam3+6.d0*Lam4))*tau1/6.d0)
    !print*, "alpha3", alpha3
    !beta1
    beta1=-4.d0-4.d0*A01*(2.d0+Lam1)-(tau2/3.d0)+((5.d0+2.d0*Lam1-Lam2+3.d0*A01)*tau1/3.d0)
    !print*, "beta1", beta1
    !beta2
    beta2=6.d0+6.d0*A01*(3.d0+2.d0*Lam1+Lam2)+4.d0*(A01**2)*((2.d0+Lam1)**2)&
    &-((9.d0*(3.d0+Lam1)+A01*(59.d0+51.d0*Lam1+3.d0*Lam2-5.d0*Lam3)&
    &+12.d0*(A01**2)*(2.d0+Lam1))*tau1/6.d0)+((31.d0+8.d0*Lam1+18.d0*Lam2&
    &-4.d0*Lam3+Lam4+12.d0*A01*(5.d0+Lam1)+9.d0*(A01**2))*tau2/36.d0)
    !print*, "beta2", beta2
    !beta3
    beta3=-4.d0-12.d0*A01*(1.d0+Lam1+Lam2)-4.d0*(A01**2)*(2.d0+Lam1)*(1.d0&
    &+2.d0*Lam1+3.d0*Lam2)+((3.d0*(4.d0+Lam1+Lam2)+A01*(29.d0&
    &+36.d0*Lam1+30.d0*Lam2+7.d0*Lam3-3.d0*Lam4)+(A01**2)&
    &*(17.d0+43.d0*Lam1+30.d0*Lam2+Lam3-Lam4))*tau1/3.d0)-((29.d0+13.d0*Lam1&
    &+27.d0*Lam2+7.d0*Lam3-4.d0*Lam4+A01*(68.d0+80.d0*Lam1+80.d0*Lam2+16.d0*Lam3&
    &+7.d0*Lam4+Lam5)+3.d0*(A01**2)*(13.d0+13.d0*Lam1+Lam2-3.d0*Lam3))*tau2/36.d0)
    !print*, "beta3", beta3
    !beta4 BE WARE, MODIFIED FROM THE ARTICLE
    dum1=(1.d0+A01*(1.d0+2.d0*Lam1+3.d0*Lam2))**2
    dum2=2.d0*A01*(7.d0+15.d0*Lam1+15.d0*Lam2+5.d0*Lam3+6.d0*Lam4)
    dum3= A01*A01*(1.d0+2.d0*Lam1+3.d0*Lam2)*(7.d0+15.d0*Lam1+3.d0*Lam2-Lam3)
    dum2=7.d0+Lam1+4.d0*Lam2+dum2+dum3
    beta4=dum1-(dum2*tau1/6.d0)
    dum1=49.d0+38.d0*Lam1+9.d0*Lam2+8.d0*Lam3+16.d0*Lam4
    dum2=2.d0*A01*(49.d0+88.d0*Lam1+112.d0*Lam2+80.d0*Lam3+35.d0*Lam4-4.d0*Lam5)
    dum3=A01*A01*(49.d0+138.d0*Lam1+219.d0*Lam2+124.d0*Lam3+27.d0*Lam4+18.d0*Lam5+Lam6)
    dum1=dum1+dum2+dum3
    dum1=dum1*tau2/144.d0
    beta4=beta4+dum1
    !print*, "beta4",beta4
    !gamma0
    gamma0=1.d0+Lam1-(tau1/6.d0)
    !print*, "gamma0",gamma0
    !gamma1
    gamma1=2.d0*(1.d0+Lam1-Lam2)-((3.d0+Lam3)*tau1/6.d0)
    !print*, "gamma1",gamma1
    !gamma2
    gamma2=Lam1*(-4.d0*Lam1+((2.d0+3.d0*Lam1+Lam2+Lam3)*tau1/3.d0))
    !print*, "gamma2",gamma2
    !L2
    L2=-1.d0+(alpha1*dumeta)+(alpha2*dumeta**2)+(alpha3*dumeta**3)+((1.d0+2.d0*dumeta)*&
    &dsqrt(1.d0+(beta1*dumeta)+(beta2*dumeta**2)+(beta3*dumeta**3)+(beta4*dumeta**4)))
    L2=L2/(12.d0*dumeta*(gamma0+gamma1*dumeta+gamma2*dumeta**2))
    !print*, "L2",L2
    !L1
    L1=1.d0+0.5d0*dumeta+2.d0*dumeta*(1.d0+Lam1+Lam2)*L2-0.5d0*dumeta*A01*(3.d0+2.d0*Lam1+Lam2)
    L1=L1/(1.d0+2.d0*dumeta)
    !print*, "L1",L1
    !S1
    S1=(-1.5d0+2.d0*(1.d0+Lam1+Lam2)*L2-0.5d0*(3.d0+2.d0*Lam1+Lam2)*A01)!Care inconsistency
    S1=S1*(dumeta/(1.d0+2.d0*dumeta))
    !print*, "S1",S1
    !S2
    S2=(-1.d0+dumeta+2.d0*(1.d0-2.d0*dumeta*Lam1*(1.d0+Lam1))*L2-(1.d0-dumeta*((1.d0+Lam1)**2))*A01)
    S2=S2/(2.d0*(1.d0+2.d0*dumeta))
    !print*, "S2",S2
    !S3
    S3=(-(((1.d0-dumeta)**2)/(12.d0*dumeta))-(0.5d0*(Lam1+1.d0)-dumeta*Lam2)*L2)
    S3=S3+(4.d0+2.d0*Lam1-dumeta*(3.d0*Lam2+2.d0*Lam1+1.d0))*(A01/12.d0)
    S3=S3/(1.d0+2.d0*dumeta)
    !print*, "S3",S3
    !F(t)->F(ik)
    Ftdum1=L1+(L2/(Lam1-1.d0))-A0*(Lam1-1.d0)
    Ftdum2=L2/(Lam1-1.d0)
    Ftdum3=Lam1-1.d0
    !Do i1=1, kpoints
      !Ftp(i1)=CMPLX(1.d0+A0,0)+CMPLX(0.D0,Ftdum1*k(i1))-(CMPLX(A0,Ftdum2*k)*exp(CMPLX(0.d0,-Ftdum3*k)))
      !Ftp(i1)=-Ftp(i1)/CMPLX((1.d0-S2*(k**2)),(S1*k-S3*(k**3)))/(12.d0*eta)
      !Ftm(i1)=CMPLX(1.d0+A0,0)+CMPLX(0.D0,-Ftdum1*k)-(CMPLX(A0,-Ftdum2*k)*exp(CMPLX(0.d0,Ftdum3*k)))
      !Ftm(i1)=-Ftp(i1)/CMPLX(12.d0*eta*(1.d0-S2*(k**2)),-12.d0*eta*(S1*k-S3*(k**3)))
      Ftp=(1.d0+A0+Ftdum1*(im*k))-(A0+Ftdum2*(im*k))*exp(-Ftdum3*(im*k))
      Ftp=Ftp/(1.d0+(S1*(im*k))+(S2*(im*k)**2)+(S3*(im*k)**3))
      Ftp=-Ftp/(12.d0*dumeta)
      Ftm=(1.d0+A0+Ftdum1*(-im*k))-(A0+Ftdum2*(-im*k))*exp(-Ftdum3*(-im*k))
      Ftm=Ftm/(1.d0+(S1*(-im*k))+(S2*(-im*k)**2)+(S3*(-im*k)**3))
      Ftm=-Ftm/(12.d0*dumeta)
      !G(t)->G(ik)
      !Gtp(i1)=CMPLX(0.d0,k)*Ftp(i1)*exp(CMPLX(0.d0,-k))/(CMPLX(1.d0,0.d0)+12.d0*eta*Ftp(i1)*exp(CMPLX(0.d0,-k)))
      Gtp=(im*k)*Ftp*exp(-im*k)/(1.d0+12.d0*dumeta*Ftp*exp(-im*k))
      !G(-t)->G(-ik)
      !Gtm(i1)=CMPLX(0.d0,-k)*Ftm(i1)*exp(CMPLX(0.d0,k))/(CMPLX(1.d0,0.d0)+12.d0*eta*Ftm(i1)*exp(CMPLX(0.d0,k)))
      Gtm=(-im*k)*Ftm*exp(im*k)/(1.d0+12.d0*dumeta*Ftm*exp(im*k))

      s=Realpart((Gtp-Gtm)/(im*k))
      s=1.d0-(12.d0*dumeta*s)
      !s=Realpart(Gtm(i1))
    !End Do
  End Function


  !Function that calculates the spinodal temperature of a Hard sphere + Square Well
  !system using Sharma-Sharma approximation, it needs the information of:
  ! "phi" -> volume fraction of the system
  ! "lambda" -> length of the square well potential
  Function Ts_hssw_gsvw(phi,lambda) result(T)
    Implicit None
    Real * 8, Intent(in) :: phi !Volume fraction value
    Real * 8, Intent(in) :: lambda !length of the square well potential
    Real * 8 :: T
    Real * 8 :: Tshvw
    Tshvw = Ts_hssw_shvw(phi,lambda)
    T = 1.d0/log(1.d0 + (1.d0/Tshvw))
  End Function

  Function u_HSSW(r,lambda) result(u)
    Implicit None
    Real * 8, Intent(in) :: r,lambda
    Real * 8 :: u
    If ( 1.d0 <= r .and. r <= lambda) Then
      u = -1.d0
    Else
      u = 0.d0
    End If
  End Function

  Function gr_tail_hssw_gsvw(k_max,lambda,T,rho,r) result(g_tail)
    Implicit None
    Real * 8, Intent(in) :: k_max, lambda, T, rho, r !
    Real * 8 :: g_tail !
    Real * 8 :: dum1,dum2,dum3,dum4,phi,Ac1,Ac2,Ac3,As1,As2,Bc1,Bs1,A1,A2
    Real * 8 :: dumAplus,dumAminus,dumBplus,dumBminus,dumA,SiAp,SiAm,SiBp,SiBm,Sirk
    Real * 8 :: sp,sp2,sp3,sp4,sm,sm2,sm3,sm4,slp,slm
    Real * 8 :: si, ci
    phi = pi * rho / 6.d0
    phi = phi*(1.d0 - (phi / 16.d0))  !Verlet-Weiss correction
    dum1 = - ( ( 1.d0 + 2.d0 * phi ) ** 2 )
    dum2 = 6.d0 * phi * ( ( 1.d0 + ( phi / 2.d0 ) ) ** 2 )
    dum3 = - 0.5d0 * phi * ( ( 1.d0 + ( 2.d0 * phi ) ) ** 2 )
    dum4 = ( 1.d0 - phi )
    dum4 = dum4 ** 4
    dum1 = dum1 / dum4 !alpha
    dum2 = dum2 / dum4  !beta
    dum3 = dum3 / dum4  !delta
    Bs1 = exp(1.d0/T) - 1.d0
    Ac1 = Bs1 - dum1 - dum2 - dum3
    Ac2 = 2.d0 * dum2 + 12.d0 * dum3
    Ac3 = -24.d0 * dum3
    As1 = -Bs1 + dum1 + 2.d0 * dum2 + 4.d0 * dum3
    As2 = - 24.d0 * dum3
    A1 = -2 * dum2
    A2 = 24.d0 * dum3
    Bc1 = - lambda * Bs1
    sp = 1.d0 + r
    sp2 = sp * sp
    sp3 = sp2 * sp
    sp4 = sp2 * sp2
    sm = 1.d0 - r
    sm2 = sm * sm
    sm3 = sm2 * sm
    sm4 = sm2 * sm2
    slp = lambda + r
    slm = lambda - r
    dumAplus = Ac1 + (As1 * sp) - (Ac2 * sp2 / 2.d0) - (As2 * sp3 / 6.d0) + (Ac3 * sp4 / 24.d0)
    dumAminus = Ac1 + (As1 * sm) - (Ac2 * sm2 / 2.d0) - (As2 * sm3 / 6.d0) + (Ac3 * sm4 / 24.d0)
    dumBplus = Bc1 + (Bs1 * slp)
    dumBminus = Bc1 + (Bs1 * slm)
    dumA = (A2 *(r ** 4)/12.d0) - (A1*(r**2))
    call cisi(sp * k_max, ci, si)
    SiAp = (pi/2.d0) - si
    call cisi(slp * k_max, ci, si)
    SiBp = (pi/2.d0) - si
    !Conditional dependendt of r values
    SiAm = 0.d0
    SiBm = - pi / 2.d0
    !if (r == 0) SiAm =-pi/2.d0
    if (r > 1.d0) Then
      SiAm = -pi/2.d0
      if (r < lambda) SiBm = -SiBm
    End If
    call cisi(sm * k_max, ci, si)
    SiAm = SiAm - si
    call cisi(slm * k_max, ci, si)
    SiBm = SiBm - si
    !
    call cisi(r * k_max, ci, si )
    Sirk = (pi/2.d0) - si
    !g_tail = A * ( r8_si( ( 1.d0 + r) * k_max ) - r8_si( (1.d0 - r) * k_max ) )
    !g_tail = g_tail + (B * (r8_si( ( lambda - r ) * k_max ) -  r8_si( ( lambda + r ) * k_max ) ))
    !g_tail = - g_tail / (pi * r) !
    !if (r > 1 .and. r < lambda ) g_tail =  g_tail - (B / r)
    !if (r > 1) g_tail = g_tail + (A/r)
    g_tail = dumAplus * SiAp + dumBplus * SiBp - dumAminus * SiAm - dumBminus * SiBm + dumA * Sirk
    g_tail = g_tail / (pi*r)
    return
  End Function

  Function gr_tail_hssw_shvw(k_max,lambda,T,rho,r) result(g_tail)
    Implicit None
    Real * 8, Intent(in) :: k_max, lambda, T, rho, r !
    Real * 8 :: g_tail !
    Real * 8 :: dum1,dum2,dum3,dum4,phi,Ac1,Ac2,Ac3,As1,As2,Bc1,Bs1,A1,A2
    Real * 8 :: dumAplus,dumAminus,dumBplus,dumBminus,dumA,SiAp,SiAm,SiBp,SiBm,Sirk
    Real * 8 :: sp,sp2,sp3,sp4,sm,sm2,sm3,sm4,slp,slm
    Real * 8 :: si, ci
    phi = pi * rho / 6.d0
    phi = phi*(1.d0 - (phi / 16.d0))  !Verlet-Weiss correction
    dum1 = - ( ( 1.d0 + 2.d0 * phi ) ** 2 )
    dum2 = 6.d0 * phi * ( ( 1.d0 + ( phi / 2.d0 ) ) ** 2 )
    dum3 = - 0.5d0 * phi * ( ( 1.d0 + ( 2.d0 * phi ) ) ** 2 )
    dum4 = ( 1.d0 - phi )
    dum4 = dum4 ** 4
    dum1 = dum1 / dum4 !alpha
    dum2 = dum2 / dum4  !beta
    dum3 = dum3 / dum4  !delta
    Bs1 = 1.d0 / T
    Ac1 = Bs1 - dum1 - dum2 - dum3
    Ac2 = 2.d0 * dum2 + 12.d0 * dum3
    Ac3 = -24.d0 * dum3
    As1 = -Bs1 + dum1 + 2.d0 * dum2 + 4.d0 * dum3
    As2 = - 24.d0 * dum3
    A1 = -2 * dum2
    A2 = 24.d0 * dum3
    Bc1 = - lambda * Bs1
    sp = 1.d0 + r
    sp2 = sp * sp
    sp3 = sp2 * sp
    sp4 = sp2 * sp2
    sm = 1.d0 - r
    sm2 = sm * sm
    sm3 = sm2 * sm
    sm4 = sm2 * sm2
    slp = lambda + r
    slm = lambda - r
    dumAplus = Ac1 + (As1 * sp) - (Ac2 * sp2 / 2.d0) - (As2 * sp3 / 6.d0) + (Ac3 * sp4 / 24.d0)
    dumAminus = Ac1 + (As1 * sm) - (Ac2 * sm2 / 2.d0) - (As2 * sm3 / 6.d0) + (Ac3 * sm4 / 24.d0)
    dumBplus = Bc1 + (Bs1 * slp)
    dumBminus = Bc1 + (Bs1 * slm)
    dumA = (A2 *(r ** 4)/12.d0) - (A1*(r**2))
    call cisi(sp * k_max, ci, si)
    SiAp = (pi/2.d0) - si
    call cisi(slp * k_max, ci, si)
    SiBp = (pi/2.d0) - si
    !Conditional dependendt of r values
    SiAm = 0.d0
    SiBm = - pi / 2.d0
    !if (r == 0) SiAm =-pi/2.d0
    if (r > 1.d0) Then
      SiAm = -pi/2.d0
      if (r < lambda) SiBm = -SiBm
    End If
    call cisi(sm * k_max, ci, si)
    SiAm = SiAm - si
    call cisi(slm * k_max, ci, si)
    SiBm = SiBm - si
    !
    call cisi(r * k_max, ci, si )
    Sirk = (pi/2.d0) - si
    !g_tail = A * ( r8_si( ( 1.d0 + r) * k_max ) - r8_si( (1.d0 - r) * k_max ) )
    !g_tail = g_tail + (B * (r8_si( ( lambda - r ) * k_max ) -  r8_si( ( lambda + r ) * k_max ) ))
    !g_tail = - g_tail / (pi * r) !
    !if (r > 1 .and. r < lambda ) g_tail =  g_tail - (B / r)
    !if (r > 1) g_tail = g_tail + (A/r)
    g_tail = dumAplus * SiAp + dumBplus * SiBp - dumAminus * SiAm - dumBminus * SiBm + dumA * Sirk
    g_tail = g_tail / (pi*r)
    return
  End Function

  Function absc(h)
    Complex * 16, intent(in) :: h
    Real * 8 :: absc
    absc = abs(real(h)) + abs(aimag(h))
  End function

  Subroutine cisi(x,ci,si)
    Implicit None
    Real * 8, intent(in) :: x
    Real * 8, intent (out) :: ci,si
    Integer, parameter :: MAXIT = 100
    Real * 8, parameter :: EPS=6.d-8, FPMIN=1.d-30, TMIN=2.d0
    Integer :: i,k
    Real * 8 :: a,err,fact,sign,sum,sumc,sums,t,term,absc
    Complex * 16 :: h,b,c,d,del
    Logical :: odd, convergence
    absc(h) = abs(real(h))+abs(aimag(h))
    t = abs(x)
    if (t.eq.0.) then
      si = 0.d0
      ci = -1.d0/FPMIN
      return
    end if
    if (t .gt. TMIN) then
      b = cmplx(1.d0,t)
      c = 1.d0/FPMIN
      d = 1.d0/b
      h = d
      i = 2
      convergence = .False.
      do while (convergence .eqv. .False.)
        a = -(i-1)**2
        b = b+2.d0
        d = 1.d0/(a*d+b)
        c = b+a/c
        del = c*d
        h = h*del
        i = i + 1
        if ( (absc(del-1.d0) .lt. EPS) .or. ( i == MAXIT) ) convergence = .True.
      end do
  !      if () pause 'cf failed in cisi'
        if (i .ge. MAXIT) Print*, 'cf failed in cisi'
        h = cmplx(cos(t),-sin(t))*h
        ci = -real(h)
        si = pih+aimag(h)
        else
          if( t .lt. sqrt(FPMIN) )then
            sumc = 0.d0
            sums = t
          else
            sum  = 0.d0
            sums = 0.d0
            sumc = 0.d0
            sign = 1.d0
            fact = 1.d0
            odd  = .true.
            k = 1
            convergence = .False.
            do while (convergence .eqv. .False.)
              fact = fact * t / k
              term = fact / k
              sum  = sum + sign * term
              err  = term / abs(sum)
              if ( odd ) then
                sign = -sign
                sums = sum
                sum  = sumc
              else
                sumc = sum
                sum  = sums
              end if
              if( (err .lt. EPS) .or. (k == MAXIT) ) convergence = .True.
              odd = .not. odd
              k = k + 1
              end do
              !pause 'maxits exceeded in cisi'
            end if
            if (k .ge. MAXIT ) print*, 'maxits exceeded in cisi'
            si = sums
            ci = sumc+log(t)+EULER
          end if
          if( x .lt. 0.d0 ) si = -si
          return
  End Subroutine

  !This subroutine writes down the structure factor, the direct correlation function and the
	!inverse of the structure factor and the associated run file info
	Subroutine structure_writting_hssw(k,sk,kw,sys_cond,approx,quad,quadp_f,quadp_i,folder,id)
    Implicit None
		Real * 8, dimension(:), intent(in) :: k, sk, kw, sys_cond,quadp_f
		Integer, dimension(:) :: quadp_i
		Character(*), intent(in) :: folder, approx, quad
		Integer, intent(in) :: id
		Real * 8 :: phi, T, lambda
		Integer :: di_phi, i_T, di_T
		Integer, parameter :: di_phi_min=4, i_T_max = 3, di_T_min = 4, id_max=2
		Character (len=di_phi_min) :: phi_char
		Character (len=i_T_max+di_T_min+1) :: T_char
		Character (len=id_max) :: id_char
		Character (len=len(folder)+len(phi_char)+len(T_char)+len(id_char)+19) :: File_name
		Character (len=len(File_name)) :: info_File_name
		Integer :: i1
		!System conditions
		phi = sys_cond(1)
		T = sys_cond(2)
		lambda = sys_cond(3)
		!Calculating the associated integers for phi and T:
		i_T = int(T)
		di_T = int( (T- dble(i_T)) * (10.d0 ** di_T_min) )
		di_phi = int(phi * (10.d0 ** di_phi_min) )
		!Writting con character variables, formats have to be consistent with char lengths
		Write(phi_char,"(I4.4)") di_phi
		Write(T_char,"(I3.3,A1,I4.4)") i_T,"_",di_T
		Write(id_char,"(I2.2)") id
		File_name = folder//"S_phi_0_"//phi_char//"_T_"//T_char//"_id_"//id_char//".dat"
		info_File_name = folder//"S_phi_0_"//phi_char//"_T_"//T_char//"_id_"//id_char//".inf"
		Open(unit=100, file=File_name, status="Replace")
		Write (100,*) "# 						k 										S(k)											Weights"
		Do i1=1, size(k)
			Write(100,*) k(i1),sk,kw(i1)
		End Do
		Close(100)
		Open (unit=100, file=info_File_name, status="replace")
		Write(100,*) "###  Run Information for the file: "//File_name
		Write(100,*) "### System parameters ###"
		Write(100,*) "Volume Fraction phi = ", phi
		Write(100,*) "Temperature      T* = ", T
		Write(100,*) "Lambda              = ", lambda
		Write(100,*) "###############################"
		Write(100,*) "### Structure Approximation ###"
		Select case (approx)
			case ("SHVW")
				Write(100,*) "				Sharma-Sharma + Verlet Weiss (HSPY)"
			case ("GSVW")
				Write(100,*) "				Generalized Sharma-Sharma + Verlet Weiss (HSPY)"
			case Default
				Write(100,*) "				"//approx
		End Select
		Write(100,*) "###############################"
		Write(100,*) "### Quadratures ###"
		Select case (quad)
			case("CLCU")
				Write(100,*) "				Clenshaw-Curtis"
				Write(100,*) "quadrature number of points = ", quadp_i(1)
				Write(100,*) "quadrature minimum value    = ", quadp_f(1)
				Write(100,*) "quadrature maximum value    = ", quadp_f(2)
			case Default
				Write(100,*) "				"//quad
				Write(100,*) "quadrature number of points = ", quadp_i(1)
				Write(100,*) "quadrature minimum value    = ", quadp_f(1)
				Write(100,*) "quadrature maximum value    = ", quadp_f(2)
		End Select
		Write(100,*) "###############################"
		Write(100,*) "### Wave-Vectors ###"
		Write(100,*) "number of points     = ", size(k)
		Write(100,*) "k minimum value      = ", k(1)
		Write(100,*) "k maximum value      = ", k(size(k))
		Close(100)
	End Subroutine
  !This subroutine writes down the structure factor, the direct correlation function and the
	!inverse of the structure factor and the associated run file info
	Subroutine writting_hssw_p(sys_cond,approx,i_unit)
    Implicit None
		Real * 8, dimension(:), intent(in) :: sys_cond
		Integer :: i_unit
		Character(*), intent(in) :: approx
		Real * 8 :: phi, T, lambda
		Integer :: i1
		!System conditions
		phi = sys_cond(1)
		T = sys_cond(2)
		lambda = sys_cond(3)
		Write(i_unit,*) "### Hard Sphere + Square Well System parameters ###"
		Write(i_unit,*) "Volume Fraction phi = ", phi
		Write(i_unit,*) "Temperature      T* = ", T
		Write(i_unit,*) "Lambda              = ", lambda
		Write(i_unit,*) "###############################"
		Write(i_unit,*) "### Structure Approximation ###"
		Select case (approx)
			case ("SHVW")
				Write(i_unit,*) "				Sharma-Sharma + Verlet Weiss (HSPY)"
			case ("GSVW")
				Write(i_unit,*) "				Generalized Sharma-Sharma + Verlet Weiss (HSPY)"
			case Default
				Write(i_unit,*) "				"//approx
		End Select
    Write(i_unit,*) "###############################"
	End Subroutine

  !Pressure for the Hard Sphere Square Well system through compressibility of Sharma-Sharma
  Function P_HSSW_SH(T,phi,lambda) result(P)
    Implicit None
    Real * 8, intent(in) :: T, phi, lambda
    Real * 8 :: P
    P = 6.d0 * phi * T * (1.d0 + phi**2 - phi **3) / ( pi * (1.d0 - phi)**3 )
    P = P - (24.d0 * (lambda**3 - 1.d0) * phi ** 2 / pi)
  End Function

  !Pressure for the Hard Sphere Square Well system through compressibility of Generalized Sharma-Sharma
  Function P_HSSW_GSH(T,phi,lambda) result(P)
    Implicit None
    Real * 8, intent(in) :: T, phi, lambda
    Real * 8 :: P
    Real * 8 :: FT, phi2, phi3
    phi2 = phi  * phi
    phi3 = phi2 * phi
    FT = exp(1.d0/T) - 1.d0
    P = 6.d0 * phi * T * (1.d0 + phi + phi2 - phi3) / ( pi * (1.d0 - phi)**3 )
    P = P - (24.d0 * (lambda ** 3 - 1.d0) * phi2 * T * FT / pi)
  End Function

  Function P_HSSW_gr(T,rho,lambda,g,r) result(P)
    Implicit None
    Real * 8, intent(in) :: T, rho, lambda
    Real * 8, Dimension(:), intent(in) :: g,r
    Real * 8 :: P, glsm,glsp,gsp,rp,rd
    integer :: i1
    rp=0.d0
    do i1=1,size(r)
      rd=r(i1)
      if ( rp <= 1.d0 .and. rd > 1.d0) gsp=g(i1)
      if ( rp < lambda .and. rd >= lambda ) glsm = g(i1-1)
      if ( rp <= lambda .and. rd > lambda ) glsp = g(i1)
      rp=rd
    end do

    P = rho * T + (2.d0*pi/3.d0) * (rho**2 * T * gsp)
    P = P - (pi/3.d0) * rho **2 * lambda**3 * (glsp + glsm)
    return
  End Function

  Function phi_c_HSSW_SH() result(phi_c)
    Implicit None
    Real * 8 :: phi_lo, phi_up, phi_test, f_lo, f_up, f_test
    Real * 8 :: phi_c
    Real * 8 :: error
    Real * 8, parameter :: tol = 1.d-12
    Integer :: i1
    phi_lo = 5.d-2
    phi_up = 3.d-1
    f_lo = F_phic_aux(phi_lo)
    f_up = F_phic_aux(phi_up)
    error = (phi_up - phi_lo) / phi_up
    Do while (error > tol)
      phi_test = (phi_lo + phi_up) / 2.d0
      f_test = F_phic_aux(phi_test)
      If (f_test*f_up<0.d0) Then
        phi_lo = phi_test
        f_lo = f_test
      Else
        phi_up = phi_test
        f_up = f_test
      End If
      error = (phi_up - phi_lo) / phi_up
    End Do
    If (abs(f_up)>abs(f_lo)) Then
      phi_c = phi_lo
    Else
      phi_c = phi_up
    End If
  End Function

  Function T_c_HSSW_SH(lambda) result(T_c)
    Implicit None
    Real * 8, intent(in) :: lambda
    Real * 8 :: T_c
    Real * 8 :: phi_c
    Real * 8 :: FT, phi_c_m5, lambda3m
    phi_c = phi_c_HSSW_SH()
    phi_c_m5 = (1.d0 - phi_c) ** 5
    lambda3m = (lambda ** 3) -1.d0
    FT = (2.d0 + 5.d0 * phi_c - phi_c ** 2) / ( 2.d0 * phi_c_m5 * lambda3m )
    T_c = 1.d0 / FT
    Return
  End Function

  Function T_c_HSSW_GSH(lambda) result(T_c)
    Implicit None
    Real * 8, intent(in) :: lambda
    Real * 8 :: T_c
    Real * 8 :: phi_c
    Real * 8 :: FT, phi_c_m5, lambda3m
    phi_c = phi_c_HSSW_SH()
    phi_c_m5 = (1.d0 - phi_c) ** 5
    lambda3m = (lambda ** 3) - 1.d0
    FT = (2.d0 + 5.d0 * phi_c - phi_c ** 2) / ( 2.d0 * phi_c_m5 * lambda3m )
    T_c = 1.d0 / log(FT + 1.d0)
    Return
  End Function

  Function F_phic_aux(phi) Result(F)
    Implicit None
    Real * 8, intent(in) :: phi
    Real * 8 :: F
    Real * 8 :: phi5, phi4, phi3, phi2
    phi2 = phi  * phi
    phi3 = phi2 * phi
    phi4 = phi3 * phi
    phi5 = phi4 * phi
    F = -phi5 + 5.d0*phi4 -4.d0*phi3 -20.d0*phi2 -5.d0*phi + 1.d0
    Return
  End Function

  Subroutine Binodal_phis_HSSW_GSH(T,lambda,phi_l,phi_u)
    Implicit None
    Real * 8, intent(in) :: T, lambda
    Real * 8, intent(out) :: phi_l, phi_u
    Real * 8 :: phi_m, error, phi_s_l, phi_s_u
    Real * 8 :: phi_c, T_c,PB1,PB2, P_test, P_up, P_lo
    Real * 8 :: FT, A1, A2, Pdum
    Integer, parameter :: qnp = 2**10
    Real * 8, dimension(qnp) :: quad_p, quad_w
    Real * 8, dimension(qnp) :: phi_q1, phi_w1, integrand1, phi_q2, phi_w2, integrand2
    Integer :: i1
    Real * 8, parameter :: tol = 1.d-10
    call Clenshaw_Curtis_quadrature(qnp,quad_p,quad_w)
    Pdum = pi / (6.d0 * T)
    phi_c = phi_c_HSSW_SH()
    T_c = T_c_HSSW_GSH(lambda)
    FT= exp(1.d0/T) - 1.d0
    call roots_phi_s_P_aux(FT,lambda,phi_c,phi_s_l,phi_s_u)
    PB1 = P_HSSW_GSH(T,phi_s_l,lambda) * Pdum
    PB2 = P_HSSW_GSH(T,phi_s_u,lambda) * Pdum
    If (PB1 > PB2) Then
      P_up = PB1
      P_lo = PB2
      If (P_lo< 0.d0) P_lo = 0.d0
    Else
      P_up = PB2
      P_lo = PB1
      If (P_lo<0.d0) P_lo = 0.d0
    End If
    P_test = (P_up + P_lo) / 2.d0
    phi_l = Pressure_root_HSSW(FT,lambda,P_test,0.d0,phi_s_l)
    phi_m = Pressure_root_HSSW(FT,lambda,P_test,phi_s_l,phi_s_u)
    phi_u = Pressure_root_HSSW(FT,lambda,P_test,phi_s_u,1.d0)
    phi_q1 = quad_p
    phi_q2 = quad_p
    phi_w1 = quad_w
    phi_w2 = quad_w
    call rescale_OW(qnp,phi_q1,phi_w1,-1.d0,1.d0,phi_l,phi_m)
    call rescale_OW(qnp,phi_q2,phi_w2,-1.d0,1.d0,phi_m,phi_u)
    !Open (unit=10, file="integrand_test.dat", status="Replace")
    Do i1 = 1, qnp
      integrand1(i1) =  (P_HSSW_GSH(T,phi_q1(i1),lambda) * Pdum - P_test)/(phi_q1(i1)**2)
      integrand2(i1) =  (P_HSSW_GSH(T,phi_q2(i1),lambda) * Pdum - P_test)/(phi_q2(i1)**2)
      !Write(10,*) integrand1(i1),phi_q1(i1),integrand2(i1),phi_q2(i1)
    End Do
    !Close (10)
    !stop
    A1 = abs(sum(integrand1 * phi_w1))
    A2 = abs(sum(integrand2 * phi_w2))
    error = Abs(A1 - A2) / A1
    Do while (error > tol)
      If (A1 > A2) Then
        P_lo  = P_test
      Else
        P_up = P_test
      End If
      P_test = (P_up + P_lo) / 2.d0
      phi_l = Pressure_root_HSSW(FT,lambda,P_test,0.d0,phi_s_l)
      phi_m = Pressure_root_HSSW(FT,lambda,P_test,phi_s_l,phi_s_u)
      phi_u = Pressure_root_HSSW(FT,lambda,P_test,phi_s_u,1.d0)
      phi_q1 = quad_p
      phi_q2 = quad_p
      phi_w1 = quad_w
      phi_w2 = quad_w
      call rescale_OW(qnp,phi_q1,phi_w1,-1.d0,1.d0,phi_l,phi_m)
      call rescale_OW(qnp,phi_q2,phi_w2,-1.d0,1.d0,phi_m,phi_u)
      Do i1 = 1, qnp
        integrand1(i1) =  (P_HSSW_GSH(T,phi_q1(i1),lambda) * Pdum - P_test)/(phi_q1(i1)**2)
        integrand2(i1) =  (P_HSSW_GSH(T,phi_q2(i1),lambda) * Pdum - P_test)/(phi_q2(i1)**2)
      End Do
      A1 = abs(sum(integrand1 * phi_w1))
      A2 = abs(sum(integrand2 * phi_w2))
      error = Abs(A1 - A2) / A1
    End Do
    !Open (unit=10, file="integrand_test.dat", status="Replace")
    !Do i1 = 1, qnp
    !  Write(10,*) integrand1(i1),phi_q1(i1),integrand2(i1),phi_q2(i1)
    !End Do
    !Close (10)
    Print*, phi_l, phi_u, T
    Return
  End Subroutine

  Function Pressure_root_HSSW(FT,lambda,Pc,phi_Bl,phi_Bu) result(phi_root)
    Implicit None
    Real * 8, intent(in) :: FT, lambda, Pc
    Real * 8, intent(in) ::  phi_Bl, phi_Bu
    Real * 8 :: phi_lo, phi_up
    Real * 8 :: phi_root
    Real * 8 :: error, P_lo, P_up, phi_test, P_test
    Real * 8, parameter :: tol = 1.d-12
    phi_lo = phi_Bl
    phi_up = phi_Bu
    error = (phi_up - phi_lo) / phi_up
    P_lo = P_root_HSSW_aux(FT,lambda,Pc,phi_lo)
    P_up = P_root_HSSW_aux(FT,lambda,Pc,phi_up)
    Do while(error > tol)
      phi_test = (phi_lo + phi_up) / 2.d0
      P_test = P_root_HSSW_aux(FT,lambda,Pc,phi_test)
      If (P_test * P_lo < 0.d0) Then
        phi_up = phi_test
        P_up = P_test
      Else
        phi_lo = phi_test
        P_lo = P_test
      End If
      error = (phi_up - phi_lo) / phi_up
    End Do
    If (abs(P_lo)<abs(P_up)) Then
      phi_root = phi_lo
    Else
      phi_root = phi_up
    End If
  End Function

  Function P_root_HSSW_aux(FT,lambda,Pc,phi) result(P)
    Implicit None
    Real * 8, intent(in) :: FT, lambda, Pc,phi
    Real * 8 :: P
    Real * 8 :: phi2,phi3,phim3
    phi2 = phi * phi
    phi3 = phi2 * phi
    phim3 = (1.d0 - phi) ** 3
    P = (phi * (1.d0 + phi + phi2 - phi3) / phim3) - (4.d0 * (lambda ** 3 - 1.d0)* FT * phi2) - Pc
  End Function

  Subroutine roots_phi_s_P_aux(FT,lambda,phi_c,phi_l,phi_u)
    Implicit None
    Real * 8, intent(in) :: FT, lambda, phi_c
    Real * 8 :: phi_l,phi_u
    Real * 8 :: error
    Real * 8, parameter :: tol = 1.d-12
    Real * 8 :: phi_test, phi_lo, phi_up
    Real * 8 :: root_lo, root_up, root_test
    phi_lo = 0.d0
    phi_up = phi_c
    error = (phi_up - phi_lo) / phi_up
    root_lo = roots_phi_s_aux(FT,lambda,phi_lo)
    root_up = roots_phi_s_aux(FT,lambda,phi_up)
    Do while (error > tol)
      phi_test  = (phi_up + phi_lo) / 2.d0
      root_test = roots_phi_s_aux(FT,lambda,phi_test)
      If (root_test * root_lo < 0.d0) Then
        phi_up  = phi_test
        root_up = root_test
      Else
        phi_lo = phi_test
        root_lo = root_test
      End If
      error = (phi_up - phi_lo) / phi_up
    End Do
    If (abs(root_lo)<abs(root_up)) Then
      phi_l = phi_lo
    Else
      phi_l = phi_up
    End If
    phi_lo = phi_c
    phi_up = 1.d0
    error  = (phi_up - phi_lo) / phi_up
    root_lo = roots_phi_s_aux(FT,lambda,phi_lo)
    root_up = roots_phi_s_aux(FT,lambda,phi_up)
    Do while (error > tol)
      phi_test  = (phi_up + phi_lo) / 2.d0
      root_test = roots_phi_s_aux(FT,lambda,phi_test)
      If (root_test * root_lo < 0.d0) Then
        phi_up  = phi_test
        root_up = root_test
      Else
        phi_lo  = phi_test
        root_lo = root_test
      End If
      error = (phi_up - phi_lo) / phi_up
    End Do
    If (abs(root_lo)<abs(root_up)) Then
      phi_u = phi_lo
    Else
      phi_u = phi_up
    End If
    return
  End Subroutine

  Function roots_phi_s_aux(FT,lambda,phi) result(F)
    Implicit None
    Real * 8, intent(in) :: FT,lambda,phi
    Real * 8 :: F
    Real * 8 :: A
    Real * 8 :: phim4, phi2
    phi2 = phi * phi
    phim4 = (1.d0 - phi) ** 4
    A = -8.d0 * FT * (lambda ** 3 - 1.d0)
    F = A * phi * phim4 + 8.d0 * phi - 2.d0 * phi2 + phim4
  End Function

  !This subroutine writes down the radial distribution function, the coordination number
	! and the associated run file info
	Subroutine g_rel_writting_hssw(r,k,gr,coor_n_r,drp,drw,sys_cond,approx,quad_k,&
            & quad_dr,quadp_f,quadp_i,quadp_drf,quadp_dri,folder,id)
		Real * 8, dimension(:), intent(in) :: r, k, gr, coor_n_r,drp, drw, sys_cond,quadp_f,quadp_drf
		Integer, dimension(:) :: quadp_i, quadp_dri
		Character(*), intent(in) :: folder, approx, quad_k, quad_dr
		Integer, intent(in) :: id
		Real * 8 :: phi, T, lambda
		Integer :: di_phi, i_T, di_T
		Integer, parameter :: di_phi_min=4, i_T_max = 3, di_T_min = 4, id_max=2
		Character (len=di_phi_min) :: phi_char
		Character (len=i_T_max+di_T_min+1) :: T_char
		Character (len=id_max) :: id_char
		Character (len=len(folder)+len(phi_char)+len(T_char)+len(id_char)+19) :: File_name
		Character (len=len(File_name)) :: info_File_name
		Integer :: i1, i2
		!System conditions
		phi = sys_cond(1)
		T = sys_cond(2)
		lambda = sys_cond(3)
		!Calculating the associated integers for phi and T:
		i_T = int(T)
		di_T = int( (T- dble(i_T)) * (10.d0 ** di_T_min) )
		di_phi = int(phi * (10.d0 ** di_phi_min) )
		!Writting con character variables, formats have to be consistent with char lengths
		Write(phi_char,"(I4.4)") di_phi
		Write(T_char,"(I3.3,A1,I4.4)") i_T,"_",di_T
		Write(id_char,"(I2.2)") id
		File_name = folder//"g_phi_0_"//phi_char//"_T_"//T_char//"_id_"//id_char//".dat"
		info_File_name = folder//"g_phi_0_"//phi_char//"_T_"//T_char//"_id_"//id_char//".inf"
		Open(unit=100, file=File_name, status="Replace")
		Write (100,*) "# 			r							g(r)							#_coor(r)"
    Do i1=1, size(r)
      Write(100,*) r(i1),gr(i1),coor_n_r(i1)
		End Do
		Close(100)
		Open (unit=100, file=info_File_name, status="replace")
		Write(100,*) "###  Run Information for the file: "//File_name
		Write(100,*) "### System parameters ###"
		Write(100,*) "Volume Fraction phi = ", phi
		Write(100,*) "Temperature      T* = ", T
		Write(100,*) "Lambda              = ", lambda
		Write(100,*) "###############################"
		Write(100,*) "### Structure Approximation ###"
		Select case (approx)
			case ("SHVW")
				Write(100,*) "				Sharma-Sharma + Verlet Weiss (HSPY)"
			case ("GSVW")
				Write(100,*) "				Generalized Sharma-Sharma + Verlet Weiss (HSPY)"
			case Default
				Write(100,*) "				"//approx
		End Select
		Write(100,*) "###############################"
		Write(100,*) "### Structure Quadratures ###"
		Select case (quad_k)
			case("CLCU")
				Write(100,*) "				Clenshaw-Curtis"
				Write(100,*) "quadrature number of points = ", quadp_i(1)
				Write(100,*) "quadrature minimum value    = ", quadp_f(1)
				Write(100,*) "quadrature maximum value    = ", quadp_f(2)
			case Default
				Write(100,*) "				"//quad_k
				Write(100,*) "quadrature number of points = ", quadp_i(1)
				Write(100,*) "quadrature minimum value    = ", quadp_f(1)
				Write(100,*) "quadrature maximum value    = ", quadp_f(2)
		End Select
		Write(100,*) "###############################"
		Write(100,*) "### Wave-Vectors ###"
		Write(100,*) "number of points     = ", size(k)
		Write(100,*) "k minimum value      = ", k(1)
		Write(100,*) "k maximum value      = ", k(size(k))
    Write(100,*) "###############################"
		Write(100,*) "### Radial Distance ###"
		Write(100,*) "number of points     = ", size(k)
		Write(100,*) "r minimum value      = ", r(1)
		Write(100,*) "r maximum value      = ", r(size(r))
    Write(100,*) "###############################"
		Write(100,*) "### Differential Radial Distance Quadratures ###"
		Select case (quad_dr)
			case("CLCU")
				Write(100,*) "				Clenshaw-Curtis"
				Write(100,*) "quadrature number of points = ", quadp_dri(1)
				Write(100,*) "quadrature minimum value    = ", quadp_drf(1)
				Write(100,*) "quadrature maximum value    = ", quadp_drf(2)
			case Default
				Write(100,*) "				"//quad_dr
				Write(100,*) "quadrature number of points = ", quadp_i(1)
				Write(100,*) "quadrature minimum value    = ", quadp_f(1)
				Write(100,*) "quadrature maximum value    = ", quadp_f(2)
		End Select
    Write(100,*)"# 						r										Weight"
    Do i1 =1 , size(drp)
      Write(100,*) drp(i1), drw(i1)
    End Do
		Close(100)
	End Subroutine

  Function P0(k,a) Result(P)
    Implicit None
    Real * 8, intent(in) :: k,a
    Real * 8 :: P
    Real * 8 :: ka, k2, ka2, ka4
    ka = k*a
    k2 = k**2
    ka2 = ka * ka
    ka4 = ka2 * ka2
    If (ka < 0.075 ) Then
      If (k > 0.d0) Then
        P = ( -1.d0 / k2) + a*a* (0.5d0 - (ka2/24.d0))
      Else
        P = - 1.d20
      End If
    Else
      P = - cos(ka) / k2
    End If
    Return
  End function

  Function P1(k,a) Result(P)
    Implicit None
    Real * 8, intent(in) :: k,a
    Real * 8 :: P
    Real * 8 :: ka,k3, ka2,a2,a3
    ka = k*a
    k3 = k**3
    a2 = a*a
    ka2 = ka*ka
    If (ka < 0.075 ) Then
        a3=a2*a
        P = ( (10.d0 - ka2)/(30.d0) )
        P = P * a3
    Else
      P = ( sin(ka) - ka * cos(ka) ) / k3
    End If
    Return
  End function

  Function P2(k,a) Result(P)
    Implicit None
    Real * 8, intent(in) :: k,a
    Real * 8 :: P
    Real * 8 :: ka,k2,k4, ka2, a4
    ka = k*a
    k2 = k**2
    ka2 = ka * ka
    k4 = k2 * k2
    a4 = a**4
    If (ka < 0.075 ) Then
      P = (2.d0/k4) + a4*(0.25d0 - (1.d0 / 36.d0) * ka2)
    Else
      P = ( 2.d0*ka*sin(ka) +( 2.d0-ka2 )*cos(ka) ) / k4
    End If
    Return
  End function

  Function P4(k,a) Result(P)
    Implicit None
    Real * 8, intent(in) :: k,a
    Real * 8 :: P
    Real * 8 :: ka,ka2,ka4,k2,k4,k6,a6,a4
    ka = k*a
    ka2 = ka * ka
    ka4 = ka2 * ka2
    k2 = k*k
    k4 = k2 * k2
    k6 = k4 * k2
    ka2 = ka * ka
    a6 = a**6
    a4 = a**4
    If (ka < 0.075 ) Then
      P = (-24.d0/k6) + a6 * ( (1.d0/6.d0) - (1.d0/48.d0) * ka2 )
    Else
      P = ( (4.d0*ka*(ka2-6.d0))* sin(ka) - (ka4 - 12.d0*ka2 + 24.d0) * cos(ka) ) / k6
    End If
    Return
  End function
  !!!!!!!!!!!!!!!!!!Delta Ps functions

  Function DP0(k,a,b) Result(DP)
    Implicit None
    Real * 8, intent(in) :: k,a,b
    Real * 8 :: DP
    Real * 8 :: ka, k2, ka2, kb, kb2,a2,b2
    ka = k*a
    kb = k*b
    If (ka < 0.075d0 .And. kb < 0.075d0) Then
      a2 = a * a
      b2 = b * b
      ka2 = a2 * ka * ka
      kb2 = b2 * kb * kb
      DP = 0.5d0 * (a2-b2) - (ka2-kb2)/24.d0
    Else
      DP = P0(k,a) - P0(k,b)
    End If
    Return
  End function

  Function DP1(k,a,b) Result(DP)
    Implicit None
    Real * 8, intent(in) :: k,a,b
    Real * 8 :: DP
    DP = P1(k,a)-P1(k,b)
    Return
  End function

  Function DP2(k,a,b) Result(DP)
    Implicit None
    Real * 8, intent(in) :: k,a,b
    Real * 8 :: DP
    Real * 8 :: ka2, kb2, ka, kb, a4,b4,a2,b2, k2
    ka = k*a
    kb = k*b
    If (ka < 0.075 .And. kb < 0.075d0 ) Then
      a2=a*a
      b2=b*b
      a4 = a2 * a2
      b4 = b2 * b2
      k2 = k * k
      ka2 = a2*k2*a4
      kb2 = b2*k2*b4
      DP = (a4-b4)*0.25d0 - (ka2-kb2)*(1.d0/36.d0)
    Else
      DP = P2(k,a)-P2(k,b)
    End If
    Return
  End function

  Function DP4(k,a,b) Result(DP)
    Implicit None
    Real * 8, intent(in) :: k,a,b
    Real * 8 :: DP
    Real * 8 :: ka, kb
    Real * 8 :: a6, b6, ka2, kb2,a2,b2,k2
    ka = k*a
    kb = k*b
    If (ka < 0.075d0 .And. kb < 0.075d0 ) Then
      k2 = k * k
      a2 = a*a
      b2 = b*b
      a6 = a**6
      b6 = b**6
      ka2 = a2 * k2 * a6
      kb2 = b2 * k2 * b6
      DP = (a6-b6)* (1.d0/6.d0) - (ka2-kb2) * (1.d0/48.d0)
    Else
      DP = P4(k,a) - P4(k,b)
    End If
    Return
  End function

  Function c1_hssw(T,lambda,k) Result(c1k)
    Implicit None
    Real * 8, intent(in) :: T,lambda,k
    Real * 8 :: c1k
    Real * 8 :: x
    Real * 8 :: A,B,C
    Real * 8 :: A1, A2, A3
    Real * 8 :: B1, B2, B3, B4, B5
    Real * 8 :: C1, C2, C3, C4, C5
    Real * 8 :: LM1
    LM1 = lambda - 1
    x = exp(1.d0/T) - 1.d0
    If ( 1.d0 < lambda .and. lambda < 2.d0 ) Then
      A1 = 4.d0/3.d0
      A2 = -1.d0
      A3 = 1.d0/12.d0
      B1 = -2.d0
      B2 = 1.d0/6.d0
      B3 = 0.5d0 * (lambda**2 - 1.d0)
      B4 = - 4.d0 * ( lambda**3 - 1.d0 ) / 3.d0
      B5 = lambda **2 - 1
      C1 = B4
      C2 = -B5
      C3 = B2
      C4 = 2.d0 * B3 * B3
      C5 = C4
      A = - ( A1 * DP1(k,1.d0,0.d0) + A2 * DP2(k,1.d0,0.d0) + A3 * DP4(k,1.d0, 0.d0) )
      A = A + x * ( A1 * DP1(k,lambda,1.d0) + A2 * DP2(k,lambda,1.d0) + A3 * DP4(k,lambda,1.d0) )
      B = - ( B1 * DP2(k,LM1,0.d0) + B2 * DP4(k,LM1,0.d0) )
      B = B - ( B3 * DP0(k,1.d0,LM1) + B4 * DP1(k,1.d0,LM1) + B5 * DP2(k,1.d0,LM1) )
      B = B + x * ( B3 * DP0(k,1.d0,lambda) + B4 * DP1(k,1.d0,lambda) + B5 * DP2(k,1.d0,lambda) )
      C = - ( C1 * DP1(k,LM1,0.d0) + C2 * DP2(k,LM1,0.d0) + C3 * DP4(k,LM1,0.d0) )
      C = C - ( C4 * DP0(k,1.d0,LM1) )
      C = C + x * ( C4 * DP0(k,lambda,1.d0) )
    Else If ( 2.d0 .le. lambda .and. lambda < 3.d0) Then
    Else If ( 3.d0 .le. lambda ) Then
    Else
    End If
    c1k = A + (x*B) + (x*x*C)
    c1k = c1k * ( (2.d0 * pi) ** 2 )
    Return
  End Function

  Function dc1_hssw(T,lambda,k) Result(c1k)
    Implicit None
    Real * 8, intent(in) :: T,lambda,k
    Real * 8 :: c1k
    Real * 8 :: x,x2
    Real * 8 :: A,B,C
    Real * 8 :: A1, A2, A3
    Real * 8 :: B1, B2, B3, B4, B5
    Real * 8 :: C1, C2, C3, C4, C5
    Real * 8 :: LM1
    LM1 = lambda - 1
    x = exp(1.d0/T) - 1.d0
    x2 = x * x
    If ( 1.d0 < lambda .and. lambda < 2.d0 ) Then
      A1 = 4.d0/3.d0
      A2 = -1.d0
      A3 = 1.d0/12.d0
      B1 = -2.d0
      B2 = 1.d0/6.d0
      B3 = 0.5d0 * (lambda**2 - 1.d0)
      B4 = - 4.d0 * ( lambda**3 - 1.d0 ) / 3.d0
      B5 = lambda **2 - 1
      C1 = B4
      C2 = -B5
      C3 = B2
      C4 = 2.d0 * B3 * B3
      C5 = C4
      A =  x * ( A1 * DP1(k,lambda,1.d0) + A2 * DP2(k,lambda,1.d0) + A3 * DP4(k,lambda,1.d0) )
      B = - ( B1 * DP2(k,LM1,0.d0) + B2 * DP4(k,LM1,0.d0) )
      B = B - ( B3 * DP0(k,1.d0,LM1) + B4 * DP1(k,1.d0,LM1) + B5 * DP2(k,1.d0,LM1) )
      B = B + x * ( B3 * DP0(k,1.d0,lambda) + B4 * DP1(k,1.d0,lambda) + B5 * DP2(k,1.d0,lambda) )
      C = - ( C1 * DP1(k,LM1,0.d0) + C2 * DP2(k,LM1,0.d0) + C3 * DP4(k,LM1,0.d0) )
      C = C - ( C4 * DP0(k,1.d0,LM1) )
      C = C + x * ( C4 * DP0(k,lambda,1.d0) )
    Else If ( 2.d0 .le. lambda .and. lambda < 3.d0) Then
    Else If ( 3.d0 .le. lambda ) Then
    Else
    End If
    c1k = A + (x*B) + (x2*C)
    c1k = c1k * ( (2.d0 * pi) ** 2 )
    Return
  End Function

  Function dc0_hssw(T,lambda,k) Result(c0k)
    Implicit None
    Real * 8, intent(in) :: T,lambda,k
    Real * 8 :: c0k
    Real * 8 :: x
    Real * 8 :: A,B,C
    Real * 8 :: A1, A2, A3
    Real * 8 :: B1, B2, B3, B4
    Real * 8 :: C1, C2, C3, C4
    Real * 8 :: LM1
    x = exp(1.d0/T) - 1.d0
    !c0k = x * ( 4.d0 * pi ) * DP1(k,lambda,1.d0)
    c0k = ( 4.d0 * pi ) * DP1(k,lambda,1.d0) / T
    Return
  End Function

  Function c1_hs_py(k) Result (c1k)
    Implicit None
    Real * 8, intent(in) :: k
    Real * 8 :: c1k
    Real * 8 :: A1, A2, A3
    A1 =  4.d0 / 3.d0
    A2 = -1.d0
    A3 = 1.d0/12.d0
    c1k = - ( A1 * DP1(k,1.d0,0.d0) + A2 * DP2(k,1.d0,0.d0) + A3 * DP4(k,1.d0, 0.d0) )
    c1k = c1k * ( (2.d0 * pi) ** 2 )
    Return
  End Function

  Function ck_hssw_py_1_order(phi,T,lambda,k) Result(ck)
    Implicit None
    Real * 8, intent(in) :: phi,T,lambda,k
    Real * 8 :: ck
    Real * 8 :: rho
    rho = 6.d0 * phi / pi
    ck = c_hs_py(phi,k) + dc0_hssw(T,lambda,k) + rho * dc1_hssw(T,lambda,k)
  End Function

  Function ck_hssw_pyvw_1_order(phi,T,lambda,k) Result(ck)
    Implicit None
    Real * 8, intent(in) :: phi,T,lambda,k
    Real * 8 :: ck
    Real * 8 :: rho
    rho = 6.d0 * phi / pi
    ck = c_hs_vw(phi,k) + dc0_hssw(T,lambda,k) + rho * dc1_hssw(T,lambda,k)
  End Function

  Function isk_hssw_py_1_order(phi,T,lambda,k) Result(isk)
    Implicit None
    Real * 8, intent(in) :: phi,T,lambda,k
    Real * 8 :: isk
    Real * 8 :: rho
    rho = 6.d0 * phi / pi
    isk = 1.d0 - (rho * ck_hssw_py_1_order(phi,T,lambda,k) )
  End Function

  Function isk_hssw_pyvw_1_order(phi,T,lambda,k) Result(isk)
    Implicit None
    Real * 8, intent(in) :: phi,T,lambda,k
    Real * 8 :: isk
    Real * 8 :: rho, rho_vw
    Real * 8 :: phi_vw
    Real * 8 :: chs, csw
    phi_vw = phi*(1.d0 - (phi / 16.d0))
    rho = 6.d0 * phi / pi
    rho_vw = 6.d0 * phi_vw / pi
    chs = c_hs_vw(phi,k)
    csw = dc0_hssw(T,lambda,k) !+ rho * dc1_hssw(T,lambda,k)
    isk = 1.d0 - ( (rho_vw * chs) + (rho * csw) )
  End Function

  Function sk_hssw_py_1_order(phi,T,lambda,k) Result(sk)
    Implicit None
    Real * 8, intent(in) :: phi,T,lambda,k
    Real * 8 :: sk
    Real * 8 :: rho
    sk = 1.d0 / isk_hssw_py_1_order(phi,T,lambda,k)
  End Function

  Function sk_hssw_pyvw_1_order(phi,T,lambda,k) Result(sk)
    Implicit None
    Real * 8, intent(in) :: phi,T,lambda,k
    Real * 8 :: sk
    Real * 8 :: rho
    sk = 1.d0 / isk_hssw_pyvw_1_order(phi,T,lambda,k)
  End Function

  Function Ts_hssw_pyvw_1_order(phi,lambda) Result(Ts)
    Implicit None
    Real * 8, intent(in) :: phi, lambda
    Real * 8 :: Ts
    Real * 8 :: T_u, T_l, T_test, isk
    Logical :: convergence
    Real * 8, Parameter :: tol=1.d-7
    Real * 8 :: error

    !Finding Temperature Bounds:
    T_test = 1.d0
    isk = isk_hssw_pyvw_1_order(phi,T_test,lambda,0.d0)
    convergence = .False.
    If (isk > 0.d0) Then
      T_u = T_test
      Do while (convergence .eqv. .False.)
        T_test = T_test * 0.5d0
        isk = isk_hssw_pyvw_1_order(phi,T_test,lambda,0.d0)
        If (isk < 0.d0) Then
          T_l = T_test
          convergence = .True.
        Else
          T_u = T_test
        End If
      End Do
    Else
      T_l = T_test
      Do while (convergence .eqv. .False.)
        T_test = T_test * 2.d0
        isk = isk_hssw_pyvw_1_order(phi,T_test,lambda,0.d0)
        If (isk < 0.d0) Then
          T_l = T_test
        Else
          T_u = T_test
          convergence = .True.
        End If
      End Do
    End If
    !!!!!!!!!!!!!!!!!!!!
    !Reducing bounds to a given tolerance
    convergence = .False.
    Do while (convergence .eqv. .False.)
      T_test = (T_u + T_l) * 0.5d0
      isk = isk_hssw_pyvw_1_order(phi,T_test,lambda,0.d0)
      If (isk < 0.d0) Then
        T_l = T_test
      Else
        T_u = T_test
      End If
      error = (T_u - T_l ) / T_u
      If (error < tol) convergence = .True.
    End Do
    Ts = T_test
  End Function

End Module
