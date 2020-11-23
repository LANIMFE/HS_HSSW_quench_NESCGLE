Module hsay_structure
  Use math_constants
  Use hs_structure
  Implicit None
  Contains

    !Function that calculates the direct correlation function of a Hard Sphere +
    !Attractive Yukawa System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    !for the hard core part
    ! c(k,ϕ,T,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function c_hsay_shvw(phi,T,z,k) result(c)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8 :: c
      Real * 8 :: k_func,ckdum,chs
      Real * 8 :: k2,z2,k4
      chs = c_hs_vw(phi,k)
      !chs = c_hs_py(phi,k)
      k2 = k * k
      z2 = z ** 2
      if (k > 0.075d0) then
        c = ( k * cos(k) ) + ( z * sin(k) )
        c = c / k
      else
        k4 = k2 * k2
        c = (1.d0 + z) - ((3.d0 + z) / 6.d0) * k2 + ((5.d0 + z) / 120.d0) * k4
      end if
      c = 4.d0 * pi * c / (T * (k2 + z2))
      c = chs + c
      return
    End Function

    !Function that calculates the inverse structure factor of a Hard Sphere +
    !Attractive Yukawa System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! 1/S(k,ϕ,T,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function is_hsay_shvw(phi,T,z,k) result(is)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent (in) :: k !Value of the wave vector
      Real * 8 :: is
      Real * 8 :: c
      Real * 8 :: chs
      Real * 8 :: k2,z2,k4, phi_vw
      phi_vw = phi*(1.d0 - (phi / 16.d0))
      chs = c_hs_vw(phi,k)
      !chs = c_hs_py(phi,k)
      k2 = k * k
      z2 = z ** 2
      if (k > 0.075d0) then
        c = ( k * cos(k) ) + ( z * sin(k) )
        c = c / k
      else
        k4 = k2 * k2
        c = (1.d0 + z) - ((3.d0 + z) / 6.d0) * k2 + ((5.d0 + z) / 120.d0) * k4
      end if
      c = 4.d0 * pi * c / (T * (k2 + z2))
      is =  (phi_vw * chs) + (phi * c)
      is = 1.d0 - 6.d0 * is / pi
      return
    End Function

    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function s_hsay_shvw(phi,T,z,k) result(s)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent (in) :: k !Value of the wave vector
      Real * 8 :: s
      s = 1.d0 / is_hsay_shvw(phi,T,z,k)
      return
    End Function

    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function c_hsay_gsvw(phi,T,z,k,qn,qw) result(c)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: c
      Real * 8 :: cay,u,chs,r,w, dumsin,kr
      Real * 8 :: k2,z2,k4
      Integer :: i1
      cay = 0.d0
      do i1=1, size(qn)
        r = qn(i1)
        kr = k * r
        if (kr > 0.075d0 ) then
          w = qw(i1) * r * sin(k * r) / k
        else
          w = qw(i1) * r**2 * (1.d0 + ( (k**2) * r/6.d0) )
        end if
        u = - exp(z * (r-1.d0) )/r
        cay = cay + ( ( exp(-u/T) - 1.d0 ) * w )
      end do
      chs = c_hs_vw(phi,k)
      c = chs + cay
      return
    End Function

    !Function that calculates the inverse structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! 1/S(k,ϕ,T,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function is_hsay_gsvw(phi,T,z,k,qn,qw) result(is)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: is
      is = c_hsay_gsvw(phi,T,z,k,qn,qw)
      is = 1.d0 - (6.d0 * phi * is / pi)
      return
    End Function

    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function s_hsay_gsvw(phi,T,z,k,qn,qw) result(s)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: s
      s = 1.d0 / is_hsay_gsvw(phi,T,z,k,qn,qw)
      return
    End Function

    Function Ts_hsay_shvw(phi,z) result(Ts)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8 :: Ts !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8 :: S_HS
      S_HS = s_hs_vw(phi,0.d0)
      Ts = 24.d0 * phi * (z + 1.d0)  * S_HS / (z*z)
      return
    End Function

End Module
