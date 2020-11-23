Module scgle
Use math_constants
Use omp_lib
Use structure_function_selector
Use quadratures
Implicit None

contains

	!Function that calculates the wave vector dependent part of the memory function
	!lambda(k,kc)
	!"k" -> wave vector
	!"kc" ->  Adjustable parameter, found that for HS systems kc ≈ 2 π * 1.305
	function lambda(k,kc) result(l)
		Implicit None
		Real * 8, intent(in) :: k, kc
		Real * 8 :: l
		l = 1.d0 + ((k/kc)*(k/kc))
		l = 1.d0 / l
		return
	end function
	!Function that calculates the wave vector dependent part of the memory function
	!lambda(k,kc)
	!"k" -> wave vector
	!"kc" ->  Adjustable parameter, found that for HS systems kc ≈ 2 π * 1.305
	function lambda_test(k,kc) result(l)
		Implicit None
		Real * 8, intent(in) :: k, kc
		Real * 8 :: l
		l = exp(-(k/kc)**10)
		return
	end function

	!Function that calculates the Intermediate scattering function
	!F(k,D0,t,S)
	!The function depends on:
	!"k" -> wave vector
	!"D0" -> Free particle diffussion coefficient
	!"t" -> correlation time
	!"s" -> Structure factor evaluated at "k"
	function int_scat_F(k,D0,t,s) result(f)
		real * 8, intent (in) :: k, D0, t, s
		real * 8 :: f
 		f = exp (- D0 * t * k**2 / s) * s
	end function

	!Function that calculates the self part of Intermediate scattering function
	!Fs(k,D0,t,S)
	!The function depends on:
	!"k" -> wave vector
	!"D0" -> Free particle diffussion coefficient
	!"t" -> correlation time
	function int_scat_Fs(k,D0,t) result(fs)
		real * 8, intent(in) :: k, D0, t
		real * 8 :: fs
		fs = exp (- D0 * t * k**2 )
	end function

	!Function that calculates the integral for the memory function for a time t given:
	! "D0" -> free particle diffussion coeffcient
	! "rho" -> system density
	! "k" -> wave vector array
	!	"kw" -> wave vector weights array
	! "S" -> structure factor array
	! "F" -> Intermediate Scattering function at a given time t
	! "Fs" -> Self part of "F" for the same given time t
	! "d" -> system dimensionality
	Function delta_z(D0,rho,k,kw,S,F,Fs,d) result(z)
		real * 8, intent(in), Dimension(:) :: k ! Wave vectors
		real * 8, intent(in), Dimension(:) :: kw ! Wave vectors weights
		real * 8, intent(in), Dimension(:) :: S ! Structure factor evaluated at k
		real * 8, intent(in), Dimension(:) :: F ! Intermediate Scattering Function
		real * 8, intent(in), Dimension(:) :: Fs !Self Part  Intermediate Scattering Function
		real * 8, intent(in) :: D0 !diffussion coefficient of free particle
		real * 8, intent(in) :: rho ! density of the system
		real * 8 :: HC_dim
		integer, intent(in) :: d !dimension of the system
		real * 8 :: z
		real * 8, dimension(size(k)) :: Integrand
		!HC_dim = (3.d0 - 0.5d0 * d )
		HC_dim = 1.d0
		!OMP PARALLEL WORKSHARE
		Integrand = (k ** (1+d)) * F * Fs * ((S - 1.d0) / S) ** 2
		z = D0 * HC_dim * sum(Integrand * kw) / ( 2 * d * rho * pi ** (d-1) )
		!OMP END PARALLEL WORKSHARE
		return
	End function

	!Subroutine that calculates the integral for the time dependent diffussion coefficient
	! "D0" -> free particle diffussion coeffcient
	! "it" -> integer for the wanted time array value
	! "dt" -> time differential
	! "Delz" -> time dependent memory array
	! "Difft" -> Diffussion coefficient array
	Subroutine DtS(it,D0,dt,Delz,Difft)
		integer, intent(in) :: it !time integer for which we want to calculate Dt
		real * 8, intent(in) :: D0 !diffussion coefficient of free particle
		real * 8, intent(in) :: dt !time differential
		real * 8, intent(in), dimension(:) :: Delz !time dependent memory function
		real * 8, dimension(:), intent(inout) :: Difft
		real * 8 :: DDz, Difftsum
		integer :: i1
		Difftsum = 0.d0
		do i1 = 1, it-1
			DDz = Delz(it-i1+1) - Delz(it-i1)
			Difftsum = Difftsum + (Difft(i1) * DDz)
		end do
		Difft(it) = -Difftsum + Difft(it-1)/dt
		Difft(it) = Difft(it) / ( (1.d0/dt) + Delz(1))
		return
	End Subroutine

	!Subroutine that calculates the mean squared displacement
	! "D0" -> free particle diffussion coeffcient
	! "it" -> integer for the wanted time array value
	! "dt" -> time differential
	! "Delz" -> time dependent memory array
	! "msdt" -> msdt array
	Subroutine MSDt_it(it,D0,dt,Delz,msdt)
		integer, intent(in) :: it !time integer for which we want to calculate Dt
		real * 8, intent(in) :: D0 !diffussion coefficient of free particle
		real * 8, intent(in) :: dt !time differential
		real * 8, intent(in), dimension(:) :: Delz !time dependent memory function
		real * 8, dimension(:), intent(inout) :: msdt
		real * 8 :: DDz, msdtsum, t
		integer :: i1
		!msdtsum = 0.d0
		!do i1 = 1, it-1
		!t = it*dt
		msdtsum = - (Delz(1) + (1.d0/dt)) * msdt(it-1)
		do i1 = 2, it-1
			!msdtsum = msdtsum + (msdt(i1) * (Delz(it-i1+1) - Delz(it-i1)) )
			msdtsum = msdtsum + (Delz(i1) * (msdt(it-i1+1) - msdt(it-i1)) )
		end do
		!msdt(it) = ((dt * D0) - (dt * msdtsum) + msdt(it-1) ) / ( 1.d0 + dt * Delz(1) )
		msdt(it) = dt * (1.d0 - msdtsum - (Delz(it) * msdt(1)) ) / (1.d0 + dt * Delz(1))
		!Print*, dt * msdtsum / (1.d0 + dt * Delz(1)), dt * (Delz(it) * msdt(1)) / (1.d0 +  dt *Delz(1))
		!If ((msdt(it) - msdt(it-1))/dt > 1.d0) stop
		!Print*, dt, dt * msdtsum, msdt(it-1)
		return
	End Subroutine

	!Subroutine that calculates the mean squared displacement
	! "D0" -> free particle diffussion coeffcient
	! "it" -> integer for the wanted time array value
	! "dt" -> time differential
	! "Delz" -> time dependent memory array
	! "msdt" -> msdt array
	Subroutine MSDt_it_old(it,D0,dt,Delz,msdt)
		integer, intent(in) :: it !time integer for which we want to calculate Dt
		real * 8, intent(in) :: D0 !diffussion coefficient of free particle
		real * 8, intent(in) :: dt !time differential
		real * 8, intent(in), dimension(:) :: Delz !time dependent memory function
		real * 8, dimension(:), intent(inout) :: msdt
		real * 8 :: DDz, msdtsum
		integer :: i1, n2
		msdtsum = 0.d0

		Do i1 = 1, n2

		End Do
		do i1 = 1, it-1
			msdtsum = msdtsum + (msdt(i1) * (Delz(it-i1+1) - Delz(it-i1)) )
		end do
		msdt(it) = ((dt * D0) - (dt * msdtsum) + msdt(it-1) ) / ( 1.d0 + dt * Delz(1) )
		!Print*, (msdt(it) - msdt(it-1))/dt
		If ((msdt(it) - msdt(it-1))/dt > 1.d0) stop
		!Print*, dt, dt * msdtsum, msdt(it-1)
		return
	End Subroutine

!Subroutine that calculates the dynamics of short times given:
! "d" -> system density
! "nst" -> number of small times in the wanted array results
! "rho" -> system density
! "D0" -> free diffussion coefficient
! "k" -> wave vector array
! "kw" -> wave vector weights array
! "t" -> times array
! "S" -> structure factor
! RESULTS IN CALCULATION OF the "nst * size(k)" array for:
! "F" -> Intermediate Scattering Function array
! "Fs" -> Self part of "F"
! AND "nst" first values of:
! "Delz" -> memory function array
! "Difft" -> time dependent diffussion coefficient
Subroutine short_times_dynamics(d,nst,rho,D0,k,kw,t,S,F,Fs,Delz,Difft,msd)
	Implicit None
	Integer, intent(in) ::  d !system dymension
	Integer, intent(in) :: nst !Number of small times that we are going to use
	Real * 8, intent(in) :: rho !density of the system
	Real * 8, intent(in) :: D0 !free particle diffussion coefficient
	Real * 8, intent(in), dimension(:) :: k, kw !Wave vector array and weights array
	Real * 8, intent(in), dimension(:) :: t !Array of times
	Real * 8, intent(in), dimension(:) :: S !Structure factor array
	Real * 8, intent(out), dimension(:,:) :: F,Fs !Arrays of Intermediate Scattering Function and its self part, as well as the writting variables
	Real * 8, intent(out), dimension(:) :: Delz !Array of the memory function
	Real * 8, intent(out), dimension(:) :: Difft !Array of the time dependent diffussion coefficient
	Real * 8, intent(out), dimension(:) :: msd !mean square displacement array
	Real * 8, dimension(size(k)) :: Fv, Fsv
	Real * 8 :: sumDelz,dt
	Integer :: i1,i2
	sumDelz = 0.d0
	dt = t(1)
	Do i2=1, nst
		Do i1=1, size(k)
			Fv(i1)  = int_scat_F(k(i1), D0, t(i2), s(i1))
			Fsv(i1) = int_scat_Fs(k(i1), D0 , t(i2))
		End Do
		F(:,i2)  = Fv
		Fs(:,i2) = Fsv
		Delz(i2) = delta_z(D0,rho,k,kw,S,Fv,Fsv,d)
		sumDelz =  sumDelz + Delz(i2)
		Difft(i2) = D0 * (1.d0 - dt*sumDelz)
		 !IntDt(i2,dt,Difft)
	End Do
	msd(1) = Difft(1) * dt
	Do i1=2, nst
		msd(i1) = msd(i1-1) + Difft(i1) * dt
	End do
End Subroutine

!
Subroutine short_times_dynamics_writting(nst,rho,D0,k,t,S,F,Fs)
	Implicit None
	Integer, intent(in) :: nst !Number of small times that we are going to use
	Real * 8, intent(in) :: rho !density of the system
	Real * 8, intent(in) :: D0 !free particle diffussion coefficient
	Real * 8, intent(in), dimension(:) :: k !Wave vector array
	Real * 8, intent(in), dimension(:) :: t !Array of times
	Real * 8, intent(in), dimension(:) :: S !Structure factor array
	Real * 8, intent(out), dimension(:,:) :: F,Fs !Arrays of Intermediate Scattering Function and its self part, as well as the writting variables
	Real * 8, dimension(size(k)) :: Fv, Fsv
	Real * 8 :: sumDelz,dt
	Integer :: i1,i2
	Do i2=1, nst
		Do i1=1, size(k)
			Fv(i1)  = int_scat_F(k(i1), D0, t(i2), s(i1))
			Fsv(i1) = int_scat_Fs(k(i1), D0 , t(i2))
		End Do
		F(:,i2)  = Fv
		Fs(:,i2) = Fsv
	End Do
End Subroutine

!Function that constructs the dummy vector α(k) used for the calculation of
!intermediate times, it needs:
! "D0" -> Free particle diffussion coefficient
! "dt" -> time differential
! "k" -> wave vector array
! "S" -> Structure factor array
! "lam" -> interpolating function λ(k) array
! "mem1" -> first value of the memory function "Delz(1)"
Function alphac(D0,dt,k,S,lam,mem1) result(ac)
	Implicit None
	Real * 8, intent(in) :: D0 !Free particle diffusion coefficient
	Real * 8, intent(in) :: dt !time differential
	Real * 8, intent(in), dimension(:) :: k,S,lam !Wave vector, Structure Factor, Lambda(k)
	Real * 8, intent(in) :: mem1 !First value of the memory function
	Real * 8, dimension(size(k)) :: ac !Dummy vector for discretization algorithm
	!OMP WORKSHARE
	ac = 1.d0 / S
	ac = (1.d0/dt) + (D0 * ((k ** 2) * ac)) + ( lam * mem1)
	ac = 1.d0 / ac
	!OMP END WORKSHARE
End Function
!Function that constructs the dummy vector αs(k) used for the calculation of
!intermediate times, it needs:
! "D0" -> Free particle diffussion coefficient
! "dt" -> time differential
! "k" -> wave vector array
! "lam" -> interpolating function λ(k) array
! "mem1" -> first value of the memory function "Delz(1)"
Function alphas(D0,dt,k,lam,mem1) result(as)
	Implicit None
	Real * 8, intent(in) :: D0 !Free particle diffusion coefficient
	Real * 8, intent(in) :: dt !time differential
	Real * 8, intent(in), dimension(:) :: k,lam !Wave vector, Structure Factor, Lambda(k)
	Real * 8, intent(in) :: mem1 !First value of the memory function
	Real * 8, dimension(size(k)) :: as	!Dummy vector for discretization algorithm
	!OMP WORKSHARE
	as = (1.d0/dt) + (D0 * k ** 2) + ( lam * mem1)
	as = 1.d0 / as
	!OMP END WORKSHARE
End Function

!Subroutine that calculates the times values of "F", "Fs" and "Delz" from "it"
!up to size(t). The subroutine needs:
! "d" -> system dimensionality
! "it" -> initial time for which we want to start to calculate the array values
! "rho" -> system density
! "D0" -> Free particle diffussion coefficient
! "dt" -> time differential
! "k" -> wave vector array
! "kw" -> wave vector weights array
! "lam" -> interpolating function λ(k) array
! "t" -> times array
! "S" -> Structure factor array
! THE FUNCTION RETURNS THE VALUES OF
! "F(it:end,:)" -> Intermediate Scattering Function
! "Fs(it:end,:)" -> Self part of "F"
! "Delz(it:end)" -> Memory function
! "Difft(it:end)" -> Diffussion Coefficient function
Subroutine medium_times_dynamics(d,it,rho,D0,dt,k,kw,lam,t,S,Fc,Fs,Delz,Difft,msd)
	Implicit None
	Integer, intent(in) :: d !Dimension of the system
	Integer, intent(in) :: it !Beggining of medium times where we are going to start
	Real * 8, intent(in) :: rho !density of the system
	Real * 8, intent(in) :: D0 !free particle diffussion coefficient
	Real * 8, intent(in) :: dt !time differential
	Real * 8, intent(in), dimension(:) :: k, kw !Wave vector array and weights array
	Real * 8, intent(in), dimension(:) :: t !Array of times
	Real * 8, intent(in), dimension(:) :: S, lam !Structure factor array, Lambda(k) array
	Real * 8, intent(out), dimension(:,:) :: Fc,Fs !Arrays of Intermediate Scattering Function and its self part
	Real * 8, intent(out), dimension(:) :: Delz !Array of the memory function
	Real * 8, intent(out), dimension(:) :: Difft !Array of the Diffussion coefficient
	Real * 8, intent(out), dimension(:) :: msd !mean square displacement array
	Real * 8, dimension(size(k)) :: Fcv, Fsv,dFc1,dFs1,Fcdum,Fsdum,Fs1,Fc1 !Views on a time for F and Fs
	Real * 8, dimension(size(k)) :: ac, as, lamc,lams!Dummies for discretization
	Real * 8, dimension(size(Delz)) :: DeltaDelz !Array of the memory function
	Real * 8 :: delztest, error, a_approx, b_approx
	Integer :: i1,i2
	Logical :: convergence
	Real * 8, parameter :: tol = 1.d-7
	ac = alphac(D0,dt,k,S,lam,Delz(1))
	as = alphas(D0,dt,k,lam,Delz(1))
	lamc = ac * lam
	lams = as * lam
	Fs1 = Fs(:,1)
	Fc1 = Fc(:,1)
	dFc1 = S - Fc1
	dFs1 = 1.d0 - Fs1
	Do i1=2, it-1
		deltadelz(i1) = Delz(i1)-Delz(i1-1)
	End Do
	Do i1=it, size(t)
		!NON ITERATIVE PART
		Fsdum = 0.d0
		Fcdum = 0.d0
		!$OMP PARALLEL DO DEFAULT(NONE) SHARED(Fc,Fs,i1,deltadelz), PRIVATE(i2) REDUCTION(+:Fsdum,Fcdum)
		Do i2=2, i1-1
			!OMP WORKSHARE
			Fsdum = Fsdum - (deltadelz(i1+1-i2) * Fs(:,i2))
			Fcdum = Fcdum - (deltadelz(i1+1-i2) * Fc(:,i2))
			!OMP END WORKSHARE
		End Do
		!$OMP END PARALLEL DO
		!Fcv = F(:,i1-1)
		!Fsv = Fs(:,i1-1)
		Fsdum = as * (lam * (Fsdum + Delz(i1-1) * Fs1) + (Fs(:,i1-1)/dt))
		Fcdum = ac * (lam * (Fcdum + Delz(i1-1) * Fc1) + (Fc(:,i1-1)/dt))
		!ITERATIVE PART
		convergence = .False.
		delztest = Delz(i1-1) !Can be interpolated for a better approximation
		!b_approx = (log(delz(i1-2)) - log(delz(i1-1))) / (dt)
		!a_approx = exp(log(delz(i1-1)) + ( b_approx * dt * (i1-1) ) )
		!delztest = a_approx * exp(-b_approx * dt * i1)
		Do while (convergence .eqv. .False.)
			Fcv =  Fcdum + lamc * delztest * dFc1
			Fsv =  Fsdum + lams * delztest * dFs1
			delz(i1) = delta_z(D0,rho,k,kw,S,Fcv,Fsv,d)
			error = 1.d0 - dabs(delztest / delz(i1))
			if (error < tol) convergence = .True.
			delztest = delz(i1)
		End Do
		Fc(:,i1)  = Fcv
		Fs(:,i1) = Fsv
		Call DtS(i1,D0,dt,Delz,Difft)
		Call MSDt_it(i1,D0,dt,Delz,msd)
		!Print*, (Delz(i1)-Delz(i1-1))/dt
		!msd(i1) = msd(i1-1) + Difft(i1) * dt
		deltadelz(i1) = Delz(i1)-Delz(i1-1)
	End Do
End Subroutine

!Subroutine that calculates the times values of "F", "Fs" and "Delz" from "it"
!up to size(t). The subroutine needs:
! "d" -> system dimensionality
! "it" -> initial time for which we want to start to calculate the array values
! "rho" -> system density
! "D0" -> Free particle diffussion coefficient
! "dt" -> time differential
! "k" -> wave vector array
! "kw" -> wave vector weights array
! "lam" -> interpolating function λ(k) array
! "t" -> times array
! "S" -> Structure factor array
! THE FUNCTION RETURNS THE VALUES OF
! "F(it:end,:)" -> Intermediate Scattering Function
! "Fs(it:end,:)" -> Self part of "F"
! "Delz(it:end)" -> Memory function
! "Difft(it:end)" -> Diffussion Coefficient function
Subroutine medium_times_dynamics_old(d,it,rho,D0,dt,k,kw,lam,t,S,Fc,Fs,Delz,Difft,msd)
	Implicit None
	Integer, intent(in) :: d !Dimension of the system
	Integer, intent(in) :: it !Beggining of medium times where we are going to start
	Real * 8, intent(in) :: rho !density of the system
	Real * 8, intent(in) :: D0 !free particle diffussion coefficient
	Real * 8, intent(in) :: dt !time differential
	Real * 8, intent(in), dimension(:) :: k, kw !Wave vector array and weights array
	Real * 8, intent(in), dimension(:) :: t !Array of times
	Real * 8, intent(in), dimension(:) :: S, lam !Structure factor array, Lambda(k) array
	Real * 8, intent(out), dimension(:,:) :: Fc,Fs !Arrays of Intermediate Scattering Function and its self part
	Real * 8, intent(out), dimension(:) :: Delz !Array of the memory function
	Real * 8, intent(out), dimension(:) :: Difft !Array of the Diffussion coefficient
	Real * 8, intent(out), dimension(:) :: msd !mean square displacement array
	Real * 8, dimension(size(k)) :: Fcv, Fsv,dFc1,dFs1,Fcdum,Fsdum,Fs1,Fc1 !Views on a time for F and Fs
	Real * 8, dimension(size(k)) :: ac, as, lamc,lams!Dummies for discretization
	Real * 8, dimension(size(Delz)) :: DeltaDelz !Array of the memory function
	Real * 8 :: delztest, error, a_approx, b_approx
	Integer :: i1,i2,n2
	Logical :: convergence
	Real * 8, parameter :: tol = 1.d-7
	ac = alphac(D0,dt,k,S,lam,Delz(1))
	as = alphas(D0,dt,k,lam,Delz(1))
	lamc = ac * lam
	lams = as * lam
	Fs1 = Fs(:,1)
	Fc1 = Fc(:,1)
	dFc1 = S - Fc1
	dFs1 = 1.d0 - Fs1
	Do i1=2, it-1
		deltadelz(i1) = Delz(i1)-Delz(i1-1)
	End Do
	n2 = it/2
	if (mod(n2,2)==1) n2= n2+1
	Do i1=it, size(t)
		!NON ITERATIVE PART
		Fsdum = 0.d0
		Fcdum = 0.d0
		!$OMP PARALLEL DO DEFAULT(NONE) SHARED(Fc,Fs,i1,deltadelz), PRIVATE(i2) REDUCTION(+:Fsdum,Fcdum)
		Do i2=2, i1-1
			!OMP WORKSHARE
			Fsdum = Fsdum - (deltadelz(i1+1-i2) * Fs(:,i2))
			Fcdum = Fcdum - (deltadelz(i1+1-i2) * Fc(:,i2))
			!OMP END WORKSHARE
		End Do
		!$OMP END PARALLEL DO
		!Fcv = F(:,i1-1)
		!Fsv = Fs(:,i1-1)
		Fsdum = as * (lam * (Fsdum + Delz(i1-1) * Fs1) + (Fs(:,i1-1)/dt))
		Fcdum = ac * (lam * (Fcdum + Delz(i1-1) * Fc1) + (Fc(:,i1-1)/dt))
		!ITERATIVE PART
		convergence = .False.
		delztest = Delz(i1-1) !Can be interpolated for a better approximation
		!b_approx = (log(delz(i1-2)) - log(delz(i1-1))) / (dt)
		!a_approx = exp(log(delz(i1-1)) + ( b_approx * dt * (i1-1) ) )
		!delztest = a_approx * exp(-b_approx * dt * i1)
		Do while (convergence .eqv. .False.)
			Fcv =  Fcdum + lamc * delztest * dFc1
			Fsv =  Fsdum + lams * delztest * dFs1
			delz(i1) = delta_z(D0,rho,k,kw,S,Fcv,Fsv,d)
			error = 1.d0 - dabs(delztest / delz(i1))
			if (error < tol) convergence = .True.
			delztest = delz(i1)
		End Do
		Fc(:,i1)  = Fcv
		Fs(:,i1) = Fsv
		Call DtS(i1,D0,dt,Delz,Difft)
		Call MSDt_it(i1,D0,dt,Delz,msd)
		!Print*, (Delz(i1)-Delz(i1-1))/dt
		!msd(i1) = msd(i1-1) + Difft(i1) * dt
		deltadelz(i1) = Delz(i1)-Delz(i1-1)
	End Do
End Subroutine

Subroutine medium_times_dynamics_writting(d,it,rho,D0,dt,k,lam,t,S,F,Fs,Delz)
	Implicit None
	Integer, intent(in) :: d !Dimension of the system
	Integer, intent(in) :: it !Beggining of medium times where we are going to start
	Real * 8, intent(in) :: rho !density of the system
	Real * 8, intent(in) :: D0 !free particle diffussion coefficient
	Real * 8, intent(in) :: dt !time differential
	Real * 8, intent(in), dimension(:) :: k !Wave vector array and weights array
	Real * 8, intent(in), dimension(:) :: t !Array of times
	Real * 8, intent(in), dimension(:) :: S, lam !Structure factor array, Lambda(k) array
	Real * 8, intent(out), dimension(:,:) :: F,Fs !Arrays of Intermediate Scattering Function and its self part
	Real * 8, intent(in), dimension(:) :: Delz !Array of the memory function
	Real * 8, dimension(size(k)) :: Fcv, Fsv,dFc1,dFs1,Fcdum,Fsdum,Fs1,Fc1 !Views on a time for F and Fs
	Real * 8, dimension(size(k)) :: ac, as, lamc,lams!Dummies for discretization
	Real * 8, dimension(size(Delz)) :: deltadelz
	Real * 8 :: delztest, error
	Integer :: i1,i2
	Logical :: convergence
	Real * 8, parameter :: tol = 1.d-8
	ac = alphac(D0,dt,k,S,lam,Delz(1))
	as = alphas(D0,dt,k,lam,Delz(1))
	lamc = ac * lam
	lams = as * lam
	Fs1 = Fs(:,1)
	Fc1 = F(:,1)
	dFc1 = S - Fc1
	dFs1 = 1.d0 - Fs1
	Do i1=2, size(Delz)
		deltadelz(i1) = Delz(i1)-Delz(i1-1)
	End Do
	Do i1=it, size(t)
		Fsdum = 0.d0
		Fcdum = 0.d0
		!OMP PARALLEL DO DEFAULT(NONE) SHARED(F,Fs,i1,deltadelz), PRIVATE(i2,Fcv,Fsv) REDUCTION(+:Fsdum,Fcdum)
		Do i2=2, i1-1
			Fcv = F(:,i2)
			Fsv = Fs(:,i2)
			Fsdum = Fsdum - (deltadelz(i1+1-i2) * Fsv)
			Fcdum = Fcdum - (deltadelz(i1+1-i2) * Fcv)
		End Do
		!OMP END PARALLEL DO
		Fsdum = as * (lam * (Fsdum + Delz(i1-1) * Fs1) + (Fsv/dt))
		Fcdum = ac * (lam * (Fcdum + Delz(i1-1) * Fc1) + (Fcv/dt ))
		Fcv =  Fcdum + lamc * Delz(i1) * dFc1
		Fsv =  Fsdum + lams * Delz(i1) * dFs1
		F(:,i1)  = Fcv
		Fs(:,i1) = Fsv
	End Do
End Subroutine

!Subroutine that saves half of the calculated times for the dynamic variables:
! "F" -> Intermediate scattering function
! "Fs" -> Self part of "F"
! "Delz" -> Memory function
! "Difft" -> Diffussion Coefficient array
! Along with the new time vector in terms of the time differential "dt"
! "t" -> time vector
! "dt" -> time differential
Subroutine half_save(dt,t,Delz,F,Fs,Difft,msd)
	Implicit None
	Real * 8, intent(in) :: dt
	Real * 8, dimension(:), intent(inout) :: t, Delz, Difft, msd
	Real * 8, dimension(:,:), intent(inout) :: F, Fs
	integer :: i1, sizet
	sizet = size(t)
	Do i1=1, sizet
		t(i1) = i1 * dt
	End Do
	!Do i1=1, sizet/2
	Do i1=1, sizet/2
	F(:,i1)  = (F(:,2*i1) + F(:,2*i1-1)) * 0.5d0
	Fs(:,i1) = (Fs(:,2*i1)+ Fs(:,2*i1-1)) * 0.5d0
	Delz(i1) = (Delz(2*i1)+ Delz(2*i1-1)) * 0.5d0
	Difft(i1) = (Difft(2*i1)+ Difft(2*i1-1)) * 0.5d0
	msd(i1) = (msd(2*i1)+ msd(2*i1-1)) * 0.5d0
	End DO
	!Do i1=sizet/4+1, sizet/2
	!F(:,i1)  = (F(:,2*i1) + 4.d0*F(:,2*i1-1)+ F(:,2*i1-2)) / 6.d0
	!Fs(:,i1) = (Fs(:,2*i1)+ 4.d0*Fs(:,2*i1-1)+ Fs(:,2*i1-2)) / 6.d0
	!Delz(i1) = (Delz(2*i1)+ 4.d0*Delz(2*i1-1)+Delz(2*i1-2)) / 6.d0
	!Difft(i1) = (Difft(2*i1)+ 4.d0*Difft(2*i1-1)+Difft(2*i1-2)) / 6.d0
	!msd(i1) = (msd(2*i1)+ 4.d0*msd(2*i1-1)+ msd(2*i1-2)) / 6.d0
	!End do
	!Do i1 = 3 , sizet / 2
	!	F(:,i1)  = (F(:,2*i1) + 2.d0*F(:,2*i1-1) + 2.d0*F(:,2*i1-2) + F(:,2*i1-3)) / 6.d0
	!	Fs(:,i1) = (Fs(:,2*i1)+ 2.d0*Fs(:,2*i1-1)+ 2.d0*Fs(:,2*i1-2)+ Fs(:,2*i1-3)) / 6.d0
	!	Delz(i1) = (Delz(2*i1)+ 2.d0*Delz(2*i1-1)+ 2.d0*Delz(2*i1-2)+ Delz(2*i1-3)) / 6.d0
	!	Difft(i1) = (Difft(2*i1)+ 2.d0*Difft(2*i1-1)+ 2.d0*Difft(2*i1-2)+ Difft(2*i1-3)) / 6.d0
	!	msd(i1) = (msd(2*i1)+ 2.d0*msd(2*i1-1)+ 2.d0*msd(2*i1-2)+ msd(2*i1-3)) / 6.d0
	!End Do
	!i1=sizet/2
	!F(:,i1)  = F(:,2*i1)
	!Fs(:,i1) = Fs(:,2*i1)
	!Delz(i1) = Delz(2*i1)
	!Difft(i1) = Difft(2*i1)
	!msd(i1) = msd(2*i1)
	!Do i1 = sizet/4+1, sizet/2
	!	F(:,2*i1) = (F(:,2*i1) +4.d0 *F(:,2*i1-1) +F(:,2*i1-2))/6.d0
	!	Fs(:,2*i1) = (Fs(:,2*i1) +4.d0 *Fs(:,2*i1-1) +Fs(:,2*i1-2))/6.d0
	!	Delz(2*i1) = (Delz(2*i1) +4.d0 *Delz(2*i1-1) +Delz(2*i1-2))/6.d0
	!	msd(2*i1) = (msd(2*i1) +4.d0 *msd(2*i1-1) +msd(2*i1-2))/6.d0
	!End Do
	!Do i1 = 1 , sizet / 2
	!	F(:,i1)  = F(:,2*i1)
	!	Fs(:,i1) = Fs(:,2*i1)
	!	Delz(i1) = Delz(2*i1)
	!	Difft(i1) = Difft(2*i1)
	!	msd(i1) = msd(2*i1)
	!End Do
End Subroutine

!Subroutine that saves half of the calculated times for the dynamic variables:
! "F" -> Intermediate scattering function
! "Fs" -> Self part of "F"
! "Delz" -> Memory function
! "Difft" -> Diffussion Coefficient array
! Along with the new time vector in terms of the time differential "dt"
! "t" -> time vector
! "dt" -> time differential
Subroutine half_save_old(dt,t,Delz,F,Fs,Difft,msd)
	Implicit None
	Real * 8, intent(in) :: dt
	Real * 8, dimension(:), intent(inout) :: t, Delz, Difft, msd
	Real * 8, dimension(:,:), intent(inout) :: F, Fs
	integer :: i1, sizet
	sizet = size(t)
	Do i1=1, sizet
		t(i1) = i1 * dt
	End Do
	Do i1 = 1 , sizet / 2
		F(:,i1)  = F(:,2*i1)
		Fs(:,i1) = Fs(:,2*i1)
		Delz(i1) = Delz(2*i1)
		Difft(i1) = Difft(2*i1)
		msd(i1) = msd(2*i1)
	End Do
End Subroutine

!Subroutine that saves the wanted information for writting later on
!The writting variables needs to be previously allocated with the
!maximum size already known.
! "iks" -> the integer assigned for the wave-vector that we want to save
! "ist" -> the initial position for the variables to save into the writting variables
! "iwt" -> the initial position for the writting variables into which we are going to save
! "t" -> time vector array
! "Delz" -> memory vector array
! "F" -> Intermediate Scattering Function array
! "Fs" -> self part of "F"
! "tw" -> writting vector array for "t"
! "Delzw" -> writting vector array for "Delz"
! "Fw" -> writting vector array for "F"
! "Fsw" -> writting vector array for "Fsw"
Subroutine writting_save(iks,ist,iwt,t,Delz,F,Fs,tw,Delzw,Fw,Fsw)
	Implicit None
	Integer, intent(in) :: ist, iwt !initial position for the entering arrays and writting arrays
	Integer, intent(in) :: iks !Position of the wanted wave-vector follow up for "F" and "Fs"
	Real * 8, Dimension(:), intent(in) :: t, Delz !initial position for
	Real * 8, Dimension(:,:), intent(in) :: F,Fs !Intermediate scattering Function and its self part
	Real * 8, Dimension(:), intent(out) :: tw, DelzW, Fw, Fsw !Writting variables
End Subroutine

Function IntDelz(it,dt,Delz)
	Implicit None
	Integer :: it !initial time integer for the integration
	Real * 8, intent(in) :: dt !time differential
	Real * 8, dimension(:), intent(in) :: Delz !Memory function
	Real * 8 :: IntDelz
	Integer :: sizet
	sizet= size(Delz)
	IntDelz = sum(Delz(it:sizet)) * dt
End Function

Function IntDt(ft,dt,Difft)
	Implicit None
	Integer :: ft !final time integer for the integration
	Real * 8, intent(in) :: dt !time differential
	Real * 8, dimension(:), intent(in) :: Difft !Diffussion Coefficient
	Real * 8 :: IntDt
	Integer :: sizet
	sizet= size(Difft)
	IntDt = sum(Difft(1:ft)) * dt
End Function

!Subroutine that calculates the dynamics from small to long times
!The subroutine NEEDS:
! "d" -> system dimensionality
! "it" -> initial time for which we want to use small times dynamics subroutine
! "rho" -> system density
! "D0" -> Free particle diffussion coefficient
! "dt" -> time differential
! "indt" -> initial time differential
! "k" -> wave vector array
! "kw" -> wave vector weights array
! "lam" -> interpolating function λ(k) array
! "t" -> times array
! "S" -> Structure factor array
! THE FUNCTION RETURNS THE LONG TIME VALUES OF:
! "F" -> Intermediate Scattering function array
! "Fs" -> self part of "F"
! "Delz" -> memory function
! AND WRITES UPON something ALL THE TIMES VALUES OF:
! "F(iktest,:)"
! "Fs(iktest,:)"
Subroutine Long_t_dynamics(dyn_p,arrays,Dl,w_op,w_arrays_in,w_arrays_out)
	Implicit None
	Real * 8, dimension(7), intent(in) :: dyn_p
	Real * 8, dimension(:,:), intent(in) :: arrays
	Real * 8, dimension(:,:), intent(in) , optional :: w_arrays_in
	Real * 8, dimension(:,:), intent(out), optional :: w_arrays_out
	Real * 8, dimension(nint(dyn_p(3))) :: t !Array of times
	Real * 8, intent(out) :: Dl !long time diffussion coefficient, MSD
	Real * 8, dimension(:), Allocatable :: k, kw !Wave vector array and weights array
	Real * 8, dimension(:), Allocatable :: S, lam !Structure factor array, Lambda(k) array
	Real * 8, dimension(size(arrays,2),nint(dyn_p(3))) :: F,Fs !Arrays of Intermediate Scattering Function and its self part
	Real * 8, dimension(nint(dyn_p(3))) :: Delz !Array of the memory function
	Real * 8, dimension(nint(dyn_p(3))) :: Difft !Diffussion coefficient array
	Real * 8, dimension(nint(dyn_p(3))) :: msd !Mean square displacement array
	Integer :: d !Dimension of the system
	Integer :: it,nt !Short times points and medium times points
	Real * 8 :: rho !density of the system
	Real * 8 :: D0 !free particle diffussion coefficient
	Real * 8 :: indt !time differential
	Integer ::	decimations !Number of decimations used
	Integer :: i1, i2, sizek, sizekw, itw, alW, itm
	real * 8 :: time1,time2
	real * 8 :: dt
	!Writting Variables:
	Logical :: w_op !Select .True. for writting down the information
	Real * 8, dimension(:), allocatable :: k_w, S_w, lam_w !writting variables needed
	Real * 8, dimension(:,:), allocatable :: F_w, Fs_w
	Real * 8, dimension(:,:), allocatable :: FW, FsW
	Real * 8, dimension(:), allocatable :: DelzW, tW, DifftW, msdW
	!Data unpacking!!!!!
	d  = nint(dyn_p(1))
	it = nint(dyn_p(2))
	nt = nint(dyn_p(3))
	decimations = nint(dyn_p(4))
	rho = dyn_p(5)
	D0  = dyn_p(6)
	indt = dyn_p(7)
	sizek  = size(arrays,2)
	!k-Functions
	Allocate(k(sizek),kw(sizek),S(sizek),lam(sizek))
	k   = arrays(1,:)
	kw  = arrays(2,:)
	S   = arrays(3,:)
	lam = arrays(4,:)
	If (w_op .eqv. .True.) Then
		!k-Functions writting-arrays inputs
		sizekw  = size(w_arrays_in,2)	!Size of writting wave vectors array
		alW = (nt + ( decimations * nt/2 )) !Size of array needed for the writting variables
		!Internal memory allocation
		Allocate(k_w(sizekw),S_w(sizekw),lam_w(sizekw))
		Allocate( DelzW(alW), FW(sizekw,alW), FsW(sizekw,alW), tW(alW), DifftW(alw), msdW(alw))
		Allocate (F_w(sizekw,nt), Fs_w(sizekw,nt))
		k_w   = w_arrays_in(1,:)
		S_w   = w_arrays_in(2,:)
		lam_w = w_arrays_in(3,:)
	End If
	!!!!!!!!!!!!!!!!!
	do i1=1, nt
		t(i1) = i1 * indt
	end do
	dt = indt
	Print*, "Long times dynamics initialization"
	call short_times_dynamics(d,it,rho,D0,k,kw,t,S,F,Fs,Delz,Difft,msd)
	Print*, "Short Times done"
	call medium_times_dynamics(d,it + 1,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,Difft,msd)
	Print*, "Medium Times done"
	Dl = IntDelz(1,dt,Delz)
	If (w_op .eqv. .True.) Then
		call short_times_dynamics_writting(it,rho,D0,k_w,t,S_w,F_w,Fs_w)
		call medium_times_dynamics_writting(d,it + 1,rho,D0,dt,k_w,lam_w,t,S_w,F_w,Fs_w,Delz)
		itw = 0
		Do i1 = 1, nt
			itw = itw + 1
			FW(:,itw) = F_w(:,i1)
			FsW(:,itw) = Fs_w(:,i1)
			DelzW(itw) = Delz(i1)
			tW(itw) = t(i1)
			DifftW(itw) = Difft(i1)
			msdW(itw) = msd(i1)
		End Do
	End If
	itm = (nt / 2) + 1
	Do i1=1, decimations
		print*, "decimation", i1, " of ", decimations
		dt = 2.d0 * dt
		call half_save(dt,t,Delz,F,Fs,Difft,msd)
		If (w_op .eqv. .True.) Then
			Do i2 = 1 , nt / 2
				F_w(:,i2)  = F_w(:,2*i2)
				Fs_w(:,i2) = Fs_w(:,2*i2)
			End Do
		End If
		call medium_times_dynamics(d,itm,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,Difft,msd)
		If (w_op .eqv. .True.) Then
			call medium_times_dynamics_writting(d,itm,rho,D0,dt,k_w,lam_w,t,S_w,F_w,Fs_w,Delz)
		 	Do i2 = 1 + (nt/2), nt
				itw = itw + 1
				tW(itw) = t(i2)
				FW(:,itw)  = F_w(:,i2)
				FsW(:,itw) = Fs_w(:,i2)
				DelzW(itw) = Delz(i2)
				DifftW(itw) = Difft(i2)
				msdW(itw) = msd(i2)
			End Do
		End If
		Dl = Dl + IntDelz(itm,dt,Delz)
	End Do
	Dl = 1.d0 / (1.d0 + Dl)
	If (w_op .eqv. .True.) Then
		w_arrays_out(1,:) = tW
		Do i1 = 1, sizekw
			w_arrays_out(1+i1,:) = FW(i1,:)
			w_arrays_out(sizekw + 1 + i1,:) = FsW(i1,:)
		End Do
		w_arrays_out( 2 + (2*sizekw),: ) = DelzW
		w_arrays_out( 3 + (2*sizekw),: ) = DifftW
		w_arrays_out( 4 + (2*sizekw),: ) = msdW
		Deallocate(FW,FsW,DelzW,tW,F_w,Fs_w,DifftW,msdW)
	End If
	Deallocate(k,kw,S,lam)
	print*, "Dl = ", Dl
End Subroutine


!Subroutine that calculates the dynamics from small to long times until convergence
!for the diffussion coefficient is met
!The subroutine NEEDS:
! "d" -> system dimensionality
! "it" -> initial time for which we want to use small times dynamics subroutine
! "rho" -> system density
! "D0" -> Free particle diffussion coefficient
! "dt" -> time differential
! "indt" -> initial time differential
! "k" -> wave vector array
! "kw" -> wave vector weights array
! "lam" -> interpolating function λ(k) array
! "t" -> times array
! "S" -> Structure factor array
! THE FUNCTION RETURNS THE LONG TIME VALUES OF:
! "F" -> Intermediate Scattering function array
! "Fs" -> self part of "F"
! "Delz" -> memory function
! AND WRITES UPON something ALL THE TIMES VALUES OF:
! "F(iktest,:)"
! "Fs(iktest,:)"
Subroutine Long_t_dynamics_Dl_convergence(dyn_p,arrays,Dl,w_op,w_arrays_in,dyn_units)
	Implicit None
	Real * 8, dimension(7), intent(in) :: dyn_p
	Real * 8, dimension(:,:), intent(in) :: arrays
	Real * 8, dimension(:,:), intent(in) , optional :: w_arrays_in
	Integer, dimension(:) :: dyn_units
	Real * 8, dimension(nint(dyn_p(3))) :: t !Array of times
	Real * 8, intent(out) :: Dl !long time diffussion coefficient, MSD
	Real * 8, dimension(:), Allocatable :: k, kw !Wave vector array and weights array
	Real * 8, dimension(:), Allocatable :: S, lam !Structure factor array, Lambda(k) array
	Real * 8, dimension(size(arrays,2),nint(dyn_p(3))) :: F,Fs !Arrays of Intermediate Scattering Function and its self part
	Real * 8, dimension(nint(dyn_p(3))) :: Delz !Array of the memory function
	Real * 8, dimension(nint(dyn_p(3))) :: Difft !Diffussion coefficient array
	Real * 8, dimension(nint(dyn_p(3))) :: msd !Mean square displacement array
	Real * 8, dimension(:), Allocatable :: tau_a
	Integer :: d !Dimension of the system
	Integer :: it,nt !Short times points and medium times points
	Real * 8 :: rho !density of the system
	Real * 8 :: D0 !free particle diffussion coefficient
	Real * 8 :: indt !time differential
	Integer ::	decimations !Number of decimations used
	Integer :: i1, i2, i3, sizek, sizekw, itw, alW, itm
	real * 8 :: time1,time2
	real * 8 :: dt
	!Writting Variables:
	Logical :: w_op !Select .True. for writting down the information
	Real * 8, dimension(:), allocatable :: k_w, S_w, lam_w !writting variables needed
	Real * 8, dimension(:,:), allocatable :: F_w, Fs_w
	!Array for the integral of the memory function
	Real * 8, dimension(nint(dyn_p(3))) :: ti,tiw, Delzi
	Real * 8, dimension(1 + (nint(dyn_p(3)))/2 ) :: t_d,t_dw, Delz_d
	!Diffussion Coefficient convergence parameters
	Real * 8, parameter :: Dl_zero = 1.d-10, Dl_tol = 1.d-5
	Logical :: Dl_convergence
	Real * 8 :: int_Delz, int_Delz_dum, Dl_prev, error_Dl
	!Data unpacking!!!!!
	d  = nint(dyn_p(1))
	it = nint(dyn_p(2))
	nt = nint(dyn_p(3))
	decimations = nint(dyn_p(4))
	rho = dyn_p(5)
	D0  = dyn_p(6)
	indt = dyn_p(7)
	sizek  = size(arrays,2)
	!k-Functions
	Allocate(k(sizek),kw(sizek),S(sizek),lam(sizek))
	k   = arrays(1,:)
	kw  = arrays(2,:)
	S   = arrays(3,:)
	lam = arrays(4,:)
	If (w_op .eqv. .True.) Then
		!k-Functions writting-arrays inputs
		sizekw  = size(w_arrays_in,2)	!Size of writting wave vectors array
	!	alW = (nt + ( decimations * nt/2 )) !Size of array needed for the writting variables
		!Internal memory allocation
		Allocate(k_w(sizekw),S_w(sizekw),lam_w(sizekw),tau_a(sizek))
		Allocate (F_w(sizekw,nt), Fs_w(sizekw,nt))
		k_w   = w_arrays_in(1,:)
		S_w   = w_arrays_in(2,:)
		lam_w = w_arrays_in(3,:)
		write (dyn_units(1),*)"# 			time			k:", k_w(:)
		write (dyn_units(2),*)"#			time				k:", k_w(:)
		write (dyn_units(3),*)"#				time", "						Delz(t)", "							D(t)", "								msd(t)"
		write (dyn_units(5),*)"# 			k			tau_alpha"
		tau_a = 0.d0
	End If
	!!!!!!!!!!!!!!!!!
	do i1=1, nt
		t(i1) = i1 * indt
	end do
	dt = indt
	Print*, "Long times dynamics initialization"
	call short_times_dynamics(d,it,rho,D0,k,kw,t,S,F,Fs,Delz,Difft,msd)
	!Print*, "Short Times done"
	call medium_times_dynamics(d,it + 1,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,Difft,msd)
!	Print*, "Medium Times done"
	!Initial Dl* calculation
	Call Simpson_3_8_quad(ti,tiw,size(ti),0.d0,t(nt))
	Delzi(:) = Delz(:)
	int_Delz = sum( Delzi(:) * tiw(:) )
	Dl = 1.d0 / ( 1.d0 + int_Delz )
	If (w_op .eqv. .True.) Then
		call short_times_dynamics_writting(it,rho,D0,k_w,t,S_w,F_w,Fs_w)
		call medium_times_dynamics_writting(d,it + 1,rho,D0,dt,k_w,lam_w,t,S_w,F_w,Fs_w,Delz)
		Do i1 =1, nt
			write(dyn_units(1),*) ti(i1), F_w(:,i1)
			write(dyn_units(2),*) ti(i1), Fs_w(:,i1)
			write(dyn_units(3),*) t(i1), Delz(i1), Difft(i1), msd(i1)
		End Do
		call Fs_tau_a_close_val(Fs,t,1,tau_a)
	End If
	itm = (nt / 2) + 1
	Dl_convergence = .False.
	i3 = 0
	Do while (Dl_convergence .eqv. .False.)
		i3 = i3 + 1
		!print*, i3, "'  Decimation"
		dt = 2.d0 * dt
		call half_save(dt,t,Delz,F,Fs,Difft,msd)
		call medium_times_dynamics(d,itm,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,Difft,msd)
		If (w_op .eqv. .True.) Then
			Do i2 = 1 , nt / 2
				F_w(:,i2)  = F_w (:,2*i2)
				Fs_w(:,i2) = Fs_w(:,2*i2)
			End Do
			call medium_times_dynamics_writting(d,itm,rho,D0,dt,k_w,lam_w,t,S_w,F_w,Fs_w,Delz)
			Do i1 =itm, nt
				write(dyn_units(1),*) t(i1), F_w(:,i1)
				write(dyn_units(2),*) t(i1), Fs_w(:,i1)
				write(dyn_units(3),*) t(i1), Delz(i1), Difft(i1), msd(i1)
			End Do
			call Fs_tau_a_close_val(Fs,t,itm-1,tau_a)
		End If
		!Dl* calculation!!!!!!!!!!!!!!!
		Call Simpson_3_8_quad(t_d,t_dw,size(t_d),t(itm-1),t(nt))
		int_Delz = int_Delz + sum(Delz(itm-1:nt) * t_dw(:))
		Dl_prev = Dl
		Dl = 1.d0 / (1.d0 + int_Delz)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		error_Dl = Abs(Dl_prev/Dl - 1.d0)
		!print*, error_Dl, Dl
		!Exit Conditions
		If (error_Dl < Dl_tol) Dl_convergence = .True.
		If (Dl < Dl_zero) Dl_convergence = .True.
	End Do
	If (w_op .eqv. .True.) Then
		Do i1=1, sizek
			write(dyn_units(5),*) k(i1), tau_a(i1)
		End Do
		Deallocate(F_w,Fs_w,tau_a)
	End If
	Deallocate(k,kw,S,lam)
	print*, "Dl = ", Dl
End Subroutine

!Subroutine that discriminates between the closest value to t=e^-1 for the matrix Fs(k,t)
Subroutine Fs_tau_a_close_val(Fs,t,start,tau_a)
	Implicit None
	Real * 8, dimension(:), intent(out) :: tau_a
	Real * 8, dimension(:,:), intent(in) :: Fs
	Real * 8, dimension(:), intent(in) :: t
	Integer, intent(in) :: start
	Real * 8, parameter :: e_val = exp(-1.d0)
	Real * 8 :: error, F0,F1, t0, t1,m,b
	Integer :: i1,i2, sizek,sizet
	sizek = size(tau_a)
	sizet = size(t)
	Do i1 = start, sizet-1
		Do i2=1, sizek
			F0 = Fs(i2,i1)
			F1 = Fs(i2,i1+1)
			If ( F0 > e_val .and. F1 < e_val ) Then
				t0 = t(i1)
				t1 = t(i1+1)
				m = (F1 - F0) / (t1-t0)
				b = F0 - (m * t0)
				tau_a(i2) = (e_val - b)/m
			Else
			End If
		End Do
	End Do
End Subroutine



!Subroutine that calculates the dynamics from small to long times until convergence
!for the diffussion coefficient is met
!The subroutine NEEDS:
! "d" -> system dimensionality
! "it" -> initial time for which we want to use small times dynamics subroutine
! "rho" -> system density
! "D0" -> Free particle diffussion coefficient
! "dt" -> time differential
! "indt" -> initial time differential
! "k" -> wave vector array
! "kw" -> wave vector weights array
! "lam" -> interpolating function λ(k) array
! "t" -> times array
! "S" -> Structure factor array
! THE FUNCTION RETURNS THE LONG TIME VALUES OF:
! "F" -> Intermediate Scattering function array
! "Fs" -> self part of "F"
! "Delz" -> memory function
! AND WRITES UPON something ALL THE TIMES VALUES OF:
! "F(iktest,:)"
! "Fs(iktest,:)"
Subroutine Long_t_dynamics_Dl_convergence_zero_del_z(dyn_p,arrays,Dl,w_op,w_arrays_in,dyn_units)
	Implicit None
	Real * 8, dimension(7), intent(in) :: dyn_p
	Real * 8, dimension(:,:), intent(in) :: arrays
	Real * 8, dimension(:,:), intent(in) , optional :: w_arrays_in
	Integer, dimension(3) :: dyn_units
	Real * 8, dimension(nint(dyn_p(3))) :: t !Array of times
	Real * 8, intent(out) :: Dl !long time diffussion coefficient, MSD
	Real * 8, dimension(:), Allocatable :: k, kw !Wave vector array and weights array
	Real * 8, dimension(:), Allocatable :: S, lam !Structure factor array, Lambda(k) array
	Real * 8, dimension(size(arrays,2),nint(dyn_p(3))) :: F,Fs !Arrays of Intermediate Scattering Function and its self part
	Real * 8, dimension(nint(dyn_p(3))) :: Delz !Array of the memory function
	Real * 8, dimension(nint(dyn_p(3))) :: Difft !Diffussion coefficient array
	Real * 8, dimension(nint(dyn_p(3))) :: msd !Mean square displacement array
	Integer :: d !Dimension of the system
	Integer :: it,nt !Short times points and medium times points
	Real * 8 :: rho !density of the system
	Real * 8 :: D0 !free particle diffussion coefficient
	Real * 8 :: indt !time differential
	Integer ::	decimations !Number of decimations used
	Integer :: i1, i2, i3, sizek, sizekw, itw, alW, itm,i4
	real * 8 :: time1,time2
	real * 8 :: dt
	!Writting Variables:
	Logical :: w_op !Select .True. for writting down the information
	Real * 8, dimension(:), allocatable :: k_w, S_w, lam_w !writting variables needed
	Real * 8, dimension(:,:), allocatable :: F_w, Fs_w
	!Array for the integral of the memory function
	Real * 8, dimension(nint(dyn_p(3))) :: ti,tiw, Delzi
	Real * 8, dimension(1 + (nint(dyn_p(3)))/2 ) :: t_d,t_dw, Delz_d
	!Diffussion Coefficient convergence parameters
	Real * 8, parameter :: Dl_zero = 1.d-14, Dl_tol = 1.d-10
	Logical :: Dl_convergence
	Real * 8 :: int_Delz, int_Delz_dum, Dl_prev, error_Dl
	!Zeros of Del Z variables:
	Real * 8, dimension(2**6) :: z_tests
	Real * 8 , dimension(2**6,size(t)):: zero_delz_tests,zero_delz_tests_prod
	!Data unpacking!!!!!
	d  = nint(dyn_p(1))
	it = nint(dyn_p(2))
	nt = nint(dyn_p(3))
	decimations = nint(dyn_p(4))
	rho = dyn_p(5)
	D0  = dyn_p(6)
	indt = dyn_p(7)
	sizek  = size(arrays,2)
	!k-Functions
	Allocate(k(sizek),kw(sizek),S(sizek),lam(sizek))
	k   = arrays(1,:)
	kw  = arrays(2,:)
	S   = arrays(3,:)
	lam = arrays(4,:)
	If (w_op .eqv. .True.) Then
		!k-Functions writting-arrays inputs
		sizekw  = size(w_arrays_in,2)	!Size of writting wave vectors array
	!	alW = (nt + ( decimations * nt/2 )) !Size of array needed for the writting variables
		!Internal memory allocation
		Allocate(k_w(sizekw),S_w(sizekw),lam_w(sizekw))
		Allocate (F_w(sizekw,nt), Fs_w(sizekw,nt))
		k_w   = w_arrays_in(1,:)
		S_w   = w_arrays_in(2,:)
		lam_w = w_arrays_in(3,:)
		write (dyn_units(1),*)"# 			time			k:", k_w(:)
		write (dyn_units(2),*)"#			time				k:", k_w(:)
		write (dyn_units(3),*)"#				time", "						Delz(t)", "							D(t)", "								msd(t)"
	End If
	!!!!!!!!!!!!!!!!!
	do i1=1, nt
		t(i1) = i1 * indt
	end do
	dt = indt
	Print*, "Long times dynamics initialization"
	call short_times_dynamics(d,it,rho,D0,k,kw,t,S,F,Fs,Delz,Difft,msd)
	Print*, "Short Times done"
	call medium_times_dynamics(d,it + 1,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,Difft,msd)
	Print*, "Medium Times done"
	!Initial Dl* calculation
	Call Simpson_3_8_quad(ti,tiw,size(ti),0.d0,t(nt))
	Delzi(:) = Delz(:)
	int_Delz = sum( Delzi(:) * tiw(:) )
	Dl = 1.d0 / ( 1.d0 + int_Delz )
	If (w_op .eqv. .True.) Then
		call short_times_dynamics_writting(it,rho,D0,k_w,t,S_w,F_w,Fs_w)
		call medium_times_dynamics_writting(d,it + 1,rho,D0,dt,k_w,lam_w,t,S_w,F_w,Fs_w,Delz)
		Do i1 =1, nt
			write(dyn_units(1),*) ti(i1), F_w(:,i1)
			write(dyn_units(2),*) ti(i1), Fs_w(:,i1)
			write(dyn_units(3),*) t(i1), Delz(i1), Difft(i1), msd(i1)
		End Do
		!zeros of Del Z
		Open (unit= 301, file="zero_delz_test.dat", status="Replace")
		call zero_Delta_zeta(d,it+1,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,z_tests,zero_delz_tests)
		Do i1=it, nt
			Write (301,*) t(i1), Delz(i1)
			Do i4=1, size(z_tests)
				write(301,*) z_tests(i4), zero_delz_tests(i4,i1)
			End Do
		End Do
	End If
	itm = (nt / 2) + 1
	Dl_convergence = .False.
	i3 = 0
	Do while (Dl_convergence .eqv. .False.)
		i3 = i3 + 1
		!print*, i3, "'  Decimation"
		dt = 2.d0 * dt
		call half_save(dt,t,Delz,F,Fs,Difft,msd)
		call medium_times_dynamics(d,itm,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,Difft,msd)
		If (w_op .eqv. .True.) Then
			Do i2 = 1 , nt / 2
				F_w(:,i2)  = F_w (:,2*i2)
				Fs_w(:,i2) = Fs_w(:,2*i2)
			End Do
			call medium_times_dynamics_writting(d,itm,rho,D0,dt,k_w,lam_w,t,S_w,F_w,Fs_w,Delz)
			Do i1 =itm, nt
				write(dyn_units(1),*) t(i1), F_w(:,i1)
				write(dyn_units(2),*) t(i1), Fs_w(:,i1)
				write(dyn_units(3),*) t(i1), Delz(i1), Difft(i1), msd(i1)
			End Do
			call zero_Delta_zeta(d,it+1,rho,D0,dt,k,kw,lam,t,S,F,Fs,Delz,z_tests,zero_delz_tests)
			Do i1=itm, nt
				Write (301,*) t(i1), Delz(i1)
				Do i4=1, size(z_tests)
					write(301,*) z_tests(i4), zero_delz_tests(i4,i1)
				End Do
			End Do
		End If
		!Dl* calculation!!!!!!!!!!!!!!!
		Call Simpson_3_8_quad(t_d,t_dw,size(t_d),t(itm-1),t(nt))
		int_Delz = int_Delz + sum(Delz(itm-1:nt) * t_dw(:))
		Dl_prev = Dl
		Dl = 1.d0 / (1.d0 + int_Delz)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		error_Dl = Abs(Dl_prev/Dl - 1.d0)
		!print*, error_Dl, Dl
		!Exit Conditions
		If (error_Dl < Dl_tol) Dl_convergence = .True.
		If (Dl < Dl_zero) Dl_convergence = .True.
	End Do
	If (w_op .eqv. .True.) Then
		Deallocate(F_w,Fs_w)
	End If
	Deallocate(k,kw,S,lam)
	print*, "Dl = ", Dl
End Subroutine

!This subroutine calculates the limits of Fc and Fs whe tau tends to infinity
Subroutine limit_t_infty_Fc_Fs_mono(D0,gamma_val,k,lamk,Sk,fc,fs)
	Implicit None
	Real * 8, intent(in) :: D0, gamma_val
	Real * 8, dimension(:), intent(in) :: k, lamk, Sk
	Real * 8, dimension(:), intent(out) :: fc, fs
	Real * 8 :: Delz_infty
	Real * 8, dimension(size(k)) :: k2, lamkDelz
	If (gamma_val < 1.d20 ) Then
		k2 = k*k
		Delz_infty = 1.d0 / gamma_val
		lamkDelz = lamk * Delz_infty
		fc = lamkDelz * Sk     / (lamkDelz + (k2 * D0 / Sk) )
		fs = lamkDelz / (lamkDelz + (k2 * D0) )
	Else
		Print*, "System is not arrested"
	End If
End Subroutine


!Subroutine that writes down the dynamic properties and info file
Subroutine dynamics_writting(k_w, array_w, units)
	Implicit None
	Integer, intent(in), dimension(:) :: units
	Real * 8, dimension(:,:), intent(in) :: array_w
	Real * 8, dimension(:), intent(in) ::  k_w
	Integer :: i1,sizek,sizet
	sizek = size(k_w)
	sizet = size(array_w,2)
	write (units(1),*)"# 			time			k:", k_w(:)
	write (units(2),*)"#			time				k:", k_w(:)
	write (units(3),*)"#				time", "						Delz(t)", "							D(t)", "								msd(t)"
	Do i1 = 1,sizet
		write(units(1),*) array_w(1,i1), array_w(2:sizek+1,i1)
		write(units(2),*) array_w(1,i1), array_w(sizek+2: 2 * sizek + 1,i1)
		write(units(3),*) array_w(1,i1), array_w(2*sizek+2,i1), array_w(2*sizek+3,i1), array_w(2*sizek+4,i1)
	End do
End Subroutine

!
!This subroutines calculates the dynamics of a quench given the initial and final structure factor
Subroutine SCGLE_dynamics(dyn_p,s_arrays,bt,writting_op,s_w_arrays,dyn_w_arrays,dyn_units)
	Implicit None
	Real * 8, dimension(:,:), Intent(in) :: s_arrays
	Real * 8, dimension(:,:), Intent(in) :: s_w_arrays
	Real * 8, dimension(:,:), Intent(out) :: dyn_w_arrays
	Real * 8, dimension(:), Intent(in) :: dyn_p
	Integer, dimension(:), Intent(in) :: dyn_units
	Real * 8, intent(out) :: bt
	Logical, intent(in) :: writting_op
	Real * 8 	:: u_val, ua_val
	call Long_t_dynamics(dyn_p,s_arrays,bt,writting_op,s_w_arrays,dyn_w_arrays)
	If (writting_op .eqv. .True. ) Then
		call dynamics_writting(s_arrays(1,:), dyn_w_arrays, dyn_units)
		call dynamics_info_writting(dyn_p, dyn_units(4))
	Else
	End If
End Subroutine

!Subroutine that calculates the alpha relaxation time from a writting array
Subroutine alpha_relax_time(w_arrays_out,tau_alphas)
	Implicit None
	Real * 8, dimension(:,:), intent(in) :: w_arrays_out
	Real * 8, dimension(:), intent(out) :: tau_alphas
	Real * 8 :: tau1,tau2, F1, F2, a, b
	Integer :: sizekw, sizetw, i1, i2, itau1, itau2
	Real * 8 :: error,error_test
	sizetw = size(w_arrays_out,2)
	sizekw = (size(w_arrays_out,1) - 4) / 2
	Print*, sizetw, sizekw
	tau_alphas = 0.d0
	Do i1 = 1, sizekw
		error = 1.d0
		Do i2 = 1, sizetw
			error_test =  abs( exp(-1.d0) - w_arrays_out(sizekw + 1 + i1,i2) )
			If (error > error_test ) then
				tau2 = w_arrays_out(1,i2)
				error = error_test
				itau2 = i2
				itau1 = itau2 - 1
				F2 = w_arrays_out(sizekw + 1 + i1,i2)
			Else
			End If
		End Do
		F1 = w_arrays_out(sizekw + 1 + i1,itau1)
		tau1 = w_arrays_out(1,itau1)
		a = (F2 - F1) / (tau2 - tau1)
		b = F1 - a * tau1
		tau_alphas(i1) = ( exp(-1.d0) - b ) / a
	End Do
	Print*, "Taus done"
	Return
End Subroutine

!Subroutine that writes down the dynamic properties and info file
Subroutine dynamics_info_writting(dyn_p, i_unit)
	Implicit None
	Integer, intent(in) :: i_unit
	Real * 8, dimension(:), intent(in) :: dyn_p
	Integer :: it, mt, decim
	Real * 8 :: dt, D0, rho, t_max
	Integer :: i1, d
	d     = nint(dyn_p(1) )
  it   = nint(dyn_p(2))
	mt   = nint(dyn_p(3))
	decim = nint(dyn_p(4) )
	rho   = dyn_p(5)
	D0    = dyn_p(6)
	t_max = mt * dt
	Do i1=1, decim
		dt = dyn_p(7) * (2 ** i1)
		t_max = t_max + (mt * dt / 2.d0)
	End Do
	dt = dyn_p(7)
	Write(i_unit,*) "### Dynamics General Info ###"
	Write(i_unit,*) "number of points						=", mt + (mt * decim / 2)
	Write(i_unit,*) "time minimum value					=", dt
	Write(i_unit,*) "time maximum value					=", t_max
	Write(i_unit,*) "number of initial times		=", it
	Write(i_unit,*) "number of medium times 		=", mt
	Write(i_unit,*) "number of decimations			=", decim
	Write(i_unit,*) "initial time differential	=", dt
	Write(i_unit,*) "system density							=", rho
	Write(i_unit,*) "Free particle Diff Coeff 	=", D0
	Write(i_unit,*) "###############################"
End Subroutine

Subroutine dynamics_Long_t_info_writting(Dl,gam,i_unit)
	Implicit None
	Integer, intent(in) :: i_unit
	Real * 8, intent(in) :: Dl, gam
	Write(i_unit,*) "### Long time dynamics information: ###"
	Write(i_unit,*) "Diffussion coefficient			=", Dl
	Write(i_unit,*) "MSD limit or gamma					=", gam
	Write(i_unit,*) "###############################"
End Subroutine

!Subroutine that writes down the dynamic properties and info file
Subroutine lambda_info_writting(kc,i_unit)
	Implicit None
	Integer, intent(in) :: i_unit
	Real * 8, intent(in) :: kc
	Integer :: i1
	Write(i_unit,*) "### Wave-Vector Interpolation Function ###"
	Write(i_unit,*) "k_c used						=", kc
	Write(i_unit,*) "###############################"
End Subroutine

!Gamma function of the theory, asociated with msd limit at t -> infinity, in terms of:
! "d" -> system dimensionality
! "rho" -> system density
! "k" -> wave vector arrays
! "kw" -> wave vector weights array
! "S" -> structure factor array
! "lam" -> lambda function array
Function gamma(d,rho,k,kw,S,lam)
	Implicit None
	Real * 8, dimension(:), intent(in) :: k,kw,S,lam
	Real * 8, intent(in) :: rho
	Real * 8, dimension(size(k)) :: kdum,lams,k2,Integrand, gammak2
	Integer, intent(in) :: d
	Real * 8 :: gamma
	Real * 8, parameter :: tol=1.d-8, max_val = 1.d20
	Real * 8 :: dimen_dum, error, gamma_test
	Logical :: convergence
	dimen_dum = d * ( (2.d0 * pi ) ** d ) * rho
	dimen_dum = dimen_dum / (2.d0 * (d-1.d0) * pi)
	k2 =  k ** 2
	kdum = ( ( (S-1.d0) * lam) ** 2 ) * kw * ( k ** (d+1) )
	lams = S * lam
	convergence = .False.
	gamma = 1.d-7
	Do while (convergence .eqv. .False.)
		gammak2 = gamma * k2
		Integrand = kdum / ((lams + gammak2) * (lam + gammak2))
		gamma_test =  dimen_dum / sum(Integrand)
		error = abs(1.d0 - (gamma/gamma_test))
		gamma = gamma_test
		if ( error < tol ) convergence = .True.
		if ( gamma > max_val ) convergence = .True.
	End do
End Function

!This subroutine calculates the arrest value of the "i_param" of a structure with
!upper limit "u_l" and low limit "l_l" while giving their respectives values "g_u"
!and "g_l", given a tolerance "tol"
Subroutine arrest_bound(k,kw,lam,sys_p,sys,sys_approx,i_param,d,u_l,l_l,g_u,g_l,tol)
	Implicit None
	Real * 8, intent(in), dimension(:) :: k, kw, sys_p,lam
	Integer, intent(in) :: i_param, d
	Real * 8, intent(inout) :: u_l, l_l
	Real * 8, intent(in) :: tol
	Real * 8, dimension(size(k)) :: S_up,S_lo,S_test
	Real * 8, intent(out) :: g_u, g_l
	Real * 8 :: g_test, k_val, error
	Real * 8 :: rho_up, rho_lo, rho_test
	Real * 8, parameter :: g_max = 1.d20
	Real * 8 :: rho_dum
	Real * 8, dimension(size(sys_p)) :: s_p_u, s_p_l, s_p_t
	Character(*), intent(in) :: sys, sys_approx
	Integer :: i1, knp, i_k
	Logical :: convergence
	rho_dum = 2.d0 * d / pi
	i_k = size(sys_p)
	knp = size(k)
	!INITIAL PARAMETERS
	s_p_u = sys_p
	s_p_l = sys_p
	s_p_t = sys_p
	s_p_u(i_param) = u_l
	s_p_l(i_param) = l_l
	!INITIAL STRUCTURES
	Do i1 = 1, knp
		s_p_u(i_k) = k(i1)
		s_p_l(i_k) = k(i1)
		S_up(i1) = structure("s0",sys,sys_approx,s_p_u)
		S_lo(i1) = structure("s0",sys,sys_approx,s_p_l)
	End Do
	rho_up = s_p_u(1) * rho_dum
	rho_lo = s_p_l(1) * rho_dum
	g_u = gamma(d,rho_up,k,kw,S_up,lam)
	Print*, g_u
	g_l = gamma(d,rho_lo,k,kw,S_lo,lam)
	Print*, g_l
	convergence = .False.
	If (g_u > g_max .and. g_l < g_max) Then
		Do while (convergence .eqv. .False.)
			s_p_t(i_param) = (u_l + l_l) / 2.d0
			!Test structure
			Do i1 = 1, knp
				s_p_t(i_k) = k(i1)
				S_test(i1) = structure("s0",sys,sys_approx,s_p_t)
			End Do
			rho_test =  s_p_t(1) * rho_dum
			g_test = gamma(d,rho_test,k,kw,S_test,lam)
			If ( g_test > g_max ) Then
				g_u = g_test
				u_l = s_p_t(i_param)
			Else
				g_l = g_test
				l_l = s_p_t(i_param)
			End If
			error = (u_l - l_l) / u_l
			If (error < tol ) convergence = .True.
			print*, u_l,l_l,error,g_test
		End Do
	Else If (g_u < g_max .and. g_l > g_max)	Then
		Do while (convergence .eqv. .False.)
			s_p_t(i_param) = (u_l + l_l) / 2.d0
			!Test structure
			Do i1 = 1, knp
				s_p_t(i_k) = k(i1)
				S_test(i1) = structure("s0",sys,sys_approx,s_p_t)
			End Do
			rho_test =  s_p_t(1) * rho_dum
			g_test = gamma(d,rho_test,k,kw,S_test,lam)
			If ( g_test < g_max ) Then
				g_u = g_test
				u_l = s_p_t(i_param)
			Else
				g_l = g_test
				l_l = s_p_t(i_param)
			End If
			error = (u_l - l_l) / u_l
			If (error < tol ) convergence = .True.
			print*, u_l,l_l,error,g_test
		End Do
	Else If (g_u > g_max .and. g_l > g_max) Then
		Print*, "No arrest found in between the limiting parameters"
	Else If (g_u < g_max .and. g_l < g_max) Then
		Print*, "No ergodic state found in between the limiting parameters"
	End If
	!open(unit=12, file= "test.dat", status="Replace")
	!Do i1 = 1 , knp
	!	Write(12,*) k(i1), s_up(i1), s_lo(i1), s_test(i1)
	!End Do
	!close(12)
	return
End Subroutine

Subroutine zero_Delta_zeta(d,it,rho,D0,dt,k,kw,lam,t,S,Fc,Fs,Delz,z_tests,zero_delz_tests)
	Implicit None
	Integer, intent(in) :: d !Dimension of the system
	Integer, intent(in) :: it !Beggining of medium times where we are going to start
	Real * 8, intent(in) :: rho !density of the system
	Real * 8, intent(in) :: D0 !free particle diffussion coefficient
	Real * 8, intent(in) :: dt !time differential
	Real * 8, intent(in), dimension(:) :: k, kw !Wave vector array and weights array
	Real * 8, intent(in), dimension(:) :: t !Array of times
	Real * 8, intent(in), dimension(:) :: S, lam !Structure factor array, Lambda(k) array
	Real * 8, intent(in), dimension(:,:) :: Fc,Fs !Arrays of Intermediate Scattering Function and its self part
	Real * 8, intent(in), dimension(:) :: Delz !Array of the memory function
	Real * 8, dimension(size(k)) :: Fcv, Fsv,dFc1,dFs1,Fcdum,Fsdum,Fs1,Fc1 !Views on a time for F and Fs
	Real * 8, dimension(size(k)) :: ac, as, lamc,lams!Dummies for discretization
	Real * 8, dimension(size(Delz)) :: DeltaDelz !Array of the memory function
	Real * 8 :: delztest, error, a_approx, b_approx
	Real * 8, intent(out), dimension(:,:) :: zero_delz_tests
	Real * 8, intent(in), dimension(:) :: z_tests
	Integer :: i1,i2
	Logical :: convergence
	Real * 8, parameter :: tol = 1.d-7
	ac = alphac(D0,dt,k,S,lam,Delz(1))
	as = alphas(D0,dt,k,lam,Delz(1))
	lamc = ac * lam
	lams = as * lam
	Fs1 = Fs(:,1)
	Fc1 = Fc(:,1)
	dFc1 = S - Fc1
	dFs1 = 1.d0 - Fs1
	Do i1=2, size(t)
		deltadelz(i1) = Delz(i1)-Delz(i1-1)
	End Do
	Do i1=it, size(t)
		!NON ITERATIVE PART
		Fsdum = 0.d0
		Fcdum = 0.d0
		!$OMP PARALLEL DO DEFAULT(NONE) SHARED(Fc,Fs,i1,deltadelz), PRIVATE(i2) REDUCTION(+:Fsdum,Fcdum)
		Do i2=2, i1-1
			!OMP WORKSHARE
			Fsdum = Fsdum - (deltadelz(i1+1-i2) * Fs(:,i2))
			Fcdum = Fcdum - (deltadelz(i1+1-i2) * Fc(:,i2))
			!OMP END WORKSHARE
		End Do
		!$OMP END PARALLEL DO
		!Fcv = F(:,i1-1)
		!Fsv = Fs(:,i1-1)
		Fsdum = as * (lam * (Fsdum + Delz(i1-1) * Fs1) + (Fs(:,i1-1)/dt))
		Fcdum = ac * (lam * (Fcdum + Delz(i1-1) * Fc1) + (Fc(:,i1-1)/dt))
		!ITERATIVE PART
		convergence = .False.
		!delztest = Delz(i1-1) !Can be interpolated for a better approximation
		b_approx = (log(delz(i1-2)) - log(delz(i1-1))) / (dt)
		a_approx = exp(log(delz(i1-1)) + ( b_approx * dt * (i1-1) ) )
		!delztest = a_approx * exp(-b_approx * dt * i1)
		Do i2=1, size(z_tests)
			delztest = z_tests(i2)
			!Do while (convergence .eqv. .False.)
			Fcv =  Fcdum + lamc * delztest * dFc1
			Fsv =  Fsdum + lams * delztest * dFs1
			zero_delz_tests(i2,i1) = delta_z(D0,rho,k,kw,S,Fcv,Fsv,d) - delztest
		End Do
	End Do
End Subroutine

End Module
