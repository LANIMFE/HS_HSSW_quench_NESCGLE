Module quadratures
  Use math_constants
	Implicit None

Contains
  !This subroutine select the appropiate subroutine for whatever quadrature in terms
  !of the parameters of each quadrature; look inside the function for the parameters
  !order for each quadrature
	! "quad" accepts:
  !                 "CLCU" -> Clenshaw Curtis quadrature
  ! "params_i" -> array with the needed integer parameters
  ! "params_i" -> array with the needed floating point parameters
  ! "qp" -> quadrature points array
  ! "qw" -> quadrature weights array
	Subroutine quad_select(quad,params,qp,qw)
		Implicit None
		character(len=4), intent(in) :: quad !quadrature selector variable
    real * 8, dimension(:), intent(in) :: params !quadrature floating parameters
		real * 8, dimension(:), intent(out) :: qp, qw !quadrature points
    !real * 8, dimension(:), optional, intent(out) :: qw !quadrature weights
    integer :: np
    real * 8 :: xmin,xmax,dr
    xmin = params(1)
    xmax = params(2)
    np   = nint( params(3) )
    select case(quad)
    !Clenshaw Curtis quadrature
    case("CLCU")
      !INTEGER parameters needed order:
      !   Number of points "np"
      !Floating point parameters needed order:
      !   Domain minimum value "xmin"
      !   Domain maximum value "xmax"
      call Clenshaw_Curtis_quadrature(np,qp,qw)
      call rescale_OW(np,qp,qw,-1.d0,1.d0,xmin,xmax)
    !Rectangles quadratures
    case("RECT")
        !INTEGER parameters needed order:
        !   Number of points "np"
        !Floating point parameters needed order:
        !   Domain minimum value "xmin"
        !   Domain maximum value "xmax"
        dr = (xmax-xmin) / dble(np)
        qw = dr
        call rectangles(qp,xmin,dr)
    end select
    return
	End Subroutine

	Subroutine Clenshaw_Curtis_quadrature(order,x,w)
		Implicit None
    Integer, intent(in) :: order
    Real * 8, dimension(order), intent(out) :: x, w
	  Real * 8 :: b
	  Integer :: i
	  Integer :: j
	  Real * 8 :: theta
	  If ( order < 1 ) Then
	    write ( *, '(a)' ) ' '
	    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
	    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
	    stop
	  End If
	  If ( order == 1 ) Then
	    x(1) = 0.0D+00
	    w(1) = 2.0D+00
	    return
	  End If
		Do i = 1, order
	    x(i) = cos ( dble ( order - i ) * pi / dble ( order - 1 ) )
	  End Do
	  x(1) = -1.0D+00
	  If ( mod ( order, 2 ) == 1 ) Then
	    x((order+1)/2) = 0.0D+00
	  End If
	  x(order) = +1.0D+00
	  Do i = 1, order
	    theta = dble( i - 1 ) * pi / dble( order - 1 )
	    w(i) = 1.0D+00
	    Do j = 1, ( order - 1 ) / 2
	      If ( 2 * j == ( order - 1 ) ) Then
	        b = 1.0D+00
	      Else
	        b = 2.0D+00
	      End If
	      w(i) = w(i) - b * cos ( 2.0D+00 * dble( j ) * theta ) / dble( 4 * j * j - 1 )
	    End Do
	  End Do
	  w(1) = w(1) / dble( order - 1 )
	  w(2:order-1) = 2.0D+00 * w(2:order-1) / dble( order - 1)
	  w(order) = w(order) / dble( order - 1 )
	  return
	End Subroutine

  Subroutine rectangles(x,x_min,dr)
		Implicit None
	  Real * 8, intent(out), dimension(:) :: x
    Real * 8, intent(in) :: dr,x_min
    integer :: i1
	  Do i1=1, size(x)
      x(i1) = x_min + i1 * dr
    End Do
	  return
	End Subroutine

  !This subroutines rescales the abscissas and weights values from a range of
  !x ϵ [xa, xb] to y ϵ [ya, yb] and obtains the new weights in the new arrays
  !y, wy
  Subroutine rescale_NA(n,x,wx,xa,xb,y,wy,ya,yb)
  	implicit none
  	integer, intent(in) :: n
  	real * 8, intent(in) :: xa,xb,ya,yb
  	real * 8, dimension(n), intent(in) :: x,wx
  	real * 8, dimension(n), intent(out) :: y,wy
  	y(1:n) = ya  + (( yb - ya ) * (x(1:n)-xa) / ( xb - xa ))
    wy(1:n) = ( yb - ya ) * wx(1:n) / ( xb - xa)
  End Subroutine

  !This subroutines rescales the abscissas and weights values from a range of
  !x ϵ [xa, xb] to y ϵ [ya, yb] and obtains the new wy in the arrays x,wx
  Subroutine rescale_OW(n,x,wx,xa,xb,ya,yb)
  	implicit none
  	integer, intent(in) :: n
  	real * 8, intent(in) :: xa,xb,ya,yb
  	real * 8, dimension(n), intent(inout) :: x,wx
  	x(1:n) = ( ( ya + yb ) + ( yb - ya ) * x(1:n) ) / ( xb - xa )
    wx(1:n) = ( yb - ya ) * wx(1:n) / ( xb - xa)
  End Subroutine

  !This subroutine writes down the information of a given quadrature used
	Subroutine quad_writting(quad,i_unit,params)
    Implicit None
    Real * 8, dimension(:), intent(in) :: params
    Character(len=4), intent(in) :: quad
    Integer, intent(in) :: i_unit
    Select case (quad)
			case("CLCU")
				Write(i_unit,*) "				Clenshaw-Curtis"
				Write(i_unit,*) "quadrature minimum value    = ", params(1)
				Write(i_unit,*) "quadrature maximum value    = ", params(2)
        Write(i_unit,*) "quadrature number of points = ", params(3)
      case("RECT")
        Write(i_unit,*) "				Rectangles"
        Write(i_unit,*) "quadrature minimum value    = ", params(1)
        Write(i_unit,*) "quadrature differential     = ", params(2)
			case Default
				Write(i_unit,*) "				"//quad
        Write(i_unit,*) "quadrature minimum value    = ", params(1)
				Write(i_unit,*) "quadrature maximum value    = ", params(2)
        Write(i_unit,*) "quadrature number of points = ", params(3)
		End Select
		Write(i_unit,*) "###############################"
  End Subroutine

  !This subroutine generates the quadrature points and weights of
  !up to Simpson 3/8 quadrature given an order greater than 0
  Subroutine Simpson_3_8_quad(qp,qw,order,qpmin,qpmax)
    Implicit None
    Real * 8, intent(out), dimension(:) :: qp, qw !Quadrature points and quadrature weights
    Real * 8, intent(in) :: qpmin, qpmax
    Integer, intent(in) :: order !Order of the quadrature
    Real * 8 :: dx
    Integer :: i1, i2, ratio, mod_ratio, imax
		dx = (qpmax-qpmin) / dble(order)
    !Calculating the quadrature points
    Do i1 = 1, order
      qp(i1) = qpmin + (i1 * dx)
    End Do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (order<1) Then
      Print *, 'Error at "Simpson_3_8_quad", not a valid order '
    Else If (order == 1) Then
      qw(:) = dx
      return
    Else If (order==2) Then
      qw(:) =  0.5d0 * dx
			return
		Else If (order==3) Then
      qw(1) = 1.d0
      qw(2) = 4.d0
      qw(3) = 1.d0
      qw(:) = qw(:) * dx / 3.d0
			return
		Else If (order==4) Then
      qw(1) = 1.d0
      qw(2) = 3.d0
      qw(3) = 3.d0
      qw(4) = 1.d0
      qw(:) = qw(:) * dx * 0.375d0
			return
		Else If (order==5) Then
      qw(1) = 1.d0
      qw(2) = 4.d0
      qw(3) = 2.d0
      qw(4) = 4.d0
      qw(5) = 1.d0
      qw(:) = qw(:) * dx / 3.d0
			return
		Else If (order>5) Then
      ratio = order/3
			mod_ratio = mod(order,3)
      qw(1) = 1.d0
      Do i1=1, ratio
        i2 = (i1 * 3) - 1
        qw(i2)   = 3.d0
        qw(i2+1) = 3.d0
        if(i2+2 < order) qw(i2+2) = 2.d0
      End Do
      qw(:) = qw(:) * 3.d0 * dx / 8.d0
			If	(mod_ratio == 0) Then
				qw(order-2) = ((3.d0/8.d0) + (1.d0/3.d0)) * dx
        qw(order-1) = 4.d0 * dx / 3.d0
        qw(order)   = 1.d0 * dx / 3.d0
        return
      Else If (mod_ratio==1) Then
        qw(order-1) = (3.d0 / 8.d0 + 0.5d0) * dx
        qw(order)   =  0.5d0 * dx
        return
      Else If (mod_ratio==2) Then
				qw(order) = (3.d0/8.d0) * dx
        return
      End If
    End If
  End Subroutine


End Module
