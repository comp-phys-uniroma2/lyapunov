module functions
  use precision
  implicit none
  private

  public :: ff
  public :: pendolo
  public :: harmonic
  public :: pendolo_auto
  public :: lpendolo_auto
  public :: lorenz
  public :: lorenz_linear

  real(dp), public :: k1 
  real(dp), public :: k2
  real(dp), public :: Q
  real(dp), public :: A
  real(dp), public :: w
  real(dp), public :: rr
  real(dp), public :: ss
  real(dp), public :: bb
  real(dp), public :: u1,u2,u3

contains

  subroutine fsub(t,u,up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)     
    real(dp), intent(out) :: up(:)       

    up(1) = u(2)
    up(2) = k1*u(1) + k2*u(2)

  end subroutine fsub

  function f1(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = k1*u(1) 

  end function f1


  function ff(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = u(2)
    up(2) = k1*u(1) + k2*u(2)

  end function ff

  function pendolo(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = u(2)
    up(2) = -sin(u(1)) - u(2)/Q + A*cos(w*t)

  end function pendolo

  function lpendolo_auto(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = u(2)
    up(2) = -cos(u(1))*u(1) - u(2)/Q - A*sin(u(3))*u(3)
    up(3) = w

  end function lpendolo_auto

  function pendolo_auto(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = u(2)
    up(2) = -sin(u(1)) - u(2)/Q + A*cos(u(3))
    up(3) = w

  end function pendolo_auto

  function harmonic(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = u(2)
    up(2) = -u(1) - u(2)/Q + A*cos(w*t)

  end function harmonic

  function lorenz(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = ss*(u(2) - u(1))
    up(2) = u(1)*( rr - u(3) ) - u(2)
    up(3) = u(1)*u(2) - bb*u(3)

  end function lorenz

  function lorenz_linear(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = ss*(u(2) - u(1))
    up(2) = (rr-u3)*u(1) -u(2) - u1*u(3)
    up(3) = u(1)*u2 +u1*u(2)- bb*u(3)

  end function lorenz_linear


end module functions


