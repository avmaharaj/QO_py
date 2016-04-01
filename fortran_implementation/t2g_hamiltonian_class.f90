module t2g_hamiltonian_class
  !The Hamiltonian class - returns the diagonal and off diagonal
  !blocks of the Hamiltonian


  implicit none
  private
  real :: pi = 3.1415926535897931d0 ! Class-wide private constant
  integer,parameter :: dp = selected_real_kind(15,307)
  real(kind=8), parameter :: ta = 1.d0
  integer, parameter :: dim_fix = 6

  !Construct the Hamiltonian class
  type, public :: Hamiltonian
     real(kind=8) :: spin_orbit 
     real(kind=8) :: mu
     real(kind=8) :: t2 
     real(kind=8) :: theta 
     real(kind=8) :: delta
   contains
     procedure :: dim => print_dim
     procedure :: diag => get_diagonal_block
     procedure :: rBlock => get_right_block
     procedure :: lBlock => get_left_block
  end type Hamiltonian



contains
  subroutine print_dim(this)
    class(Hamiltonian), intent(in) :: this
    print *, "Dimension of Hamiltonian is", this%spin_orbit
  end subroutine print_dim

    
  function get_diagonal_block(this, x, b, ky, kz) result(block)
    !Gets a diagonal block of the Hamiltonian
    class(Hamiltonian), intent(in) :: this
    integer, intent(in) :: x
    real(kind=8), intent(in) :: b,ky,kz

    complex(kind=8), dimension(6,6) :: block

    block(:,:) = (0.d0,0.d0)



    block(1,1) = dcmplx(2.d0*ta*cos(2.d0*pi*b*x*cos(this%theta*pi/180.d0) - ky), this%delta)
    block(1,1) = block(1,1) + dcmplx(2.d0*ta*cos(-2.d0*pi*b*x*sin(this%theta*pi/180.d0) - kz), 0.d0)
    block(1,1) = block(1,1) + dcmplx(this%mu, 0.d0)
    block(4,4) = dcmplx(2.d0*ta*cos(2.d0*pi*b*x*cos(this%theta*pi/180.d0) - ky) &
                + 2.d0*ta*cos(-2.d0*pi*b*x*sin(this%theta*pi/180.d0) - kz) &
                + this%mu, this%delta)


    block(2,2) = dcmplx(2.d0*this%t2*cos(2.d0*pi*b*x*cos(this%theta*pi/180.d0) - ky) &
                + 2.d0*ta*cos(-2.d0*pi*b*x*sin(this%theta*pi/180.d0) - kz) &
                + this%mu, this%delta)
    block(5,5) = dcmplx(2.d0*this%t2*cos(2.d0*pi*b*x*cos(this%theta*pi/180.d0) - ky) &
                + 2.d0*ta*cos(-2.d0*pi*b*x*sin(this%theta*pi/180.d0) - kz) &
                + this%mu, this%delta)


    block(3,3) = dcmplx(2.d0*ta*cos(2.d0*pi*b*x*cos(this%theta*pi/180.d0) - ky) &
                + 2.d0*this%t2*cos(-2.d0*pi*b*x*sin(this%theta*pi/180.d0) - kz) &
                + this%mu, this%delta)
    block(6,6) = dcmplx(2.d0*ta*cos(2.d0*pi*b*x*cos(this%theta*pi/180.d0) - ky) &
                + 2.d0*this%t2*cos(-2.d0*pi*b*x*sin(this%theta*pi/180.d0) - kz) &
                + this%mu, this%delta)

    !Now add the (complex) spin orbit terms
    block(1,2) = dcmplx(0.d0, this%spin_orbit)
    block(2,1) = dcmplx(0.d0, -1.d0*this%spin_orbit)

    block(1,6) = dcmplx(-1.d0*this%spin_orbit,0.d0)
    block(6,1) = dcmplx(-1.d0*this%spin_orbit,0.d0)

    block(2,6) = dcmplx(0.d0, this%spin_orbit)
    block(6,2) = dcmplx(0.d0, -1.d0*this%spin_orbit)

    block(3,4) = dcmplx(this%spin_orbit,0.d0)
    block(4,3) = dcmplx(this%spin_orbit,0.d0)

    block(3,5) = dcmplx(0.d0, -1.d0*this%spin_orbit)
    block(5,3) = dcmplx(0.d0, this%spin_orbit)

    block(4,5) = dcmplx(0.d0,-1.d0*this%spin_orbit)
    block(5,4) = dcmplx(0.d0, this%spin_orbit)

!     write(6,"(8F20.10)") x, b, block(1,1),block(2,2),block(3,3)

  end function get_diagonal_block



  function get_right_block(this, x, b, ky, kz) result(block)
    !Gets a diagonal block of the Hamiltonian
    class(Hamiltonian), intent(in) :: this
    integer, intent(in) :: x
    real(kind=8), intent(in) :: b,ky,kz

    complex(kind=8), dimension(6,6) :: block

    block(:,:) = (0.d0,0.d0)

    !Now add the (complex) spin orbit terms
    block(1,1) = dcmplx(this%t2,0.d0)
    block(4,4) = dcmplx(this%t2,0.d0)

    block(2,2) = dcmplx(ta,0.d0)
    block(5,5) = dcmplx(ta,0.d0)

    block(3,3) = dcmplx(ta,0.d0)
    block(6,6) = dcmplx(ta,0.d0)

  end function get_right_block


  function get_left_block(this, x, b, ky, kz) result(block)
    !Gets a diagonal block of the Hamiltonian
    class(Hamiltonian), intent(in) :: this
    integer, intent(in) :: x
    real(kind=8), intent(in) :: b,ky,kz

    complex(kind=8), dimension(6,6) :: block

    block(:,:) = (0.d0,0.d0)

    !Now add the (complex) spin orbit terms
    block(1,1) = dcmplx(this%t2,0.d0)
    block(4,4) = dcmplx(this%t2,0.d0)

    block(2,2) = dcmplx(ta,0.d0)
    block(5,5) = dcmplx(ta,0.d0)

    block(3,3) = dcmplx(ta,0.d0)
    block(6,6) = dcmplx(ta,0.d0)

  end function get_left_block

end module t2g_hamiltonian_class


