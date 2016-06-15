module solver_class
  
  use bilayer_hamiltonian_class
  use linalg

  !The Hamiltonian class - returns the diagonal and off diagonal
  !blocks of the Hamiltonian

  implicit none
  private
  real :: pi = 3.1415926535897931d0 ! Class-wide private constant
  integer,parameter :: dp = selected_real_kind(15,307)
  integer :: dim_fix
  integer :: system_size
  complex(kind=8), dimension(:,:,:), allocatable :: dLeft,cLeft,dRight,cRight

  !Construct the Solver class
  type, public :: Solver
    type(Hamiltonian) :: Ham
    integer :: num_steps = 1500
    integer :: Nsize = 8192
    integer :: dim = 2
     real(kind=8) :: min_inv_field = 30.d0
     real(kind=8) :: max_inv_field = 6000.d0
     logical :: verbose = .true.
   contains
     procedure :: run => run_iter
!     procedure :: dos => get_dos
  end type Solver



contains
 
  subroutine initialize(this)
    class(Solver), intent(in) :: this
 
    system_size = this%Nsize
    dim_fix = this%dim

    allocate(dLeft(system_size,dim_fix,dim_fix),cLeft(system_size,dim_fix,dim_fix))
    allocate(dRight(system_size,dim_fix,dim_fix),cRight(system_size,dim_fix,dim_fix))

  end subroutine initialize



  subroutine run_iter(this)
    class(Solver), intent(in) :: this
    integer :: b
    real(kind=8) :: dos, dfield, bfield, invbfield

    


    
    !system_size = this%Nsize
    do b=1,this%num_steps
      call initialize(this)
      dos = 0.d0
      dfield = (this%max_inv_field - this%min_inv_field)/(1.d0*this%num_steps)
      bfield = 1.d0/(1.d0*this%min_inv_field + 1.d0*(b-1)*dfield)
      invbfield = 1.d0/bfield


      dos = get_dos(this,bfield)

    

      deallocate(dLeft,cLeft,dRight,cRight)

    enddo

    

  end subroutine run_iter

    
  



  function get_dos(this, bfield) result(dos2)
    !Gets a diagonal block of the Hamiltonian
    class(Solver), intent(in) :: this
    real(kind=8), intent(in) :: bfield
    real(kind=8) :: dos2
    complex(kind=8), dimension(dim_fix,dim_fix) :: block

    dLeft(:,:,:) = dcmplx(0.d0,0.d0)
    cLeft(:,:,:) = dcmplx(0.d0,0.d0)
    dRight(:,:,:) = dcmplx(0.d0,0.d0)
    cRight(:,:,:) = dcmplx(0.d0,0.d0)


    call downsweep(this, bfield)
    call upsweep(this, bfield)
    dos2 = calc_dos(this, bfield)
   
  end function get_dos



  subroutine downsweep(this, bfield)
    class(Solver), intent(in) :: this
    real(kind=8), intent(in) :: bfield
    complex(kind=8), dimension(dim_fix,dim_fix) :: extra_term
    complex(kind=8), dimension(dim_fix,dim_fix) :: dinv
    integer :: i

    extra_term(:,:)  = dcmplx(0.d0,0.d0)
    dinv(:,:)  = dcmplx(0.d0,0.d0)
    dLeft(:,:,:) = dcmplx(0.d0,0.d0)
    cLeft(:,:,:) = dcmplx(0.d0,0.d0)

    do i=1,system_size-1
      if (i.eq.1) then
        dLeft(i,:,:) = this%Ham%diag(i,bfield,0.d0,0.d0)
      else
        extra_term = mymatmul(cLeft(i-1,:,:),this%Ham%rBlock(i-1,bfield,0.d0,0.d0),dim_fix)
        dLeft(i,:,:) = this%Ham%diag(i,bfield,0.d0,0.d0) + extra_term
      endif

      dinv = findinv(dLeft(i,:,:), dim_fix)
      cLeft(i,:,:) = dcmplx(-1.d0,0.d0)*(mymatmul(this%Ham%lBlock(i+1,bfield,0.d0,0.d0),dinv,dim_fix))
    enddo

    i = system_size
    dLeft(i,:,:) = this%Ham%diag(i, bfield,0.d0,0.d0) &
                   + mymatmul(cLeft(i-1,:,:),this%Ham%rBlock(i-1,bfield,0.d0,0.d0),dim_fix)

  end subroutine downsweep


  subroutine upsweep(this, bfield)
    class(Solver), intent(in) :: this
    real(kind=8), intent(in) :: bfield
    complex(kind=8), dimension(dim_fix,dim_fix) :: extra_term
    complex(kind=8), dimension(dim_fix,dim_fix) :: dinv
    integer :: i,j

    extra_term(:,:)  = dcmplx(0.d0,0.d0)
    dinv(:,:)  = dcmplx(0.d0,0.d0)
    dRight(:,:,:) = dcmplx(0.d0,0.d0)
    cRight(:,:,:) = dcmplx(0.d0,0.d0)

    do j=0,system_size-2
      i = system_size - j
      if (i.eq.system_size) then
        dRight(i,:,:) = this%Ham%diag(i,bfield,0.d0,0.d0)
      else
        extra_term = mymatmul(cRight(i+1,:,:),this%Ham%lBlock(i+1,bfield,0.d0,0.d0),dim_fix)
        dRight(i,:,:) = this%Ham%diag(i,bfield,0.d0,0.d0) + extra_term
      endif

      dinv = findinv(dRight(i,:,:), dim_fix)
      cRight(i,:,:) = dcmplx(-1.d0,0.d0)*(mymatmul(this%Ham%rBlock(i-1,bfield,0.d0,0.d0),dinv,dim_fix))
    enddo

    i = 1
    dRight(i,:,:) = this%Ham%diag(i, bfield,0.d0,0.d0) &
                   + mymatmul(cRight(i+1,:,:),this%Ham%lBlock(i+1,bfield,0.d0,0.d0),dim_fix)

  end subroutine upsweep




  function calc_dos(this, bfield) result(dos1)
    class(Solver), intent(in) :: this
    real(kind=8), intent(in) :: bfield
    real(kind=8) :: dos1
    complex(kind=8) :: trace

    complex(kind=8), dimension(dim_fix,dim_fix) :: ginv, g
    integer :: i,j

    dos1 = 0.d0

    do i=1,system_size
      ginv = dcmplx(-1.d0,0.d0)*this%Ham%diag(i, bfield,0.d0,0.d0) + dLeft(i,:,:) + dRight(i,:,:)
      g = findinv(ginv, dim_fix)

      trace = dcmplx(0.d0,0.d0)

      do j=1,dim_fix
        trace = trace + g(j,j)
      enddo

      dos1 = dos1 - ( 1.d0/ (1.d0*pi * system_size) )  * dimag(trace)
      !write(6,*) i,dos1
    enddo

    if (this%verbose) then
      write(6,"(3F30.15,1I10)"), bfield, 1.d0/bfield, dos1, system_size
    endif

  end function calc_dos

   
  
end module solver_class


