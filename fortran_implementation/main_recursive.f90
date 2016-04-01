program recursive_dos
  use solver_class
  use t2g_hamiltonian_class
  use get_command_line_input
  implicit none

  type(InputParameters) :: params !Declare an instance of the inputvariables
  type(Hamiltonian) :: h1     ! Declare an instance of the Hamiltonian
  type(Solver) :: solver_instance       ! Declare the solver instance
  
  !We will use double precision throughout
  integer,parameter :: dp = selected_real_kind(15,307)
  !Variables we need from the command line
  integer :: length, numsteps, polarAngle


  !1. Grab some input from the command line
  params = InputParameters()  
  call params%grabinput()


  !2. Construct the Hamiltonian
  h1 = Hamiltonian(spin_orbit=params%spinOrbit, mu=params%mu, t2=params%t2, &
                    theta=params%theta, delta=params%delta )       

  !3. Now construct the solver instance
  solver_instance = Solver(Ham=h1, Nsize=params%length, num_steps=params%nsteps, &
                          min_inv_field = params%minInvField, &
                          max_inv_field = params%maxInvField)

  !4. And now, run the solver!!
  call solver_instance%run()



 


end program recursive_dos


 