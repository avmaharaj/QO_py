module get_command_line_input
  implicit none
  private
  integer,parameter :: dp = selected_real_kind(15,307)
  character(32), dimension(8) :: inchars



  type, public :: InputParameters
     integer :: length = 16394
     integer :: nsteps = 2000
     real(kind=8) :: theta = 0.d0
     real(kind=8) :: delta = 0.001d0
     real(kind=8) :: t2 = 0.5d0
     real(kind=8) :: mu = -4.96d0
     real(kind=8) :: spinOrbit = 0.005d0
     real(kind=8) :: minInvField = 30.d0
     real(kind=8) :: maxInvField = 6000.d0
  contains
      procedure :: grabinput => get_input1
  end type InputParameters

contains

  subroutine get_input1(this)
    class(InputParameters), intent(in) :: this
    integer :: i
    character(len=32) :: arg
            
    i = 0
    do
      call get_command_argument(i, arg)
      if (len_trim(arg) == 0) EXIT
        inchars(i+1) = trim(arg) 
      i = i+1
    enddo

    call initialize_vars(this,i)

  end subroutine get_input1



  subroutine initialize_vars(this,imax)
    class(InputParameters) :: this
    integer, intent(in) :: imax
    integer :: i,intvalue
    real(kind=8) :: realvalue
    character(len=32) :: argtot,name,val

    do i=2,imax
      call split_string(inchars(i),name,val,'=')
      
      select case (name)
        case ('length')
          read(val,'(i10)') intvalue
          this%length = intvalue
        case ('mu')
          read(val, '(F7.5)') realvalue
          this%mu = realvalue
        case ('theta')
          read(val, '(F10.5)') realvalue
          this%theta = realvalue
        case ('nsteps')
          read(val,'(i10)') intvalue
          this%nsteps = intvalue
        case ('delta')
          read(val,'(F10.5)') realvalue
          this%delta = realvalue
        case ('t2')
          read(val, '(F4.2)') realvalue
          this%t2 = realvalue
        case ('spinOrbit')
          read(val, '(F7.5)') realvalue
          this%spinOrbit = realvalue
        case ('minInvField')
          read(val,'(F12.5)') realvalue
          this%minInvField = realvalue
        case ('maxInvField')
          read(val,'(F12.5)') realvalue
          this%maxInvField = realvalue
        case default
          write(6,*) trim(inchars(i))," is an invalid Input argument!!"
          write(6,*) "Input must be of the form arg=xxx"
          write(6,*) "With arg= length, theta, nsteps, delta,t2, spinOrbit minInvField, maxInvField"
          write(6,*) "Exiting the program ....."
          call exit 
        end select
    enddo

  end subroutine initialize_vars




  subroutine split_string(instring, string1, string2, delim)
  ! split a string into 2 either side of a delimiter token
  ! taken from https://gist.github.com/ponderomotion/
    character(32) :: instring
    character :: delim
    character(32),INTENT(OUT):: string1,string2
    integer :: index

    instring = TRIM(instring)

    index = SCAN(instring,delim)
    string1 = instring(1:index-1)
    string2 = instring(index+1:)

  end subroutine split_string

end module get_command_line_input
