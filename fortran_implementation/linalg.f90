 module linalg

  Implicit none
  private
  integer,parameter :: dp = selected_real_kind(15,307)
  complex(kind=8),allocatable,dimension(:,:)::A
  complex(kind=8),allocatable,dimension(:,:)::res
  complex(kind=8),allocatable,dimension(:)::WORK
  integer,allocatable,dimension(:)::IPIV
  integer i,j,info,error, M
  public :: findinv

  contains

  function findinv(Ain, dim) result(res)
    integer :: dim 
    complex(kind=8), dimension(dim,dim)::Ain
    complex(kind=8),dimension(dim,dim)::res
    !definition of the test matrix A
    M= dim

    allocate(A(M,M),WORK(M),IPIV(M),stat=error)

    A = Ain

    
    if (error.ne.0)then
      print *,"error:not enough memory"
      stop
    end if
     
    call ZGETRF(M,M,A,M,IPIV,info)
    if(info .eq. 0) then
     ! write(*,*)"succeded"
    else
     write(*,*)"failed"
    end if

    call ZGETRI(M,A,M,IPIV,WORK,M,info)
    if(info .eq. 0) then
     ! write(*,*)"succeded"
    else
     write(*,*)"failed"
    end if

    res = A

    deallocate(A,IPIV,WORK)
    
    

 end function findinv

end module linalg
