 module linalg

  Implicit none
  private
  integer,parameter :: dp = selected_real_kind(15,307)
  complex(kind=8),allocatable,dimension(:,:)::A
  complex(kind=8),allocatable,dimension(:,:)::res
  complex(kind=8),allocatable,dimension(:)::WORK
  integer,allocatable,dimension(:)::IPIV
  integer i,j,info,error, M
  public :: findinv, sncndn, mymatmul

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



 function mymatmul(A1,B1,dim) result(prod)
    integer :: dim
    complex(kind=8), dimension(dim,dim) :: A1,B1
    complex(kind=8), dimension(dim,dim) :: prod
    complex(kind=8) :: const

    const = dcmplx(1.d0,0.d0)

    call zgemm('N','N',dim,dim,dim,const,A1,dim,B1,dim,const,prod,dim)

  end function






  SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
  REAL(kind=8) :: cn,dn,emmc,sn,uu,CA
  PARAMETER (CA=.0003) 
  INTEGER i,ii,l
  REAL(kind=8):: a,b,c,d,emc,u,em(13),en(13) 
  LOGICAL bo
  emc=emmc
  u=uu
  if (emc.ne.0.)then
      bo=(emc.lt.0.) 
      if(bo)then
        d=1.-emc 
        emc=-emc/d 
        d=sqrt(d) 
        u=d*u
      endif
        
    a=1.
    dn=1.
    
    do i=1,13 
      l=i
      em(i)=a
      emc=sqrt(emc)
      en(i)=emc
      c=0.5*(a+emc) 
      if(abs(a-emc).le.CA*a) goto 1 
      emc=a*emc
      a=c 
    enddo 
    1 u=c*u 
    sn=sin(u)
    cn=cos(u) 
    if(sn.eq.0.) goto 2 
    a=cn/sn
    c=a*c
    do ii=l,1,-1
      b=em(ii)
      a=c*a
      c=dn*c 
      dn=(en(ii)+a)/(b+a) 
      a=c/b
    enddo
    a=1./sqrt(c**2+1.) 
    if(sn.lt.0.)then
      sn=-a 
    else
      sn=a 
    endif
    cn=c*sn
    2 if(bo)then
       a=dn
       dn=cn
       cn=a
       sn=sn/d
    endif 
  else
    cn=1./cosh(u) 
    dn=cn 
    sn=tanh(u)
  endif
  return
END



    
    

end module linalg
