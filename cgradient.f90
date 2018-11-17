subroutine show_matrix(B,m,k)
    implicit none
    integer m,k,j
    real(kind=16), dimension(m,k) :: B
	do j=1,size(B,1) 
      write(*,'(20G12.4)') B(j,:) 
    end do      
    write(*,*) 
end subroutine show_matrix

program cgradient
      implicit none
      integer i, it,n
      real*16,allocatable ::A(:,:),b(:,:)    
      real*16 ::  error, lambda,q
      real*16,allocatable :: v(:,:),v0(:,:), r(:,:),r0(:,:),p(:,:)
      real*16, allocatable :: pt1(:,:)

      !define data 
        !A*x = b
        n = 3
        
       allocate(A(n,n))
       allocate(r(n,1))
       allocate(r0(n,1))    
       allocate(v(n,1))
       allocate(v0(n,1))
       allocate(p(n,1))
       allocate(pt1(n,1))

        A(1,1) = 10
        A(1,2) = 1
        A(1,3) = 0
        
        A(2,1) = 1
        A(2,2) = 10
        A(2,3) = 1
        
        A(3,1) = 0
        A(3,2) = 1
        A(3,3) = 10
        
        r0(1,1) = 11
        r0(2,1) = 11
        r0(3,1) = 1
        
        write(*,*) ' type  it'
        read (*,*) it
        !first guess: x = 0
        r0 = -r0
        v(1:n,1) = 0
        
        lambda = sum(r0*r0)/sum(matmul(A,r0)*r0)
        v =v -lambda*r0
        r = r0 -lambda*matmul(A,r0)
        p = -r
        error = norm2(v - v0)/norm2(v)

        !write(*,*) lambda
        
       ! write (*,*) 'solution'
        !call show_matrix(v,n,1)
        
        !write (*,*) 'error' 
        !write (*,'(20G12.4)') error
        
         do i=2,it
           v0 = v
           !necessita de um ponteiro para a variável r
           pt1=r0
           r0 = r
           r = pt1
           lambda = sum(r*r)/sum(r0*r0)
           p = -r+lambda*p
           q = sum(r*r)/sum(matmul(A,p)*p)
           v = v+q*p
           r = r0+q*matmul(A,p)
          enddo

        error = norm2(v - v0)/norm2(v)
        
        write (*,*) 'solution'
        call show_matrix(v,n,1)
        
        write (*,*) 'error' 
        write (*,'(20G12.4)') error
        deallocate(A,v,v0,r,r0,p)
endprogram cgradient
 
