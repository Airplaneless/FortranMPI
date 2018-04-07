program main
    implicit none

    interface
        subroutine READ_FORMULATION(B,g)
            real, allocatable, intent(out) :: B(:,:)
            real, allocatable, intent(out) :: g(:)
        end subroutine
    end interface

    integer i,j
    real, allocatable :: B(:,:), g(:), X(:)

    call READ_FORMULATION(B,g)
    allocate(X(SIZE(g)))

    X = B*g

    do i=1,SIZE(B,2)
        do j=1,SIZE(B,1)
            X(i) = X(i) + B(i,j)*g(j)
        enddo
    enddo

    do i=1,SIZE(X)
        write(6,*)X(i)
    enddo

end program main

subroutine READ_FORMULATION(B, g)
    real, allocatable, intent(out) :: B(:,:)
    real, allocatable, intent(out) :: g(:)

    integer :: n, i, j
    real :: div

    open(10,file='Btest',form='formatted',status='old')
    read(10,*)n
    allocate (B(n,n))
    allocate (g(n))
    
    do i=1,n
        read(10,*)(B(j,i), j=1,n),g(i) 
    enddo
    
    do i=1,n
        div = B(i,i)
        g(i) = g(i) / div
        B(i,i) = 0
        do j=2,n
            B(i,j) = -B(i,j) / div
        enddo
    enddo
    
end subroutine READ_FORMULATION