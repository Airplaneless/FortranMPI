program main
    implicit none

    interface
        subroutine readFormulation(B,g)
            real, allocatable, intent(out) :: B(:,:)
            real, allocatable, intent(out) :: g(:)
        end subroutine
    end interface

    integer i,j
    real, allocatable :: B(:,:), g(:)

    call readFormulation(B,g)
    do i=1,SIZE(B, 1)
        write(6,*)(B(i,j), j=1,SIZE(B, 1))
    enddo

end program main

subroutine readFormulation(B, g)
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
    
end subroutine readFormulation