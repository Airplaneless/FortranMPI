program main
implicit none

    interface
        subroutine VCONCAT(A, v, bottom)
            integer, allocatable :: A(:,:), B_lin(:), B(:,:)
            integer, allocatable :: v(:)
            logical :: bottom
        end subroutine VCONCAT
    end interface

    integer, allocatable :: A(:,:), B(:,:), c(:)
    integer :: i,j

    allocate(A(5,5))
    do i = 1, 5
        do j = 1, 5
            A(i,j) = 10*i+j
        enddo
    enddo
    
    allocate(c(5))
    do i = 1, 5
        c(i) = i*10
    enddo

    write(6,*)'A: '
    do i=lbound(A,1),ubound(A,1)
        write(6,*)(A(i,j), j=lbound(A,2),ubound(A,2))
    enddo

    write(6,*)'c = ',c

    call VCONCAT(A, c, .false.)
    call VCONCAT(A, c, .true.)

    write(6,*)'new A: '
    write(6,*)'bounds: ',lbound(A),ubound(A)
    do i=lbound(A,1),ubound(A,1)
        write(6,*)(A(i,j), j=lbound(A,2),ubound(A,2))
    enddo

    allocate(B(3,5))
    B = A(0:2, :)
    write(6,*)'B: '
    do i=lbound(B,1),ubound(B,1)
        write(6,*)(A(i,j), j=lbound(B,2),ubound(B,2))
    enddo

end program main

    subroutine VCONCAT(A, v, bottom)
        integer, allocatable :: A(:,:), B_lin(:), B(:,:)
        integer, allocatable :: v(:)
        logical :: bottom
        
        allocate(B_lin(SIZE(A,1)*(SIZE(A,2)+1)))
        if (bottom) then
            allocate(B(lbound(A,1):ubound(A,1)+1, lbound(A,2):ubound(A,2)))
            B_lin = [RESHAPE(TRANSPOSE(A), [SIZE(A,1)*SIZE(A,2)]), v]
        else
            allocate(B(lbound(A,1)-1:ubound(A,1), lbound(A,2):ubound(A,2)))
            B_lin = [v, RESHAPE(TRANSPOSE(A), [SIZE(A,1)*SIZE(A,2)])]
        endif
        B = TRANSPOSE(RESHAPE(B_lin, [SIZE(A,2), SIZE(A,1)+1]))
        call MOVE_ALLOC(B, A)
    
    end subroutine VCONCAT