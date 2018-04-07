program main
    implicit none

    interface
        subroutine READ_FORMULATION(B,g)
            real, allocatable, intent(out) :: B(:,:)
            real, allocatable, intent(out) :: g(:,:)
        end subroutine
        subroutine EVAL_STEP(X0, X1, Bsendcounts, g_sendcounts, g_part, B_part, B_size)
            real, allocatable,intent(in)  :: X0(:)
            integer, intent(in) :: B_size
            integer, allocatable,intent(in) :: Bsendcounts(:), g_sendcounts(:)
            real, allocatable,intent(in) :: B_part(:,:), g_part(:,:)
            real, allocatable,intent(out) :: X1(:)
        end subroutine
        subroutine SPLIT_MATRIX(A, dims, nproc, sendcounts, displs)
            real, allocatable, intent(in)  :: A(:,:)
            integer, intent(in)  :: nproc
            integer, intent(in)  :: dims
            integer, allocatable, intent(out) :: sendcounts(:)
            integer, allocatable, intent(out) :: displs(:)
        end subroutine
        subroutine PRINT_2DMATRIX(A, id)
            integer,intent(in) :: id
            real,allocatable,intent(in)  :: A(:,:)
        end subroutine
        subroutine PRINT_1DMATRIX(A, id)
            integer,intent(in) :: id
            real,allocatable,intent(in)  :: A(:)
        end subroutine
        subroutine PRINT_1DMATRIX_INT(A, id)
            integer,intent(in) :: id
            integer,allocatable,intent(in)  :: A(:)
        end subroutine
    end interface   


    integer i,j
    real, allocatable :: B(:,:), g(:,:), X(:)
    integer, allocatable :: sendcounts(:), displs(:)

    call READ_FORMULATION(B,g)
    B = RESHAPE(B, (/SIZE(B,1)*SIZE(B,2), 1/))
    call SPLIT_MATRIX(B, 6, 4, sendcounts, displs)
    !call SPLIT_MATRIX(g, 1, 4, sendcounts, displs)

    write(6,*)'B:'
    write(6,*)(B(i,1), i=1,SIZE(B))

    !write(6,*)'g:'
    !write(6,*)(g(i,1), i=1,SIZE(g))

    write(6,*)'sendcounts: '
    do i=1, SIZE(sendcounts)
        write(6,*)sendcounts(i)
    enddo

    write(6,*)'displs: '
    do i=1, SIZE(displs)
        write(6,*)displs(i)
    enddo

end program main


    subroutine EVAL_STEP(X0, X1, Bsendcounts, g_sendcounts, g_part, B_part, B_size)
        real, allocatable,intent(in)  :: X0(:)
        integer, intent(in) :: B_size
        integer, allocatable,intent(in) :: Bsendcounts(:), g_sendcounts(:)
        real, allocatable,intent(in) :: B_part(:,:), g_part(:,:)
        real, allocatable,intent(out) :: X1(:)
    
        allocate(X1(Bsendcounts(myid + 1)))
        do i=1,g_sendcounts(myid + 1)
            X1(i) = g_part(i,1)
            do j=1,B_size
                X1(i) = X1(i) + B_part(i,j)*X0(j)
            enddo
        enddo    
    end subroutine EVAL_STEP
    
    subroutine SPLIT_MATRIX(A, dims, nproc, sendcounts, displs)
        real, allocatable, intent(in)  :: A(:,:)
        integer, intent(in)  :: nproc
        integer, intent(in)  :: dims
        integer, allocatable, intent(out) :: sendcounts(:)
        integer, allocatable, intent(out) :: displs(:)
    
        integer :: i
        allocate(sendcounts(nproc))
        allocate(displs(nproc))
    
        sendcounts = SIZE(A,1) / (nproc * dims)
        do i=1,mod(SIZE(A,1),nproc*dims)/dims
            sendcounts(i) = sendcounts(i) + 1
        enddo

        sendcounts = sendcounts * dims
        
        displs(1) = 0
        do i=2,nproc
            displs(i) = displs(i-1) + sendcounts(i)*SIZE(A,2)
        enddo
            
    end subroutine SPLIT_MATRIX
    
    subroutine READ_FORMULATION(B, g)
        real, allocatable, intent(out) :: B(:,:)
        real, allocatable, intent(out) :: g(:,:)
    
        integer :: n, i, j
        real :: div
    
        open(10,file='Btest',form='formatted',status='old')
        read(10,*)n
        allocate (B(n,n))
        allocate (g(n,1))
        
        do i=1,n
            read(10,*)(B(i,j), j=1,n),g(i,1) 
        enddo
        
        do i=1,n
            do j=1,n
                if(i /= j) then
                    B(i,j) = -B(i,j) / B(i,i)
                else
                    B(i,i) = B(i,i)
                endif
            enddo
            g(i,1) = g(i,1) / B(i,i)
            B(i,i) = 0
        enddo
        
    end subroutine READ_FORMULATION
    
    subroutine PRINT_1DMATRIX(A, id)
        integer,intent(in) :: id
        real,allocatable,intent(in)  :: A(:)
        write(6,*)''
        write(6,*)'id:',id,'DEBUG 1D MATRIX' 
        do i=1,SIZE(A)
            write(6,*)A(i)
        enddo
    end subroutine PRINT_1DMATRIX
    
    subroutine PRINT_1DMATRIX_INT(A, id)
        integer,intent(in) :: id
        integer,allocatable,intent(in)  :: A(:)
        write(6,*)''
        write(6,*)'id:',id,'DEBUG 1D MATRIX' 
        do i=1,SIZE(A)
            write(6,*)A(i)
        enddo
    end subroutine PRINT_1DMATRIX_INT
    
    subroutine PRINT_2DMATRIX(A, id)
        integer,intent(in) :: id
        real,allocatable,intent(in)  :: A(:,:)
        write(6,*)''
        write(6,*)'id:',id,'DEBUG 2D MATRIX' 
        do i=1,SIZE(A,1)
            write(6,*)(A(i,j), j=1,SIZE(A,2))
        enddo
    end subroutine PRINT_2DMATRIX