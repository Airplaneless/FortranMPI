program main
    implicit none
    include 'mpif.h'

    interface
        subroutine READ_FORMULATION(B,g)
            real, allocatable, intent(out) :: B(:,:)
            real, allocatable, intent(out) :: g(:,:)
        end subroutine
        subroutine EVAL_STEP(X0, X1, g_sendcounts, g_part, B_part, dims)
            real, allocatable,intent(in)  :: X0(:)
            integer, intent(in) :: dims
            integer, allocatable,intent(in) :: g_sendcounts(:)
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

    integer :: myid,nproc,stat(MPI_STATUS_SIZE),ierr
    integer :: B_size, dims, steps
    real :: diff
    real, allocatable :: B(:,:), g(:,:), B_part(:,:), g_part(:,:), X(:), X_next(:), X_old(:) 
    integer, allocatable :: Bsendcounts(:), Bdispls(:), g_sendcounts(:), g_displs(:)
    integer :: i,j

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr) 
    
    if(myid == 0) then 
        call READ_FORMULATION(B, g)
        dims = SIZE(B,1)
        B = RESHAPE(TRANSPOSE(B), (/dims*dims, 1/))
        call SPLIT_MATRIX(B, dims, nproc, Bsendcounts, Bdispls)
        call SPLIT_MATRIX(g, 1, nproc, g_sendcounts, g_displs)
    endif

    !call PRINT_2DMATRIX(B, myid)

    call MPI_BCAST(dims, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if(myid /= 0) then 
        allocate(g(dims,1))
        allocate(Bsendcounts(nproc))
        allocate(g_sendcounts(nproc))
        allocate(g_displs(nproc))
    endif

    call MPI_BCAST(Bsendcounts, SIZE(Bsendcounts), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(g_sendcounts, SIZE(g_sendcounts), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(g_displs, SIZE(g_displs), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    allocate(B_part(Bsendcounts(myid + 1),1))

    if(myid==0) then
        write(6,*)'INFO: B matrix',(B(i,1), i = 1,SIZE(B))
        write(6,*)'INFO: g matrix:',(g(i,1), i = 1,SIZE(g,1))
    endif

    call MPI_SCATTERV(B, Bsendcounts, Bdispls, MPI_REAL, B_part, Bsendcounts(myid+1), MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
    allocate(g_part(g_sendcounts(myid + 1), 1))
    call MPI_SCATTERV(g, g_sendcounts, g_displs, MPI_REAL, g_part, g_sendcounts(myid+1), MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
    
    allocate(X(dims))
    allocate(X_old(dims))
    do i=1,dims
        X(i) = 0.0
    enddo
    write(6,*)'dims:',dims
    B_part = TRANSPOSE(RESHAPE(TRANSPOSE(B_part), (/ dims, Bsendcounts(myid+1)/dims /)))

    call PRINT_2DMATRIX(B_part, myid)
    write(6,*)'INFO: id: ',myid,'g_part: ',(g_part(i,1), i=1,SIZE(g_part))

    allocate(X_next(g_sendcounts(myid + 1)))
    diff = 1000
    steps = 0
    do while(abs(diff) > 0.0001)
        X_old = X
        steps = steps + 1
        do i=1,SIZE(B_part,1)
            X_next(i) = g_part(i,1)
            do j=1,SIZE(B_part,2)
                X_next(i) = X_next(i) + B_part(i,j)*X(j)
            enddo
        enddo
        !X = 0
        call MPI_ALLGATHERV(X_next, g_sendcounts(myid+1), MPI_REAL, X, g_sendcounts, g_displs, MPI_REAL, MPI_COMM_WORLD, ierr)

        diff = SUM(X - X_old)
        write(6,*)'DEBUG: Step :',steps,'x = ',(X(i), i=1,SIZE(X))
        write(6,*)'DEBUG: diff: ',diff
    enddo

    if(myid == 0) then
        write(6,*)'INFO: Result:',(X(i), i=1,SIZE(X))
    endif
    
    call MPI_FINALIZE(ierr)

end program main

subroutine EVAL_STEP(X0, X1, g_sendcounts, g_part, B_part, dims)
    real, allocatable,intent(in)  :: X0(:)
    integer, intent(in) :: dims
    integer, allocatable,intent(in) :: g_sendcounts(:)
    real, allocatable,intent(in) :: B_part(:,:), g_part(:,:)
    real, allocatable,intent(out) :: X1(:)

    do i=1,g_sendcounts(myid + 1)
        X1(i) = g_part(i,1)
        do j=1,dims
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