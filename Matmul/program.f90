program main
implicit none
include 'mpif.h'

    interface
        subroutine GET_FORMULATION(A, B, a_x_size, a_y_size, b_x_size, b_y_size)
            real, allocatable :: A(:,:), B(:,:)
            integer :: b_x_size, b_y_size, a_x_size, a_y_size
            integer :: i,j
        end subroutine GET_FORMULATION
        subroutine PRINT_2DARRAY(A, thread)
            real, allocatable :: A(:,:)
            integer :: thread
            integer :: i, j
        end subroutine PRINT_2DARRAY
        subroutine SLICE_MATRIX(A, num, map)
            real, allocatable :: A(:,:)
            integer :: num, i
            integer, allocatable :: map(:,:)
        end subroutine SLICE_MATRIX
    end interface

    integer :: ierr,rank,size,stat(MPI_STATUS_SIZE), subarray_datatype
    real, allocatable :: a(:,:), b(:,:), c(:,:), a_slice(:,:), c_slice(:,:), c_lin(:), c_slice_lin(:)
    integer, allocatable :: a_map(:,:), a_map_lin(:), rank_map(:), c_displs(:), c_counts(:)
    integer :: sizes(2), subsizes(2), start(2)
    integer :: nax, nay, nbx, nby
    integer :: i,j,k

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

    if (rank == 0) then
        call GET_FORMULATION(a,b,nax,nay,nbx,nby)
        call SLICE_MATRIX(a, size, a_map)
        a_map_lin = RESHAPE(a_map, [size*2])
        a_slice = a(1:a_map(1,1),1:nay)
        allocate(c(nax, nby))
        allocate(c_lin(nay*nbx))
    endif
    ! Send map
    allocate(rank_map(2))
    call MPI_Scatter(a_map_lin, 2, MPI_INTEGER, rank_map, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! Send B matrix to all process
    call MPI_BCAST(nbx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nby, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (rank/=0) then
        allocate(b(nbx, nby))
        allocate(a_slice(rank_map(1),nbx))
    endif
    call MPI_BCAST(b, nbx*nby, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    ! Send rows of A matrix   
    do i = 1, size-1
        if (rank == 0) then 
            sizes = [nby, nbx]
            subsizes = [a_map(1,i+1), nbx]
            start = [a_map(2,i+1)-1, 0]
            call MPI_Type_create_subarray(2, sizes, subsizes, start, MPI_ORDER_FORTRAN, MPI_REAL, subarray_datatype, ierr)
            call MPI_Type_commit(subarray_datatype, ierr)
            call MPI_Send(a, 1, subarray_datatype, i, i, MPI_COMM_WORLD, ierr)
            call MPI_Type_Free(subarray_datatype, ierr)
        endif
        if (rank == i) then
            call MPI_Recv(a_slice, rank_map(1)*nbx, MPI_REAL, 0, i, MPI_COMM_WORLD, stat, ierr)
        endif
    enddo
    ! Eval rows of C = A x B
    allocate(c_slice(rank_map(1), nby))
    allocate(c_slice_lin(rank_map(1)*nby))
    c_slice = 0.d0
    do i=lbound(c_slice,1), ubound(c_slice,1)
        do j=lbound(c_slice,2), ubound(c_slice,2)
            do k=lbound(a_slice,2),ubound(a_slice,2)
                !write(6,*)'i',i,'j',j,'k',k
                c_slice(i,j) = c_slice(i,j) + a_slice(i,k)*b(k,j)
            enddo
        enddo
    enddo
    c_slice_lin = RESHAPE(TRANSPOSE(c_slice), [rank_map(1)*nby])
    ! Collect C rows from proc on 0 rank
    if (rank==0) then
        allocate(c_displs(size))
        allocate(c_counts(size))
        c_displs = 0
        c_counts = a_map(1,:)*nby
        do i=2,size
            c_displs(i) = c_displs(i-1) + nby*a_map(1,i)
        enddo
    endif
    call MPI_Gatherv(c_slice_lin, rank_map(1)*nby, MPI_REAL, c_lin, c_counts, c_displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    ! Reshape C array and write in file
    if (rank==0) then
        c = TRANSPOSE(RESHAPE(c_lin, [nax, nby]))
        open(unit = 9, file = 'C', action="write")
        call PRINT_2DARRAY(c, 9)
    endif

    call MPI_FINALIZE(ierr)

end program main

subroutine PRINT_2DARRAY(A, thread)
    real, allocatable :: A(:,:)
    integer :: thread
    integer :: i, j

    do i=lbound(A,1) ,ubound(A,1)
        write(thread,*)(A(i,j), j=lbound(A,2),ubound(A,2))
    enddo
end subroutine PRINT_2DARRAY


subroutine GET_FORMULATION(A, B, a_x_size, a_y_size, b_x_size, b_y_size)
    real, allocatable :: A(:,:), B(:,:)
    integer :: b_x_size, b_y_size, a_x_size, a_y_size
    integer :: i,j

    open(unit = 10, file = 'A')
    read(10, *)a_x_size
    read(10, *)a_y_size
    allocate(A(a_x_size, a_y_size))
    do i=1,a_x_size
        read(10,*)(A(i,j), j=1,a_y_size)
    enddo
    close(10)
    open(unit = 10, file = 'B')
    read(10, *)b_x_size
    read(10, *)b_y_size
    allocate(B(b_x_size, b_y_size))
    do i=1,b_x_size
        read(10,*)(B(i,j), j=1,b_y_size)
    enddo
    close(10)
end subroutine GET_FORMULATION

subroutine SLICE_MATRIX(A, num, map)
    real, allocatable :: A(:,:)
    integer :: num, i
    integer, allocatable :: map(:,:)
    ! Build map array for cutting matrix A(m x n) on A_i(m x k_i), where i = 1
    ! map = [
    !       [*k_i*]
    !       [*indexes of start A_i matrix in A*]
    !]
    allocate(map(2,num))
    map(1,:) = SIZE(A,1) / num
    do i = 1, MOD(SIZE(A,1), num)
        map(1,i) = map(1,i) + 1
    enddo
    map(2,1) = 1
    do i = 2, num
        map(2,i) =  map(2,i-1) + map(1,i-1)
    enddo 
end subroutine SLICE_MATRIX
    