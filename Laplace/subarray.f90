program main
implicit none
include 'mpif.h'

    interface
        subroutine PRINT_2DARRAY(A, thread)
            real, allocatable :: A(:,:)
            integer :: thread
        end subroutine PRINT_2DARRAY
    end interface

    integer :: ierr,myid,nproc,stat(MPI_STATUS_SIZE)
    real, allocatable :: A(:,:), B(:,:)
    integer :: n, i, j
    integer :: subarray_datatype

    integer :: sizes(2)
    integer :: subsizes(2)
    integer :: start(2)
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

    sizes = [6, 6]
    subsizes = [3, 3]
    start = [0, 0]

    if (myid == 0) then
        open(7,file='A',form='formatted',status='old')
        open(8,file='debug.log',form='formatted',status='old')
        read(7,*)n
        allocate(A(n,n))
        do i=1,n
            read(7,*)(A(i,j), j=1,n)
        enddo
        call PRINT_2DARRAY(A, 8)

        call MPI_Type_create_subarray(2, sizes, subsizes, start, MPI_ORDER_FORTRAN, MPI_REAL, subarray_datatype, ierr)
        call MPI_Type_commit(subarray_datatype, ierr)
        call MPI_Send(A, 1, subarray_datatype, 1, 1, MPI_COMM_WORLD, ierr)
    endif

    if (myid == 1) then
        allocate(B(subsizes(1),subsizes(2)))
        call MPI_Recv(B, subsizes(1)*subsizes(2), MPI_REAL, 0, 1, MPI_COMM_WORLD, stat, ierr)
        call PRINT_2DARRAY(B,6)
    endif

    call MPI_FINALIZE(ierr)

end program main

subroutine PRINT_2DARRAY(A, thread)
    real, allocatable :: A(:,:)
    integer :: thread
    do i=lbound(A,1),ubound(A,1)
        write(thread,*)(A(i,j), j=lbound(A,2),ubound(A,2))
    enddo        
end subroutine PRINT_2DARRAY

    do while(proceed)
        ! Neigbours exchange
        do i = 0, nproc-1
            if (myid == i) then
                if (id_map(3) /= -1) then
                    call MPI_Send(phi_part(id_map(3)+1,:), phi_size, MPI_REAL, i-1, i+nproc+42, MPI_COMM_WORLD, ierr)
                endif
            endif
            if (myid == i-1) then
                call MPI_Recv(lower, phi_size, MPI_REAL, i, i+nproc+42, MPI_COMM_WORLD, stat, ierr)
            endif
        enddo

        do i = 0, nproc-1
            if (myid == i) then
                if (id_map(4) /= -1) then
                    call MPI_Send(phi_part(id_map(4)-1,:), phi_size, MPI_REAL, i+1, i+nproc+84, MPI_COMM_WORLD, ierr)
                endif
            endif
            if (myid == i+1) then
                call MPI_Recv(upper, phi_size, MPI_REAL, i, i+nproc+84, MPI_COMM_WORLD, stat, ierr)
            endif
        enddo

        ! Concatanate neigbours
        if (id_map(3) == -1) then
            call VCONCAT(phi_part, lower, .true.)
        elseif (id_map(4) == -1) then
            call VCONCAT(phi_part, upper, .false.)
        else
            call VCONCAT(phi_part, upper, .false.)
            call VCONCAT(phi_part, lower, .true.)
        endif

        ! Eval relaxation
        call EVAL_RELAXATION(phi_part, new_phi_part)
        call L2_NORM(phi_part, new_phi_part, error)
        call MOVE_ALLOC(new_phi_part, phi_part)
        ! Eval error
        all_error = 0.d0
        call MPI_Reduce(error, all_error, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (myid == 0) then
            write(6,*)'error = ',all_error
            if (all_error < 1000) then
                proceed = .false.
            endif
        endif
        call MPI_BCAST(proceed, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    enddo