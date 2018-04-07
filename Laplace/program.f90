program main
    implicit none
    include 'mpif.h'
    
        interface
            subroutine PRINT_2DARRAY(A, thread)
                real, allocatable :: A(:,:)
                integer :: thread
                integer :: i, j
            end subroutine PRINT_2DARRAY
            subroutine GET_FORMULATION(A, n)
                real, allocatable :: A(:,:)
                integer :: n    
                real :: lphi, rphi, uphi, bphi
            end subroutine GET_FORMULATION
            subroutine SLICE_MATRIX(A, num, map)
                real, allocatable :: A(:,:)
                integer :: num, i
                integer, allocatable :: map(:,:)
            end subroutine SLICE_MATRIX
            subroutine EVAL_RELAXATION(A, B)
                real, allocatable :: A(:,:), B(:,:)
                integer :: i, j
            end subroutine EVAL_RELAXATION
            subroutine VCONCAT(A, v, bottom)
                real, allocatable :: A(:,:), B_lin(:), B(:,:)
                real, allocatable :: v(:)
                logical :: bottom
            end subroutine VCONCAT
            subroutine L2_NORM(A,B, diff)
                real, allocatable :: A(:,:), B(:,:)
                real :: diff
            end subroutine L2_NORM
        end interface

        integer :: ierr,myid,nproc,stat(MPI_STATUS_SIZE),subarray_datatype, tag
        real, allocatable :: phi(:,:), phi_part(:,:), new_phi_part(:,:), upper(:), lower(:), temp(:,:)
        real :: lphi, rphi, uphi, bphi, error, all_error
        integer :: phi_size, sizes(2), subsizes(2), start(2)
        integer, allocatable :: phi_map(:,:), phi_map_lin(:), id_map(:)
        integer :: i, j, step
        real :: tol
        logical :: proceed
        integer :: proc_print

        tol = 0.1
        proceed = .true.

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        ! Initialization of problem and splitting matrix on processes
        if (myid == 0) then
            open(unit = 8, file = 'debug.log')
            call GET_FORMULATION(phi, phi_size)
            call SLICE_MATRIX(phi, nproc, phi_map)
            allocate(phi_map_lin(nproc*4))
            phi_map_lin = RESHAPE(phi_map, [nproc*4])
        endif

        allocate(id_map(4))
        call MPI_Bcast(phi_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(phi_map_lin, 4, MPI_INTEGER, id_map, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        allocate(phi_part(id_map(1), phi_size))
        allocate(temp(id_map(1), phi_size))
        allocate(upper(phi_size))
        allocate(lower(phi_size))

        if (myid == 0) then
            phi_part = phi(1:id_map(1), 1:phi_size)
        endif
        ! Split phi
        do i = 1, nproc-1
            if (myid == 0) then 
                sizes = [phi_size, phi_size]
                subsizes = [phi_map(1,i+1), phi_size]
                start = [phi_map(2,i+1)-1, 0]
                call MPI_Type_create_subarray(2, sizes, subsizes, start, MPI_ORDER_FORTRAN, MPI_REAL, subarray_datatype, ierr)
                call MPI_Type_commit(subarray_datatype, ierr)
                call MPI_Send(phi, 1, subarray_datatype, i, i, MPI_COMM_WORLD, ierr)
                call MPI_Type_Free(subarray_datatype, ierr)
            endif
            if (myid == i) then
                call MPI_Recv(phi_part, id_map(1)*phi_size, MPI_REAL, 0, i, MPI_COMM_WORLD, stat, ierr)
            endif
        enddo

        ! Main iterations
        step = 1
        tag = nproc
        do while(proceed)
            ! Neigbours exchange
            ! upper
            do i = 0, nproc-1
                if (myid == i) then
                    if (id_map(3) /= -1) then
                        call MPI_Send(phi_part(1,:), phi_size, MPI_REAL, i-1, tag, MPI_COMM_WORLD, ierr)
                        !call MPI_Send(phi_part(ubound(phi_part,1)+1,:), phi_size, MPI_REAL, i-1, tag, MPI_COMM_WORLD, ierr)
                    endif
                endif
                if (myid == i-1) then
                    call MPI_Recv(lower, phi_size, MPI_REAL, i, tag, MPI_COMM_WORLD, stat, ierr)
                    !write(6,*)'lower, id:',myid,'array = ',lower
                endif
                tag = tag + 1
            enddo
            ! and lower
            do i = 0, nproc-1
                if (myid == i) then
                    if (id_map(4) /= -1) then
                        call MPI_Send(phi_part(id_map(1),:), phi_size, MPI_REAL, i+1, tag, MPI_COMM_WORLD, ierr)
                        !call MPI_Send(phi_part(lbound(phi_part,1)-1,:), phi_size, MPI_REAL, i+1, tag, MPI_COMM_WORLD, ierr)
                    endif
                endif
                if (myid == i+1) then
                    call MPI_Recv(upper, phi_size, MPI_REAL, i, tag, MPI_COMM_WORLD, stat, ierr)
                    !write(6,*)'upper, id:',myid,'array = ',upper
                endif
                tag = tag + 1
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
            !write(6,*)'id: ',myid
            !call PRINT_2DARRAY(phi_part, 6)
            call L2_NORM(phi_part, new_phi_part, error)
            ! Realloc phi
            temp = new_phi_part(lbound(temp,1):ubound(temp,1), 1:phi_size)
            call MOVE_ALLOC(temp, phi_part)
            deallocate(new_phi_part)
            ! Eval error
            all_error = 0.d0
            call MPI_Reduce(error, all_error, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            if (myid == 0) then
                if (MOD(step, 20) == 0) then
                    write(6,*)'step: ',step,'error = ',all_error
                endif
                if (all_error < tol) then
                !if (step > 1) then
                    proceed = .false.
                endif
                !proceed = .false.
            endif
            call MPI_BCAST(proceed, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            step = step + 1
        enddo

        !write(6,*)'id: ',myid
        !call PRINT_2DARRAY(phi_part, 6)
        proc_print = 0
        do i = 0, nproc-1
            if (myid == i .and. myid == proc_print) then
                if (i == 0) then
                    open(unit = 9, file = 'result', status="new", action="write")
                    write(9,*)phi_size
                    write(9,*)phi_size
                    call PRINT_2DARRAY(phi_part, 9)
                    close(9)
                else
                    open(unit = 9, file = 'result', status="old", position="append", action="write")
                    call PRINT_2DARRAY(phi_part, 9)
                    close(9)
                endif
                proc_print = proc_print + 1
            endif
            call MPI_BCAST(proc_print, 1, MPI_INTEGER, i, MPI_COMM_WORLD, ierr)
        enddo

        !if (myid == nproc-1) then
        !    call PRINT_2DARRAY(phi_part, 6)
        !endif

        call MPI_FINALIZE(ierr)

end program main


subroutine L2_NORM(A,B, diff)
    real, allocatable :: A(:,:), B(:,:)
    real :: diff

    diff = 0.d0
    do i = lbound(B,1), ubound(B,1)
        do j = lbound(B,2), ubound(B,2)
            diff = diff + (A(i,j) - B(i,j))**2
        enddo
    enddo

end subroutine L2_NORM

subroutine VCONCAT(A, v, bottom)
    real, allocatable :: A(:,:), B_lin(:), B(:,:)
    real, allocatable :: v(:)
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

subroutine EVAL_RELAXATION(A, B)
    real, allocatable :: A(:,:), B(:,:)
    integer :: i, j

    allocate(B(lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2)))
    
    !B = A
    do i = lbound(B,1), ubound(B,1)
        do j = lbound(B,2), ubound(B,2)
            if (i == lbound(B,1) .or. i == ubound(B,1) .or. j == lbound(B,2) .or. j == ubound(B,2)) then
                B(i,j) = A(i,j)
            else
                B(i,j) = ( A(i+1,j) + A(i-1,j) + A(i,j+1) + A(i,j-1) ) / 4.d0
            endif
        enddo
    enddo

end subroutine EVAL_RELAXATION

subroutine SLICE_MATRIX(A, num, map)
    real, allocatable :: A(:,:)
    integer :: num, i
    integer, allocatable :: map(:,:)
    ! Build map array for cuttin matrix A(m x n) on A_i(m x k_i), where i = 1
    ! map = [
    !       [*k_i*]
    !       [*indexes of start A_i matrix in A*]
    !       [*indexes of upper neighbours for A_i*]
    !       [*indexes of lower neighbours for A_i*]
    !]
    allocate(map(4,num))
    map(1,:) = SIZE(A,1) / num
    do i = 1, MOD(SIZE(A,1), num)
        map(1,i) = map(1,i) + 1
    enddo
    map(2,1) = 1
    do i = 2, num
        map(2,i) =  map(2,i-1) + map(1,i-1)
    enddo 
    map(3,:) = -1
    map(3,:) = -1
    do i = 1, num
        if (i == 1) then
            map(3,i) = -1
            map(4,i) = map(2,i) + map(1,i)
        elseif (i==num) then
            map(4,i) = -1
            map(3,i) = map(2,i) - 1
        else    
            map(4,i) = map(2,i) + map(1,i)
            map(3,i) = map(2,i) - 1
        endif
    enddo
end subroutine SLICE_MATRIX
    
subroutine PRINT_2DARRAY(A, thread)
    real, allocatable :: A(:,:)
    integer :: thread
    integer :: i, j

    do i=lbound(A,1) ,ubound(A,1)
        write(thread,*)(A(i,j), j=lbound(A,2),ubound(A,2))
    enddo
end subroutine PRINT_2DARRAY

subroutine GET_FORMULATION(A, n)
    real, allocatable :: A(:,:)
    integer :: n    
    real :: lphi, rphi, uphi, bphi

    open(unit = 10, file = 'formulation')

    !write(6,*)'Input number of points'
    read(10, *)n
    allocate(A(n, n))
    A = 0.d0
    !write(6,*)'Input lbound: '
    read(10,*)lphi
    A(:,1) = lphi
    !write(6,*)'Input rbound: '
    read(10,*)rphi
    A(:,n) = rphi
    !write(6,*)'Input ubound: '
    read(10,*)uphi
    A(n,:) = uphi
    !write(6,*)'Input bbound: '
    read(10,*)bphi
    A(1,:) = bphi
    close(10)
end subroutine GET_FORMULATION