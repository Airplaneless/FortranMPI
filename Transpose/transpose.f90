program main
implicit none
include 'mpif.h'
integer :: ierr,myid,nproc,stat(MPI_STATUS_SIZE),i,j,n
real :: tmp
real, allocatable :: A(:,:), A0(:,:), A1(:,:), A2(:,:), A3(:,:)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)


if(myid.eq.0) then
    write(6,*)'Read A'
    open(10,file='A',form='formatted',status='old')
    read(10,*)n
    read(10,*)
    allocate (A(n,n))
    do i=1,n
        read(10,*)(A(i,j), j=1,n)
    enddo
    allocate (A1(n/2,n/2), A2(n/2,n/2), A3(n/2, n/2), A0(n/2, n/2))
    do i=1,n/2
        do j=1,n/2
            A0(i,j) = A(i,j)
            A1(i,j) = A(i+n/2, j)
            A2(i,j) = A(i, j+n/2)
            A3(i,j) = A(i+n/2, j+n/2)
        enddo
    enddo
    write(6,*)'A: '
    do i=1,n
        write(6,*)(A(i,j), j=1,n)
    enddo
endif !myid=0
    
call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

if(myid.eq.0) then
    call MPI_SEND(A1,n*n/4,MPI_REAL,1,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(A2,n*n/4,MPI_REAL,2,2,MPI_COMM_WORLD,ierr)
    call MPI_SEND(A3,n*n/4,MPI_REAL,3,3,MPI_COMM_WORLD,ierr)
elseif(myid.ne.0) then
    allocate (A0(n/2,n/2))
    call MPI_RECV(A0,n*n/4,MPI_REAL,0,myid,MPI_COMM_WORLD,stat,ierr)
endif

do i=1,n/2
    do j=i+1,n/2
        tmp = A0(i,j)
        A0(i,j) = A0(j,i)
        A0(j,i) = tmp
    enddo
enddo
        
if(myid.ne.0) then
    call MPI_SEND(A0, n*n/4,MPI_REAL,0,myid+nproc,MPI_COMM_WORLD,ierr)
elseif(myid.eq.0) then
    call MPI_RECV(A1,n*n/4,MPI_REAL,1,1+nproc,MPI_COMM_WORLD,stat,ierr)
    call MPI_RECV(A2,n*n/4,MPI_REAL,2,2+nproc,MPI_COMM_WORLD,stat,ierr)
    call MPI_RECV(A3,n*n/4,MPI_REAL,3,3+nproc,MPI_COMM_WORLD,stat,ierr)
    do i=1,n/2
        do j=1,n/2
            A(i,j) = A0(i,j)
            A(i+n/2,j) = A1(i,j)
            A(i,j+n/2) = A2(i,j)
            A(i+n/2,j+n/2) = A3(i,j)
        enddo
    enddo
    write(6,*)'Transpose A: '
    do i=1,n
        write(6,*)(A(i,j), j=1,n)
    enddo
endif
call MPI_FINALIZE(ierr)
end
