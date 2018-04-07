program main
implicit none
include 'mpif.h'
integer :: ierr,myid,nproc,stat(MPI_STATUS_SIZE),n1,n2,i,j
real, allocatable :: d(:,:), c(:)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

if(myid.eq.0) then
    write(6,*)'Read A'
    open(10,file='A',form='formatted',status='old')
    read(10,*)n1
    read(10,*)
    allocate (d(n1,n1))
    do i=1,n1
        read(10,*)(d(i,j), j=1,n1)
    enddo
    call MPI_SEND(n1,1,MPI_INTEGER,nproc-1,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(d(1,:),n1,MPI_REAL,nproc-1,2,MPI_COMM_WORLD,ierr)
elseif(myid.eq.nproc-1) then
    call MPI_RECV(n2,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,stat,ierr)
    allocate (c(n2))
    call MPI_RECV(c,n2,MPI_REAL,0,2,MPI_COMM_WORLD,stat,ierr)
    write(6,*)'A: '
    do i=1,n2
        write(6,*)(c(i))
    enddo
endif

call MPI_FINALIZE(ierr)
end
 
