program example1
implicit none
include 'mpif.h'
integer :: ierr,myid,nproc,stat(MPI_STATUS_SIZE)
real :: a,b 

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

if(myid.eq.0) then
    write(6,*)'Input A'
    read(5,*)a
    call MPI_SEND(a,1,MPI_REAL,nproc-1,1,MPI_COMM_WORLD,ierr)
elseif(myid.eq.nproc-1) then
    call MPI_RECV(b,1,MPI_REAL,0,1,MPI_COMM_WORLD,stat,ierr)
    write(6,*)'A = ',b
endif

call MPI_FINALIZE(ierr)
end
 
