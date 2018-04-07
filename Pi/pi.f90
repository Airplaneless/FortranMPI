program main
implicit none
include 'mpif.h'
integer :: N,myid,nproc,stat(MPI_STATUS_SIZE),ierr
real :: psum, Pi
integer :: npart, nstart, nstop, i

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

if(myid.eq.0) then
    write(6,*)'Number:'
    read(5,*)N
    write(6,*)'Evaluate on ',nproc,'processes'
endif

call MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

psum = 0.d0
npart = N/nproc
nstart = myid*npart
nstop = nstart + npart - 1
do i=nstart, nstop
    psum = psum + (-1.d0)**i/(2.d0*i + 1.d0)
    !psum = psum + i
enddo

call MPI_REDUCE(psum,Pi,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

if(myid.eq.0) then
    write(6,*)'Pi = ',Pi*4
endif

call MPI_FINALIZE(ierr)
end
