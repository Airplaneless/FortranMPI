program main
implicit none
include 'mpif.h'
integer  myid,nproc,stat(MPI_STATUS_SIZE),ierr

integer :: n,napp,i,j,nsize
real :: div
real, allocatable :: B(:,:), d(:), Bi(:,:), X(:)

integer, allocatable :: ncounts(:)
integer, allocatable :: nlines(:)
integer, allocatable :: displ(:)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr) 

allocate (nlines(nproc))
allocate (displ(nproc))
allocate (ncounts(nproc))

if(myid == 0) then
    
    write(6,*)'Read A'
    open(10,file='B',form='formatted',status='old')
    read(10,*)n
    read(10,*)napp
    allocate (B(n,n))
    allocate (d(n))
    do i=1,n
        read(10,*)(B(j,i), j=1,n),d(i) 
    enddo
    
    do i=1,n
        div = B(i,i)
        B(i,i) = 0
        do j=2,n
            B(i,j) = -B(i,j)/div
            !B(i,j) = -B(i,j)
        enddo
    enddo

    nlines(:) = n/nproc
    do i=1,mod(n,nproc)
        nlines(i) = nlines(i) + 1
    enddo

    do i=1,nproc
        if(i.eq.1) then
            displ(i) = 0
        elseif (i.ne.1) then
            displ(i) = displ(i-1) + nlines(i)*n
        endif
    enddo
    ncounts = nlines*n
    write(6,*)'nlines:',(nlines(i), i=1,nproc)
    write(6,*)'ncounts:',(ncounts(i), i=1,nproc)
    write(6,*)'displ:',(displ(i), i=1,nproc)
    write(6,*)'b vec:',(d(i), i=1,n)
    write(6,*)'B: '
    do i=1,n
        write(6,*)(B(i,j), j=1,n)
    enddo
endif

call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (myid.ne.0) allocate(d(n))
call MPI_BCAST(d, n, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(nlines, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(displ, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(ncounts, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

allocate(X(n))
X(:) = 0.0

allocate(Bi(n,nlines(myid+1)))
nsize = nlines(myid+1)*n
call MPI_SCATTERV(B, ncounts, displ, MPI_REAL, Bi, nsize, MPI_REAL, 0, MPI_COMM_WORLD, ierr) 


if (myid.eq.0) then
    do i=1,n
        write(6,*)'My id: ',myid,'Bi:',(Bi(i,j), j=1,nlines(1))
    enddo
endif

if (myid.eq.1) then
    do i=1,n
        write(6,*)'My id: ',myid,'Bi:',(Bi(i,j), j=1,nlines(2))
    enddo
endif



call MPI_FINALIZE(ierr)
end
