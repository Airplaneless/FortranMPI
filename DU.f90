program prg1
		implicit none
		integer ::  nsh,ierr,myid,nproc,n,nstart,nstop,i
		double precision:: mysumel,mysumte,S
		integer sts(MPI_STATUS_SIZE)
		MPI_INIT(IERR)
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD,myid,ierr)
		CALL MPI_COMM_RANK(MPI_COMM_WORLD,nproc,ierr)
		n=15
		myid=
		nsh=
		nproc=
		   if(myid.eq.0)then
		    nsh=n/nproc
		    nsh=n-(nsh*nproc)
		    n=n-nsh
		    nsh=nsh+1
		     else then nsh=0
		   endif
		MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		
		nstart=n-((n/nproc)*(myid+1))
		write(6,*)myid,nstart
		nstop=n-((n/nproc)*myid)+nsh-1
		write(6,*)myid,nstop
		do i=nstart,nstop
		  mysumte=mysumte+(-1.d0)**i/(2.d0*i+1.d0)
		enddo
		MPI_REDUCE(mysumte,S,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		if(myid.eq.0)then
		  S=4*S
		  Write(6,*)'pi = ',S 
		MPI_FINALIZE(IERR)
		end
	
