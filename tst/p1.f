        program p1
!        include 'mpif.h'
        integer :: i,j,n,m,l,h,k
        real ::a
        double precision ::b,c(10,10)
!        double precision, allocatable ::d(:,:)
        
        write(6,*)'input a'
        
        read(5,*)a
        b = 10.0d0
        write(6,*)'a =', a
        write(6,*)'b =', b
        if (a.eq.0) then
          goto 111
        else
        
          open(10,file='A',form='formatted',status='old')
           read(10,*)n
        call wr(n)
        
!        deallocate (d)        
         endif    
111     continue

        end
    
       subroutine wr(n)
       integer :: n
        double precision, allocatable ::d(:,:)
       
           allocate (d(n,n))
           do i=1,n
             read(10,*)(d(i,j),j=1,n)
             write(6,100)(d(i,j),j=1,n)
           enddo
 100    format(100f5.2)
       end
