        program p1
!        include 'mpif.h'
        integer :: i,j,n,m,l,h,k
        real ::a
        double precision ::b,c(10,10)
        double precision, allocatable ::d(:,:)
        
        write(6,*)'input a'
        
        read(5,*)a
        write(6,*)'a =', a*a
        if (a.eq.0) then
          goto 111
        else
        
          open(10,file='A',form='formatted',status='old')
           read(10,*)n
           allocate (d(n,n))
           do i=1,n
             read(10,*)(d(i,j),j=1,n)
             write(6,100)(d(i,j),j=1,n)
           enddo
         endif    
        
        
        deallocate (d)        
111     continue
 100    format(100f5.2)
        end
       
