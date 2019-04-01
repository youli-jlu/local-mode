       program main
         implicit none

         integer,parameter :: natoms=12
         real*8 hess(natoms*3,natoms*3)
         integer i,j
         integer k
         character(len=1024) line

         integer*4,parameter :: readfileunit=770
         integer*4,parameter :: writefileunit=771

         open (unit=readfileunit,file="frequencies.out",action="READ")
         open (unit=writefileunit,file="hessian-dimer.dat",action="WRITE")

         line="guess"
         do while( line .NE. "####")
         read (readfileunit,*) line
         enddo
         k=1
         do while(k <= natoms*3)
           read (readfileunit,100) 
100        format (BN,A100)
           call readhessian(k,hess,readfileunit,natoms*3)
         end do

         do j=1,natoms*3
          do i=j,natoms*3
           write (writefileunit,*) hess(i,j)
          enddo
           write (writefileunit,*) " "
         enddo

         end program main

         subroutine readhessian(k,hess,readfileunit,natoms)
           integer i,j
           integer k
           integer readfileunit
           integer natoms
           real*8 hess(natoms,natoms)
           character(len=24)  ch
         
           do i=k,natoms
             read(readfileunit,'(A22)',advance="no") ch
             write(*,*) ch
           do j=k,natoms
             if (i.EQ.j.OR.mod(j,5).EQ.0) then
               read(readfileunit,*) hess(i,j)
               write(*,*) hess(i,j)
               exit
             else
               read(readfileunit,'(2x,F10.7)',advance="no") hess(i,j)
               write(*,*) hess(i,j)
             end if 
           enddo 
           enddo
           k=j+1
           write(*,*) k

           end subroutine readhessian
      

