      include 'jacobi.f90'
      program PH
        implicit none
        integer,parameter :: n=4 !N=mobile atoms' number
        integer,parameter :: tn=12 !tn=total number of atoms

        real*8 h(3*n,3*n),v(3*n,3*n),d(3*n),fh(3*tn,3*tn)
        real*8 mass(n)
        integer itmax,itnum,rotnum,atom(3*n),l(n)
        integer i,j,k

        open(unit=778,file='hessian-dimer.dat',action='READ')
        open(unit=779,file='select.dat',action='READ')
        open(unit=780,file='lm.dat',action='write')

        !read atoms' index and mass in select.dat
        do i=1,n
          read(779,*) k,mass(i)
          !for convenience
          mass(i)=sqrt(mass(i))
          write(*,*) mass(i)
          do j=1,3
            !atom() means index of lm atoms in f-hessian
            atom(3*i-3+j)=3*k-3+j
          enddo
        enddo
        close(779)


      !  do i=1,n
      !  do j=1,3
      !    atom(3*i-3+j)=3*l(i)-3+j
      !  enddo
      !  enddo

        do i=1,3*tn
          do j=i,3*tn
          !you should put full-mass-weighted Hessian in 778
          read(778,*) fh(i,j)
          fh(j,i)=fh(i,j)
          enddo
        enddo
        close(778)
        
        do i=1,3*n
         do j=1,3*n
           !get ph in f-hessian
           h(i,j)=fh(atom(i),atom(j))
         enddo
        enddo 
        
        write(*,*) 'iteration steps:'
        read(*,*)  itmax 
        call jacobi_eigenvalue(3*n,h,itmax,v,d,itnum,rotnum)

        
        write(780,*) "vibrational mode from patial hessian method"
        do i=1,n
          write(780,'(3x,F10.7,3x,F10.7,3x,F10.7)') d(3*i-2),d(3*i-1),d(3*i)
          write(780,'(3x,F10.7,3x,F10.7,3x,F10.7)') sqrt(d(3*i-2)),sqrt(d(3*i-1)),sqrt(d(3*i))
          d(3*i-2)=sqrt(d(3*i-2))*5140.46  !!transfer unit
          d(3*i-1)=sqrt(d(3*i-1))*5140.46
          d(3*i)=sqrt(d(3*i))*5140.46
          write(780,'(3x,F10.3,3x,F10.3,3x,F10.3)') d(3*i-2),d(3*i-1),d(3*i)
          write(780,*) '  '
          do j=1,3*n
          !v(j,i)=i-mode j-line
          k=(j-1)/3+1
          v(j,3*i-2)=v(j,3*i-2)/mass(k)
          v(j,3*i-1)=v(j,3*i-1)/mass(k)
          v(j,3*i)=v(j,3*i)/mass(k)
          write(780,'(3x,F10.7,3x,F10.7,3x,F10.7)') v(j,3*i-2),v(j,3*i-1),v(j,3*i)
          enddo
          write(780,*) '-------'
        enddo

        write(780,*) itnum

        close(780)

        stop
        end program PH



        
            
