      PROGRAM linear_function
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      n=462

      do i=1,n
         read(*,*) v1, v2, v3
         write(*,25) -v1/v2, '*x+', -v3/v2, ','
      end do

25    format(f7.4,a3,f7.4,a1)
      
c---- for the vertical lines
      m=113

      do i=1,m
         read(70,*) v1, v2, v3
         if (v1.ne.0.d0) then
         write(71,*) -v3/v1
         end if
      end do

      END PROGRAM

