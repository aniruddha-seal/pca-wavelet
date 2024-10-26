      program vertical_lines
      implicit double precision (a-h,o-z)
      
      do i=1,98
         read(71,*) x
         write(72,*) x, 0
         write(72,*) x, 7
      end do

      end program
