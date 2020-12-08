program problem4
   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: O = real64
   real(O) :: rx, ry, rz, vx, vy, vz, t, h
   integer :: n

   open(21, file="id000002-XV.csv", status='old')
   open(31, file="id000003-XV.csv", status='old')
   open(41, file="id000004-XV.csv", status='old')
   open(51, file="id000005-XV.csv", status='old')
   open(61, file="id000006-XV.csv", status='old')
   open(71, file="id000007-XV.csv", status='old')
   open(81, file="id000008-XV.csv", status='old')
   open(91, file="id000009-XV.csv", status='old')
   open(1, file="id000010-XV.csv", status='old')

   read(21,*)
   ! read(31,*)
   ! read(41,*)
   ! read(51,*)
   ! read(61,*)
   ! read(71,*)
   ! read(81,*)
   ! read(91,*)
   ! read(1,*)

   read(21,*) t, rx, ry, rz, vx, vy, vz
   call leapfrog(t,rx,ry,rz,vx,vy,vz,h)


   close(21)
   ! close(31)
   ! close(41)
   ! close(51)
   ! close(61)
   ! close(71)
   ! close(81)
   ! close(91)
   ! close(1)

   contains
   subroutine leapfrog(t, rx, ry, rz, vx, vy, vz, h)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz, t, h, mu
      real(O), dimension(3) :: x, y, f
      real(O), dimension(3) :: xn12, yn1, xn1


      mu = 1.0_O
      x(1) = rx
      x(2) = ry
      x(3) = rz
      y(1) = vx
      y(2) = vy
      y(3) = vz
      write(*,*) x

      ! f(1) = -mu * (x(1) / ((x(1)**2) + (x(2)**2) + (x(3)**2))**(3.0_O/2.0_O))
      ! f(2) = -mu * (x(2) / ((x(1)**2) + (x(2)**2) + (x(3)**2))**(3.0_O/2.0_O))
      ! f(3) = -mu * (x(3) / ((x(1)**2) + (x(2)**2) + (x(3)**2))**(3.0_O/2.0_O))

      xn12(:) = x(:) + (0.5*h*y(:))
      f(1) = -mu * (xn12(1) / ((xn12(1)**2) + (xn12(2)**2) + (xn12(3)**2))**(3.0_O/2.0_O))
      f(2) = -mu * (xn12(2) / ((xn12(1)**2) + (xn12(2)**2) + (xn12(3)**2))**(3.0_O/2.0_O))
      f(3) = -mu * (xn12(3) / ((xn12(1)**2) + (xn12(2)**2) + (xn12(3)**2))**(3.0_O/2.0_O))
      yn1(:) = y(:) + (h*f(:))
      xn1(:) = x(:) + (0.5_O*h*yn1(:))

      !convert back to rx, ry, rz, etc. to use next integration

      rx = xn1(1)
      ry = xn1(2)
      rz = xn1(3)
      vx = yn1(1)
      vy = yn1(2)
      vz = yn1(3)

   end subroutine leapfrog

end program problem4