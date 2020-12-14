program problem4
   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: O = real64
   real(O) :: rx, ry, rz, vx, vy, vz, t, h!, energy, am!, n
   integer :: i

   !h = 365250.0_O
   h = 0.008_O
   !h = 3652500.0_O

   open(21, file="id000002-XV.csv", status='old')
   open(31, file="id000003-XV.csv", status='old')
   open(41, file="id000004-XV.csv", status='old')
   open(51, file="id000005-XV.csv", status='old')
   open(61, file="id000006-XV.csv", status='old')
   open(71, file="id000007-XV.csv", status='old')
   open(81, file="id000008-XV.csv", status='old')
   open(91, file="id000009-XV.csv", status='old')
   open(1, file="id000010-XV.csv", status='old')

   open(24, file="id000002-LF.csv", status='replace')
   open(34, file="id000003-LF.csv", status='replace')
   open(44, file="id000004-LF.csv", status='replace')
   open(54, file="id000005-LF.csv", status='replace')
   open(64, file="id000006-LF.csv", status='replace')
   open(74, file="id000007-LF.csv", status='replace')
   open(84, file="id000008-LF.csv", status='replace')
   open(94, file="id000009-LF.csv", status='replace')
   open(4, file="id000010-LF.csv", status='replace')

   read(21,*)
   read(31,*)
   read(41,*)
   read(51,*)
   read(61,*)
   read(71,*)
   read(81,*)
   read(91,*)
   read(1,*)

   read(21,*) t, rx, ry, rz, vx, vy, vz
   i = 2
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(21)
   
   read(31,*) t, rx, ry, rz, vx, vy, vz
   i = 3
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(31)

   read(41,*) t, rx, ry, rz, vx, vy, vz
   i = 4
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(41)

   read(51,*) t, rx, ry, rz, vx, vy, vz
   i = 5
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(51)

   read(61,*) t, rx, ry, rz, vx, vy, vz
   i = 6
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(61)

   read(71,*) t, rx, ry, rz, vx, vy, vz
   i = 7
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(71)

   read(81,*) t, rx, ry, rz, vx, vy, vz
   i = 8
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(81)

   read(91,*) t, rx, ry, rz, vx, vy, vz
   i = 9
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(91)

   read(1,*) t, rx, ry, rz, vx, vy, vz
   i = 10
   call go(t,rx,ry,rz,vx,vy,vz,h,i)
   close(1)

   close(24)
   close(34)
   close(44)
   close(54)
   close(64)
   close(74)
   close(84)
   close(94)
   close(4)

   contains
   subroutine leapfrog(rx, ry, rz, vx, vy, vz)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz, t, mu
      real(O), dimension(3) :: x, y, f
      real(O), dimension(3) :: xn12, yn1, xn1


      !mu = 1.0_O
      mu = (6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3) ! G (au/yr) * Msun where Msun=1
      x(1) = rx
      x(2) = ry
      x(3) = rz
      y(1) = vx
      y(2) = vy
      y(3) = vz
      !write(*,*) x

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

   subroutine go(t, rx, ry, rz, vx, vy, vz, h, i)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz, t, h, mu, e0, energy, de, am, am0, dam, r, v
      real(O), dimension(3) :: x, y, f
      real(O), dimension(3) :: xn12, yn1, xn1
      integer :: n, nf, i

      !nf = int((1e6*365)/h)
      nf = int(1e6/h)
      n = 0

      !do n=1,1000 !while(n<[1My in d])
      do while(n<nf)
         t = n*h

         call calc(rx, ry, rz, vx, vy, vz, energy, am)
         e0 = energy
         am0 = am

         call leapfrog(rx, ry, rz, vx, vy, vz)
         call calc(rx, ry, rz, vx, vy, vz, energy, am)
         dam = am - am0
         de = energy - e0

         if (mod(t, 10000.0_O) == 0) then
            if (i==2) then
               write(24, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==3) then
               write(34, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==4) then
               write(44, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==5) then
               write(54, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==6) then
               write(64, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==7) then
               write(74, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==8) then
               write(84, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==9) then
               write(94, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            else if (i==10) then
               write(4, *) t, rx, ry, rz, vx, vy, vz, energy, am, de, dam
            end if
         end if

         n = n + 1
      end do


   end subroutine go 

   subroutine calc(rx, ry, rz, vx, vy, vz, energy, am)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz, t, h, mu, energy, am
      real(O) :: hx, hy, hz, r, v
      real(O), dimension(3) :: x, y, f
      real(O), dimension(3) :: xn12, yn1, xn1

      r = sqrt((rx**2) + (ry**2) + (rz**2))
      v = sqrt((vx**2) + (vy**2) + (vz**2))
      energy = (0.5_O*v**2) - (1.0_O/r)
      hx = ((ry*vz) - (rz*vy))
      hy = ((rz*vx) - (rx*vz))
      hz = ((rx*vy) - (ry*vx))
      am = sqrt((hx**2) + (hy**2) + (hz**2)) !this is a bit sus.. could be one issue
   
   end subroutine calc

end program problem4