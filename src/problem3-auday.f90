program problem3

   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: O = real64
   real(O) :: rx, ry, rz, vx, vy, vz, t
   real(O) :: a, e, i, omega, varpi, f
   real(O) :: mmer, mven, mear, mmar, mjup, msat, mura, mnep, mplu, msun
   integer :: j

   msun = 1.98911e30_O !kg
   mmer = 0.3302e24_O / msun
   mven = 4.8685e24_O / msun
   mear = 5.9736e24_O / msun
   mmar = 0.64185e24_O / msun
   mjup = 1898.6e24_O / msun
   msat = 568.46e24_O / msun
   mura = 86.832e24_O / msun
   mnep = 102.43e24_O / msun
   mplu = 0.0127e24_O / msun

   rx = 0.0_O
   ry = 0.0_O
   rz = 0.0_O
   vx = 0.0_O
   vy = 0.0_O
   vz = 0.0_O

   open(23, file="id000002-XV-new.csv", status='replace')
   open(22, file="id000002-EL.csv", status='old')
   open(33, file="id000003-XV-new.csv", status='replace')
   open(32, file="id000003-EL.csv", status='old')
   open(43, file="id000004-XV-new.csv", status='replace')
   open(42, file="id000004-EL.csv", status='old')
   open(53, file="id000005-XV-new.csv", status='replace')
   open(52, file="id000005-EL.csv", status='old')
   open(63, file="id000006-XV-new.csv", status='replace')
   open(62, file="id000006-EL.csv", status='old')
   open(73, file="id000007-XV-new.csv", status='replace')
   open(72, file="id000007-EL.csv", status='old')
   open(83, file="id000008-XV-new.csv", status='replace')
   open(82, file="id000008-EL.csv", status='old')
   open(93, file="id000009-XV-new.csv", status='replace')
   open(92, file="id000009-EL.csv", status='old')
   open(3, file="id000010-XV-new.csv", status='replace')
   open(2, file="id000010-EL.csv", status='old')

   do j=1,1001
      
      read(22,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mmer)

      write(23,*) t, rx, ry, rz, vx, vy, vz

      read(32,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mven)

      write(33,*) t, rx, ry, rz, vx, vy, vz

      read(42,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mear)

      write(43,*) t, rx, ry, rz, vx, vy, vz

      read(52,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mmar)

      write(53,*) t, rx, ry, rz, vx, vy, vz

      read(62,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mjup)

      write(63,*) t, rx, ry, rz, vx, vy, vz

      read(72,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, msat)

      write(73,*) t, rx, ry, rz, vx, vy, vz

      read(82,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mura)

      write(83,*) t, rx, ry, rz, vx, vy, vz

      read(92,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mnep)

      write(93,*) t, rx, ry, rz, vx, vy, vz

      read(2,*) t, a, e, i, omega, varpi, f

      call el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mplu)

      write(3,*) t, rx, ry, rz, vx, vy, vz

   end do

   close(23)
   close(22)
   close(33)
   close(32)
   close(43)
   close(42)
   close(53)
   close(52)
   close(63)
   close(62)
   close(73)
   close(72)
   close(83)
   close(82)
   close(93)
   close(92)
   close(3)
   close(2)

   contains
   subroutine el2xv(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, m)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz, t, m
      real(O) :: a, e, i, omega, varpi, f
      real(O) :: R, V, h, hx, hy, hz, w, wplusf, Rdot
      real(O) :: rfactor, vfactor, G, msun, mu, n, period, pi

      pi = 4.0_O*atan(1.0_O)

      G = (6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3)
      msun = 1.0_O
      mu = G*(msun + m)

      period = sqrt(((4._O*pi**2)/(mu))*(a**3))

      n = ((2*pi)/(period))

      w = varpi - omega

      r = (a*(1._O-e**2)/(1._O+e*cos(f)))
      rx = r*((cos(omega)*(cos(w+f)))-(sin(omega)*(sin(w+f))*cos(i)))
      ry = r*((sin(omega)*(cos(w+f)))+(cos(omega)*(sin(w+f))*cos(i)))
      rz = r*(sin(w+f)*(sin(i)))

      v = sqrt(mu*((2._O/r)-(1._O/a)))
      rdot = ((n*a)/(sqrt(1-e**2)))*(e*(sin(f)))
      vx = rdot*((cos(omega)*(cos(w+f)))-(sin(omega)*(sin(w+f))*cos(i)))
      vy = rdot*((sin(omega)*(cos(w+f)))+(cos(omega)*(sin(w+f))*cos(i)))
      vz = rdot*(sin(w+f)*(sin(i)))

      vfactor = 365.25_O

      vx = vx * vfactor
      vy = vy * vfactor
      vz = vz * vfactor

      return
   end subroutine el2xv
end program problem3