program problem2

   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: O = real64
   real(O) :: rx, ry, rz, vx, vy, vz, t
   real(O) :: a, e, i, omega, varpi, f, rp, ra
   real(O) :: mmer, mven, mear, mmar, mjup, msat, mura, mnep, mplu, msun
   integer :: n

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

   !write(*,*) msun, mear

   a = 0.0_O
   e = 0.0_O
   i = 0.0_O
   omega = 0.0_O
   varpi = 0.0_O
   f = 0.0_O
   rp = 0.
   ra = 0.

   open(21, file="id000002-XV.csv", status='old')
   open(22, file="id000002-EL.csv", status='replace')
   open(31, file="id000003-XV.csv", status='old')
   open(32, file="id000003-EL.csv", status='replace')
   open(41, file="id000004-XV.csv", status='old')
   open(42, file="id000004-EL.csv", status='replace')
   open(51, file="id000005-XV.csv", status='old')
   open(52, file="id000005-EL.csv", status='replace')
   open(61, file="id000006-XV.csv", status='old')
   open(62, file="id000006-EL.csv", status='replace')
   open(71, file="id000007-XV.csv", status='old')
   open(72, file="id000007-EL.csv", status='replace')
   open(81, file="id000008-XV.csv", status='old')
   open(82, file="id000008-EL.csv", status='replace')
   open(91, file="id000009-XV.csv", status='old')
   open(92, file="id000009-EL.csv", status='replace')
   open(1, file="id000010-XV.csv", status='old')
   open(2, file="id000010-EL.csv", status='replace')

   read(21,*)
   read(31,*)
   read(41,*)
   read(51,*)
   read(61,*)
   read(71,*)
   read(81,*)
   read(91,*)
   read(1,*)

   do n=1,1001
      
      read(21,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mmer, rp, ra)

      write(22,*) t, a, e, i, omega, varpi, f, rp, ra

      read(31,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mven, rp, ra)

      write(32,*) t, a, e, i, omega, varpi, f, rp, ra

      read(41,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mear, rp, ra)

      write(42,*) t, a, e, i, omega, varpi, f, rp, ra

      read(51,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mmar, rp, ra)

      write(52,*) t, a, e, i, omega, varpi, f, rp, ra

      read(61,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mjup, rp, ra)

      write(62,*) t, a, e, i, omega, varpi, f, rp, ra

      read(71,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, msat, rp, ra)

      write(72,*) t, a, e, i, omega, varpi, f, rp, ra

      read(81,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mura, rp, ra)

      write(82,*) t, a, e, i, omega, varpi, f, rp, ra

      read(91,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mnep, rp, ra)

      write(92,*) t, a, e, i, omega, varpi, f, rp, ra

      read(1,*) t, rx, ry, rz, vx, vy, vz

      call xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, mplu, rp, ra)

      write(2,*) t, a, e, i, omega, varpi, f, rp, ra

   end do

   close(21)
   close(22)
   close(31)
   close(32)
   close(41)
   close(42)
   close(51)
   close(52)
   close(61)
   close(62)
   close(71)
   close(72)
   close(81)
   close(82)
   close(91)
   close(92)
   close(1)
   close(2)

   contains
   subroutine xv2el(t, rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, m, rp, ra)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz, t
      real(O) :: a, e, i, omega, varpi, f, rp, ra
      real(O) :: R, V, h, hx, hy, hz, w, wplusf, Rdot
      real(O) :: rfactor, vfactor, G, msun, mu, m
      real(O) :: sinf, cosf, sino, coso, sinwf, coswf

      G = (6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3)
      msun = 1.0_O
      mu = G*(msun + m)

      vfactor = 365.25_O

      vx = vx * vfactor
      vy = vy * vfactor
      vz = vz * vfactor

      R = sqrt((rx**2) + (ry**2) + (rz**2))
      V = sqrt((vx**2) + (vy**2) + (vz**2))

      Rdot = sign(sqrt(v**2-((h**2)/(r**2))), ((rx*vx)+(ry*vy)+(rz*vz)))

      hx = ((ry*vz) - (rz*vy))
      hy = ((rz*vx) - (rx*vz))
      hz = ((rx*vy) - (ry*vx))
      h = sqrt((hx**2) + (hy**2) + (hz**2))

      a = (((2.0_O/R) - ((V**2)/mu)))**(-1)
      e = sqrt(1.0_O-((h**2)/(mu*a)))
      i = acos(hz/h)

      if (hz > 0) then
         sino = (hx/(h*sin(i)))
         coso = (-hy/(h*sin(i)))
      else
         sino = (-hx/(h*sin(i)))
         coso = (hy/(h*sin(i)))
      end if

      omega = atan2(sino, coso)

      sinwf = (rz/(R*sin(i)))
      coswf = ((1/(cos(omega)))*((rx/R)+(sin(omega)*sinwf*cos(i))))

      wplusf = atan2(sinwf, coswf)

      sinf = (((a*(1-e**2))/(h*e))*Rdot)
      cosf = ((1/e)*((a*(1-e**2))/(R)-1))
      f = atan2(sinf,cosf)
      w = wplusf - f

      varpi = omega + w

      rp = a*(1-e)
      ra = a*(1+e)

   end subroutine xv2el
end program problem2