program problem5
   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: O = real64
   real(O) :: t, rxn, ryn, rzn, vxn, vyn, vzn, rxp, ryp, rzp, vxp, vyp, vzp, mn, mp, M
   real(O), dimension(3) :: rp, vp, rn, vn
   real(O) :: an, ap, eccn, eccp, in, ip, omegap, omegan, varpin, varpip, fn, fp, h
   integer :: n, nf
   real(O) :: HKepler, HInteraction, HTotal, Man, Map

   mn = 102.43e24_O / 1.98911e30_O
   mp = 0.0127e24_O / 1.98911e30_O

   open(91, file="id000009-XV.csv", status='old')
   open(1, file="id000010-XV.csv", status='old')
   open(99, file="id000009-INT.csv", status='replace')
   open(9, file="id000010-INT.csv", status='replace')
   open(10, file="eandphi.csv", status='replace')

   read(91,*)
   read(1,*)

   read(91,*) t, rxn, ryn, rzn, vxn, vyn, vzn
   read(1,*) t, rxp, ryp, rzp, vxp, vyp, vzp

   rn(1) = rxn
   rn(2) = ryn
   rn(3) = rzn
   vn(1) = vxn
   vn(2) = vyn
   vn(3) = vzn
   rp(1) = rxp
   rp(2) = ryp
   rp(3) = rzp
   vp(1) = vxp
   vp(2) = vyp
   vp(3) = vzp

   h = 0.008_O !timestep in years-- who knows what it should be..
   nf = int(1e5_O / h)
   n = 0
   Man = 0
   Map = 0

   do while(n < nf)

      t = h*n

      call c2b(vn, vp)
      call xv2el(rn(1), rn(2), rn(3), vn(1), vn(2), vn(3), an, eccn, in, omegan, varpin, fn, mn)
      call x2vel(rp(1), rp(2), rp(3), vp(1), vp(2), vp(3), ap, eccp, ip, omegap, varpip, fp, mp)
      call danby(t, rn(:), vn(:), an, eccn, in, omegan, varpin, fn, Man, mn, h, n)
      call danby(t, rp(:), vp(:), ap, eccp, ip, omegap, varpip, fp, Map, mp, h, n) !how to get M and E in these oordinates? Swifter? <--remember 1/2 timestep
      call el2xv(rn(1), rn(2), rn(3), vn(1), vn(2), vn(3), an, eccn, in, omegan, varpin, fn, mn) !currently wrong
      call el2xv(rp(1), rp(2), rp(3), vp(1), vp(2), vp(3), ap, eccp, ip, omegap, varpip, fp, mp) !currently wrong
      call kick(rn, vn, rp, vp, h)
      call xv2el(rn(1), rn(2), rn(3), vn(1), vn(2), vn(3), an, eccn, in, omegan, varpin, fn, mn)
      call x2vel(rp(1), rp(2), rp(3), vp(1), vp(2), vp(3), ap, eccp, ip, omegap, varpip, fp, mp)
      call danby(t, rn(:), vn(:), an, eccn, in, omegan, varpin, fn, Man, mn, h, n)
      call danby(t, rp(:), vp(:), ap, eccp, ip, omegap, varpip, fp, Map, mp, h, n) !how to get M and E in these oordinates? Swifter? <--remember 1/2 timestep
      call el2xv(rn(1), rn(2), rn(3), vn(1), vn(2), vn(3), an, eccn, in, omegan, varpin, fn, mn) !currently wrong
      call el2xv(rp(1), rp(2), rp(3), vp(1), vp(2), vp(3), ap, eccp, ip, omegap, varpip, fp, mp) !currently wrong
      call b2c(vn, vp) !is this even necessary? Also the math is probably wrong
      if(mod(t, 1000.0_O) == 0) then !output every 1 Kyr for now..
         call calcandwrite(t, rn, rp, vn, vp, varpin, varpip, Man, Map) ! delta e will have to be calculated in Pandas
      end if

      n = n + 1
   
   end do

   !H = Hkepler + Hint

   !Hkepler = 

   close(91)
   close(1)
   close(99)
   close(9)
   close(10)

   contains
   subroutine c2b(vn, vp)
      implicit none
      real(O), dimension(3) :: vn, vp, vs
      real(O) :: ms, mn, mp

      ms = 1.0_O
      mn = 102.43e24_O / 1.98911e30_O
      mp = 0.0127e24_O / 1.98911e30_O

      vs(:) = -(((mn * vn(:)) + (mp * vp(:))) / (ms + mn + mp))
      vn(:) = vn(:) + vs(:) !might be vn - vs; not sure
      vp(:) = vp(:) + vs(:)

   end subroutine c2b

   subroutine xv2el(rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, m)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz
      real(O) :: a, e, i, omega, varpi, f
      real(O) :: R, V, h, hx, hy, hz, w, wplusf, Rdot
      real(O) :: rfactor, vfactor, G, msun, mu, m
      real(O) :: sinf, cosf, sino, coso, sinwf, coswf

      G = (6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3)
      msun = 1.0_O
      mu = G*(msun + m)

      rfactor = 1.459787e11_O
      vfactor = rfactor / (60._O * 60._O * 24._O)

      rx = rx * rfactor
      ry = ry * rfactor
      rz = rz * rfactor
      !rx = rx * vfactor
      !ry = ry * vfactor
      !rz = rz * vfactor
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

   end subroutine xv2el

   subroutine danby(t, x, y, a, ecc, i, omega, varpi, f, Ma, m, h, n)
      implicit none
      real(O) :: t, Ma, h, step, r, v, E0, r0, v0
      real(O), dimension(3) :: x, y, x0, y0
      real(O) :: a, ecc, i, omega, varpi, f, M, E, G, mu, pi, period, ms
      real(O) :: ft, gt, fdot, gdot
      integer :: n, j

      step = h / 2.0_O

      G = (6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3)
      ms = 1.0_O
   
      mu = G*(ms + m)

      period = sqrt(((4._O*pi**2)/(mu))*(a**3))
      
      pi = 4.0_O * atan(1.0_O)

      M = M + (n*step) * ((2.0_O*pi)/(period))

      if (M > (2.0_O*pi)) then
         M = mod(M, 2.0_O*pi)
      end if

      r = sqrt((x(1)**2)+(x(2)**2)+(x(3)**2))
      v = sqrt((y(1)**2)+(y(2)**2)+(y(3)**2))

      E = acos(r/((a*(1-ecc))))

      !check if it's going to or away from periapsis

      if(v>0) then
         E = E + pi
      end if

      E = E0

      do j=1,20
         E = solve(E, M, ecc)
      end do

      !f and g functions-- wrong for now (old code needs to be converted)

      r0 = r
      v0 = v
      x0(:) = x(:)
      y0 = y(:)
      !tt = h !not sure if this is the correct way to get the right value of tau
      ft = (a/r0) * (cos(E-E0)-1.0_O) +1.0_O
      !gt = (t-t0) + (1/n) * (sin(Ei-E0) - (Ei-E0))
      gt = step + 1.0_O /n * (sin(E-E0) - (E-E0)) !step = h = t-t0 <--I think??
      x(:) = (ft*x(:)) + (gt*y(:))
      fdot = -(a**2 * n / (r * r0)) * sin(E-E0)
      gdot = (a/r) * (cos(E-E0) - 1.0_O) + 1.0_O
      y(:) = (fdot*x0(:)) + (gdot*y0(:))

   end subroutine danby

   subroutine el2xv(rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, m)
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

   subroutine kick(rn, vn, rp, vp, h)
      implicit none
      real(O) :: h, mu, t, HKepler, HInt, Hsun, magrn, magrp, magvp, magvn, ms, mn, mp, G
      real(O), dimension(3) :: rn, vn, rp, vp, a!, f
      integer :: n

      G = (6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3)

      !Hkepler = 

      ms = 1.0_O
      mn = 102.43e24_O / 1.98911e30_O
      mp = 0.0127e24_O / 1.98911e30_O

      magrn = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))
      magvn = sqrt((vn(1)**2) + (vn(2)**2) + (vn(3)**2))
      magrp = sqrt((rp(1)**2) + (rp(2)**2) + (rp(3)**2))
      magvp = sqrt((vp(1)**2) + (vp(2)**2) + (vp(3)**2))
      
      Hsun = (1.0_O/(2.0_O*ms))* abs((mn*magvn + mp*magvp))**2
      HInt = ((G*mn*mp) / (abs(magrn - magrp)))

      a(:) = ???

      vn(:) = vn(:) + a(:)*h
      vp(:) = vp(:) + a(:)*h

      ! f(1) = -mu * (x(1) / ((x(1)**2) + (x(2)**2) + (x(3)**2))**(3.0_O/2.0_O))
      ! f(2) = -mu * (x(2) / ((x(1)**2) + (x(2)**2) + (x(3)**2))**(3.0_O/2.0_O))
      ! f(3) = -mu * (x(3) / ((x(1)**2) + (x(2)**2) + (x(3)**2))**(3.0_O/2.0_O))
      ! y(:) = y(:) + (h*f(:))

   end subroutine kick

   subroutine b2c(vn, vp)
      implicit none
      real(O), dimension(3) :: vn, vp, vs
      real(O) :: ms, mn, mp

      ms = 1.0_O
      mn = 102.43e24_O / 1.98911e30_O
      mp = 0.0127e24_O / 1.98911e30_O

      vs(:) = -(((mn * vn(:)) + (mp * vp(:))) / (ms + mn + mp))
      vn(:) = vn(:) - vs(:)! might be vn + vs
      vp(:) = vp(:) - vs(:)

   end subroutine b2c

   subroutine calcandwrite(t, rn, rp, vn, vp, varpin, varpip, Man, Map)
      implicit none
      real(O), dimension(3) :: rn, rp, vn, vp
      real(O) :: t, varpin, varpip, Man, Map, lambdan, lambdap, pi, phi, energy, vs
      real(O) :: rnscalar, rpscalar, vnscalar, vpscalar, pn, pp, ms, mn, mp, rnp, G

      ms = 1.0_O
      mn = 102.43e24_O / 1.98911e30_O
      mp = 0.0127e24_O / 1.98911e30_O
      vs = 0._O !I think..
      G = (6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3)
      
      lambdan = Man + varpin
      lambdap = Map + varpip

      phi = (3*lambdap) - (2 * lambdan) - (varpip)

      rnscalar = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))
      rpscalar = sqrt((rp(1)**2) + (rp(2)**2) + (rp(3)**2))
      vnscalar = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))
      vpscalar = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))

      !calculate total energy (could be wrong-- might be - instead of +, or something completely different...)
      energy = (((ms*vs)**2)/(2*ms) + ((mn*vnscalar)**2)/(2*mn) + ((mp*vpscalar)**2)/(2*mp)) + ((G*mn*mp)/(abs(rnp)))

      write(99,*) t, rn(:), vn(:)
      write(9,*) t, rp(:), vp(:)
      write(10, *) t, energy, phi

   end subroutine calcandwrite

   function solve(E0rq, Mrq, ecc) result(Eirq)
      implicit none
      real(O), intent(in) :: E0rq, Mrq
      real(O) :: fq, f1q, f2q, f3q !f1 is f', f2 is f'', etc.
      real(O) :: i1q, i2q, i3q !delta i1, 2, 3, etc.
      real(O) :: Eirq, ecc
      
      fq = E0rq - ecc * sin(E0rq) - Mrq
      f1q = 1.0_O - ecc*cos(E0rq)
      f2q = ecc*sin(E0rq)
      f3q = ecc*cos(E0rq)
      i1q = -(fq/f1q)
      i2q = -(fq / (f1q + i1q * f2q / 2._O))
      i3q = -(fq/(f1q + i2q * f2q / 2._O + i2q**2 * f3q / 6._O))
      
      Eirq = E0rq + i3q
      
      return
  end function solve

end program problem5