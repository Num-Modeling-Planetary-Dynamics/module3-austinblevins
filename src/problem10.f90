module symplectic
   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: O = real64
   real(O) :: pi = 4.0_O * atan(1.0_O)
   real(O) :: ms = 1.0_O
   real(O) :: mn = 102.43e24_O / 1.98911e30_O
   real(O) :: mp =  0.0127e24_O / 1.98911e30_O
   real(O) :: G = ((6.672e-11_O * 1.98911e30_O * ((60._O*60._O*24._O*365._O)**2)) / ((1.495978707e11_O)**3))
   !real(O) :: G = 4.0_O * pi

   contains
   subroutine c2b(vn, vp, vs)
      implicit none
      real(O), dimension(3) :: vn, vp, vs

      vs(:) = -(((mn * vn(:)) + (mp * vp(:))) / (ms + mn + mp))
      vn(:) = vn(:) + vs(:) !might be vn - vs; not sure
      vp(:) = vp(:) + vs(:)

      !do we need vs for anything else?

   end subroutine c2b

   subroutine xv2el(rx, ry, rz, vx, vy, vz, a, e, i, omega, varpi, f, m)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz
      real(O) :: a, e, i, omega, varpi, f
      real(O) :: R, V, h, hx, hy, hz, w, wplusf, Rdot
      real(O) :: rfactor, vfactor, m, mu
      real(O) :: sinf, cosf, sino, coso, sinwf, coswf

      mu = G*(ms + m)

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
      real(O) :: a, ecc, i, omega, varpi, f, M, Ea, mu, period
      real(O) :: ft, gt, fdot, gdot
      integer :: n, j

      write(*,*) "rx, ry, rz, vx, vy, vz [before solve]"
      write(*,*) x, y

      step = h
   
      mu = G*(ms + m)

      period = sqrt(((4._O*pi**2)/(mu))*(a**3))
      
      pi = 4.0_O * atan(1.0_O)

      Ma = Ma + (n*step) * ((2.0_O*pi)/(period))



      if (Ma > (2.0_O*pi)) then
         Ma = mod(Ma, 2.0_O*pi)
      end if

      r = sqrt((x(1)**2)+(x(2)**2)+(x(3)**2))
      v = sqrt((y(1)**2)+(y(2)**2)+(y(3)**2))

      !write(*,*) (((1-(r/a)))/ecc), r, a, ecc, Ea

      Ea = acos((((1-(r/a)))/ecc))
      
      write(*,*) "cos(E), r, a, ecc, Ea, Ma [before solve]"
      write(*,*) (((1-(r/a)))/ecc), r, a, ecc, Ea, Ma

      !check if it's going to or away from periapsis

      !if(v>0) then
        ! E = E + pi
      !end if

      if ((Ea + pi) > (2*pi)) then
         Ea = Ea - pi
      end if

      E0 = Ea

      do j=1,20
         Ea = solve(Ea, Ma, ecc)
         !Ea = test(Ea, Ma, ecc)
      end do

      write(*,*) "Ea, Ma, ecc [after solve]"
      write(*,*) Ea, Ma, ecc

      !f and g functions

      r0 = r
      v0 = v
      x0(:) = x(:)
      y0 = y(:)
      !tt = h !not sure if this is the correct way to get the right value of tau
      ft = (a/r0) * (cos(Ea-E0)-1.0_O) +1.0_O
      !gt = (t-t0) + (1/n) * (sin(Ei-E0) - (Ei-E0))
      gt = step + 1.0_O /n * (sin(Ea-E0) - (Ea-E0)) !step = h = t-t0 <--I think??
      x(:) = (ft*x(:)) + (gt*y(:))
      !x(:) = (ft*x(:)) + (gy*y(:))
      fdot = -(a**2 * n / (r * r0)) * sin(Ea-E0)
      gdot = (a/r) * (cos(Ea-E0) - 1.0_O) + 1.0_O
      y(:) = (fdot*x0(:)) + (gdot*y0(:))

      write(*,*) "rx, ry, rz, vx, vy, vz [after solve]"
      write(*,*) x(:), y(:)

      !E0 = Ea

   end subroutine danby

   !not copying subroutine el2xv from problem5.f90 because I'm not sure if it's needed

   subroutine kick(rn, vn, rp, vp, dt)
      implicit none
      real(O) :: dt, mu, t, HInt, Hsun, magrn, magrp, magvp, magvn, ms, mn
      real(O), dimension(3) :: rn, vn, rp, vp, a
      integer :: n

      magrn = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))
      magvn = sqrt((vn(1)**2) + (vn(2)**2) + (vn(3)**2))
      magrp = sqrt((rp(1)**2) + (rp(2)**2) + (rp(3)**2))
      magvp = sqrt((vp(1)**2) + (vp(2)**2) + (vp(3)**2))
      
      !Hsun = (1.0_O/(2.0_O*ms))* abs((mn*magvn + mp*magvp))**2
      !HInt = ((G*mn*mp) / (abs(magrn - magrp)))

      a(:) = ((G*mn*mp)/(mn+mp)) / (abs(rn(:)-rp(:))**2)

      vn(:) = vn(:) + a(:)*dt
      vp(:) = vp(:) + a(:)*dt

   end subroutine kick

   subroutine b2c(vn, vp,  vs)
      implicit none
      real(O), dimension(3) :: vn, vp, vs

      vs(:) = -(((mn * vn(:)) + (mp * vp(:))) / (ms + mn + mp))
      vn(:) = vn(:) - vs(:)! might be vn + vs
      vp(:) = vp(:) - vs(:)

   end subroutine b2c

   subroutine calcandwrite(t, rn, rp, vn, vp, varpin, varpip, Man, Map, vs)
      implicit none
      real(O), dimension(3) :: rn, rp, vn, vp, vs
      real(O) :: t, varpin, varpip, Man, Map, lambdan, lambdap, pi, phi, energy, vsscalar
      real(O) :: rnscalar, rpscalar, vnscalar, vpscalar, pn, pp, ms, mn, mp, rnp, G

      vsscalar = sqrt((vs(1)**2) + (vs(2)**2) + (vs(3)**2))
      
      lambdan = Man + varpin
      lambdap = Map + varpip

      phi = (3*lambdap) - (2 * lambdan) - (varpip)

      rnscalar = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))
      rpscalar = sqrt((rp(1)**2) + (rp(2)**2) + (rp(3)**2))
      vnscalar = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))
      vpscalar = sqrt((rn(1)**2) + (rn(2)**2) + (rn(3)**2))

      !calculate total energy (could be wrong-- might be - instead of +, or something completely different...)
      energy = (((ms*vsscalar)**2)/(2*ms) + ((mn*vnscalar)**2)/(2*mn) + ((mp*vpscalar)**2)/(2*mp)) + ((G*mn*mp)/(abs(rnp)))

      write(99,*) t, rn(:), vn(:)
      write(9,*) t, rp(:), vp(:)
      write(10, *) t, energy, phi

   end subroutine calcandwrite

   subroutine lindrift(rn, vn, rp, vp, dt)
      implicit none
      real(O), dimension(3) :: rp, vp, rn, vn, vs, Ptot, sigmap
      real(O) :: dt, mtot

      !"Each body takes a linear drift in position of the amount (h/2*ms)*Ptot" (p. 5, Duncan et al.)

      mtot = ms + mp + mn

      !sigmap = (ms*vs(:)) + (mn*vn(:)) + (mp*vp(:))

      !Ptot(:) = (mn*vn(:) - (mn/mtot)*(sigmap)) + (mp*vp(:) - (mp/mtot)*(sigmap)) + (sigmap)

      !rn(:) = rn(:) + (dt/(2._O * ms)*(Ptot(:)))
      !rp(:) = rp(:) + (dt/(2._O * ms)*(Ptot(:)))

      vs(:) = (((mn * vn(:)) + (mp*vp(:))) / mtot)

      rn(:) = rn(:) + (dt*vs(:))
      rp(:) = rp(:) + (dt*vs(:))

   end subroutine lindrift

   subroutine init(rx, ry, rz, vx, vy, vz, r, v)
      implicit none
      real(O) :: rx, ry, rz, vx, vy, vz
      real(O), dimension(3) :: r, v

      r(1) = rx
      r(2) = ry
      r(3) = rz
      v(1) = vx * 365._O !to get into years
      v(2) = vy * 365._O
      v(3) = vz * 365._O

   end subroutine init

   function test(E0r, Mr, ecc) result(Eir)
      implicit none
      real(O), intent(in) :: E0r, Mr, ecc
      real(O) :: Eir

      Eir = Mr + ecc*sin(E0r)

      write (*,*) Eir

   end function test

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
      i2q = -(fq / (f1q + (0.5_O * i1q * f2q)))
      i3q = -(fq/(f1q + (0.5 * i2q * f2q) + (1._O/6._O)*i2q**2 * f3q))
      
      Eirq = E0rq + i3q

      if(Eirq > 2.0_O*pi) then
         Eirq = mod(Eirq, 2.0_O*pi)
      end if

      !write(*,*) Eirq
      
      return
  end function solve

end module symplectic

program integrator
   use symplectic
   implicit none
   real(O) :: t, rxn, ryn, rzn, vxn, vyn, vzn, rxp, ryp, rzp, vxp, vyp, vzp, Man, Map
   real(O), dimension(3) :: rp, vp, rn, vn, vs
   real(O) :: an, ap, eccn, eccp, in, ip, omegap, omegan, varpin, varpip, fn, fp, h
   integer :: n, nf

   open(91, file="id000009-XV.csv", status='old')
   open(1, file="id000010-XV.csv", status='old')
   open(99, file="id000009-INT.csv", status='replace')
   open(9, file="id000010-INT.csv", status='replace')
   open(10, file="eandphi.csv", status='replace')

   read(91,*)
   read(1,*)

   read(91,*) t, rxn, ryn, rzn, vxn, vyn, vzn
   read(1,*) t, rxp, ryp, rzp, vxp, vyp, vzp

   call init(rxn, ryn, rzn, vxn, vyn, vzn, rn, vn)
   call init(rxp, ryp, rzp, vxp, vyp, vzp, rp, vp)

   h = 0.008_O !timestep in years-- who knows what it should be..
   nf = int(1e5_O / h)
   n = 1
   Man = 0
   Map = 0

   do while(n < nf)

      t = h*n

      call c2b(vn, vp, vs)
      call xv2el(rn(1), rn(2), rn(3), vn(1), vn(2), vn(3), an, eccn, in, omegan, varpin, fn, mn)
      call xv2el(rp(1), rp(2), rp(3), vp(1), vp(2), vp(3), ap, eccp, ip, omegap, varpip, fp, mp)
      call lindrift(rn, vn, rp, vp, h/2._O)
      call kick(rn, vn, rp, vp, h/2.0_O)
      call xv2el(rn(1), rn(2), rn(3), vn(1), vn(2), vn(3), an, eccn, in, omegan, varpin, fn, mn)!not in structure but I think it should be there
      call xv2el(rp(1), rp(2), rp(3), vp(1), vp(2), vp(3), ap, eccp, ip, omegap, varpip, fp, mp)!not in structure but I think it should be there
      call danby(t, rn(:), vn(:), an, eccn, in, omegan, varpin, fn, Man, mn, h, n)
      call danby(t, rp(:), vp(:), ap, eccp, ip, omegap, varpip, fp, Map, mp, h, n)
      call kick(rn, vn, rp, vp, h/2.0_O)
      call lindrift(rn, vn, rp, vp, h/2._O)
      call b2c(vn, vp, vs)
      if(mod(t, 1000.0_O) == 0) then !output every 1 Kyr for now..
         call xv2el(rn(1), rn(2), rn(3), vn(1), vn(2), vn(3), an, eccn, in, omegan, varpin, fn, mn)
         call xv2el(rp(1), rp(2), rp(3), vp(1), vp(2), vp(3), ap, eccp, ip, omegap, varpip, fp, mp)
         call calcandwrite(t, rn, rp, vn, vp, varpin, varpip, Man, Map, vs) ! delta e will have to be calculated in Pandas
      end if

      n = n + 1
   
   end do

   close(91)
   close(1)
   close(99)
   close(9)
   close(10)

end program integrator