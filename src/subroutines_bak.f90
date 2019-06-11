!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          flux computes the mutual inductance/2/pi between        **
!**          two circulars of rectangular cross section.             **
!**                                                                  **
!**                                                                  **
!**     calling arguments:                                           **
!**       r1,r2...........radius of first and second coil            **
!**       z1,z2...........elevation                                  **
!**       w1,w2...........width                                      **
!**       h1,h2...........height                                     **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
subroutine flux(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,fuxx)
  use exparm,only:mgaus1,mgaus2
  implicit integer*4 (i-n), real*8 (a-h, o-z)
  dimension post1(mgaus1),wght1(mgaus1),post2(mgaus2) &
           ,wght2(mgaus2)
  data init/0/

  if (init > 0) go to 100
!vas introduced to make it ok in hydra
!org      ngaus1=mgaus1
!org      ngaus2=mgaus2
!org      call lgauss(post1,wght1,ngaus1,ner)
!org      call lgauss(post2,wght2,ngaus2,ner)
!      print*,'initializing',init,mgaus1,mgaus2
  call lgauss(post1,wght1,mgaus1,ner)
  call lgauss(post2,wght2,mgaus2,ner)
  init=1
  100 continue
!
  x = 0.
  fuxx = 0.
  zf = z1
  zs = z2
  hf2 = h1*.5
  hs2 = h2*.5

  if (t12==0) go to 200
  xl1 = r1-0.5*w1-0.5*h1/abs(t12)
  xr1 = r1+0.5*w1+0.5*h1/abs(t12)
  xbm1 = xl1+w1
  xtm1 = xr1-w1
  if (t12 < 0.) xbm1 = xr1 - w1
  if (t12 < 0.) xtm1 = xl1 + w1
!
  200 if (t22==0) go to 300
  xl2 = r2-0.5*w2-0.5*h2/abs(t22)
  xr2 = r2+0.5*w2+0.5*h2/abs(t22)
  xbm2 = xl2+w2
  xtm2 = xr2-w2
  if (t22 < 0.) xbm2 = xr2 - w2
  if (t22 < 0.) xtm2 = xl2 + w2
!
  300 continue
!      print*,'vasan : ',init,mgaus1,mgaus2
  do i = 1,mgaus1
!org      do i = 1,ngaus1
    rf = r1+.5*w1*post1(i)
    if (t12 ~= 0) rf = r1+(0.5*w1+0.5*h1/abs(t12))*post1(i)
    do j = 1,mgaus2
!org         do j = 1,ngaus2
      rs = r2+0.5*w2*post2(j)
      if (t22 ~= 0) rs = r2+(0.5*w2+0.5*h2/abs(t22))*post2(j)
      call soleno(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2,&
                  xtm1,xtm2,hfa,hsa,rf,rs,solx)
      fuxx = fuxx+wght1(i)*wght2(j)/hfa/hsa*solx
    enddo 
  enddo 
!
return
end subroutine flux


!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          lgauss computes the zeroes of the legendre polynomial   **
!**          and their associated weights for a gaussian quadrature. **
!**                                                                  **
!**     calling arguments:                                           **
!**       x...............zeroes of legendre polynomial              **
!**       w...............weights                                    **
!**       n...............order of legendre polynomial               **
!**       nn..............error flag                                 **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine lgauss(x,w,n,nn)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(n) ::  x,w
!
      nn = 0
      if (n<1) then
        nn = 1
        return
      elseif(n==1) then
!----------------------------------------------------------------------
!-- request for a zero point formula is meaningless                  --
!----------------------------------------------------------------------
        x(1) = 0.
        w(1) = 2.
        return
      endif
!----------------------------------------------------------------------
!-- for a one point formula, send back                               --
!-- results without computing.                                       --
!----------------------------------------------------------------------
      r = n
      g = -1.
!----------------------------------------------------------------------
!-- the initial guess for the smallest root                          --
!-- of p(n) is taken as -1.                                          --
!----------------------------------------------------------------------
      do i = 1,n
        test = -2.
        ic = n+1-i
!----------------------------------------------------------------------
!-- whenever we find a root of the                                   --
!-- polynomial, its negative is also a root.                         --
!-- the index ic tells WHERE to store the other root                 --
!----------------------------------------------------------------------
        if (ic < i) go to 950
     40 continue 
        s = g
        t = 1.
        u = 1.
        v = 0.
!----------------------------------------------------------------------
!-- evaluation of the n-th legendre polynomial                       --
!-- and its first derivative.                                        --
!-- WHERE   u = ds/dx                                                --
!--         v = dt/dx                                                --
!--         dp=dp/dx                                                 --
!----------------------------------------------------------------------
      do k = 2,n
      a = k
      p = ((2.0*a-1.0)*s*g-(a-1.0)*t)/a
      dp = ((2.0*a-1.0)*(s+g*u)-(a-1.0)*v)/a
      v = u
      u = dp
      t = s
      s = p
      enddo 
      if (abs((test-g)/(test+g)) < 0.0000005) go to 100
      sum = 0.
      if (i==1) go to 70
!----------------------------------------------------------------------
!-- the following computes the reduced                               --
!-- legendre polynomial and its derivative.                          --
!----------------------------------------------------------------------
      do k = 2,i
         sum = sum+1./(g-x(k-1))
      enddo 
   70 test = g
      g = g-p/(dp-p*sum)
      go to 40
  100 x(ic) = -g
      x(i) = g
      w(i) = 2./(r*t*dp)
      w(ic) = w(i)
  900 g = g-r*t/((r+2.)*g*dp+r*v-2.*r*t*sum)
	    enddo
  950 return

      end subroutine lgauss

!**********************************************************************
!**                                                                  **
!**     main program:  MHD fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          soleno computes the inductance for a solenoid           **
!**                                                                  **
!**     calling arguments:                                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine soleno(ra,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2, &
                        xtm1,xtm2,hfa,hsa,rf,rs,sol)
      use consta,only:pi
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension z(2,2)
      real*8 ksq,kpsq,msl,mut
      data init/0/
!
      if (init > 0) go to 5
      rpi=pi
      rpi2=rpi*0.5
      rh=rpi*2.0e-07
      ut=rh/3.0
      err=1.0e-05
      init=1
    5 continue
!
      dr = rf-rs
      drsq = dr*dr
      tr = rf+rs
      trsq = tr*tr
      fr = 4.*rf*rs
      sol = 0.
      csq = fr/trsq
      cpsq = 1.-csq
!
   10 if (t12 ~= 0.) go to 20
      r = rf-ra+0.5*w1
      zc1 = z1-0.5*h1-0.5*w1*t1
      zb1 = zc1+t1*r
      zt1 = zb1+h1
      go to 30
!
   20 zb1 = z1-h1/2.
      if (t12 < 0.) go to 25
      if (rf > xbm1) zb1 = zb1+t12*(rf-xbm1)
      zt1 = z1+h1/2.
      if (rf < xtm1) zt1 = zt1-t12*(xtm1-rf)
      go to 30
!
   25 if (rf < xbm1) zb1 = zb1 + t12*(rf-xbm1)
      zt1 = z1 + h1/2.
      if (rf > xtm1) zt1 = zt1 - t12*(xtm1-rf)
!
   30 if (t22 ~= 0.) go to 40
      r = rs-r2+0.5*w2
      zc2 = z2-0.5*h2-0.5*w2*t2
      zb2 = zc2+t2*r
      zt2 = zb2+h2
      go to 50
!
   40 zb2 = z2-h2/2.
      if (t22 < 0.) go to 45
      if (rs > xbm2) zb2 = zb2+t22*(rs-xbm2)
      zt2 = z2+h2/2.
      if (rs < xtm2) zt2 = zt2-t22*(xtm2-rs)
      go to 50
!
   45 if (rs < xbm2) zb2 = zb2 + t22*(rs-xbm2)
      zt2 = z2 + h2/2.
      if (rs > xtm2) zt2 = zt2 - t22*(xtm2-rs)
!
   50 z(1,1) = zb1
      z(2,1) = zb2
      z(1,2) = zt1
      z(2,2) = zt2
      hfa = zt1-zb1
      hsa = zt2-zb2
!
      do i = 1,2
      do j = 1,2
      sign = -.25
      if (i ~= j) sign = .25
      dz = z(1,i)-z(2,j)
      dzsq = dz*dz
      r2sq = dzsq+drsq
      r1sq = dzsq+trsq
      r1 = sqrt(r1sq)
      ksq = fr/r1sq
      t = 2./ksq-1.
      kpsq = 1.-ksq
      alpha = 1.
!--------------------------------------------------------------------------
!--  To avoid numerical truncation                                       --
!--------------------------------------------------------------------------
      if (kpsq < 1.0e-30) kpsq = 1.0e-30
      beta = sqrt(kpsq)
      if (beta < 1.0e-30) beta = 1.0e-10
      if (cpsq < 1.0e-30) cpsq = 1.0e-10
      delta = cpsq/beta
      epsi = csq/cpsq
      zeta = 0.
      sinf = 0.
      sa = .25
!
  100 continue
      sa = 2.*sa
      ambsq = (alpha-beta)*(alpha-beta)
      sinf = sinf+sa*ambsq
      alphat = alpha
      epsit = epsi
      alpha = .5*(alpha+beta)
      beta = sqrt(alphat*beta)
      epsi = (delta*epsi+zeta)/(1.+delta)
      delta = beta/4./alpha*(2.+delta+1./delta)
      zeta = .5*(epsit+zeta)
      if (abs(delta-1.) > err) go to 100
      if (ambsq > 1.e-14) go to 100
      cay = rpi2/alpha
      pik = cay*zeta
      ek = .5*cay*(ksq+sinf)
      msl = rh*dzsq*(r1*ek-drsq*pik/r1)
      if (csq-1.) 290,190,290
  190 msl = rh*dzsq*(r1*ek-cay*fr/r1*.5)
  290 continue
!
      mut = msl+ut*fr*r1*(cay-t*ek)
      sol = sol+sign*mut
      enddo 
      enddo 
!
      return
      end subroutine soleno

