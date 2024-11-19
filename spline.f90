      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      use constants,only:p_
      implicit none
      INTEGER n,NMAX
      REAL(p_) yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=800)
      INTEGER i,k
      REAL(p_) p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+&
     &1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
     &u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      use constants,only:p_
      implicit none
      INTEGER n
      REAL(p_) x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL(p_) a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** &
     &2)/6.
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


            SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      use constants,only:p_
      implicit none
      INTEGER m,n,NN
      REAL(p_) x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=200)
!CU    USES spline
      INTEGER j,k
      REAL(p_) y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      use constants,only:p_
      implicit none
      INTEGER m,n,NN
      REAL(p_) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=200)
!CU    USES spline,splint
      INTEGER j,k
      REAL(p_) y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
