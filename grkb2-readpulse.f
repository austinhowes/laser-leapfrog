      PROGRAM HYD1
C**********************************************************************
C  Nonstationary processes with hydrogen atom in pulsed laser fields  *
C  version 1.1  Drake, April-May 2005                                 *
C**********************************************************************
c      USE PORTLIB
      implicit double precision (a-h,o-z)
      integer xx, tt
      double precision ss, norm
      parameter(nxmax=381000,ntmax=100000,ncmax=61,nprmax=20000,
     1  nfmax=10)
      complex*16 czero,chalf,cone,ctwo,cthree,cfour,ci 
      complex*16 ch,cst,cst1,cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      complex*16 a,b,c,r,gam,u,cmn,cpl,a1,b1,c1,r1,gam1,u1,q,v,ff,bb
      complex*16 bet,xj
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1      ngob1,nerg
      common/ccnst1/czero,chalf,cone,ctwo,cthree,cfour,ci,ch,cst,cst1,
     1    cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/carry1/a(1:nxmax),b(1:nxmax),c(1:nxmax),r(1:nxmax),
     1         gam(1:nxmax),u(1:nxmax),cmn(0:ncmax),cpl(0:ncmax)
      common/carry2/a1(1:ncmax+1),b1(1:ncmax+1),c1(1:ncmax+1),
     1        r1(1:ncmax+1),gam1(1:ncmax+1),u1(1:ncmax+1)
      common/carry3/q(0:nxmax,0:ncmax),v(0:nxmax,0:ncmax),
     1               ff(0:nxmax,0:ncmax),bb(0:nxmax,0:ncmax)
      common/arry1/x(0:nxmax),g(0:nxmax),f(0:ntmax),fd(0:nxmax,1:nfmax),
     2             cek3(0:ncmax,0:ncmax,0:2*ncmax)      
      common/keys/ key1,key2,key3,key4,key5
      common/in1/w1,e1,s1,ton1,toff1,del1,w2,e2,s2,ton2,toff2,del2,agbr 
      common/discr/nf,nn(nfmax),ll(nfmax)
      common/in2/emin,de,enelec

      open(10,file='hyd.out',status='unknown')
      open(99,file='betas.out',status='unknown')
      open(40,file='read.me',status='unknown')
      open(50,file='hyd7.inp',status='old')
      open(30,file='mathem.scr',status='unknown')
      open(91,file='norm.out', status='unknown')
      open(92,file='field.out',status='old')
      open(93,file='wavefn.out',status='unknown')
      open(94,file='cross.out',status='unknown')

      ! leapfrog testing files
      open(95, file="leapfrogtesting.out", status="unknown")
      open(96, file="initial.out", status="unknown")

C**** READ INITIAL DATA
      read(50,*) dt,h
      read(50,*) nx, nprint, nc, key1, key2, key3, key4, key5
         if(key1.eq.0) goto 1
      read(50,*) nf
      read(50,*) (nn(i),ll(i), i=1,nf)
1        continue
      if(key4.eq.1) read(50,*) emin,de,emax
      nerg = idint((emax-emin)/de) + 1
      read(50,*) npulse
      if(npulse.eq.0) goto 2
      read(50,*) nshape
      read(50,*) e1,w1,fwhm1,plato1,del1
      if(npulse.eq.2) read(50,*) e2,w2,fwhm2,plato2,del2,delay
2     continue
      read(50,*) gbr, agbr
      close(50)
C**************************************************
      call cnstnt   ! definition of constants
      call arrays   ! initiate arrays
      call pulse    ! laser pulse
      call radial   ! radial: functions, potentials, gobbler
      call inprint  ! short printout of initial data

C**** SET COEFFICIENTS a,b,c AND INPUT FOR SOLVING DE
      do 15 i=1,nx-1
           a(i) = -cst2
           c(i) = -cst2  
15    continue
      a(1) = czero
      c(nx-1) = czero
      do 19 i=0,nx
         q(i,0) = dcmplx(g(i),0.d0)
19    continue
      do 20 i=1,nx-1
         do 16 j=0,nc
             bb(i,j) = cone + cst3*v(i,j) + cst1
             ff(i,j) = cone - cst3*v(i,j) - cst1 
16       continue
20    continue
      do 26 jj=0,nc
          q(nx,jj) = czero
          q(0,jj) = czero
26    continue


      a1(1) = czero
      c1(nc2) = czero
      do 52 jj=1,nc2
          b1(jj) = cone
52    continue

      do k = 1, nx
         write(96,*) k * h, real(q(k,0)), imag(q(k,0)),
     >        real(q(k,0))**2 + imag(q(k,0))**2
      end do

      kkk=0

C  MAIN TIME LOOP
*******************************************************
      if(nc.eq.0)  ntfin=-1
      if(ntfin.eq.0) ntfin=-1

!     calculate second time step
      q = q * exp(ci * dt / 2)
      ss = dt / h ** 2

      do tt = 3, ntmax

!     leapfrog propagation
         do xx = 1, nx - 1
            psi_r = real(q(xx,0))
            psi_im = imag(q(xx,0))
            if (mod(tt, 2) .eq. 1) then
               ! update the real part
               psi_r = psi_r - ss * imag(q(xx+1,0))
     >              -  ss * imag(q(xx-1,0))
     >              + 2 * (ss + real(v(xx,0))*dt) * psi_im
            else
               ! update the imaginary part
               psi_im = psi_im + ss * real(q(xx+1,0))
     >              + ss * real(q(xx-1,0))
     >              - 2 * (ss + real(v(xx,0))*dt) * psi_r
            end if
            q(xx,0) = cmplx(psi_r, psi_im)
         end do
         norm = 0.d0
         do k = 1, nx
            norm = norm + h * real(q(k,0))**2
         end do
         write(91,*) tt, norm
      end do

      ! output final wavefunction
      do k = 1, nx
         write(95,*) k * h, real(q(k,0)), imag(q(k,0)),
     >        real(q(k,0))**2 + imag(q(k,0))**2
      end do
     
      stop

C*****  PROPAGATION FROM DIAGONAL TERM BY tau = dt/2
      do 17 jj=0,nc
          do 18 i=1,nx-1  
             r(i) = q(i,jj)*ff(i,jj) + cst2*(q(i+1,jj)+q(i-1,jj))
             b(i) = bb(i,jj)
18        continue
C         if(b(1).eq.0.)pause 'tridag: rewrite equations'
          bet=b(1)
          u(1)=r(1)/bet
          do 111 j=2,nx-1
             gam(j)=c(j-1)/bet
             bet=b(j)-a(j)*gam(j)
C            if(bet.eq.0.)pause 'tridag failed'
             u(j)=(r(j)-a(j)*u(j-1))/bet
111       continue
          do 112 j=nx-2,1,-1
             u(j)=u(j)-gam(j+1)*u(j+1)
112       continue
          do 113 j=1,nx-1
             q(j,jj) = u(j)
113       continue
17    continue
C*****end of propagation by dt/2

C*****PROPAGATION FROM PULSE BY dt
      do 51 j=1,nx-1
         xj = dcmplx(x(j)*f(k),0.d0) 
         r1(1) = q(j,0) - cpl(0)*xj*q(j,1)
         r1(nc2) = q(j,nc)  - cmn(nc)*xj*q(j,nc1)
         do 58 jj=1,nc1
            r1(jj+1) = q(j,jj)
     1               - xj*(cmn(jj)*q(j,jj-1)+cpl(jj)*q(j,jj+1))
58       continue
         do 56 jj=2,nc2
            a1(jj) = cmn(jj-1)*xj
56       continue
         do 55 jj=1,nc
            c1(jj) = cpl(jj-1)*xj
55       continue

C         if(b1(1).eq.0.)pause 'tridag: rewrite equations'
         bet=b1(1)
         u1(1)=r1(1)/bet
         do 211 jj=2,nc2
            gam1(jj)=c1(jj-1)/bet
            bet=b1(jj)-a1(jj)*gam1(jj)
C            if(bet.eq.0.)pause 'tridag failed'
            u1(jj)=(r1(jj)-a1(jj)*u1(jj-1))/bet
211      continue
         do 212 jj=nc,1,-1
            u1(jj)=u1(jj)-gam1(jj+1)*u1(jj+1)
212      continue
         do 213 jj=0,nc
            q(j,jj) = u1(jj+1)
213      continue
51    continue
C*****end propagation by dt
c      do 3111 i=0,nx
c       write(10,1112) i,x(i),g(i),q(i,0),q(i,1),q(i,2)
c3111  continue

C***** SECOND PROPAGATION FROM DIAGONAL TERM BY tau = dt/2
      do 67 jj=0,nc
          do 68 i=1,nx-1  
              r(i) = q(i,jj)*ff(i,jj) + cst2*(q(i+1,jj)+q(i-1,jj))
              b(i) = bb(i,jj)
68        continue
C         if(b(1).eq.0.) pause 'tridag: rewrite equations'
          bet=b(1)
          u(1)=r(1)/bet
          do 311 j=2,nx-1
              gam(j)=c(j-1)/bet
              bet=b(j)-a(j)*gam(j)
C             if(bet.eq.0.)pause 'tridag failed'
              u(j)=(r(j)-a(j)*u(j-1))/bet
311       continue
          do 312 j=nx-2,1,-1
              u(j)=u(j)-gam(j+1)*u(j+1)
312       continue
          do 313 j=1,nx-1
              q(j,jj) = u(j)
313       continue
67    continue
C*****end of second diagonal propagation

C***PRINTOUT FOR TIME LOOP WITH NUMBER nprint*N (N=1,2,...)
      if(kkk.eq.nprint) then
          call output(k)
          kkk=0
      endif
100   format(10e12.5)
C      write(*,*) k
5     continue
C*****end of propagation
*kb   save the wavefunction
      do 137 jj=0,nc
       write(93,9393) (q(i,jj),i=0,nx)
 137  continue
 9393 format(1p10e14.6)
*kb
      if(key4.eq.1) call dstrm  ! photoelectron spectral distribution
      end

C********************************************************
C********************************************************
 
      subroutine output(k)
      implicit double precision (a-h,o-z)
      character*3 zchar
      character*8 probnm,curren
      character*19 genr 
      character*27 flpath,flpat1
      parameter
     1(nxmax=381000,ncmax=61,nprmax=20000,ntmax=100000,
     1 nfmax=10)
      complex*16 czero,chalf,cone,ctwo,cthree,cfour,ci 
      complex*16 ch,cst,cst1,cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      complex*16 q,v,ff,bb,psi
      common/ccnst1/czero,chalf,cone,ctwo,cthree,cfour,ci,ch,cst,cst1,
     1    cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1      ngob1,nerg
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/carry3/q(0:nxmax,0:ncmax),v(0:nxmax,0:ncmax),
     1               ff(0:nxmax,0:ncmax),bb(0:nxmax,0:ncmax)
      common/arry1/x(0:nxmax),g(0:nxmax),f(0:ntmax),fd(0:nxmax,1:nfmax),
     2             cek3(0:ncmax,0:ncmax,0:2*ncmax)      
      common/keys/ key1,key2,key3,key4,key5
      common/discr/nf,nn(nfmax),ll(nfmax)
      dimension anorm(1:nxmax+1),
     1     auto(1:nxmax+1),conv(0:nxmax),
     1     vrlp(0:nfmax)

      k1=k+1
      write(10,200) k1

C**** AUTOCORRELATION FUNCTION
      if(key2.eq.0) goto 231
      do 232 i=1,ngob1
            auto(i) = real(q(i-1,0),8)*g(i-1)
232       continue
      call arsimd(ngob1,h,auto,are)
      do 233 i=1,ngob1
            auto(i) = dimag(q(i-1,0))*g(i-1)
233       continue
      call arsimd(ngob1,h,auto,aie)
      vrlp(0) = are**2 + aie**2
      write(10,300) are,aie,vrlp(0)
300   format(/,' Autocorr. function: Re =',
     1   f11.8 ,' Im =',f11.8,'  Square = ',f10.8)
231   continue

C**** POPULATION OF DISCRETE STATES FOR key1 neq 0
      if(key1.eq.0) goto 301
      do 249 kd=1,nf
          do 240 j=1,ngob1
              auto(j)=0.d0
240       continue
          do 248 i=1,ngob1                    
             auto(i) = real(q(i-1,ll(kd)),8)*fd(i-1,kd)
248       continue
          call arsimd(ngob1,h,auto,are)
          do 243 i=1,ngob1                    
             auto(i) = dimag(q(i-1,ll(kd)))*fd(i-1,kd)
243       continue
          call arsimd(ngob1,h,auto,aie)
          vrlp(kd) = are**2 + aie**2
          write(10,500) kd,nn(kd),ll(kd),are,aie,vrlp(kd)
249   continue           
500   format(/,' Overlap with state',i2,' : n =',i3,'  l =',i2,
     1 ' Re =',d12.4 ,' Im =',d12.4,'  Square = ',d12.4)
301   continue

C*** MONITOR CONVERGENCY WITH RESPECT TO PARTIAL WAVES AND NORMALIZATION
      if(key3.eq.0) goto 241
      cnorm=0.d0
      do 51 j=0,nc
         do 52 i=1,ngob1
           conv(i)=abs(q(i-1,j))**2
52       continue
      call arsimd(ngob1,h,conv,awh)
      write(10,700) j,awh
      cnorm = cnorm + awh
51    continue
700   format(/,' channel = ',i3,'   contribution = ',e13.4)      
241   continue

C*** CONTROL OF NORMALIZATION
      if(key3.ne.0) goto 242
      do 21 i=1,ngob1
         anorm(i) = 0.d0
         do 22 j=0,nc
            anorm(i) = anorm(i) + abs(q(i-1,j))**2
22       continue
21    continue
      call arsimd(ngob1,h,anorm,cnorm)
242   continue
      write(10,400) cnorm
      write(91,401) dble(k)*dt, cnorm
401   format(f10.3,f11.8)
400   format(/,' norm = ',f10.8)

100        format(8e13.5)
200        format(//,'after the time step number = ',I7)

      return
      end

C****************************************************
      subroutine inprint
      implicit double precision (a-h,o-z)
      parameter(nxmax=381000,ntmax=100000,ncmax=61,nprmax=20000,
     1  nfmax=10)
      complex*16 q,v,ff,bb
      common/arry1/x(0:nxmax),g(0:nxmax),f(0:ntmax),fd(0:nxmax,1:nfmax),
     2             cek3(0:ncmax,0:ncmax,0:2*ncmax)
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1     ngob1,nerg
      common/carry3/q(0:nxmax,0:ncmax),v(0:nxmax,0:ncmax),
     1               ff(0:nxmax,0:ncmax),bb(0:nxmax,0:ncmax)
      common/keys/ key1,key2,key3,key4,key5
      common/in1/w1,e1,s1,ton1,toff1,del1,w2,e2,s2,ton2,toff2,del2,agbr 
      dimension anrm(1:nxmax)

      write(40,901)
901   format(' Calculation parameters',/,' Pulse parameters:',/)
      write(40,902) w1,e1,s1,ton1,toff1,del1
902   format(' w1   = ',f6.4,' au',/,' e1   = ',e10.4,' au',/,' s1   =',
     1 f10.3,' au',/,' ton1   =',f10.3,' au',/,' toff1  =',f10.3,' au',
     2  /,' del1 = ',e10.4,/) 
      write(40,9002) w2,e2,s2,ton2,toff2,del2
9002  format(' w2   = ',f6.4,' au',/,' e2   = ',e10.4,' au',/,' s2   =',
     1 f10.3,' au',/,' ton2   =',f10.3,' au',/,' toff2  =',f10.3,' au',
     2  /,' del2 = ',e10.4,/) 
      write(40,903) nc
903   format(' Number of partial waves (number of els): nc = ',i4)
      write(40,904) nt, nx, dt, h
904   format(/,' Mesh parameters:',/, ' nt  =',i8,/,' nx  =',i8,/,
     1       ' dt  =',f8.4,' au',/,' h   =',f8.4,' au',/)
      write(40,905) gbr,agbr,ngob1
905   format(/,' Gobbler parameters',//,
     1   ' gbr   = ',f7.2,/,' agbr  =',f7.2,' au',/,' ngob1 =',i7,/)
      write(40,906) nprint,ntfin
906   format(/,' nprint = ',i6,' ntfin = ',i6,/)
      do 21 i=1,ngob1
         anrm(i) = g(i-1)**2
21    continue
      call arsimd(ngob1,h,anrm,cn)
      write(40,907) cn
907   format(/,' Initial norm = ', e15.8,/)
      return
      end

C****************************************************
 
      SUBROUTINE FACTOR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IFAK=500)
      COMMON /FACT / FACT(IFAK)
      FACT(1) = 0.d0
      DO 10 I=1,IFAK-1
       FACT(I+1) = FACT(I) + DLOG(real(I,8))
10    CONTINUE
      RETURN
      END

************************************************************************
      DOUBLE PRECISION FUNCTION CLB(J1,J2,J3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0,IFAK=500)
      COMMON /FACT / FACT(IFAK)
      N = ITRI(J1,J2,J3)
      IF(N.EQ.0) GOTO 10
      L = J1+J2+J3
      M = L/2*2
      IF(M.NE.L) GO TO 10
      CLB = ONE
      B = ONE
      K = (J1+J2-J3)/2
      IF(K.NE.K/2*2) B=-ONE
      AA3 = DBLE(J3)
      KG = (J1+J2+J3)/2
      KG1 = KG+1
      KG2 = 2*KG + 1
      CLB = B * DEXP ( HALF*DLOG(TWO*AA3+ONE) + FACT(KG1) -
     *  FACT(KG1-J1) - FACT(KG1-J2) - FACT(KG1-J3) +
     *  HALF*(FACT(KG2-2*J1) + FACT(KG2-2*J2) + FACT(KG2-2*J3) -
     *  FACT(KG2+1)) )
      RETURN
10    CLB=ZERO
      RETURN
      END

*********************************************
      FUNCTION ITRI(J1,J2,J3)
      ITRI=0
      IF(J1+J2.LT.J3) GOTO 10
      IF(J1+J3.LT.J2) GOTO 10
      IF(J2+J3.LT.J1) GOTO 10
      ITRI=1
10    RETURN
      END

*************************************
      SUBROUTINE ARSIMD(N,DEL,A,R)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N)
      L=N
      SUM=A(1)-A(L)
      DO 1 I=2,L,2
1     SUM=SUM+4.D0*A(I)+2.D0*A(I+1)
      R=(DEL*SUM)/3.D0
      RETURN
      END

**************************************
      subroutine cnstnt
      implicit double precision (a-h,o-z)
      complex*16 czero,chalf,cone,ctwo,cthree,cfour,ci 
      complex*16 ch,cst,cst1,cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/ccnst1/czero,chalf,cone,ctwo,cthree,cfour,ci,ch,cst,cst1,
     1    cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1          ngob1,nerg

      pi    = dacos(-1.d0)
      iz    = 100
      nbf   = nx-1
      nx1   = nx+1
      nc1   = nc-1
      nc2   = nc+1
      czero = (0.d0,0.d0)
      cone  = (1.d0,0.d0)
      ctwo  = (2.d0,0.d0)
      cthree= (3.d0,0.d0)
      cfour = (4.d0,0.d0)
      chalf = (.5d0,0.d0)
      ci    = (0.d0,1.d0)
      c27   = (27.d0,0.d0)
      tau   = dt/2.d0
      zet   = 1.d0
      ch    = dcmplx(h,0.d0)
      cst   = ci*dcmplx(tau,0.d0)
      cst1  = cst/(ctwo*ch*ch)
      cst2  = cst1/ctwo
      cst3  = cst/ctwo
      cst4  = ci*dcmplx(dt,0.d0)
      cst5  = cst/(ch*ch)
      h2    = ctwo*ch
      h3    = cthree*ch
      h24   = cone/(dcmplx(24.d0,0.d0)*ch)
      h1d5  = (1.5d0,0.d0)*ch
      hd25  = (.25d0,0.d0)*ch

      return
      end      

*************************************
      subroutine arrays
      implicit double precision (a-h,o-z)
      parameter(nxmax=381000,ntmax=100000,ncmax=61,nprmax=20000,
     1 nfmax=10)
      complex*16 czero,chalf,cone,ctwo,cthree,cfour,ci,ck
      complex*16 ch,cst,cst1,cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      complex*16 a,b,c,r,gam,u,cmn,cpl,a1,b1,c1,r1,gam1,u1,q,v,ff,bb
      common/ccnst1/czero,chalf,cone,ctwo,cthree,cfour,ci,ch,cst,cst1,
     1    cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1   ngob1,nerg
      common/carry1/a(1:nxmax),b(1:nxmax),c(1:nxmax),r(1:nxmax),
     1         gam(1:nxmax),u(1:nxmax),cmn(0:ncmax),cpl(0:ncmax)
      common/carry2/a1(1:ncmax+1),b1(1:ncmax+1),c1(1:ncmax+1),
     1        r1(1:ncmax+1),gam1(1:ncmax+1),u1(1:ncmax+1)
      common/carry3/q(0:nxmax,0:ncmax),v(0:nxmax,0:ncmax),
     1              ff(0:nxmax,0:ncmax),bb(0:nxmax,0:ncmax)
      common/arry1/x(0:nxmax),g(0:nxmax),f(0:ntmax),fd(0:nxmax,1:nfmax),
     2             cek3(0:ncmax,0:ncmax,0:2*ncmax)

C  SET INITIAL ZEROS
      do 24 i=1,nxmax
          a(i)=czero
          b(i)=czero
          c(i)=czero
          r(i)=czero
          gam(i)=czero
          u(i)=czero
24    continue
      do 22 i=0,nxmax
          x(i)=0.d0
          g(i)=0.d0
          do 21 j=0,nc
             q(i,j)  = czero
             v(i,j)  = czero
             ff(i,j) = czero
             bb(i,j) = czero
21        continue
22    continue
      do 34 jj=1,ncmax+1
          a1(jj)=czero
          b1(jj)=czero
          c1(jj)=czero
          r1(jj)=czero
          gam1(jj)=czero
          u1(jj)=czero
34    continue
      do 35 k=0,nt
          f(k)=0.d0
35    continue

C***RADIAL MESH
      do 401 i=0,nx
        x(i)=h*real(i,8)
401   continue
C*** find ngob
      ngob=nx
      if(gbr.lt.1.d-10) goto 4
        do 3 i=0,nx
          if(abs(x(i)-gbr).lt.1.d-10) then
            ngob=i
            goto 4
          endif
3     continue
4     continue
      if(ngob.eq.0) ngob=nx
      write(10,*) 'ngob =  ', ngob
      ngob1=ngob
      if(mod(ngob,2).eq.0) ngob1=ngob+1  ! makes ngob1 odd

C  SET MIXING COEFFICIENTS (multiplied by i*tau)
      cmn(0) = czero
      cpl(0) = cone/cdsqrt(cthree) * cst
      do 14 i=1,nc
        ck=dcmplx(real(i,8),0.d0)
        cmn(i) = ck/cdsqrt((ctwo*ck+cone)*(ctwo*ck-cone))*cst
        cpl(i) = (ck+cone)/cdsqrt((ctwo*ck+cone)*(ctwo*ck+cthree))*cst
14    continue
      cpl(nc) = czero

C log OF FACTORIALS
      call factor

C**** GEOMETRICAL PART OF MULTIPOLE COEFFICIENTS AND nbmax
      do 1 n1=0,nc                 
        cek1 = dsqrt(2.d0*real(n1,8)+1.d0)
        do 2 n2 = 0,nc
          cek2 = dsqrt(2.d0*real(n2,8)+1.d0)*cek1
          do 13 nbig = iabs(n1-n2), n1+n2, 2
            if(nbig.gt.nbmax) nbmax=nbig
            cek3(n1,n2,nbig) = cek2*clb(n1,n2,nbig)**2/(4.d0*pi)
13         continue
2       continue
1     continue


      return
      end

      subroutine dstrm
C***energy distribution in the outgoing state
C*********************************************
      implicit double precision (a-h,o-z)
      complex*16 tint,q,v,ff,bb,cdst,be,ct1
      parameter(nxmax=381000,ntmax=100000,ncmax=61,nprmax=20000,
     1 nemax=15000,nfmax=10)
      common/in2/emin,de,enelec
      common/carry3/q(0:nxmax,0:ncmax),v(0:nxmax,0:ncmax),
     1              ff(0:nxmax,0:ncmax),bb(0:nxmax,0:ncmax)
      common/arry1/x(0:nxmax),g(0:nxmax),f(0:ntmax),fd(0:nxmax,1:nfmax),
     2             cek3(0:ncmax,0:ncmax,0:2*ncmax)
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1   ngob1,nerg
      common/adis/cdst(0:ncmax)
      dimension tint(0:ncmax,0:nemax),ph(0:ncmax),bet(0:2*ncmax)
      dimension dm(0:nemax),ener(0:nemax),wr(nxmax),wi(nxmax),
     1   ur(0:nxmax,0:ncmax),ui(0:nxmax,0:ncmax),scpart(0:ncmax),
     2   fcp(0:200),gc(0:200),gcp(0:200),fc(0:200)
      dimension fcstore(0:200)
      ncu=nc
      do 9 nen=0,nerg
         dm(nen)=0.d0
9     continue
      ener(0) = emin
         do 8 jn=0,nc
            ur(0,jn)=0.d0
            ui(0,jn)=0.d0
8        continue

      write(10,300)
300   format(/,3x,'E(au)',6x,'Spectr    bet0    bet1     bet2     bet3',
     1 '     bet4     bet5     bet6     bet7     bet8     bet9',
     2 '    bet10    bet11    bet12    bet13    bet14')
      write(94,200)
200   format(/,3x,'# E(au)',6x,'Spectr    crs0    crs1     crs2    crs3',
     1 '     crs4     crs5     crs6     crs7     crs8     crs9',
     2 '    crs10    crs11    crs12    crs13    crs14')
      do 3 nen = 1,nerg
         ener(nen)=ener(nen-1) + de
         p = dsqrt(2.d0 * ener(nen))
         eta = -1.d0/p
         do 5 i=ngob,1,-1
             pr = p * x(i) 
                do 4 kn=0,ncu
                   fc(kn)  = 0.d0
                   fcp(kn) = 0.d0
                   gc(kn)  = 0.d0
                   gcp(kn) = 0.d0
4               continue
                call coul90(pr,ETA,0.d0,ncu,FC,FCP,GC,GCP,0,ifail)
c             call RCWFN(pr,eta,0,ncu+1,fc,fcp,gc,gcp,1.d-15,1.d3)
            if (ifail.eq.0) then
             do 7 jj=0,nc
                ur(i,jj)=fc(jj)*real(q(i,jj),8)/p
                ui(i,jj)=fc(jj)*dimag(q(i,jj))/p
7            continue
             prstore = pr
             fcstore(jj) = fc(jj)/prstore**(jj+1)
            else
              do 787 jj=0,nc
                fc(jj) = fcstore(jj)*pr**(jj+1)
                if (abs(fc(jj)).le.1.0d-80) fc(jj) = 0.0d0 
                ur(i,jj)=fc(jj)*real(q(i,jj),8)/p
                ui(i,jj)=fc(jj)*dimag(q(i,jj))/p
787           continue            
            endif
5        continue
         do 17 jj=0,nc
            ph(jj) = facouz(2.d0*ener(nen),jj,zet)
            do 18 i=0,ngob
              wr(i+1)=ur(i,jj)
              wi(i+1)=ui(i,jj)
18          continue
            call arsimd(ngob1,h,wr,rtr)
            call arsimd(ngob1,h,wi,rti)
            tint(jj,nen) = dcmplx(rtr,rti)
            dm(nen) = dm(nen) + cdabs(tint(jj,nen))**2
17       continue
c---- beta parameters
      do 13 lam=0,2*nc
        be = (0.d0,0.d0)
        do 11 j1 = 0,nc
        do 12 j2 = 0,nc         
         aw = dsqrt((2.d0*dble(j1)+1.d0)/(2.d0*dble(j2)+1.d0))
     1  *clb(j1,lam,j2)**2    
         rph=ph(j1)-ph(j2)
         ct1 = (0.d0,1.d0)**(j2-j1) * cdexp(dcmplx(0.d0,rph)) 
     1    * tint(j1,nen)*dconjg(tint(j2,nen))*dcmplx(aw,0.d0)
         be=be+ct1
12      continue
11      continue
        bet(lam) = real(be,8)*(2.d0*dble(lam)+1.d0)/dm(nen)  
13    continue 
      do 14 j3=0,nc
        scpart(j3)=cdabs(tint(j3,nen))**2*2.d0/pi*dsqrt(2.d0*ener(nen))
14    continue
        dcr = dm(nen)*2.d0/pi * dsqrt(2.d0*ener(nen))
        write(10,100) ener(nen),dcr,(bet(nn),nn=0,min(2*nc,100))
        write(99,100) ener(nen),dcr,(bet(nn),nn=0,min(2*nc,200))
        write(94,100) ener(nen),dcr,(scpart(nn),nn=0,min(nc,200))
3     continue
100   format(f9.5,e10.3,201f9.5)
      return
      end


C***********************************************************
      DOUBLE PRECISION FUNCTION FACOUZ(E,L,Z)
      implicit double precision (a-h,o-z)
      GAM=-Z/DSQRT(E)
      M=201-L
      DO 19 K=1,M
      AL=DATAN(GAM/(202-K))
      IF(K.EQ.1) GO TO 18
      FACOUZ=FACOUZ-AL
      GO TO 19
18    BE=DSQRT(GAM*GAM+(202-K)**2)
      FACOUZ=AL*200.5d0+GAM*(DLOG(BE)-1.d0)
     *+(-SIN(AL)/12.0d0+SIN(3.d0*AL)/(360.d0*BE*BE))/BE
19    CONTINUE
      RETURN
      END

      SUBROUTINE PARINV(X,A,F,N,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),F(N)
      IF(X.LT.A(1)) GO TO 11
      IF(X.GT.A(N)) GO TO 4
      K1=1
      K2=N
2     K3=K2-K1
      IF(K3.LE.1) GO TO 6
      K3=K1+K3/2
      IF(A(K3)-X) 7,8,9
7     K1=K3
      GO TO 2
9     K2=K3
      GO TO 2
8     R=F(K3)
      RETURN
3     B1=A(K1)
      B2=A(K1+1)
      B3=A(K1+2)
      B4=F(K1)
      B5=F(K1+1)
      B6=F(K1+2)
      R=B4*((X-B2)*(X-B3))/((B1-B2)*(B1-B3))+B5*((X-B1)*(X-B3))/
     *((B2-B1)*(B2-B3))+B6*((X-B1)*(X-B2))/((B3-B1)*(B3-B2))
      RETURN
6     IF(K2.NE.N) GO TO 3
      K1=N-2
      GO TO 3
4     C=DABS(X-A(N))
      IF(C.LT.0.1D-7) GO TO 5
      K1=N-2
13    WRITE(*,41) X
41    FORMAT(25H X IS OUT OF THE INTERVAL,3H X=,F13.7)
      GO TO 3
5     R=F(N)
      RETURN
11    C=DABS(X-A(1))
      IF(C.LT.0.1D-7) GO TO 12
      K1=1
      GO TO 13
12    R=F(1)
      RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE COUL90(X, ETA, XLMIN,LRANGE, FC,GC,FCP,GCP, KFN,IFAIL)
C----------------------------------------------------------------------
C
C  COULOMB & BESSEL FUNCTION PROGRAM-- COUL90 -- USING STEED'S METHOD
C
C  COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = (D/DX) F, GCP = (D/DX) G
C   FOR REAL X .GT. 0. ,REAL ETA (INCLUDING 0.), AND REAL XLMIN .GT.-1.
C   FOR (LRANGE+1) INTEGER-SPACED LAMBDA VALUES.
C   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
C   EQUATION, TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
C   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
C   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
C----------------------------------------------------------------------
C   CALLING VARIABLES; ALL REALS ARE DOUBLE PRECISION (REAL*8)
C
C   X       - REAL ARGUMENT FOR COULOMB FUNCTIONS > 0.0 
C             [ X > SQRT(ACCUR) : ACCUR IS TARGET ACCURACY 1.0D-14 ]
C   ETA     - REAL SOMMERFELD PARAMETER, UNRESTRICTED > = < 0.0
C   XLMIN   - REAL MINIMUM LAMBDA-VALUE (L-VALUE OR ORDER),
C             GENERALLY IN RANGE 0.0 - 1.0 AND MOST USUALLY 0.0
C   LRANGE  - INTEGER NUMBER OF ADDITIONAL L-VALUES : RESULTS RETURNED
C             FOR L-VALUES XLMIN TO XLMIN + LRANGE INCLUSIVE
C   FC ,GC  - REAL VECTORS F,G OF REGULAR, IRREGULAR COULOMB FUNCTIONS
C   FCP,GCP - REAL VECTORS FOR THE X-DERIVATIVES OF  F,G
C             THESE VECTORS TO BE OF LENGTH AT LEAST MINL + LRANGE
C             STARTING ELEMENT MINL = MAX0( IDINT(XLMIN+ACCUR),0 )
C   KFN     - INTEGER CHOICE OF FUNCTIONS TO BE COMPUTED :
C           = 0         REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
C           = 1    SPHERICAL BESSEL      "      "     "        j & y
C           = 2  CYLINDRICAL BESSEL      "      "     "        J & Y
C
C   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
C   IN OSCILLATING REGION X .GE. [ETA + SQRT{ETA**2 + XLM*(XLM+1)}]
C   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8 IS
C   THE SMALLEST NUMBER WITH 1.+ACC8.NE.1. FOR OUR WORKING PRECISION.
C   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
C   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
C   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
C   THE VARIABLE PACCQ IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
C----------------------------------------------------------------------
C   ERROR RETURNS                THE USER SHOULD TEST IFAIL ON EXIT
C
C   IFAIL ON INPUT IS SET TO 0                        LIMIT = 20000
C   IFAIL IN OUTPUT =  0 : CALCULATIONS SATISFACTORY
C                   =  1 : CF1 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C                   =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C                   = -1 : X < 1D-7 = SQRT(ACCUR)
C                   = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES) 
C----------------------------------------------------------------------
C  MACHINE-DEPENDENT PARAMETERS:    ACCUR - SEE ABOVE
C           SMALL - OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
C           IE 1D-30 FOR IBM REAL*8,    1D-150 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
C  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
C  ORIGINAL PROGRAM  RCWFN       IN    CPC  8 (1974) 377-395
C                 +  RCWFF       IN    CPC 11 (1976) 141-142
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
C  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
C  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188         
C  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
C  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
C  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
C  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509     1.4.94
C----------------------------------------------------------------------
C  AUTHOR: A. R. BARNETT           MANCHESTER  MARCH   1981/95
C                                  AUCKLAND    MARCH   1991
C----------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          LRANGE, KFN, IFAIL
      DOUBLE PRECISION X, ETA, XLMIN
      DOUBLE PRECISION FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
C----- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM MINL
      DOUBLE PRECISION ACCUR,ACCH,SMALL, ONE,ZERO,HALF,TWO,TEN2, RT2DPI
      DOUBLE PRECISION XINV,PK,CF1,C,D,PK1,ETAK,RK2,TK,DCF1,DEN,XLM,XLL
      DOUBLE PRECISION EL,XL,RL,SL, F,FCMAXL,FCMINL,GCMINL,OMEGA,WRONSK
      DOUBLE PRECISION WI, A,B, AR,AI,BR,BI,DR,DI,DP,DQ, ALPHA,BETA
      DOUBLE PRECISION E2MM1, FJWKB,GJWKB, P,Q,PACCQ, GAMMA,GAMMAI
      INTEGER          IEXP, NFP, NPQ, L, MINL,MAXL, LIMIT
      LOGICAL          ETANE0, XLTURN
      PARAMETER      ( LIMIT = 30000, SMALL = 1.0D-150 )
      COMMON  /STEED/  PACCQ,NFP,NPQ,IEXP,MINL    !not required in code
      COMMON  /DESET/  CF1,P,Q,F,GAMMA,WRONSK     !information only
C----------------------------------------------------------------------
C     COUL90 HAS CALLS TO: DSQRT,DABS,MAX0,IDINT,DSIGN,DFLOAT,DMIN1
C----------------------------------------------------------------------
      DATA ZERO,ONE,TWO,TEN2,HALF /0.0D0, 1.0D0, 2.0D0, 1.0D2, 0.5D0/
      DATA RT2DPI /0.79788 45608 02865  D0/ 
CQ    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
C-----THIS CONSTANT IS  DSQRT(TWO / PI):
C-----USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND DOUBLE PRECISION
C----------------CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
                        ACCUR = 1.0D-18
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO                                
      ACCH  = DSQRT(ACCUR)
C-----   TEST RANGE OF X, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
                IF( X .LE. ACCH )                GO TO 100
      IF( KFN.EQ.2 )   THEN
         XLM = XLMIN - HALF                                  
        ELSE
         XLM = XLMIN                                                     
        ENDIF
      IF( XLM.LE.-ONE .OR. LRANGE.LT.0 )         GO TO 105 
      E2MM1  = XLM * XLM + XLM
      XLTURN = X * (X -  TWO * ETA) .LT. E2MM1
      E2MM1  = E2MM1  +  ETA * ETA
      XLL    = XLM + DFLOAT(LRANGE)
C-----  LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C-----  XLL  IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
C-----  DETERMINE STARTING ARRAY ELEMENT (MINL) FROM XLMIN
      MINL  = MAX0( IDINT(XLMIN + ACCUR),0 )     ! index from 0
      MAXL  = MINL + LRANGE
C-----   EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
      XINV = ONE / X
      DEN  = ONE                       ! unnormalised F(MAXL,ETA,X)
      PK   = XLL + ONE
      CF1  = ETA / PK  +  PK * XINV                                             
           IF( DABS(CF1).LT.SMALL )    CF1 = SMALL
      RK2  = ONE
         D = ZERO
         C = CF1
C----- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA + 1: LENTZ-THOMPSON
      DO 10 L =  1 , LIMIT             ! abort if reach LIMIT (20000)    
          PK1 = PK + ONE
          IF( ETANE0 ) THEN
                ETAK = ETA / PK
                RK2  = ONE + ETAK * ETAK
                 TK  = (PK + PK1) * (XINV + ETAK / PK1)
             ELSE
                 TK  = (PK + PK1) * XINV
             ENDIF
          D   =  TK - RK2 * D          ! direct  ratio of B convergents    
          C   =  TK - RK2 / C          ! inverse ratio of A convergents
            IF( DABS(C).LT.SMALL ) C = SMALL
            IF( DABS(D).LT.SMALL ) D = SMALL
          D   = ONE / D
          DCF1=   D * C
          CF1 = CF1 * DCF1
              IF( D.LT.ZERO )    DEN = -DEN
          PK  = PK1
          IF( DABS(DCF1-ONE).LT.ACCUR )     GO TO  20 ! proper exit
   10 CONTINUE
                                            GO TO 110 ! error exit 
   20       NFP = PK - XLL - 1                        ! number of steps
              F = CF1                                 ! need DEN later
C----DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
      IF( LRANGE.GT.0 )       THEN
          FCMAXL    = SMALL  * DEN 
          FCP(MAXL) = FCMAXL * CF1
          FC (MAXL) = FCMAXL
                    XL = XLL                   
                    RL = ONE
          DO 30 L =  MAXL, MINL+1, -1
             IF( ETANE0 )  THEN
                    EL = ETA / XL                
                    RL = DSQRT( ONE + EL * EL )
                    SL = XL * XINV  + EL
                    GC (L) = RL                  ! storage
                    GCP(L) = SL
                ELSE
                    SL = XL * XINV
                ENDIF
             FC (L-1)  = ( FC(L)   * SL  +  FCP(L) ) / RL
             FCP(L-1)  =   FC(L-1) * SL  -  FC (L) * RL
             XL    =  XL - ONE                   ! end value is XLM
   30     CONTINUE
         IF( DABS(FC(MINL)).LT.ACCUR*SMALL )  FC(MINL) = ACCUR * SMALL
          F   = FCP(MINL) / FC(MINL)             ! F'/F at min L-value
          DEN = FC (MINL)                        ! normalisation
      ENDIF
C---------------------------------------------------------------------
C-----   NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C-----   EVALUATE CF2 = P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)
C---------------------------------------------------------------------
      IF( XLTURN ) CALL JWKB( X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP )
      IF( IEXP.GT.1 .OR. GJWKB.GT.(ONE / (ACCH*TEN2)) ) THEN
          OMEGA = FJWKB
          GAMMA = GJWKB * OMEGA
          P     = F
          Q     = ONE
        ELSE                                     ! find cf2     
          XLTURN = .FALSE.
          PK =  ZERO
          WI =  ETA + ETA
          P  =  ZERO
          Q  =  ONE - ETA * XINV
          AR = -E2MM1
          AI =  ETA
          BR =  TWO * (X - ETA)
          BI =  TWO
          DR =  BR / (BR * BR + BI * BI)
          DI = -BI / (BR * BR + BI * BI)
          DP = -XINV * (AR * DI + AI * DR)
          DQ =  XINV * (AR * DR - AI * DI)
          DO 40 L = 1, LIMIT
             P  = P  + DP
             Q  = Q  + DQ
             PK = PK + TWO
             AR = AR + PK
             AI = AI + WI                                                   
             BI = BI + TWO                                                  
             D  = AR * DR - AI * DI + BR                                        
             DI = AI * DR + AR * DI + BI                                        
             C  = ONE / (D * D + DI * DI)                  
             DR =  C * D                                                      
             DI = -C * DI                                                     
             A  = BR * DR - BI * DI - ONE                                       
             B  = BI * DR + BR * DI                                             
             C  = DP * A  - DQ * B
             DQ = DP * B  + DQ * A                                              
             DP = C
      IF( DABS(DP)+DABS(DQ).LT.(DABS(P)+DABS(Q)) * ACCUR ) GO TO 50
   40     CONTINUE
                                              GO TO 120 ! error exit
   50     NPQ   = PK / TWO                              ! proper exit
          PACCQ = HALF * ACCUR / DMIN1( DABS(Q),ONE )
          IF( DABS(P).GT.DABS(Q) ) PACCQ = PACCQ * DABS(P)
C---------------------------------------------------------------------
C    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
C---------------------------------------------------------------------
          GAMMA   = (F - P) / Q
          GAMMAI  = ONE / GAMMA
          IF( DABS(GAMMA) .LE. ONE )  THEN 
                 OMEGA  = DSQRT( ONE  +  GAMMA * GAMMA )
            ELSE
                 OMEGA  = DSQRT( ONE  +  GAMMAI* GAMMAI) * DABS(GAMMA)
            ENDIF 
          OMEGA  = ONE / ( OMEGA * DSQRT(Q) )
          WRONSK = OMEGA
        ENDIF   
C--------------------------------------------------------------------- 
C    RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C---------------------------------------------------------------------
      IF( KFN.EQ.1 )       THEN         !   spherical Bessel functions
                 ALPHA = XINV
                 BETA  = XINV
        ELSEIF( KFN.EQ.2 ) THEN         ! cylindrical Bessel functions
                 ALPHA = HALF * XINV
                 BETA  = DSQRT( XINV ) * RT2DPI
        ELSE                            ! kfn = 0,   Coulomb functions
                 ALPHA = ZERO     
                 BETA  = ONE
        ENDIF
      FCMINL = DSIGN( OMEGA,DEN ) * BETA
      IF( XLTURN )   THEN
                        GCMINL =   GJWKB * BETA
        ELSE
                        GCMINL =  FCMINL * GAMMA
        ENDIF
      IF( KFN.NE.0 )    GCMINL = -GCMINL         ! Bessel sign differs
      FC (MINL) = FCMINL
      GC (MINL) = GCMINL
      GCP(MINL) = GCMINL * (P - Q * GAMMAI - ALPHA) 
      FCP(MINL) = FCMINL * (F - ALPHA)
      IF( LRANGE.EQ.0 )                          RETURN
C---------------------------------------------------------------------
C    UPWARD RECURRENCE FROM GC(MINL),GCP(MINL) STORED VALUES ARE RL,SL
C    RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C      XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
C---------------------------------------------------------------------
      OMEGA = BETA * OMEGA / DABS(DEN)
                 XL = XLM
                 RL = ONE 
      DO 60  L = MINL+1 , MAXL                   ! indexed from 0
                 XL = XL + ONE
          IF( ETANE0 ) THEN
                 RL = GC (L)
                 SL = GCP(L)
            ELSE 
                 SL =  XL * XINV
            ENDIF
          GC (L)  = ( (SL - ALPHA) * GC(L-1) - GCP(L-1) ) / RL
          GCP(L)  =    RL *  GC(L-1)  -  (SL + ALPHA) * GC(L)
          FCP(L)  = OMEGA * ( FCP(L)  -  ALPHA * FC(L) )
          FC (L)  = OMEGA *   FC (L)
   60 CONTINUE
      RETURN
C------------------   ERROR MESSAGES
  100 IFAIL = -1
      WRITE(6,1000) X,ACCH
 1000 FORMAT(' FOR X = ',1PD12.3,'     TRY SMALL-X  SOLUTIONS,',
     *' OR X IS NEGATIVE'/ ,' SQUARE ROOT (ACCURACY) =  ',D12.3/)
                     RETURN
  105 IFAIL = -2                                                        
      WRITE (6,1005) LRANGE,XLMIN,XLM                                    
 1005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',    
     *I10,1P2D15.6/)                                                        
                     RETURN                                   
  110 IFAIL =  1                                                        
      WRITE (6,1010) LIMIT, CF1,DCF1, PK,ACCUR                              
 1010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',I10,' ITERATIONS',/ 
     *' CF1,DCF1,PK,ACCUR =  ',1P4D12.3/)                               
                     RETURN                                       
  120 IFAIL =  2                                                        
      WRITE (6,1020) LIMIT,P,Q,DP,DQ,ACCUR
 1020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',I7,' ITERATIONS',/  
     *' P,Q,DP,DQ,ACCUR =  ',1P4D17.7,D12.3/)
                     RETURN                                              
      END                                                               
C----------------------------------------------------------------------  
      SUBROUTINE  JWKB   (X,ETA,XL, FJWKB,GJWKB, IEXP)            
      DOUBLE PRECISION    X,ETA,XL, FJWKB,GJWKB, DZERO                      
C----------------------------------------------------------------------
C-----COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
C-----AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C-----CALCULATED IN SINGLE, RETURNED IN DOUBLE PRECISION VARIABLES
C-----CALLS DMAX1, SQRT, ALOG, EXP, ATAN2, FLOAT, INT     
C     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
C----------------------------------------------------------------------
      REAL    ZERO,HALF,ONE,SIX,TEN,RL35,ALOGE
      REAL    GH2,XLL1,HLL,HL,SL,RL2,GH,PHI,PHI10
      INTEGER IEXP, MAXEXP
      PARAMETER  ( MAXEXP = 300 )
      DATA  ZERO,HALF,ONE,SIX,TEN  /0.0E0, 0.5E0, 1.0E0, 6.0E0, 1.0E1/
      DATA DZERO,RL35,ALOGE /0.0D0, 35.0E0, 0.43429 45 E0 /  
C----------------------------------------------------------------------
CHOOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
      GH2   =  X * (ETA + ETA - X)                                         
      XLL1  = DMAX1( XL * XL + XL, DZERO )                                   
      IF( GH2 + XLL1 .LE. ZERO )                 RETURN
      HLL  = XLL1 + SIX / RL35                                           
      HL   = SQRT(HLL)                                                 
      SL   = ETA / HL + HL / X                                             
      RL2  = ONE + ETA * ETA / HLL                                         
      GH   = SQRT(GH2 + HLL) / X                                         
      PHI  = X*GH - HALF*( HL*ALOG((GH + SL)**2 / RL2) - ALOG(GH) )      
      IF ( ETA.NE.ZERO ) PHI = PHI - ETA * ATAN2(X*GH,X - ETA)         
      PHI10 = -PHI * ALOGE                                                
      IEXP  =  INT(PHI10)                                               
      IF ( IEXP.GT.MAXEXP ) THEN
           GJWKB = TEN**(PHI10 - FLOAT(IEXP))               
      ELSE
           GJWKB = EXP(-PHI)                                
           IEXP  = 0                                        
      ENDIF
      FJWKB = HALF / (GH * GJWKB)                                           
C---------------------------------------------------------------------
C     END OF CONTINUED-FRACTION COULOMB & BESSEL PROGRAM  COUL90
C---------------------------------------------------------------------
      RETURN                                                            
      END                                                               

*   ------------------------------------------------------------------
*               H N O R M
*   ------------------------------------------------------------------
*
*       Returns the value of the normalization constant for an (nl)
*   hydrogenic function with nuclear charge ZZ.
*
      DOUBLE PRECISION FUNCTION HNORM(N,L,ZZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      M = L + L + 1
      A = N + L
      B = M
      T = A
      D = B
      M = M - 1
      IF (M .EQ. 0) GO TO 2
      DO 1 I = 1,M
      A = A - 1.d0
      B = B - 1.d0
      T = T*A
1     D = D*B
2     HNORM = DSQRT(ZZ*T)/( N*D)
      RETURN
      END

*   ------------------------------------------------------------------
*               H W F
*   ------------------------------------------------------------------
*
*       Returns the value of an unnormalized (nl) hydrogenic function
*   with nuclear charge ZZ and radius r.
*
      DOUBLE PRECISION FUNCTION HWF(N,L,ZZ,R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      K = N-L-1
      P = 1.d0
      A = 1.d0
      B = K
      C = N+ L
      X = -2.d0*ZZ*R/N
*
*  *****  TEST IF UNDERFLOW MAY OCCUR, IF SO SET HWF = 0
*
      IF ( X .LT. -150.D0 ) GO TO 5
      IF (K) 1,2,3
3     DO 4 I = 1,K
      P = 1.d0 + A/B*P/C*X
      A = A + 1.d0
      B = B - 1.d0
4     C = C - 1.d0
2     HWF = P*DEXP(X/2.d0)*(-X)**(L+1)
      RETURN
1     WRITE(*,7) N,L,ZZ,R
7     FORMAT(' FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM',/,
     1    'N = ',I4,'   L = ',I4,'   Z = ',F6.1,'   R = ',F8.4)
       STOP
5     HWF = 0.d0
      RETURN
      END

C**************************************************************
      subroutine pulse
      implicit double precision (a-h,o-z)
      parameter(nxmax=381000,ntmax=100000,ncmax=61,nprmax=20000,
     1 nfmax=10)
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1      ngob1,nerg
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/in1/w1,e1,s1,ton1,toff1,del1,w2,e2,s2,ton2,toff2,del2,agbr 
      common/arry1/x(0:nxmax),g(0:nxmax),f(0:ntmax),fd(0:nxmax,1:nfmax),
     2             cek3(0:ncmax,0:ncmax,0:2*ncmax)      
      read(92,*)
10    read(92,*,end=20) ntfin,rdum1,rdum2,f(ntfin)
      go to 10
20    ntfin = ntfin+1
      f(ntfin) = 0.0d0
      print *,'read the field; ntfin = ',ntfin
*debug
      do 30 i=0,ntfin
       write(294,*) i,f(i)
30    continue
*debug
      return
      end

C****************************************************************
      subroutine radial
      implicit double precision (a-h,o-z)
      complex*16 czero,chalf,cone,ctwo,cthree,cfour,ci 
      complex*16 ch,cst,cst1,cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      complex*16 q,v,ff,bb
      parameter(nxmax=381000,ntmax=100000,ncmax=61,nprmax=20000,
     1  nfmax=10)
      common/cnst1/iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1      ngob1,nerg
      common/arry1/x(0:nxmax),g(0:nxmax),f(0:ntmax),fd(0:nxmax,1:nfmax),
     2 cek3(0:ncmax,0:ncmax,0:2*ncmax)      
      common/cnst2/pi,h,dt,tau,zet,gbr,xtime
      common/discr/nf,nn(nfmax),ll(nfmax)
      common/carry3/q(0:nxmax,0:ncmax),v(0:nxmax,0:ncmax),
     1               ff(0:nxmax,0:ncmax),bb(0:nxmax,0:ncmax)
      common/keys/ key1,key2,key3,key4,key5
      common/in1/w1,e1,s1,ton1,toff1,del1,w2,e2,s2,ton2,toff2,del2,agbr 
      common/ccnst1/czero,chalf,cone,ctwo,cthree,cfour,ci,ch,cst,cst1,
     1    cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27

C*** INITIAL RADIAL FUNCTION AND CHANNEL POTENTIALS
      do 11 i=0,nx
         g(i) = 2.d0*x(i)*exp(-x(i))      ! initial radial function
         xx=x(i)
         if(i.eq.0) xx=1.d-40            
         do 13 j=0,nc                     ! channel potentials
             cm = real(j,8)
             v(i,j) = dcmplx( cm*(cm+1.d0)/(2.d0*xx*xx) - zet/xx, 0.d0) 
C--- Gobbler:
             if(x(i).gt.gbr)
     1       v(i,j) = v(i,j) - ci*dcmplx(agbr*(x(i)/gbr - 1.d0)**3,0.d0)            
13       continue
11    continue

C**** RADIAL FUNCTIONS FOR DISCRETE STATES IF KEY1 = 1
      if(key1.eq.0) goto 9
      do 3 kd = 1, nf
          un = hnorm(nn(kd),ll(kd),1.d0)
          do 4 i = 0,nx
             fd(i,kd)   = un * hwf(nn(kd),ll(kd),1.d0,x(i))
4         continue
3     continue
9     return
      end
******
