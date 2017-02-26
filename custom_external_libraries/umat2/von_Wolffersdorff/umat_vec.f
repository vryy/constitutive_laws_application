c------------------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c------------------------------------------------------------------------------
c user subroutine for Abaqus 5.8-10
c------------------------------------------------------------------------------
c
c Implemented constitutive law:
c -----------------------------
c Hypoplasticity in the version of P.-A. von Wolffersdorff (1996):
c A hypoplastic relation for granular materials with a predefined limit
c state surface.
c Mechanics of Cohesive-Frictional Materials 1:251-271
c see also http://geotechnik.uibk.ac.at/res/hypopl.html
c
c The implementation is based on our article
c Fellin, W. and Ostermann, A. (2002):
c Consistent tangent operators for constitutive rate equations.
c International Journal for Numerical and Analytical Methods in Geomechanics
c
c ----------------------------------------------------------------------------
c Independent of the constitutive equation a cohesion can be prescribed.
c The string for the material name may contain 9 characters.
c ----------------------------------------------------------------------------
c Material constants: n = index of material constants in abaqus input file
c   	
c        ---------------------------------------------------------------------
c        n      
c        ---------------------------------------------------------------------
c        1      phi_c
c        2      h_s
c        3      n
c        4      e_d0
c        5      e_c0
c        6      e_i0
c        7      alpha
c        8      beta 
c        9      Tc = c/tan phi_c 
c               (isotropic stress subtracted to simulate cohesion)
c        ----------------------------------------------------------------------
c
c Solution dependent state variables (statev):
c definition via sdvini
c
c        1 ... time integration method
c              = 1 .. forward Euler, constant substepping, no error estimation
c              = 2 .. forward Euler, variable substepping, error control
c        2 ... used for estimating a practical step size [Huegel 1995]
c              used only when statev(1)=1
c              a smaller value leads to a large number of steps
c        3 ... desired accuracy of calculated stresses 
c              used only when statev(1)=2
c        4 ... maximum number of time substeps, if the limit is exceeded 
c              abaqus is forced to reduce the overall time step size 
c              (cut-back)
c        5 ... suggested size of first time substep (used only if statev(1)=2)
c        6 ... mobilized friction angle (output)
c        7 ... actual void ratio (essential input and output)
c
c Authors: 
c     W. Fellin, wolfgang.fellin@uibk.ac.at
c     Institute of Geotechnical and Tunnel Engineering
c
c     A. Ostermann, alexander.ostermann@uibk.ac.at
c     Department of Engineering Mathematics, Geometry and Computer Science
c
c     University of Innsbruck 
c
c This implementation uses ideas from:
c     H. Huegel (1995): estimation of a practical step size for 
c                       forward Euler with constant step size 
c                       the subroutine for evaluation of constitutive 
c                       law based on his umat
c     D. Roddeman (1997): idea of reducing time substep when constitutive 
C                         law is not defined
c
c Last change: 8/2002  
c 2/2006: Bug fix in subroutine evolut() 
c         pertubation in y --> symmetric pertubation in D
c 5/2006: integration method 3: Euler forward with Richardson extraploation,
c         ddsdde with central differences
c         recommended only for testing purposes
c----------------------------------------------------------------------------
      implicit none
      character*80 cmname
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,
     & layer, kspt, kstep, kinc
      double precision stress(ntens), statev(nstatv),
     &  ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),
     &  stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),
     &  props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      double precision sse, spd, scd, rpl, drpldt, dtime, temp, 
     &  dtemp, pnewdt, celent

c declaration of local variables
      integer i, error, maxnint, nintmeth, nasv, nasvdim, 
     &  nfasv, nydim, nyact
      double precision D(3,3), valD, phim, theta, tolsubt, 
     &  tolintT, tolabsT, dtsub, hmin 
      double precision valu33
c nasvdim denotes the maximum number of additional state variables
c if more than 18 are used change nasvdim and dimension of Q,
c asvr, asvrh in subroutine evolut
      parameter (nasvdim = 18)
      parameter (nydim = 42+7*nasvdim)
c additional state variables
      double precision  asv(nasvdim)
c number of first additional state variable in statev field 
      parameter (nfasv = 7)
c solution vector (stresses, Jacobian information, additional state variables)
      double precision  y(nydim)
c dummy variables for integrator
      double precision yp(nydim), v(nydim), y2(nydim), yh(nydim)
      external evolut
c switch for printing information
      logical prsw, elprsw
      parameter (prsw=.false.)
c print informations about time integration, useful when problems occur
      elprsw = .false.
      if (prsw) then
c print only in some defined elements
        if ((noel.eq.101).and.(npt.eq.1)) elprsw = .true.
c        if (noel.eq.101) elprsw = .true.
      endif

c define number of additional state variables
      call define(nasv)
      nyact = 42 + 7*nasv
      if (nyact.gt.nydim) then
         write(6,*) 'UMAT: nasvdim too small, program terminated'
         call XIT
      endif
c
c setting integer flags (sdvini can only set real statev ABAQUS ;-( )
c
c integration method  
      nintmeth = 2
      if (statev(1).gt.0.5d0.and.statev(1).lt.1.5d0) nintmeth = 1
      if (statev(1).gt.2.5d0.and.statev(1).lt.3.5d0) nintmeth = 3
c tolerance for substepping (forward Euler with constant substepping)
      tolsubt =  statev(2)
c tolerance for stress error (forward Euler with error control)
      tolintT =  statev(3)
c absolute tolerance for stress error (forward Euler with error control)
      tolabsT =  1.d-3
c suggested time substep size
      dtsub = statev(5)
c maximum number of time substeps
      maxnint =  int(statev(4))
c minimal time substep size
      hmin = 1.d-10

c vector of additional state variables
      do i=1,nasv
        asv(i) = statev(i-1+nfasv)
      enddo

c Error-Management:
c ----------------
c error =  0 ... no problem in time integration
c error =  1 ... problems in evaluation of the time rate, (e.g. undefined 
c                stress state), reduce time integration substeps
c error =  3 ... problems in time integration, reduce abaqus load increment 
c                (cut-back)
c error = 10 ... severe error in evaluating of time rate (e.g. wrong material 
c                constants), terminate calculation
      error = 0
c
c ----------------
c Time integration
c ----------------
c
      call getD(D,ndi,nshr,ntens,dstran,dtime)
      valD = valu33(D)         
      call iniy(y,nydim,asv,nasv,ndi,nshr,ntens,stress)

      if (elprsw) then
        write(6,*) '==================================================='
        write(6,*) 'Call of umat:'
        write(6,*) '==================================================='
        call wrista(3,y,nydim,D,dtime,coords,statev,nstatv,
     &              props,nprops,noel,npt,ndi,nshr,kstep,kinc)
      endif

c
c parameter of the numerical differentiation: sqrt(macheps)*||D||
c
c double precision
      theta = 1.0d-7 * max(valD,1.d0)
c quadruple precision
c      theta = 1.0d-14 * max(valD,1.d0)

      if (valD.le.0.d0) then 
c for ||D|| = 0 return with crude approximation of the Jacobian, 
c no stress and statev update      
        call getDZjac(ddsdde,ntens,ndi,nshr,
     &                stress,yp,nydim,asv,nasv,props,nprops,
     &                theta,dtime,dtime,tolabsT)
        return
      endif

      if ((nintmeth.eq.2).or.(nintmeth.eq.3)) then
c local extrapolation based on forward Euler, variable substeps, 
c consistent Jacobian and error estimation
        if ((dtsub.le.0.d0).or.(dtsub.gt.dtime)) then
          dtsub = dtime
        endif
        call eulexp(y,nyact,evolut,asv,nasv,dtsub,dtime,tolintT,tolabsT,
     &              maxnint,hmin,theta,D,valD,props,nprops,error,
     &              yp,v,y2,yh,elprsw)
      elseif (nintmeth.eq.1) then
c forward Euler with constant substeps, consistent Jacobian
c and no error estimation
        call euler(y,nyact,evolut,dtime,maxnint,tolsubt,tolabsT,theta,
     &             D,valD,asv,nasv,props,nprops,error,yp,elprsw)
      endif

      if (error.eq.3) then
c reduce abaqus load increment
        pnewdt = 0.25d0
        write(6,*) 'umat: ask ABAQUS for reduced step size (error=3)'
        call wrista(1,y,nydim,D,dtime,coords,statev,nstatv,
     &              props,nprops,noel,npt,ndi,nshr,kstep,kinc)
        return
      elseif (error.eq.10) then
         call wrista(2,y,nydim,D,dtime,coords,statev,nstatv,
     &               props,nprops,noel,npt,ndi,nshr,kstep,kinc)
c Don´t use STOP, use abaqus-subroutine XIT to terminate the program correctly
         call XIT
      endif
      statev(5) = dtsub
c updated solution for abaqus 
      call solout(stress,ntens,ndi,nshr,
     &            asv,nasv,ddsdde,dtime,y,nydim)
c updated vector of additional state variables to abaqus statev vector
      do i=1,nasv
        statev(i-1+nfasv) = asv(i) 
      enddo
c -----------------------
c End of time integration
c -----------------------

      if (nintmeth.eq.3) then
c ddssddee with central differences
        call getDZjac(ddsdde,ntens,ndi,nshr,
     &                stress,yp,nydim,asv,nasv,props,nprops,
     &                theta,dtime,dtime,tolabsT)
      endif

      call phimob(ndi,nshr,ntens,stress,phim,props(9))
      statev(6) = phim
      return
      end

c-----------------------------------------------------------------------------
      subroutine getDZjac(ddsdde,ntens,ndi,nshr,
     &                    stress,yp,nydim,asv,nasv,props,nprops,
     &                    theta,dtime,h,tolabs)
c-----------------------------------------------------------------------------
c generation of a crude approximation of the Jacobian matrix for D=0 
c
      implicit none
      integer ntens, ndi, nshr, nasv, nprops, nydim
      double precision ddsdde(ntens,ntens), stress(ntens),
     &                 asv(nasv), props(nprops), theta,
     &                 dtime, yp(nydim)

      integer i, j, k, error
      double precision T(3,3), Dh1(3,3), Dh2(3,3), TR1(3,3),
     &                 TR2(3,3), B(6,6), valDh1, valDh2, h, tolabs
      
      call getT(T,ndi,nshr,ntens,stress)
      do j=1,6
        do i=1,3
           do k=1,3
             Dh1(i,k) = 0
             Dh2(i,k) = 0
           enddo
        enddo
c perturbation of D --> Dh
        if (j.le.3) then
          Dh1(j,j) = Dh1(j,j) - theta
          Dh2(j,j) = Dh2(j,j) + theta
        else
          i = j/3
          k = j-i-1
          Dh1(i,k) = Dh1(i,k) - theta
          Dh2(i,k) = Dh2(i,k) + theta
        endif
        valDh1 = theta
        valDh2 = theta
c calculate perturbed objective time derivative of stresses           
        call getobjtr(TR1,yp(43),error,T,Dh1,valDh1,
     &                asv,nasv,props,nprops,h,tolabs)
        call getobjtr(TR2,yp(43),error,T,Dh2,valDh2,
     &                asv,nasv,props,nprops,h,tolabs)
        do i=1,3
          do k=1,3
            TR2(i,k) = ( TR2(i,k) - TR1(i,k) ) / theta / 2
          enddo
        enddo
        do i=1,3
          B(i,j) = TR2(i,i)
        enddo
        B(4,j) = TR2(1,2)
        B(5,j) = TR2(1,3)
        B(6,j) = TR2(2,3)
      enddo
      do i=1,ndi
        do j=1,ndi
          ddsdde(i,j) = B(i,j)
        enddo
      enddo
  
      do i=ndi+1,ndi+nshr
        do j=1,ndi
          ddsdde(i,j) = B(i+3-ndi,j)
        enddo
      enddo
   
      do i=1,ndi
        do j=ndi+1,ndi+nshr
          ddsdde(i,j) = B(i,j+3-ndi)
        enddo
      enddo
   
      do i=ndi+1,ndi+nshr
        do j=ndi+1,ndi+nshr
          ddsdde(i,j) = B(i+3-ndi,j+3-ndi)
        enddo
      enddo

      return
      end

c------------------------------------------------------------------------------
      double precision function valu33(x)
c------------------------------------------------------------------------------
c value of a matrix
c------------------------------------------------------------------------------
      implicit none
      double precision x(3,3)

      valu33 = sqrt(x(1,1)**2+x(2,2)**2+x(3,3)**2+
     &  2.0d0*(x(1,2)**2+x(1,3)**2+x(2,3)**2))

      return
      end

c-----------------------------------------------------------------------------
      subroutine iniy(y,nydim,asv,nasv,ndi,nshr,ntens,stress)
c-----------------------------------------------------------------------------
c initializes the vector of state variables
c-----------------------------------------------------------------------------
      implicit none
      integer nydim, nasv, ndi, nshr, ntens
      double precision y(nydim), asv(nasv), stress(ntens)

      integer i

      do i=1,nydim
        y(i) = 0
      enddo

      do i=1,ndi
        y(i) = stress(i)
      enddo
      if (nshr.ge.1) then
        y(4) = stress(ndi+1)
      endif
      if (nshr.ge.2) then
        y(5) = stress(ndi+2)
      endif
      if (nshr.ge.3) then
        y(6) = stress(ndi+3)
      endif
c additional state variables
      do i=1,nasv
        y(42+i) = asv(i)
      enddo

      return
      end

c-----------------------------------------------------------------------------
      subroutine solout(stress,ntens,ndi,nshr,
     &                  asv,nasv,ddsdde,dtime,y,nydim)
c-----------------------------------------------------------------------------
c copy the vector of state variables to umat output
c-----------------------------------------------------------------------------
      implicit none
      integer nydim, nasv, ndi, nshr, ntens
      double precision y(nydim), asv(nasv), stress(ntens),
     &                 ddsdde(ntens,ntens), dtime 

      integer i,j

c updated stresses
      do i=1,ndi
        stress(i) = y(i)
      enddo
      if (nshr.ge.1) stress(ndi+1) = y(4) 
      if (nshr.ge.2) stress(ndi+2) = y(5) 
      if (nshr.ge.3) stress(ndi+3) = y(6) 

c additional state variables
      do i=1,nasv
        asv(i) = y(42+i) 
      enddo

c Jacobian
      do i=1,ndi
        do j=1,ndi
          ddsdde(i,j) = y(i+6*j)/dtime        
        enddo
      enddo

      do i=ndi+1,ndi+nshr
        do j=1,ndi
          ddsdde(i,j) = y((i+3-ndi)+6*j)/dtime     
        enddo
      enddo

      do i=1,ndi
        do j=ndi+1,ndi+nshr
          ddsdde(i,j) = y(i+6*(j+3-ndi))/dtime    
        enddo
      enddo

      do i=ndi+1,ndi+nshr
        do j=ndi+1,ndi+nshr
          ddsdde(i,j) = y((i+3-ndi)+6*(j+3-ndi))/dtime 
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------------
      subroutine getT(T,ndi,nshr,ntens,stress)
c-----------------------------------------------------------------------------
c transform stresses from vector into matrix form 
c-----------------------------------------------------------------------------
      implicit none
      integer ndi, nshr, ntens
      double precision stress(ntens), T(3,3)

      integer i,j

      do i=1,3
        do j=1,3
        T(i,j) = 0
        enddo
      enddo

      do i=1,ndi
        T(i,i) = stress(i)
      enddo

      if (nshr.ge.1) then
        T(1,2) = stress(ndi+1)
        T(2,1) = T(1,2)
      endif
      if (nshr.ge.2) then
        T(1,3) = stress(ndi+2)
        T(3,1) = T(1,3)
      endif
      if (nshr.ge.3) then
        T(2,3) = stress(ndi+3)
        T(3,2) = T(2,3)
      endif

      return
      end

c-----------------------------------------------------------------------------
      subroutine getD(D,ndi,nshr,ntens,dstran,dtime)
c-----------------------------------------------------------------------------
c strain increment into strain rate (from vector into matrix form) 
c-----------------------------------------------------------------------------
      implicit none
      integer ndi, nshr, ntens
      double precision D(3,3), dstran(ntens), dtime 

      integer i,j

      do i=1,3
        do j=1,3
          D(i,j)=0
        enddo
      enddo

      do i=1,ndi
        D(i,i) = dstran(i)/dtime
      enddo
      if (nshr.ge.1) then
        D(1,2) = 0.5d0*dstran(ndi+1)/dtime
        D(2,1) = D(1,2)
      endif
      if (nshr.ge.2) then
        D(1,3) = 0.5d0*dstran(ndi+2)/dtime
        D(3,1) = D(1,3)
      endif
      if (nshr.ge.3) then
        D(2,3) = 0.5d0*dstran(ndi+3)/dtime
        D(3,2) = D(2,3)
      endif

      return
      end

c-----------------------------------------------------------------------------
      subroutine evolut(n,y,yp,asv,nasv,theta,D,valD,props,nprops,
     &                  h,tolabs,error)
c-----------------------------------------------------------------------------
c evolution equations of the state variables
c
      implicit none
      integer n, nprops, nasv
      double precision y(n), yp(n), D(3,3), valD,
     &                 props(nprops), theta, asv(nasv)
      integer error

      integer i,ii,j,k
      double precision T(3,3), TR(3,3), Tvec(6), TRh(3,3), Dh(3,3), 
     &                 Q(18), asvr(18), asvrh(18), valDh, valu33, h, 
     &                 tolabs

      error = 0

      do i=1,3
         T(i,i)=y(i)
      enddo
      T(1,2) = y(4)
      T(1,3) = y(5)
      T(2,3) = y(6)
      T(2,1) = y(4)
      T(3,1) = y(5)
      T(3,2) = y(6)
      do i=1,nasv
        asv(i) = y(i+42)
      enddo

c objective time derivative of stresses TR
c and additional state variables yp(43..43-1+nasv)
      call getobjtr(TR,asvr,error,T,D,valD,asv,nasv,props,nprops,
     &              h,tolabs)
      if (error.ge.1) then
        return
      endif
      do i=1,3
         yp(i)=TR(i,i)
      enddo
      yp(4) = TR(1,2)
      yp(5) = TR(1,3)
      yp(6) = TR(2,3)
      do i=1,nasv
         yp(42+i) = asvr(i) 
      enddo

c time derivative of consistent Jacobian
      do j=1,6
c differentiation of the i'th stress y(i) 
c with respect to the j´th stretching
c B_ij = y(i+6*j) = stress(i) / (dstran(j) / dtime) etc.
c B_11 = y(7) = partial T_11 / partial D_11 etc.
        do i=1,6         
          Tvec(i) = y(i)+theta*y(i+6*j) 
        enddo
        do i=1,3
           T(i,i) = Tvec(i)
        enddo
        T(1,2) = Tvec(4)
        T(1,3) = Tvec(5)
        T(2,3) = Tvec(6)
        T(2,1) = Tvec(4)
        T(3,1) = Tvec(5)
        T(3,2) = Tvec(6)
c  additional state variables
c  G_1 = partial e / partial D_11
       do i=1,nasv
         Q(i) = y(42+i)+theta*y(42+nasv+i+nasv*(j-1))                  
       enddo
        do i=1,3
          do k=1,3
            Dh(i,k) = D(i,k)
          enddo
        enddo
c perturbation of D --> Dh
        if (j.le.3) then
          Dh(j,j) = Dh(j,j) + theta
        else
          i = j/3
          k = j-i-1
c pertubation in y --> symmetric pertubation in D
          Dh(i,k) = Dh(i,k) + theta/2.0d0
        endif
        valDh = valu33(Dh)
c calculate perturbed objective time derivative of stresses           
        call getobjtr(TRh,asvrh,error,T,Dh,valDh,Q,nasv,props,nprops,
     &                h,tolabs)
        if (error.ge.1) then
          return
        endif
c dB = 1/theta ( TR(T+theta B, Dh, Q) - TR(T,D,e) )
        do i=1,3
          yp(i+6*j) = ( TRh(i,i) - TR(i,i) )/theta
        enddo
        yp(4+6*j) = ( TRh(1,2) - TR(1,2) )/theta
        yp(5+6*j) = ( TRh(1,3) - TR(1,3) )/theta
        yp(6+6*j) = ( TRh(2,3) - TR(2,3) )/theta
c dG
        do i=1,nasv
          ii = 42+nasv+i
          yp(ii+nasv*(j-1)) = ( asvrh(i)-asvr(i) ) / theta
        enddo

      enddo

      return 
      end

c------------------------------------------------------------------------------
      subroutine phimob(ndi,nshr,ntens,stress,phim,Tc)
c------------------------------------------------------------------------------
c calculate principal stresses with the abaqus internal routine sprinc and
c out from it the mobilized friction angle
c------------------------------------------------------------------------------
      implicit none
      integer ndi, nshr, ntens
      double precision stress(ntens), phim, Tc

      integer lstr
      double precision tmin, tmax, ps(3), PI, zero
      parameter(PI=3.141592653589793d0, zero=1d-10)

      lstr = 1
      call sprinc(stress,ps,lstr,ndi,nshr)
      tmax = -max(abs(ps(1)),abs(ps(2)),abs(ps(3)))
      tmin = -min(abs(ps(1)),abs(ps(2)),abs(ps(3)))
      if (abs(tmax+tmin-2.0d0*Tc).le.zero) then
        phim = 0.d0
      else
        phim = asin((tmax-tmin)/(tmax+tmin-2.0d0*Tc)) * 180.d0/PI
      endif

      return
      end

c-----------------------------------------------------------------------------
      subroutine euler(y,n,fcn,time,maxstep,tol,tolabs,theta,D,valD,
     &                 asv,nasv,props,nprops,error,yp,elprsw)
c-----------------------------------------------------------------------------
c
c  numerical solution of y'=f(y)
c  forward Euler method with constant step size 
c
      implicit none
      integer n, nprops, maxstep, nasv, error
      external fcn
      double precision y(n), D(3,3), valD, asv(nasv), yp(n),
     &                 props(nprops), theta, time, tol, tolabs
      logical elprsw

      integer i,j,nsub
      double precision h

      error = 0

      if (tol.le.0.d0) then
        nsub = 1
      else
        nsub = max(int(valD*time/tol),1)
      endif
      if (nsub.gt.maxstep) then
        if (elprsw) write(6,*) 'number of time substeps ', 
     &                         'too big, reject step'
        error = 3
        return      
      endif

      h = time/nsub
      do j=1,nsub
c
c  Euler step
c
        call fcn(n,y,yp,asv,nasv,theta,D,valD,props,nprops,
     &           h,tolabs,error)
        if (error.ge.3) then
          return
        elseif (error.eq.1) then
          if (elprsw) write(6,*) 'euler: stress rate ',
     &                           'undefined in sub step', j
          error = 3
          return 
        endif
        do i=1,n
          y(i) = y(i) + h*yp(i)
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------------
      subroutine eulexp(y,n,fcn,asv,nasv,h,time,tol,tolabs,maxnint,hmin,
     &                theta,D,valD,props,nprops,error,yp,v,y2,yh,elprsw)
c-----------------------------------------------------------------------------
c
c  numerical solution of y'=f(y)
c  forward Euler with local extrapolation
c
      implicit none
      integer n, nasv, nprops, maxnint, error
      external fcn
      double precision y(n), h, D(3,3), valD, asv(nasv),
     &                 yp(n), v(n), y2(n), yh(n),
     &                 props(nprops), theta, time, tol, tolabs
      logical elprsw 

      integer i, naccst, nrejst
      double precision sci, errt, hh, actt, fhnew, hsav, hmin
      logical final

      if (h.eq.time) then 
        final = .true.
      else
        final = .false.
      endif
      naccst = 0
      nrejst = 0
      error = 0
      actt = 0
      hsav = h
c
c  two Euler steps with size h/2
c
      if (elprsw) write(6,*) 'time integration starts with h = ',h
      if (elprsw) write(6,*) 'end time', time
c first call of time rate
      call fcn(n,y,yp,asv,nasv,theta,D,valD,props,nprops,h,tolabs,error)
        if (error.ge.3) then
          if (elprsw) then 
            write(6,*) 'eulexp: call fcn for h/2-step --> error=3'
          endif
        return
      elseif (error.eq.1) then
c undefined stress state in fcn
        write(6,*) 'UMAT - first call of constitutive law:'
        write(6,*) 'state variables (stress state) given by ABAQUS'
        write(6,*) 'can not be handled by constitutive law!'
        error = 10
        return
      endif     
 20   continue
      if ((h.lt.hmin).and.(.not. final)) then
          if (elprsw) write(6,*) 'time sub step too small, reject step'
          error = 3
          return
      endif
      hh = h/2
      do i=1,n
        y2(i) = y(i) + hh*yp(i)
      enddo
      call fcn(n,y2,v,asv,nasv,theta,D,valD,props,nprops,h,tolabs,error)
      if (error.eq.1) then
c undefined stress state in fcn (tr T <= 0), h = h/2
        h = h/2
        if (elprsw) then
          write(6,*) 'eulexp: stress rate undefined at half sub step'
          write(6,*) 'stress(i): ',(y2(i),i=1,6)
          write(6,*) 'reduce sub step size'
        endif
        goto 20
      endif     
      do i=1,n
        y2(i) = y2(i) + hh*v(i)
      enddo
c
c  difference of step with size h to y2
c
      do i=1,n
        v(i) = y2(i) - y(i) - h*yp(i)
      enddo
c
c  error estimate (of the first 6 components, stresses)
c
      errt = 0
      do i=1,6
        sci = max( abs(y2(i)), abs(y(i)) ) + tolabs
        errt = errt + ( v(i)/sci )**2
      enddo
      errt = sqrt(errt)
      errt = max(errt,1.d-16)
c
c calculate factor of step change (times a safety factor)
c
      fhnew = 0.90d0*sqrt(tol/errt)
      if (errt.gt.tol) then
c
c do not accept step size, reject step
c
        final = .false.
        nrejst = nrejst + 1
        if (elprsw) write(6,*) 'reject: h, actt', h,actt
        h = max(0.2d0, fhnew)*h
        goto 20
      else
c
c update and try to calculate time rate at the end of the substep
c second oder update
c
        do i=1,n
          yh(i) = y2(i) + v(i)
        enddo
        call fcn(n,yh,v,asv,nasv,theta,D,valD,props,nprops,
     &           h,tolabs,error)
        if (error.ge.3) then
          if (elprsw) then 
            write(6,*) 'eulexp: test call fcn end of step --> error=3'
          endif
          return
        elseif (error.eq.1) then
c updated state variables cannot be handled by constitutive law:
c reduce step size, reject step
          h = h/2
          if (elprsw) then
            write(6,*)'eulexp: stress rate undefined at end of sub step'
            write(6,*)'stress(i): ',(yh(i),i=1,6)
            write(6,*)'reduce sub step size'
          endif
          goto 20
        else
c updated state variables can be handled by constitutive law: accept
c step size and update state variables, reuse the above objective time
c rate of state variables for new step
          actt = actt + h
          naccst = naccst + 1
          if (naccst.gt.maxnint) then
             if (elprsw) write(6,*) 'number of time substeps ',naccst,
     &                              ' too big, reject step'
             error = 3
             return
          endif
          do i=1,n
            y(i) = yh(i)
            yp(i) = v(i)
          enddo  
          if (elprsw) write(6,*) 'accept: h, actt', h,actt
        endif
        if (final) then 
          h = hsav
          if (elprsw) then 
            write(6,*) 'end of time integration'
            write(6,*) 'number of accepted substeps: ', naccst
            write(6,*) 'number of rejected substeps: ', nrejst
            write(6,*) 'time integration proposes ',
     &                           'new start h = ',h
          endif
          return
        endif
c suggested new step size limited by a factor of 5
        h = min(5.d0,fhnew)*h
        if (actt+h.ge.time) then
          final = .true.
c substep size for next time step 
c last suggested substep, reduced by a safety factor
          hsav = 0.5d0*h
          h = time - actt
        endif
        goto 20
      endif

      return
      end

c-----------------------------------------------------------------------------
      subroutine wrista(mode,y,nydim,D,dtime,coords,statev,nstatv,props,
     &  nprops,noel,npt,ndi,nshr,kstep,kinc)
c-----------------------------------------------------------------------------
      implicit none
      integer mode, nydim, nstatv, nprops, noel, 
     &  npt, ndi, nshr, kstep, kinc    
      double precision y(nydim), D(3,3), coords(3), statev(nstatv),
     &  props(nprops), dtime
 
      integer i
      if (mode.eq.2) then
        write(6,*) '==================================================='
        write(6,*) 'ERROR: abaqus job failed during call of umat'
        write(6,*) '==================================================='
        write(6,*) 'state dumb:'
        write(6,*) 
      endif
      write(6,111) 'Step: ',kstep, 'increment: ',kinc,
     & 'element: ', noel, 'Integration point: ',npt
      write(6,*) 
      if (mode.eq.2) then
        write(6,*) 'Co-ordinates of material point:'
        write(6,104) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',
     &    coords(3)
        write(6,*) 
        write(6,*) 'Material parameters:'
        write(6,*) 
        do i=1,nprops
          write(6,105) 'prop(',i,') = ',props(i)
        enddo 
        write(6,*)
        write(6,102) 'No. of mean components:  ',ndi
        write(6,102) 'No. of shear components: ',nshr
        write(6,*)
      endif
      if ((mode.eq.2).or.(mode.eq.3)) then
        write(6,*) 'Stresses:'
        write(6,*) 
        write(6,106) 'T(1,1) = ',y(1),'T(1,2) = ',y(4),'T(1,3) = ',
     &    y(5)
        write(6,106) 'T(2,1) = ',y(4),'T(2,2) = ',y(2),'T(2,3) = ',
     &    y(6)
        write(6,106) 'y(3,1) = ',y(5),'T(3,2) = ',y(6),'T(3,3) = ',
     &    y(3)
        write(6,*) 
        write(6,*) 'Stretching rate:'
        write(6,*) 
        write(6,101) 'D(1,1) = ',D(1,1),'D(1,2) = ',D(1,2),'D(1,3) = ',
     &    D(1,3)
        write(6,101) 'D(2,1) = ',D(2,1),'D(2,2) = ',D(2,2),'D(2,3) = ',
     &    D(2,3)
        write(6,101) 'D(3,1) = ',D(3,1),'D(3,2) = ',D(3,2),'D(3,3) = ',
     &    D(3,3)
        write(6,*) 
        write(6,*) 'Time increment:'
        write(6,*) 
        write(6,108) 'dtime = ',dtime
        write(6,*) 
        write(6,*) 'Void ratio:'
        write(6,*) 
        write(6,109) 'e = ',statev(7)
        write(6,*) 
        write(6,*) '(Capillary) cohesion:'
        write(6,*) 
        write(6,110) 'Tc = ',props(9)
        write(6,*) 
        write(6,*) '==================================================='
      endif
  
101   format(1X,3(a9,e10.4,2X))
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a5,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,f12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,f12.4)
109   format(1X,a4,f10.4)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)
       
      return
      end

c ============================================================================
c Following subroutines have to be adapted when changing the constitutive law
c ============================================================================

c-----------------------------------------------------------------------------
      subroutine define(nasv)
c-----------------------------------------------------------------------------
      implicit none 
      integer nasv
c number of additional state variables 
c must be less than  18 (otherwise change nasvdim in umat_vec and 
c dimension of Q, asvr, asvrh in subroutine evolut)
      nasv = 1
      return
      end

c------------------------------------------------------------------------------
      subroutine getobjtr(TR,asvr,error,T,D,valD,asv,nasv,props,nprops,
     &                    dt,tolabs)
c------------------------------------------------------------------------------
c calculate objective time rate of stresses and additional state variables
c abaqus requires Green-McInnis-Naghdi stress rate
c------------------------------------------------------------------------------
      implicit none 
      integer nasv, nprops, error
      double precision TR(3,3), asvr(nasv), T(3,3), D(3,3), valD
      double precision asv(nasv), props(nprops), dt, tolabs

      integer i, j, m, n, delta(3,3)
      double precision Ts(3,3), Ts2(3,3), Tsv(3,3), Tsv2(3,3), 
     &  Tsv3(3,3), LL(3,3,3,3), NN(3,3)
      double precision void, phic, hs, en, ed0, ec0, ei0, alpha, beta, 
     &  Tc, trT, sinphi, sq2, sq3, sq6, c3t, valTsv, trTsv2, trTsv3,
     &  tpsi, tpsi2, Fm, Fm2, term1, term2, term3, term4, term5,
     &  a, a2, ed, ec, ei, fb, fe, fd, fss, trTs2, trD, trT2, voidnew
      double precision valu33
      data delta/1,0,0,0,1,0,0,0,1/ 
      double precision PI, zero
      parameter(PI=3.141592653589793d0, zero=1d-10)

      error = 0

c ---------------
c time derivative of additional state parameter void ratio     
c
        trD = D(1,1) + D(2,2) + D(3,3)
        asvr(1) = ( 1 + asv(1) )*trD

c ---------------
c time derivative of stresses: hypoplastic law
c calculate L and N as functions of T, e
c
      void = asv(1)
      if (nprops.ge.9) then
        phic = props(1)*PI/180 
        hs   = props(2)
        en   = props(3) 
        ed0  = props(4)
        ec0  = props(5)
        ei0  = props(6)
        alpha = props(7)
        beta  = props(8)
        Tc    = props(9)
      else
        write(6,*) 'UMAT: To less material constants.'
        write(6,*) 'nprops = ', nprops
        write(6,*) 'You have to define nprops=9 constants:'
        write(6,*) 
     &    'phi_c [°], h_s [F/A], n, e_d0, e_c0, e_i0, alpha, beta, Tc'
        write(6,*) 'Program terminated.'
        error = 10
        return
      endif
c check user input on severe errors
c which would lead to a crash in evaluation of constitutive law 
      if (phic.le.zero) then
        write(6,*) 'UMAT: phic = props(1) must be > 0'
        write(6,*) 'Program terminated.'
        error = 10
      endif
      if (hs.le.zero) then
        write(6,*) 'UMAT: hs = props(2) must be > 0'
        write(6,*) 'Program terminated.'
        error = 10
      endif
      if (ed0.le.zero) then
        write(6,*) 'UMAT: ed0 = props(4) must be > 0'
        write(6,*) 'Program terminated.'
        error = 10
      endif
      if (ec0.le.ed0) then
        write(6,*) 'UMAT: ec0 = props(5) must be > ed0 = props(4)'
        write(6,*) 'Program terminated.'
        error = 10
      endif
      if (ei0.le.ec0) then
        write(6,*) 'UMAT: ei0 = props(6) must be > ec0 = props(5)'
        write(6,*) 'Program terminated.'
        error = 10
      endif
c check actual void ratio
      if (void.le.0) then
        write(6,*) 'UMAT: illegal (negative) actual void ratio detected'
        write(6,*) 'Program terminated.'
        error = 10
      endif
c
c consider (capillary) cohesion as isotropic stress Tc = c / tan(phi_c)
c subtract this stress to get hypoplastic law for cohesion  
c
      do i=1,3
         T(i,i) = T(i,i) - Tc
      enddo

      trT = T(1,1) + T(2,2) + T(3,3)

c For zero stress --> hypoplastic compression law
      if (abs(trT).lt.tolabs) then
        do i=1,3
          do j=1,3
            TR(i,j) = 0
          enddo
        enddo
        if ((trD.lt.0).and.(void.le.ei0)) then
c derivative of compression law
          TR(1,1) = asvr(1)/void/en/3 * tolabs / ( tolabs/hs )**en 
          TR(2,2) = TR(1,1)
          TR(3,3) = TR(1,1)
        endif
        return
      endif

      if (trT.ge.0) then
c undefined stress state
        error = 1
        return
      endif

      sinphi = sin(phic)
      sq2    = sqrt(2.0d0)
      sq3    = sqrt(3.0d0)
      sq6    = sqrt(6.0d0)

c calculate relative stress: Ts = T / trac(T)
      term1 = 1.d0/trT
      do i=1,3
        do j=1,3
          Ts(i,j) = T(i,j)*term1
        enddo
      enddo
c calculate square of relative stress: Ts2 = Ts * Ts
      call aikbkj(Ts,Ts,Ts2)
      trTs2  = Ts2(1,1) + Ts2(2,2) + Ts2(3,3)
c calculate deviator: Tsv = Ts - 1/3 I
      do i=1,3
        do j=1,3
          Tsv(i,j) = Ts(i,j)
        enddo
        Tsv(i,i) = Tsv(i,i) - 1.0d0/3.0d0
      enddo

      call aikbkj(Tsv,Tsv,Tsv2)
      call aikbkj(Tsv2,Tsv,Tsv3)
c calculate norm of relative deviatoric stress, ||Tsv|| 
      valTsv = valu33(Tsv)
      trTsv2 = Tsv2(1,1) + Tsv2(2,2) + Tsv2(3,3)
      trTsv3 = Tsv3(1,1) + Tsv3(2,2) + Tsv3(3,3)

      if (trTsv2.le.1.d-10) then
        c3t = 1.0d0
      else
        c3t = -sq6*trTsv3/trTsv2**1.5d0
        if (c3t.gt.1.0d0)  c3t =  1.0d0
        if (c3t.lt.-1.0d0) c3t = -1.0d0 
      endif

      tpsi  = sq3*valTsv 
      tpsi2 = tpsi**2

      if ( tpsi2/8.d0+(2.d0-tpsi2)/(2.d0+sq2*tpsi*c3t ).lt.0.d0) then
        Fm = 1.d-10
      else 
        term1 = tpsi2/8.d0 + ( 2.d0-tpsi2 )/( 2.d0+sq2*tpsi*c3t )
        if (term1.lt.0) then
          write(6,*) 'UMAT - WARNING in subroutine getobjtr'
          write(6,*) 'cannot evaluate Fm'
          error = 1
          return
        else
          term1 = sqrt(term1)
        endif
        term2 = tpsi/2.d0/sq2
        Fm  = term1-term2
      endif

      Fm2 = Fm**2

      a   = sq3*( 3.0d0-sinphi )/( 2.0d0*sq2*sinphi )
      a2  = a**2

      ed = ed0*exp(-(-trT/hs)**en)
      ec = ec0*exp(-(-trT/hs)**en)
      ei = ei0*exp(-(-trT/hs)**en)

c check actual void ratio
      if (void.lt.ed) then
        write(6,*) 'UMAT: actual void ratio e:', void
        write(6,*) 'less than ed', ed
        write(6,*) 'Maybe reducing ABAQUS step size will help'
        error = 3
      endif
      if (void.gt.ei) then
        write(6,*) 'UMAT: actual void ratio e:', void
        write(6,*) 'greater than ei', ei
      endif

      term3 = ( (ei0-ed0)/(ec0-ed0) )**alpha
      term4 = ( 3 + a2 - a*sq3*term3 )**(-1.d0)
      term5 = hs/en*( ei0/ec0 )**beta*( 1+ei )/ei*( -trT/hs )**(1.d0-en)
       
      fb  = term4*term5 
      fe  = ( ec/void )**beta

      fd = ( void-ed )/( ec-ed )
      if (fd.gt.0.0d0) then
        fd = fd**alpha
      else
        fd = 0.0d0
      endif
   
      fss = fb*fe/trTs2

      do i=1,3
        do j=1,3
          do m=1,3
            do n=1,3
              LL(i,j,m,n) = fss*( Fm2*delta(i,m)*delta(j,n)+
     &                            a2*Ts(m,n)*Ts(i,j) )
            enddo
          enddo
          NN(i,j) = fss*(  fd*a*Fm*( Ts(i,j) + Tsv(i,j) )  )
        enddo
      enddo

c calculate  TR_ij = L_ijkl*D_kl+N_ij*||D||
      do i=1,3
        do j=1,3
          TR(i,j) = 0.0d0
          do m=1,3
            do n=1,3
              TR(i,j) = TR(i,j) + LL(i,j,m,n)*D(m,n)
            enddo
          enddo
          TR(i,j) = TR(i,j) + NN(i,j)*valD
        enddo
      enddo
c ---------------
c
c add capillary cohesion stress so that T stays unchanged outside 
c
      do i=1,3
         T(i,i) = T(i,i) + Tc
      enddo

      return
      end

c------------------------------------------------------------------------------
      subroutine aikbkj(a,b,c)
c------------------------------------------------------------------------------
c matrix multiplication
c------------------------------------------------------------------------------
      implicit none
      integer i, j, k
      double precision a(3,3), b(3,3), c(3,3)

      do i=1,3
        do j=1,3
          c(i,j) = 0.0d0
          do k=1,3
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          enddo
        enddo
      enddo

      return
      end



