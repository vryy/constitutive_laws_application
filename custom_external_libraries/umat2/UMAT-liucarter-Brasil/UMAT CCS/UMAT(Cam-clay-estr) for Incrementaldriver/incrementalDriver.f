! Copyright (C)  2007  Andrzej Niemunis
!
! incrementalDriver is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! incrementalDriver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
!  USA.

! August 2007: utility routines added
! Februar 2008 :   *Repetition repaired again
! March 2008 : *ObeyRestrictions added 

       PROGRAM that_calls_umat   ! written by  A.Niemunis  2007
       use unsymmetric_module
       implicit none
       character*80  cmname,rebarn
       integer ndi,nshr,ntens,nstatv,nprops,ncrds
       integer noel,npt,layer,kspt,lrebar,kstep,kinc,i
       real(8),parameter,dimension(3,3):: delta =
     &                          reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))

      parameter(ntens=6,ndi=3,nshr=3,ncrds=3) ! same ntens as in SOLVER
      parameter(noel=1,npt=1,layer=1,kspt=1,lrebar=1)
      parameter( rebarn ='xxx')
      REAL time_begin, time_end
      real*8 dtime,temp,dtemp,sse,spd,scd,rpl,drpldt,pnewdt,celent
      real*8 stress(ntens),
     &  ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &  stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     &  coords(ncrds),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character(20) :: keywords(10), outputfilename
      character(len=260) ::  inputline(6), aLine,heading
      real(8), dimension(6,6)  :: cMt , cMe
      real(8), dimension(6)  :: mb, mbinc

      real(8), allocatable :: props(:), statev(:), r_statev(:)

      real(8),dimension(3,3):: Qb33,eps33,T33

      integer:: ifstress(ntens), maxiter, ninc,kiter,last_step,ikeyword,
     &          iRepetition, nRepetitions, iStep,nSteps,ntens_in

      real(8):: r_stress(ntens),a_dstress(ntens),u_dstress(ntens),
     &         stress_Rosc(ntens),r_stress_Rosc(ntens),
     &      ddstress(ntens), c_dstran(ntens) ,
     &      trT, ed,xmax,xmin,deltaLoadCirc(6),phase0(6),deltaLoad(9),
     &      dstran_Cart(6), ddsdde_bar(6,6), deltaTime
      real(8),parameter :: sq3=1.7320508075688772935d0,
     &                     sq6=2.4494897427831780982d0,
     &                     sq2=1.4142135623730950488d0,
     &                     Pi =3.1415926535897932385d0
      real(8),parameter ::
     &                     i3=0.3333333333333333333d0,
     &                     i2=0.5d0,
     &                     isq2=1/sq2,
     &                     isq3=1.0d0/sq3,
     &                     isq6=1.0d0/sq6

      real(8), parameter,dimension(1:6,1:6)::MRoscI=reshape              !  M for isomorphic Roscoe variables P,Q,Z,....
     &  ((/-isq3,-2.0d0*isq6,0.0d0,  0.0d0, 0.0d0, 0.0d0,
     &     -isq3, isq6,      -isq2,  0.0d0, 0.0d0, 0.0d0,
     &     -isq3, isq6,       isq2,  0.0d0, 0.0d0, 0.0d0,
     &      0.0d0, 0.0d0,     0.0d0, 1.0d0, 0.0d0, 0.0d0,
     &      0.0d0, 0.0d0,     0.0d0, 0.0d0, 1.0d0, 0.0d0,
     &      0.0d0, 0.0d0,     0.0d0, 0.0d0, 0.0d0, 1.0d0
     &  /),(/6,6/))

      real(8), parameter,dimension(1:6,1:6)::MRoscImT=MRoscI            !  latest M^-T (is orthogonal)

       real(8), parameter,dimension(1:6,1:6)::MRendul=reshape              !  M for isomorphic Rendulic sigma_11 = -T_11,  sigma_22 = -(T_22 + T_33) / sqrt(2)  Z, ....
     &  ((/ -1.0d0, 0.0d0,  0.0d0,  0.0d0,0.0d0,0.0d0,
     &      0.0d0, -isq2,  -isq2,  0.0d0,0.0d0,0.0d0,
     &      0.0d0, -isq2,   isq2,  0.0d0,0.0d0,0.0d0,
     &      0.0d0,  0.0d0, 0.0d0,  1.0d0,0.0d0,0.0d0,
     &      0.0d0,  0.0d0, 0.0d0,  0.0d0,1.0d0,0.0d0,
     &      0.0d0,  0.0d0, 0.0d0,  0.0d0,0.0d0,1.0d0
     &  /),(/6,6/))

      real(8), parameter,dimension(1:6,1:6)::MRendulmT=MRendul          !  latest M^-T  (is orthogonal)

      real(8), parameter,dimension(1:6,1:6)::MRosc=reshape              !  M for Roscoe variables p,q,z,....
     &  ((/-i3,-1.0d0, 0.0d0,    0.0d0,0.0d0,0.0d0,
     &     -i3, i2, -1.0d0,      0.0d0,0.0d0,0.0d0,
     &      -i3,i2, 1.0d0,       0.0d0,0.0d0,0.0d0,
     &      0.0d0, 0.0d0,0.0d0,  1.0d0,0.0d0,0.0d0,
     &      0.0d0, 0.0d0, 0.0d0, 0.0d0,1.0d0,0.0d0,
     &      0.0d0, 0.0d0, 0.0d0, 0.0d0,0.0d0,1.0d0
     &  /),(/6,6/))
      real(8), parameter,dimension(1:6,1:6)::MRoscmT=reshape            !  latest M^-T  (is not orthogonal)
     & ((/-1.0d0, -2.0d0*i3, 0.0d0,  0.0d0, 0.0d0,0.0d0,
     &    -1.0d0,   i3,      -i2,    0.0d0,0.0d0,0.0d0,
     &    -1.0d0,   i3,       i2,    0.0d0,0.0d0,0.0d0,
     &    0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0,
     &    0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,
     &    0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0
     &  /),(/6,6/))

      real(8), parameter,dimension(1:6,1:6)::MCart=reshape              !  M for Cartesian coords T_11, T_22, T_33, T_12,.....
     & ((/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &    0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &    0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,
     &    0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0,
     &    0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,
     &    0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0
     &  /),(/6,6/))
      real(8), parameter,dimension(1:6,1:6)::MCartmT=MCart              !  latest M^-T (is orthogonal)


      real(8),dimension(1:6,1:6)::M,MmT                                 !  currrent M and M^-T for a given iStep


      type descriptionOfStep
        integer :: ninc, maxiter, ifstress(ntens)
        real(8) :: deltaLoadCirc(ntens),phase0(ntens),deltaLoad(9),
     &             dfgrd0(3,3), dfgrd1(3,3),deltaTime
        character(20) :: keyword2, keyword3
      end type  descriptionOfStep

      type (descriptionOfStep) :: ofStep(30)                            !  stores descriptions of up to 30 steps which are repeated
      CALL CPU_TIME ( time_begin )
!  read the material parameters
      open(1,err=901,file='parameters.inp',status='old')
      read(1,*) cmname
      cmname = trim(cmname)
      read(1,*) nprops
      allocate( props(nprops) )
      do i=1,nprops
         read(1,*) props(i)
      enddo
      close(1)
!  initialize everything
      open(1,err=902,file='initialconditions.inp')
      read(1,*) ntens_in
      stress(:) = 0.0d0
      time(:) = 0.0d0
      stran(:)=0.0d0
      dtime = 0.0d0
      do i=1,ntens_in
         read(1,*) stress(i)
      enddo
      read(1,*) nstatv
      if(nstatv < 1) goto 500
      allocate( statev(nstatv) , r_statev(nstatv))
      statev = 0.0d0
      do i=1,nstatv
         read(1,*,end=500) statev(i)
      enddo
  500 continue
      close(1)

 ! read the loading path
      open(1,err=903,file='test.inp')
      read(1,'(a)') aLine 
        i = index(aLine,'#')
        if(i==0) then 
          outputfilename=trim(aLine)
          heading = '#'
        else 
          outputfilename=trim(aLine(:i-1))
          heading = trim( aLine(i+1:) )
        endif 
        
      open(2,err=904,file=outputfilename)
      
      write(2,'(a10,50a17)') 'time(1)', 'time(2)',
     &'stran(1)','stran(2)','stran(3)','stran(4)','stran(5)','stran(6)',
     & 'stress(1)','stress(2)','stress(3)','stress(4)','stress(5)',
     & 'stress(6)','statev(1)','statev(2)','statev(3)','statev(4)',
     & 'statev(5)','statev(6)','statev(7)','statev(8)','statev(9)',
     & 'statev(10)','statev(11)','statev(12)','statev(13)..'
      if(heading(1:1) /= '#') write(2,*) trim(heading)
      write(2,'(500(g14.8,3h    ))') time+(/dtime,dtime/),
     &                                stran, stress, statev

      do 200 ikeyword=1,10000  !  LOOP OVER KEYWORDS
        read(1,*,end=999) keywords(1)
        keywords(1) = trim( keywords(1) )
        if(keywords(1) == '*Repetition') then
           read(1,*) nSteps, nRepetitions
        else
           nRepetitions=1
           nSteps=1
           keywords(2) = keywords(1)
        endif
      do 130  iRepetition  = 1,nRepetitions


      do 120 iStep = 1,nSteps
        continue
        if(iRepetition > 1) then
            ninc             = ofStep(istep)%ninc
            maxiter          = ofStep(istep)%maxiter
            ifstress         = ofStep(istep)%ifstress
            deltaLoadCirc    = ofStep(istep)%deltaLoadCirc
            phase0           = ofStep(istep)%phase0
            deltaLoad        = ofStep(istep)%deltaLoad
            dfgrd0           = ofStep(istep)%dfgrd0
            dfgrd1           = ofStep(istep)%dfgrd1
            deltaTime        = ofStep(istep)%deltaTime
            keywords(2)      = ofStep(istep)%keyword2
            keywords(3)      = ofStep(istep)%keyword3
            goto 10  ! reading steps only if iRepetition==1
         endif

        if(keywords(1) == '*Repetition') read(1,*) keywords(2)          ! = LinearLoad  or CirculatingLoad
        keywords(2)  = trim(keywords(2))

        ifstress(:)=0                                                      ! default strain control
        deltaLoadCirc(:)=0.0d0                                             ! default zero step increment
        phase0(:)=0.0d0                                                    ! default no phase shift
        deltaLoad(:) = 0.0d0
        dfgrd0 = delta
        dfgrd1 = delta


         if(keywords(2) == '*DeformationGradient') then
           read(1,*) ninc, maxiter, deltaTime
           keywords(3) = '*Cartesian'
           do i=1,9
           read(1,*)  deltaLoad(i)                                      !  dload means total change in the whole step here
           enddo
           goto 10
        endif
        if (keywords(2) == '*CirculatingLoad') then
           read(1,*) ninc, maxiter, deltaTime
           read(1,*) keywords(3)                                           !  = Cartesian or Roscoe  or RoscoeIsomorph or Rendulic
            keywords(3)  = trim(keywords(3))
           do i=1,6
           read(1,*) ifstress(i),deltaLoadCirc(i),phase0(i),deltaLoad(i)   !  dload means amplitude here
           enddo
           goto 10
        endif
        if(keywords(2) == '*LinearLoad') then
           read(1,*) ninc, maxiter, deltaTime
           read(1,*) keywords(3)
           do i=1,6
           read(1,*) ifstress(i), deltaLoad(i)                             !  dload means total change in the whole step here
           enddo
           goto 10
        endif
        if(keywords(2) == '*OedometricE1') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Cartesian'
          read(1,*) ninc, maxiter, deltaTime
          read(1,*)   deltaLoad(1)
          goto 10
        endif
        if(keywords(2) == '*OedometricS1') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Cartesian'
          read(1,*) ninc, maxiter, deltaTime
          ifstress(1) = 1
          read(1,*)   deltaLoad(1)
          goto 10
        endif
        if(keywords(2) == '*TriaxialE1') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Cartesian'
          read(1,*) ninc, maxiter, deltaTime
          read(1,*) deltaLoad(1)
          ifstress(2:3) = 1
          goto 10
        endif
        if(keywords(2) == '*TriaxialS1') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Cartesian'
          read(1,*) ninc, maxiter, deltaTime
          read(1,*)   deltaLoad(1)
          ifstress(1:3) = 1
          goto 10
        endif
        if(keywords(2) == '*TriaxialUEq') then
          keywords(2) = '*LinearLoad'
          read(1,*) ninc, maxiter, deltaTime
          keywords(3) ='*Roscoe'
          read(1,*)   deltaLoad(2)                                      ! = deviatoric strain
          goto 10
        endif
        if(keywords(2) == '*TriaxialUq') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Roscoe'
          read(1,*) ninc, maxiter, deltaTime
          read(1,*)  deltaLoad(2)                                       ! = deviatoric stress
          ifstress(2) = 1
        goto 10
        endif
       if(keywords(2) == '*PureRelaxation') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Cartesian'
          read(1,*) ninc, maxiter, deltaTime
          goto 10
        endif
       if(keywords(2) == '*PureCreep') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Cartesian'
          read(1,*) ninc, maxiter, deltaTime
          ifstress(:) =  1
          goto 10
        endif
       if(keywords(2) == '*UndrainedCreep') then
          keywords(2) = '*LinearLoad'
          keywords(3) ='*Roscoe'
          read(1,*) ninc, maxiter, deltaTime
          ifstress(2:6) =  1
          goto 10
       endif
       if(keywords(2) == '*ObeyRestrictions') then  ! ======================= *ObeyRestrictions ==================================
           read(1,*) ninc, maxiter, deltaTime
           do i=1,6
           read(1,'(a)')  inputline(i)     !  \com a line of form  ''-sd1 + sd2 + 3.0*sd3 = -10  ! a comment '' is expected
           if(index(inputline(i),'=')== 0)stop 'I expect"=" '
           enddo
           call parser(inputline, cMt,cMe,mb )
           mbinc = mb/ninc
           keywords(3) ='*Cartesian'
           ifstress(1:6) = 1            !  \com  because  we solve (cMt.ddsdde + cMe).dstran = mbinc for dstran
          goto 10
       endif

       if(keywords(2) == '*PerturbationsS') then
          read(1,*) ninc, maxiter, deltaTime
          read(1,*) keywords(3)  ! = *Rendulic  or *RoscoeIsomorph
           keywords(3) = trim( keywords(3) )
          if(keywords(3) .ne. '*Rendulic' .and.
     &       keywords(3) .ne. '*RoscoeIsomorph')
     &       write(*,*) 'warning: non-Isomorphic perturburbation'
          read(1,*)  deltaLoad(1)
          ifstress(1:6) =  1
          goto 10
       endif
       if(keywords(2) == '*PerturbationsE') then
          read(1,*) ninc, maxiter, deltaTime
          read(1,*) keywords(3)
           keywords(3) = trim( keywords(3) )
          if(keywords(3) .ne. '*Rendulic' .and.
     &       keywords(3) .ne. '*RoscoeIsomorph')
     &       write(*,*) 'warning: non-Isomorphic perturburbation'
          read(1,*) deltaLoad(1)
          goto 10
       endif
       if(keywords(2) == '*End') then
         stop 'stopped because *End encountered in test.inp'
       endif

       write(*,*) 'unknown keywords(2)=',keywords(2)
       stop 'stopped by unknown keyword in test.inp'
  10  keywords(3) = trim(keywords(3))

       if(keywords(1) == '*Repetition' .and. iRepetition == 1) then      !  remember the description of step for the next repetition
       ofStep(istep)%ninc          =    ninc
       ofStep(istep)%maxiter       =    maxiter
       ofStep(istep)%ifstress      =    ifstress
       ofStep(istep)%deltaLoadCirc =    deltaLoadCirc
       ofStep(istep)%phase0        =    phase0
       ofStep(istep)%deltaLoad     =    deltaLoad
       ofStep(istep)%dfgrd0        =    dfgrd0
       ofStep(istep)%dfgrd1        =    dfgrd1
       ofStep(istep)%deltaTime     =    deltaTime
       ofStep(istep)%keyword2      =    keywords(2)
       ofStep(istep)%keyword3      =    keywords(3)
       endif

      if(any(ifstress==1)) maxiter = max(maxiter,5)! at least 5 iterations
      if(all(ifstress==0) .and. keywords(2) .ne. '*ObeyRestrictions')
     &                                                     maxiter = 1  ! no iterations are necessary

c     CALCULATE EACH STEP AS A SEQUENCE OF PRESCRIBED INCREMENTS  with the preliminary zero-load  call to umat to get the 1st stiffness
        dstran=0
        dtime=0
        dtemp=0
        call  UMAT(stress,statev,ddsdde,sse,spd,scd,                    ! \com  first call  with dstrain=0 dtime=0 just to get the stiffness
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)


       if(keywords(3) == '*Cartesian') then
        M =  MCart
        MmT = MCartmT
       endif

       if(keywords(3) == '*Roscoe') then
         M = MRosc
         MmT = MRoscmT
       endif

       if(keywords(3) == '*RoscoeIsomorph') then
         M = MRoscI
         MmT = MRoscImT
       endif

       if(keywords(3) == '*Rendulic') then
         M = MRendul
         MmT = MRendulmT
       endif


      do 100 kinc=1,ninc

        call get_increment(keywords, time, deltaTime, ifstress, ninc,   ! get inc. in terms of Rosc. variables
     &                         deltaLoadCirc,phase0,deltaLoad,
     &                         dtime, ddstress,  dstran, Qb33,
     &                         dfgrd0, dfgrd1,drot )


        a_dstress(:)= 0.0d0    ! approximated Roscoe's dstress
        r_statev(:)=statev(:)  ! remember the initial state & stress till the iteration is completed
        r_stress(:)= stress    ! remembered Cartesian stress

       do 95 kiter=1, maxiter  !--------''equilibrium iteration'' --------

       if(keywords(2)== '*ObeyRestrictions'  ) then  ! ======================= ObeyRestrictions ==================================
          ddsdde_bar = matmul(cMt,ddsdde) + cMe
          u_dstress = - matmul(cMt,a_dstress)-matmul(cMe,dstran)+ mbinc
          call  USOLVER(ddsdde_bar,c_dstran,u_dstress,ifstress,ntens)
          dstran = dstran + c_dstran
          call  UMAT(stress,statev,ddsdde,sse,spd,scd,
     &    rpl,ddsddt,drplde,drpldt,
     &    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,iStep,kinc)
          if (kiter.lt.maxiter) then                                    ! continue iteration
            statev(:)=r_statev(:)                                        ! 1) forget the changes of state done in umat
            a_dstress  = stress  - r_stress                             ! 2) compute the new approximation of stress
            stress(:)=r_stress(:)                                       ! 3) forget the changes of stress done in umat
          else
            stran(:)=stran(:)+dstran(:)                            !  accept  the updated state and stress (Cartesian)
          endif
       endif  ! ==== obey-restrictions

       if(keywords(2) /= '*ObeyRestrictions'  ) then   ! ======================= DISObeyRestrictions ==================================
         where (ifstress == 1)  u_dstress =ddstress -a_dstress           !  undesired Roscoe stress
         ddsdde_bar = matmul(matmul(M,ddsdde),transpose(M))  ! Roscoe-Roscoe stiffness
         call  USOLVER(ddsdde_bar,c_dstran,u_dstress,ifstress,ntens)     ! get Rosc. correction  c_dstran() caused by undesired Rosc. dstress
         where (ifstress == 1) dstran = dstran + c_dstran                ! corrected Rosc. dstran
         dstran_Cart = matmul( transpose(M),dstran )                     ! transsform Rosc. to Cartesian dstran
         call  UMAT(stress,statev,ddsdde,sse,spd,scd,
     &   rpl,ddsddt,drplde,drpldt,
     &   stran,dstran_Cart,time,dtime,temp,dtemp,predef,dpred,cmname,
     &   ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &   celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,iStep,kinc)

        if (kiter.lt.maxiter) then                                        ! continue iteration
           statev(:)=r_statev(:)                                          ! 1) forget the changes of state done in umat
           stress_Rosc = matmul(M,stress)                                !    output from umat transform to Roscoe ?
           r_stress_Rosc = matmul(M,r_stress)
           where (ifstress ==1) a_dstress = stress_Rosc - r_stress_Rosc ! 2) compute the new approximation of stress
           stress(:)=r_stress(:)                                         ! 3) forget the changes of stress done in umat
        else
          stran(:)=stran(:)+dstran_Cart(:)                              !  accept  the updated state and stress (Cartesian)
        endif
       endif  ! ==== DISObey-restrictions

94    continue
      if((kiter==maxiter) .and. mod(kinc,10)==0 ) then    ! write to screen
      write(*,'(12H ikeyword = ,i3, 8H kstep = ,i3,7H kinc = ,i5,
     &          9H kiter = ,i2)')  ikeyword, iStep, kinc, kiter
      endif     ! write to screen


 95   continue !--------------------end of iteration

      if(keywords(2) =='*DeformationGradient' ) then                    !  rigid rotation of stress
      T33 = map2T(stress,6)
      T33 = matmul( matmul(Qb33,T33),transpose(Qb33))
      stress=map2stress(T33,6)                                          ! rigid rotation of strain
      eps33 = map2D(stran,6)
      eps33 = matmul( matmul(Qb33,eps33),transpose(Qb33))
      stran=map2stran(eps33,6)
      endif
      
      ! chop because of fortran error writing 1.3E-391
      where(abs(time) < 1.0d-99)   time = 0.0d0
      where(abs(stran) < 1.0d-99)  stran = 0.0d0
      where(abs(stress) < 1.0d-99) stress = 0.0d0
      where(abs(statev) < 1.0d-99) statev = 0.0d0

      write(2,'(500(g14.8,3h    ))') time+(/dtime,dtime/),
     &                                stran, stress, statev

      if(keywords(2) =='*PerturbationsS' .or.
     &   keywords(2) =='*PerturbationsE' ) then ! having plotted everything undo the increment
       stran(:)=stran(:) - dstran_Cart(:)
       statev(:)=r_statev(:)
       stress(:)=r_stress(:)
      endif




      time(1)=time(1)+dtime     !  step time at the beginning of the increment
      time(2)=time(2)+dtime     !  total time at the beginning of the increment


  100 continue  ! next kinc
  120 continue  ! next iStep
  130 continue  ! next iRepetition
  200 continue  ! next keyword




 998  stop 'End or record encountered in test.inp'

 999  close(1)
      CALL CPU_TIME ( time_end )

      PRINT *, 'Time of operation was ',  
     & time_end - time_begin, ' seconds'
      PAUSE
      close(2)
      stop 'I have reached end of file test.inp'
 901  stop 'I cannot open file parameters.inp'
 902  stop 'I cannot open file initialconditions.inp'
 903  stop 'I cannot open file test.inp'
 904  stop 'I cannot open outputfile'





      end program that_calls_umat



!==========================================================================
      subroutine get_increment(keywords, time, deltaTime,ifstress,ninc,
     &                         deltaLoadCirc,phase0,deltaLoad,
     &                         dtime, ddstress,  dstran , Qb33,
     &                         dfgrd0, dfgrd1,drot )
      use unsymmetric_module
      implicit none
      character(20):: keywords(10)
      integer, intent(in)  :: ifstress(6),ninc
      real(8), intent(in) :: time(2), deltaTime,
     &                       deltaLoadCirc(6),phase0(6),
     &                       deltaLoad(9)
      real(8), intent(out) ::  dtime, ddstress(6), dstran(6), Qb33(3,3)
      real(8), intent(in out) ::  dfgrd0(3,3), dfgrd1(3,3), drot(3,3)


      real(8), parameter :: Pi = 3.1415926535897932385d0
      real(8),parameter,dimension(3,3):: delta =
     &                          reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      real(8),dimension(3,3):: Fb,Fbb, dFb,aux33,dLb,depsb,dOmegab
      real(8):: wd(6),  ! angular velocity (in future individual for each component)
     &          w0(6),  ! initial phase shift for a component
     &          t       ! step time
      integer(4) :: i
      logical :: ok

      dtime =  deltaTime/ ninc
      dstran= 0
      ddstress=0
      Qb33 = delta
      drot = delta
      dfgrd0=delta
      dfgrd1=delta

      !------------------------------------------------------
      if(keywords(2) == '*LinearLoad') then                               !  proportional loading
       do i=1,6
        if (ifstress(i)==1)   ddstress(i) = deltaLoad(i)/ ninc
        if (ifstress(i)==0)    dstran(i) = deltaLoad(i)/ ninc              ! log strain -> corresp. displac. inc. not constant
       enddo
       ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
       ! for dfgrd0 use stran
       ! for dfgrd1 use stran-dstran
      endif
      !--------------------------------------------------
      if(keywords(2) == '*DeformationGradient') then                     ! full deformation gradient.
      Fb = reshape((/deltaLoad(1), deltaLoad(5), deltaLoad(7),           ! finite rotations calculated after Hughes+Winget 1980
     &               deltaLoad(4), deltaLoad(2), deltaLoad(9),
     &               deltaLoad(6), deltaLoad(8), deltaLoad(3)
     &  /) ,  (/3,3/))
      Fbb = delta + (Fb-delta)*(time(1)/deltaTime)
      dfgrd0  = Fbb
      dFb = (Fb-delta)/ninc
      aux33 =  Fbb + dFb/2.0d0
      dfgrd1   = Fbb  + dFb

      call matrix('inverse', aux33, 3, ok )
      dLb =  matmul(dFb,aux33)
      depsb = 0.5d0*(dLb + transpose(dLb))
      dstran=(/depsb(1,1), depsb(2,2),depsb(3,3),
     &         2.0d0*depsb(1,2),2.0d0*depsb(1,3),2.0d0*depsb(2,3)/)
      dOmegab =    0.5d0*(dLb - transpose(dLb))
      aux33 =  delta - 0.5d0*dOmegab
      call matrix('inverse', aux33, 3, ok )
      Qb33 = matmul(aux33, (delta+0.5d0*dOmegab))
      drot=Qb33
      endif
      !------------------------------------------------------
      if(keywords(2) == '*CirculatingLoad' )then                         !  harmonic oscillation
      wd(:) = 2*Pi/deltaTime
      w0 = phase0
      t= time(1)  + dtime/2   ! step time in the middle of the increment
      do i=1,6
      if(ifstress(i)==1)
     &  ddstress(i)=dtime*deltaLoadCirc(i)*wd(i)*Cos(wd(i)*t+w0(i))+
     &              deltaLoad(i)/ ninc
      if(ifstress(i)==0)dstran(i)=
     &               dtime*deltaLoadCirc(i)*wd(i)*Cos(wd(i)*t+w0(i))+
     &               deltaLoad(i)/ ninc
      enddo
       ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
       ! for dfgrd0 use stran
       ! for dfgrd1 use stran-dstran
      endif

      !--------------------------------------------------------
      if(keywords(2) == '*PerturbationsS' )then
      ddstress(1)= deltaLoad(1)*cos( time(1)*2*Pi/deltaTime )
      ddstress(2)= deltaLoad(1)*sin( time(1)*2*Pi/deltaTime  )
       ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
       ! for dfgrd0 use stran
       ! for dfgrd1 use stran-dstran
      endif

       !--------------------------------------------------------
      if(keywords(2) == '*PerturbationsE' )then
      dstran(1)= deltaLoad(1)*cos( time(1)*2*Pi/deltaTime )
      dstran(2)= deltaLoad(1)*sin( time(1)*2*Pi/deltaTime  )
       ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
       ! for dfgrd0 use stran
       ! for dfgrd1 use stran-dstran
      endif

      return
      end subroutine get_increment

      ! SOME  OF THE UTILITY ROUTINES OFTEN USED BY UMATS


      SUBROUTINE ROTSIG(S,R,SPRIME,LSTR,NDI,NSHR)                       !  some umats depend on the utility subroutine  ROTSIG  provided in abaqus.
      use unsymmetric_module                                            !  This is  an imitation of this subroutine.
      integer, intent(in) ::  LSTR,NDI,NSHR
      integer :: ntens
      real(8), dimension(3,3),intent(in) ::  R
      real(8), dimension(1:NDI+NSHR), intent(in) :: S
      real(8), dimension(1:NDI+NSHR) , intent(out):: SPRIME
      real(8), dimension(3,3) :: T33

      ntens = ndi+nshr
      if(LSTR==1) T33 = map2T(S,ntens)
      if(LSTR==0) T33 = map2D(S,ntens)

      T33 = matmul( matmul(R,T33),transpose(R))

      if(LSTR==1) SPRIME = map2stress(T33,ntens)
      if(LSTR==0) SPRIME = map2stran(T33,ntens)
      return
      END  SUBROUTINE ROTSIG



                                                                        !  some umats depend on the utility subroutine SINV provided in abaqus.
      subroutine SINV(STRESS,SINV1,SINV2,NDI,NSHR)                      !  This is  an imitation of this subroutine.
      real(8),intent(in) :: STRESS(NDI+NSHR)                            !  SINV returns  two invariants  of  stress
      real(8),intent(out) ::  SINV1,SINV2
      integer, intent(in) ::  NDI,NSHR
      real(8) :: devia(NDI+NSHR)
      real(8), parameter :: sq2 = 1.4142135623730950488d0
      if(NDI /= 3) stop 'stopped because ndi/=3 in sinv'
      sinv1 = (stress(1) + stress(2) + stress(3) )/3.0d0
      devia(1:3) = stress(1:3) - sinv1
      devia(3+1:3+nshr) = stress(3+1:3+nshr) * sq2
      sinv2 = sqrt(1.5d0 *  dot_product(devia, devia)  )
      end subroutine SINV

                                                                        !  some umats depend on the utility subroutine SPRINC  provided in abaqus.
      subroutine SPRINC(S,PS,LSTR,NDI,NSHR)                             !  This is  an imitation of this subroutine.
      real(8),intent(in) :: S(NDI+NSHR)
      real(8),intent(out) :: PS(NDI+NSHR)
      integer, intent(in) :: LSTR,NDI,NSHR                               !  LSTR == 1 ->   stress  or    LSTR == 2 ->   strain
      real(8):: A(3,3),AN(3,3)
      real(8) :: r(6)
      if(NDI /= 3) stop 'stopped because ndi/=3 in sprinc'
      r(1:3) = s(1:3)
      if(LSTR == 1 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)
      if(LSTR == 2 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)/2
      A= reshape((/r(1),r(4),r(5),r(4),r(2),r(6),r(5),r(6),r(3)/),
     &            (/3,3/))
      call spectral_decomposition_of_symmetric(A, PS, AN, 3)
      return
      end subroutine SPRINC


                                                                        !  some umats depend on the utility subroutine SPRIND  provided in abaqus.
      subroutine SPRIND(S,PS,AN,LSTR,NDI,NSHR)                         !  This is  an imitation of this subroutine.
      real(8),intent(in) :: S(NDI+NSHR)
      real(8),intent(out) :: PS(3),AN(3,3)
      integer, intent(in) :: LSTR,NDI,NSHR                               !  LSTR == 1 ->   stress  or    LSTR == 2 ->   strain
      real(8):: A(3,3)
      real(8) :: r(6)
      if(NDI /= 3) stop 'stopped because ndi/=3 in sprind'
      r(1:3) = s(1:3)
      if(LSTR == 1 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)
      if(LSTR == 2 .and. nshr > 0) r(4:3+nshr) = s(4:3+nshr)/2
      A= reshape((/r(1),r(4),r(5),r(4),r(2),r(6),r(5),r(6),r(3)/),
     &            (/3,3/))
      call spectral_decomposition_of_symmetric(A, PS, AN, 3)
      return
      end subroutine SPRIND



                                                                        !  some umats depend on the utility subroutine SPRIND  provided in abaqus.
      subroutine XIT                                                    !  This is  an imitation of this subroutine.
      stop 'stopped because umat called XIT'
      end subroutine XIT





       SUBROUTINE  spectral_decomposition_of_symmetric(A, Lam, G, n)
       implicit none
       integer, intent(in) :: n                                         ! size of the matrix
       real(8), INTENT(in)  :: A(n,n)                                   ! symmetric input matrix  n x n   (not destroyed in this routine)
       real(8), INTENT(out)  :: Lam(n)                                  ! eigenvalues
       real(8), INTENT(out)  :: G(n,n)                                  ! corresponding eigenvectors in columns of G
       integer ::  iter,i,j, p,q
       real(8) ::   cosine, sine
       real(8), dimension(:), allocatable :: pcol ,qcol
       real(8), dimension(:,:), allocatable :: x

       allocate(pcol(n) ,qcol(n), x(n,n) )
       x = A
       G=0.0d0
       do i=1,n
       G(i,i) = 1.0d0
       enddo

      do  iter = 1,30
        call  get_jacobian_rot(x, p ,q, cosine, sine, n)                !  find how to apply  optimal similarity  mapping
        call  app_jacobian_similarity(x, p,q, cosine, sine, n)          !   perform mapping

        pcol = G(:,p)                                                   !  collect rotations  to global similarity matrix
        qcol = G(:,q)
        G(:,p) =   pcol*cosine - qcol*sine
        G(:,q) =   pcol* sine + qcol *cosine

        ! here write a problem-oriented accuracy test max_off_diagonal < something
        ! but 30 iterations are usually ok for 3x3 stress or 6x6 stiffness matrix
      enddo

      do i=1,n
       Lam(i) = x(i,i)                                                  !  eigenvalues
      enddo
      deallocate( pcol ,qcol, x )
      return
      end


      SUBROUTINE  app_jacobian_similarity(A, p,q, c, s, n)              !  jacobian similarity tranformation of a square symmetric matrix A
      implicit none                                                     !  ( $ A : =  G^T .A . G $ with    Givens   rotation  G_pq = {{c,s},{-s,c}}  )
      INTEGER, INTENT(IN)        :: p,q                                 !   G is an identity n x n matrix overridden with  values  {{c,s},{-s,c}}  )
      real(8), INTENT(IN)        :: c ,s                                !   in cells {{pp, pq},{qp,qq}}  algorithm according to Kielbasinski  p.385
      integer , INTENT(IN)       :: n
      real(8), dimension(n,n),intent(inout) :: A
      real(8), dimension(n)  :: prow ,qrow
      real(8) :: App, Apq, Aqq



      if(p == q)  stop 'error: jacobian_similarity  p == q'
      if(p<1 .or. p>n) stop 'error: jacobian_similarity p out of range'
      if(q<1 .or. q>n) stop 'error: jacobian_similarity q out of range'

      prow(1:n) = c*A(1:n,p) - s*A(1:n,q)
      qrow(1:n) = s*A(1:n,p) + c*A(1:n,q)
      App = c*c*A(p,p) -2*c*s*A(p,q) + s*s*A(q,q)
      Aqq =  s*s*A(p,p) +2*c*s*A(p,q) + c*c*A(q,q)
      Apq = c*s*(A(p,p) - A(q,q)) + (c*c - s*s)* A(p,q)
      A(p,1:n) =   prow(1:n)
      A(1:n,p) =   prow(1:n)
      A(q,1:n) =   qrow(1:n)
      A(1:n,q) =   qrow(1:n)
      A(p,p) =   App
      A(q,q) =   Aqq
      A(p,q) =   Apq
      A(q,p) =   Apq

      END SUBROUTINE  app_jacobian_similarity



      SUBROUTINE  get_jacobian_rot(A, p,q, c, s, n)          !   \com returns jacobian similarity  tranformation param.
      implicit none                                          !  \com  for iterative diagonalization of  a square symm.  A
      integer , INTENT(IN)                :: n               !  \com  algorithm according to Kielbasinski 385-386
      real(8), dimension(n,n),intent(in) :: A
      INTEGER, INTENT(OUT)                :: p,q
      real(8), INTENT(OUT)                :: c ,s
      real(8) :: App, Apq, Aqq, d, t,maxoff
      integer ::   i,j

      p = 0
      maxoff  = tiny(maxoff)
      do i=1,n-1
      do j=i+1,n
       if( abs(A(i,j)) > maxoff ) then
         maxoff = abs(A(i,j))
         p=i
         q=j
       endif
      enddo
      enddo
      if (p > 0) then
        App = A(p,p)
        Apq = A(p,q)
        Aqq = A(q,q)
        d = (Aqq - App)/ (2.0d0*Apq)
        t = 1.0d0/ sign(abs(d) + sqrt(1.0d0 + d*d) , d )
        c = 1.0d0/sqrt(1.0d0 + t*t)
        s = t*c
      else                                                              ! \com  no rotation
        p=1
        q=2
        c=1
        s=0
      endif
      end subroutine get_jacobian_rot


      subroutine  parser(inputline, Mt,Me,mb)

      character(260), intent(in) ::  inputline(6)
      real(8), dimension(6,6), intent(out) :: Mt , Me
      real(8), dimension(6),intent(out) :: mb

      character(len=260) ::  inp, aux,aux3
      character(40) ::  summand(13)
      integer :: iis,i,iplus,iminus,iequal,imin,iex,itimes,Irestr
      real(8) :: factor(13)


       Mt = 0; Me= 0; mb= 0

           Do Irestr = 1,6! Irestr loop over restriction lines

                 inp = trim(adjustl(inputline(Irestr)))
                 iis = index(inp,'=')
                 if(iis==0) stop 'error: I expect ='

                 factor(1) =1;
                 if(inp(1:1) == '-') then ! do not treat the first minus as a separator
                   factor(1) = -1
                   inp=inp(2:)            ! remove the first character = '-' from inp
                 endif
              do i=1,13 ! loop over possible summands
                  iplus = index(inp,'+');  if(iplus==0) iplus=200
                  iminus = index(inp,'-');  if(iminus==0) iminus=200
                  iequal = index(inp,'=');  if(iequal==0) iequal=200
                  imin = min(iplus,iminus,iequal) ! choose the first separator
                  if(imin==200)  stop 'error: I expect +,- or ='
                  if(imin==iplus) then ! separator= '+'  everything left from + save as summand
                    summand(i) = inp(:imin-1) ; factor(i+1) = 1
                    inp = inp(imin+1:)
                  endif
                  if(imin==iminus) then ! separator= '+'  everything left from + save as summand
                    summand(i) = inp(:imin-1);  factor(i+1) = -1
                    inp = inp(imin+1:)
                  endif
                  if(imin==iequal) then ! separator= '='  everything left from + save as summand
                    summand(i) = inp(:imin-1);
                    inp = inp(imin+1:)  !  rhs possibly with sign
                    iminus = index(inp,'-');  if(iminus==0) iminus=200
                    iex = index(inp,'!');   if(iex==0) iex=len(inp)+1  ! right limit = comment or EOL
                    if(iminus == 200) then     ! '=' is not followed by  '-'
                      factor(i+1) = 1
                      summand(i+1) = inp(:iex-1)
                    else                       ! double separator: '=' followed by '-'
                      factor(i+1) = -1
                      summand(i+1) = inp(iminus+1:iex-1)
                    endif
                    exit   ! reading a single summand after '=' ends reading of the line
                  endif
              enddo ! i-loop
       Nsummands=i+1  ! summand()=LHS, summand(Nsummands)=RHS, signs in factor()

          Do i=1,Nsummands-1  ! for summands on the LHS
           aux =   adjustl( summand(i) )
           itimes = index(aux,'*')
           if(itimes /= 0) then ! if exists '*' then split the summand into factor and component
           read(aux(:itimes-1),*) fac   ! numeric factor  of the summand
           factor(i) = factor(i)*fac   ! the signed numeric factor  of the summand
           aux=adjustl(aux(itimes+1:))
           endif
           aux3 = aux(1:3)
           select case(aux3)
            case ('sd1')
            Mt(Irestr,1) = factor(i)
            case ('sd2')
            Mt(Irestr,2) = factor(i)
            case ('sd3')
            Mt(Irestr,3) = factor(i)
            case ('sd4')
            Mt(Irestr,4) = factor(i)
            case ('sd5')
            Mt(Irestr,5) = factor(i)
            case ('sd6')
            Mt(Irestr,6) = factor(i)
            case ('ed1')
            Me(Irestr,1) = factor(i)
            case ('ed2')
            Me(Irestr,2) = factor(i)
            case ('ed3')
            Me(Irestr,3) = factor(i)
            case ('ed4')
            Me(Irestr,4) = factor(i)
            case ('ed5')
            Me(Irestr,5) = factor(i)
            case ('ed6')
            Me(Irestr,6) = factor(i)
           end select
        enddo
        read(summand(Nsummands) ,*) mb(Irestr)     ! \com RHS numeric without sign
        mb(Irestr) =  mb(Irestr)*factor(Nsummands)  ! \com RHS numeric with sign
      enddo ! Irestr
      end subroutine parser




