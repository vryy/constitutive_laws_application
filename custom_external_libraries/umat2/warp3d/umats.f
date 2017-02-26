c
c           Updated:  4/17/12
c
c           Example UMATs from Abaqus implemented using the 
c           inteface provided by WARP3D. 
c
c           Exisiting UMAT codes should work unchanged.
c
c           Examples included here:
c
c           1  - linear-elastic, isotropic, isothermal
c           2  - mises plasticity, bilinear stress-strain curve.
c                kinematic hardening. temperature effects handled
c                inside umat
c           3  - mises model with isotropic hardening via user 
c                supplied points on segmental curve in inp.
c                umat let's WARP3D handle temperature effects.
c           4  - simple neo-hookean material. uses deformation 
c                gradients rather than strain increment passed in.
c
c           Example:
c           Comments at the end of this file have an
c           example model with a user material defined with each 
c           umat above. 
c
c
c           UMAT names have _1, _2 ... appended here. 
c           Make the one to be used by WARP3D have the name "umat".
c
c           The "service" routines that Abaqus provides for use by UMATs
c           are included at the bottom of this file. Right now that is
c           only "rotsig" and "xit". Others to be added shortly.
c
c           The first routine here, umat_set_sizes, is called by WARP3D
c           at various times to obtain the number of entries 
c           in "statev".
c
c           Set the value required by your umat. Add 21 for use by
c           WARP3D as shown below.
c
c           This file is compiled and included in the normal executable
c           build of WARP3D.
c
c           WARP3D calls umat routines in:
c
c              gplns1.f    -- linear stiffness computation
c              rstgp1.f    -- stress update and new [D] computation
c
c           These two code fies have lots of comments about setting
c           up data arrays, values for the umat from WARP3D data 
c           structures.
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine umat_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 6/17/12                     *
c     *                                                              *
c     *  called by warp3d for the umat to obtain statev size and     *
c     *  other characteristic information about the umat             *
c     *                                                              *
c     ****************************************************************
c
      subroutine umat_set_sizes( info_vector )
      implicit none
      integer info_vector(*)
c
      integer current_umat_example
c
c        set infor_data
c
c         1        number of history values per integration 
c                  point. Abaqus calles these "statev". Values
c                  double or single precsion based on hardware.
c    
c         2        number of values in the symmetric part of the 
c                  [D] for each integration point. for solid
c                  elements this is 21, for cohesive elements this 6.
c
c         3        = 0, the material model returns "unrotated"
c                       Cauchy stresses at n+1
c                  = 1, the material model returns the standard
c                       Cauchy stresses at n+1
c
c
      current_umat_example = 2 
c
      select case( current_umat_example ) 
c
       case( 1 )  !  linear-elastic, isotropic
        info_vector(1) = 1
        info_vector(2) = 21
        info_vector(3) = 0
c
       case( 2 )  !  mises plasticity, bilinear kinematic hardening
        info_vector(1) = 20  !!!!!!
        info_vector(2) = 21
        info_vector(3) = 0
c
       case( 3 ) ! mises plasticity, segmental isotropic hardening
        info_vector(1) = 13  !!!!!
        info_vector(2) = 21
        info_vector(3) = 0
c
       case( 4 ) ! neo-hookean material
        info_vector(1) = 6   !!!!!
        info_vector(2) = 21
        info_vector(3) = 1
c
       case default
        write(*,9000) 
        call die_abort
c
      end select
c
      return
c
 9000 format(".. Fatal error:  umat_set_sizes, ",
     & /,    "                 invalid current_umat_example",
     & /,    "                 execution terminated" )
c
      end

c     ****************************************************************
c     *                                                              *
c     *           subroutine umat: linear elasticity                 *
c     *                                                              *
c     *   Abaqus example umat for isothermal, linear elasticity      *
c     *                                                              *
c     ****************************************************************
c
c      
      subroutine umat_1( stress, statev, ddsdde, sse, spd, scd, rpl,
     1 ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,
     2 dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, 
     3 props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, 
     4 noel, npt, layer, kspt, kstep, kinc, kiter, kout )
c
      implicit real*8(a-h,o-z)
      parameter (nprecd=2)
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     1          enumax=.4999d0, newton=10, toler=1.0d-6)
      character*8 cmname
      dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 predef(1), dpred(1), props(nprops), coords(3), drot(3,3),
     3 dfgrd0(3,3), dfgrd1(3,3)

c ----------------------------------------------------------------
c    umat for isotropic elasticity
c    cannot be used for plane stress
c ----------------------------------------------------------------
c    props(1) - e
c    props(2) - nu
c ----------------------------------------------------------------
c
c              elastic properties
c 
      emod   = props(1)
      enu    = props(2)
      ebulk3 = emod/(one - two*enu)
      eg2    = emod/(one + enu)
      eg     = eg2/two
      eg3    = three*eg
      elam   = (ebulk3 - eg2)/three
c 
c               elastic stiffness
c
      do k1 = 1, 3
        do k2 = 1, 3
          ddsdde(k2,k1) = elam
        end do
        ddsdde(k1,k1) = eg2 + elam
      end do
c
      ddsdde(4,4) = eg
      ddsdde(5,5) = eg
      ddsdde(6,6) = eg
c
c               calculate stress
c
      do k1 = 1, 6
        do k2 = 1, 6
          stress(k2) = stress(k2) + ddsdde(k2,k1)*dstran(k1)
        end do
      end do

c      write(kout,*) '.... leaving UMAT....'
c      write(kout,9120) stress(1:6)
c 9120 format(5x,'stresses: ',6e14.6)
c
      return   
      end
c     ****************************************************************
c     *                                                              *
c     *   subroutine umat -- bilinear kinematic hardening            *
c     *                                                              *
c     *  Abaqus example umat for elastic-plastic response with       *
c     *  linear kinematic hardening                                  *
c     *                                                              *
c     ****************************************************************
      subroutine umat( stress, statev, ddsdde, sse, spd, scd, rpl,
     1 ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,
     2 dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, 
     3 props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, 
     4 noel, npt, layer, kspt, kstep, kinc, kiter, kout )
c
      implicit real*8(a-h,o-z)
      parameter (nprecd=2)
c
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     &     enumax=.4999d0, toler=1.0d-6)
c
      character*8 cmname
      dimension stress(*), statev(nstatv), ddsdde(ntens, ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 predef(1), dpred(1), props(nprops), coords(3), drot(3,3),
     3 dfgrd0(3,3), dfgrd1(3,3)

      
c
c  
c           local arrays
c
c              eelas  - elastic strains
c              eplas  - plastic strains
c              alpha  - shift tensor
c              flow   - plastic flow directions
c              olds   - stress at start of increment
c              oldpl  - plastic strains at start of increment
c
      dimension eelas(6), eplas(6), alpha(6), flow(6), olds(6), 
     &          oldpl(6)
      logical debug 
c
c ----------------------------------------------------------------
c           umat for isotropic elasticity and mises plasticity
c           with kinematic hardening
c
c                 props(1) - e
c                 props(2) - nu
c                 props(3) - syield
c                 props(4) - plastic hardening modulus (H')
c                 props(5) - coefficient of thermal expansion
c
c                 statev requries 20 entries per material point
c                        for current model
c
c ----------------------------------------------------------------
c
c 
      if( debug ) write(kout,*) ' .... '
      cte = props(5)
      dstran(1) = dstran(1) - cte * dtemp
      dstran(2) = dstran(2) - cte * dtemp
      dstran(3) = dstran(3) - cte * dtemp
      if( kiter .eq. 0 ) return
c
c           elastic material constants
c
      emod = props(1)
      enu = min(props(2), enumax)
      ebulk3 = emod/(one-two*enu)
      eg2 = emod/(one+enu)
      eg = eg2/two
      eg3 = three*eg
      elam = (ebulk3-eg2)/three
c
c           linear-elastic material stiffness
c
      do k1 = 1, 3
        do k2 = 1, 3
          ddsdde(k2,k1) = elam
        end do
        ddsdde(k1,k1) = eg2 + elam
      end do
      ddsdde(4,4) = eg
      ddsdde(5,5) = eg
      ddsdde(6,6) = eg
c 
c           recover elastic strain, plastic strain and shift tensor and rotate note: 
c           use code 1 for (tensor) stress, code 2 for (engineering) strain
c           recover plastic strain at n. set current yielding state.
c           0=never yielded, 1=active plasticity, 2=prior plasticity by
c           currently linear elastic
c
      call rotsig( statev(1),  drot, eelas, 2, 3, 3)
      call rotsig( statev(7),  drot, eplas, 2, 3, 3)
      call rotsig( statev(13), drot, alpha, 1, 3, 3)
      old_eqpl   = statev(19)
      deqpl      = zero
      state_np1  = zero
      if( old_eqpl .gt. zero ) state_np1 = two
c
c           save stress and plastic strains and calculate predictor 
c           (trial elastic) stress and elastic strain. the original
c           umat coding from Abaqus PPT slides for next few loops
c           had an error.
c
      do k1 = 1, 6
         olds(k1)  = stress(k1)
         oldpl(k1) = eplas(k1)
      end do
c
      do k1 = 1, 6
         eelas(k1) = eelas(k1) + dstran(k1)
         do k2 = 1, 6
           stress(k2) = stress(k2) + ddsdde(k2,k1)*dstran(k1)
         end do
      end do
c
c           calculate equivalent von mises stress
c
      smises = (stress(1)-alpha(1)-stress(2)+alpha(2))**2
     1      + (stress(2)-alpha(2)-stress(3)+alpha(3))**2
     2      + (stress(3)-alpha(3)-stress(1)+alpha(1))**2
      smises = smises+six*(stress(4)-alpha(4))**2
      smises = smises+six*(stress(5)-alpha(5))**2
      smises = smises+six*(stress(6)-alpha(6))**2
      smises = sqrt(smises/two)
c 
c           get yield stress and plastic hardening modulus
c
      syield = props(3)
      hard   = props(4)
c 
c           determine if actively yielding
c
      if( smises .gt. (one+toler)*syield ) then
c
c              actively yielding
c                separate the hydrostatic from the deviatoric stress
c                calculate the flow direction
c
         shydro  = (stress(1)+stress(2)+stress(3))/three
         flow(1) = (stress(1)-alpha(1)-shydro)/smises
         flow(2) = (stress(2)-alpha(2)-shydro)/smises
         flow(3) = (stress(3)-alpha(3)-shydro)/smises
         flow(4) = (stress(4)-alpha(4))/smises
         flow(5) = (stress(5)-alpha(5))/smises
         flow(6) = (stress(6)-alpha(6))/smises
c
c                 solve for equivalent plastic strain increment
c
         deqpl = (smises-syield)/(eg3+hard)
         state_np1 = one
c
c                 update shift tensor, elastic and plastic strains 
c                 and stress
c
         do k1 = 1, 3
           alpha(k1)  = alpha(k1)+hard*flow(k1)*deqpl
           eplas(k1)  = eplas(k1)+three/two*flow(k1)*deqpl
           eelas(k1)  = eelas(k1)-three/two*flow(k1)*deqpl
           stress(k1) = alpha(k1)+flow(k1)*syield+shydro
         end do
         do k1 =  4, 6
           alpha(k1)  = alpha(k1)+hard*flow(k1)*deqpl
           eplas(k1)  = eplas(k1)+three*flow(k1)*deqpl
           eelas(k1)  = eelas(k1)-three*flow(k1)*deqpl
           stress(k1) = alpha(k1)+flow(k1)*syield
         end do
c
c                 calculate increment of plastic dissipation
c
         dspd = zero
         do k1 = 1, 6
           avg_stress = ( stress(k1) + olds(k1) ) / two
           dspd = dspd +  avg_stress * ( eplas(k1) - oldpl(k1) )
         end do
c
c                 formulate the jacobian (material tangent)
c                 first calculate effective moduli
c
         effg  = eg*(syield+hard*deqpl)/smises
         effg2  = two*effg
         effg3  = three*effg
         efflam = (ebulk3-effg2)/three
         effhrd = eg3*hard/(eg3+hard)-effg3
         do k1 = 1, 3
           do k2 = 1, 3
             ddsdde(k2,k1) = efflam
           end do
           ddsdde(k1,k1) = effg2 + efflam
         end do
         ddsdde(4,4) = effg
         ddsdde(5,5) = effg
         ddsdde(6,6) = effg
         do k1 = 1, 6
           do k2 = 1, 6
             ddsdde(k2,k1) = ddsdde(k2,k1) + effhrd*flow(k2)*flow(k1)
           end do
         end do
      end if
c
c          store elastic strains, plastic strains, shift tensor
c          and updated (scalar) plastic strain in state variable array. 
c
      dsst = zero
      do k1 = 1, 6
         statev(k1)    = eelas(k1)
         statev(k1+6)  = eplas(k1)
         statev(k1+12) = alpha(k1)
         avg_stress =  ( stress(k1) + olds(k1) ) / two
         dsst = dsst + avg_stress * dstran(k1) 
      end do
c
c          store updated scalar plastic strain. material state
c          = 0 (never yielded), = 1 actively yielding, = 2
c          (previously yielded but currently linear-elastic)
c          these values can be output by WARP3D. see
c          optional ou_umat routine.
c
      statev(19) = old_eqpl + deqpl
      statev(20) = state_np1 
c
c          update specific elastic energy and dissipatation.
c
      sse = sse + (dsst - dspd)
      spd = spd + dspd
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *        subroutine umat - mises isotropic hardening           *
c     *        segmental stress-strain curve                         *
c     *                                                              *
c     ****************************************************************
      subroutine umat_3( stress, statev, ddsdde, sse, spd, scd, rpl,
     1 ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,
     2 dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, 
     3 props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, 
     4 noel, npt, layer, kspt, kstep, kinc, kiter, kout )
c
      implicit real*8(a-h,o-z)
      parameter (nprecd=2)
c
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     1      enumax=.4999d0, newton=10, toler=1.0d-6)
c
      character*8 cmname
      dimension stress(*), statev(nstatv), ddsdde(ntens,ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 predef(1), dpred(1), props(nprops), coords(3), drot(3,3),
     3 dfgrd0(3,3), dfgrd1(3,3)
c
c    local arrays
c ----------------------------------------------------------------
c    eelas  - elastic strains
c    eplas  - plastic strains
c    flow   - direction of plastic flow
c ----------------------------------------------------------------
c

      dimension eelas(6), eplas(6), flow(6), hard(3)

c
c ----------------------------------------------------------------
c    umat for isotropic elasticity and isotropic mises plasticity
c    cannot be used for plane stress
c ----------------------------------------------------------------
c    props(1) - e
c    props(2) - nu 
c    propos(3) - num point pairs on curve
c    props(4..) - syield an hardening data
c    calls uhard for curve of yield stress vs. plastic strain
c ----------------------------------------------------------------
c
c              WARP3D handles thermal strains for this umat.
c              No adjustment needed for dstran. Just return.
c
      if( kiter .eq. 0 ) return
c
c              elastic properties
c
      emod = props(1)
      enu = min(props(2), enumax)
      ebulk3 = emod/(one-two*enu)
      eg2 = emod/(one+enu)
      eg = eg2/two
      eg3 = three*eg
      elam = (ebulk3-eg2)/three
c 
c              elastic stiffness
c
      do k1 = 1, 3
        do k2 = 1, 3
           ddsdde(k2,k1) = elam
        end do
        ddsdde(k1,k1) = eg2 + elam
      end do
      ddsdde(4,4) = eg
      ddsdde(5,5) = eg
      ddsdde(6,6) = eg
c
c             recover elastic and plastic strains and rotate forward
c             also recover equivalent plastic strain
c
      call rotsig( statev(1), drot, eelas, 2, 3, 3 )
      call rotsig( statev(7), drot, eplas, 2, 3, 3 )
      eqplas = statev(13)
c 
c              calculate predictor stress and elastic strain
c
      do k1 = 1, 6
        do k2 = 1, 6
         stress(k2) = stress(k2) + ddsdde(k2,k1)*dstran(k1)
        end do
        eelas(k1) = eelas(k1) + dstran(k1)
      end do
c 
c              calculate equivalent von mises stress
c
      smises = (stress(1)-stress(2))**2
     1      + (stress(2)-stress(3))**2
     2      + (stress(3)-stress(1))**2
      smises = smises+six*stress(4)**2
      smises = smises+six*stress(5)**2
      smises = smises+six*stress(6)**2
      smises = sqrt(smises/two)
c
c               get yield stress from the specified hardening curve
c
      nvalue = int(props(3))
      call uhard( syiel0, hard, eqplas, eqplasrt, time,
     &     dtime, temp,
     &     dtemp, noel, npt, layer, kspt, kstep, kinc, cmname, nstatv,
     &     statev, numfieldv, predef, dpred, nvalue, props(4), kout )
c
c               determine if actively yielding
c
c      write(kout,*) '... syiel0: ',  syiel0
c      write(kout,*) '... smises: ', smises
      if( smises .le. (one+toler)*syiel0 ) go to 1000
c
c               actively yielding
c
c               separate the hydrostatic from the deviatoric stress
c               calculate the flow direction
c
      shydro = (stress(1)+stress(2)+stress(3))/three
      do k1 = 1, 3
        flow(k1) = (stress(k1)-shydro)/smises
      end do
      do k1 = 4, 6
        flow(k1) = stress(k1)/smises
      end do
c
c               solve for equivalent von mises stress
c               and equivalent plastic strain increment using newton iteration
c
      syield = syiel0
      deqpl = zero
      do kewton = 1, newton
         rhs = smises - eg3*deqpl-syield
         deqpl = deqpl + rhs/(eg3+hard(1))
         call uhard( syield, hard, eqplas+deqpl, eqplasrt, time,
     &     dtime, temp, 
     &     dtemp, noel, npt, layer, kspt, kstep, kinc, cmname, nstatv, 
     &     statev, numfieldv, predef, dpred, nvalue, props(4), kout )
         if( abs(rhs) .lt. toler*syiel0 ) goto 10
      end do
      write(kout,2) newton
 2    format(//,30x,'***warning - plasticity algorithm did not ',
     &  'converge after ',i3,' iterations' )
10    continue
c
c               update stress, elastic and plastic strains and
c               equivalent plastic strain
c
      do k1 = 1, 3
        stress(k1) = flow(k1)*syield+shydro
        eplas(k1) = eplas(k1)+three/two*flow(k1)*deqpl
        eelas(k1) = eelas(k1)-three/two*flow(k1)*deqpl
      end do
      do k1 = 4, 6
        stress(k1) = flow(k1)*syield
        eplas(k1) = eplas(k1)+three*flow(k1)*deqpl
        eelas(k1) = eelas(k1)-three*flow(k1)*deqpl
      end do
      eqplas = eqplas + deqpl
c
c               calculate plastic dissipation
c
      spd = deqpl*(syiel0+syield)/two 
c
c----------------
c               formulate the jacobian (material tangent)
c               first calculate effective moduli
c
      effg = eg*syield/smises
      effg2 = two*effg
      effg3 = three/two*effg2
      efflam = (ebulk3-effg2)/three
      effhrd = eg3*hard(1)/(eg3+hard(1))-effg3
      do k1 = 1, 3
        do k2 = 1, 3
           ddsdde(k2, k1) = efflam
        end do
        ddsdde(k1, k1) = effg2+efflam
      end do
      do k1 = 4, 6
        ddsdde(k1,k1) = effg
      end do
      do k1 = 1, 6
        do k2 = 1, 6
           ddsdde(k2,k1) = ddsdde(k2,k1)+effhrd*flow(k2)*flow(k1)
        end do
      end do
c
 1000 continue
c
c                store elastic and (equivalent) plastic strains
c                in state variable array
c
      do k1 = 1, 6
        statev(k1) = eelas(k1)
        statev(k1+6) = eplas(k1)
      end do
      statev(13) = eqplas 
c
      return
      end

      subroutine uhard( syield, hard, eqplas, eqplasrt, time, dtime,
     &     temp, dtemp, noel, npt, layer, kspt, kstep, kinc, 
     &     cmname, nstatv, statev, numfieldv, 
     &     predef, dpred, nvalue, table, kout ) 
      implicit real*8(a-h,o-z)
      parameter (nprecd=2)
c
      character*80 cmname
      dimension hard(3), statev(nstatv), time(*), 
     1          predef(numfieldv), dpred(*)
c 
      dimension table(2,nvalue)
c
      parameter(zero = 0.d0)
c
c            set yield stress to last value of table, 
c            hardening to zero
c
      syield  = table(1,nvalue)
      hard(1) = zero
c
c            if more than one entry, search table
c
      if( nvalue .gt. 1 ) then
         do k1 = 1, nvalue-1
           eqpl1 = table(2,k1+1)
           if( eqplas .lt. eqpl1 ) then
              eqpl0 = table(2,k1)
              if( eqpl1 .le .eqpl0 ) then
                  write(kout,100)
                  call xit
              endif
c
c            current yield stress and hardening
c
              deqpl = eqpl1 - eqpl0
              syiel0 = table(1,k1)
              syiel1 = table(1,k1+1)
              dsyiel = syiel1 - syiel0
              hard(1) = dsyiel / deqpl
              syield = syiel0 + ( eqplas-eqpl0 ) * hard(1)
              goto 10
            endif 
           end do
10         continue
      endif
      return 
c
 100  format(//,  30x,  '***error - plastic strain must be', 
     1         ' entered in ascending order' )
c
      end
c     ****************************************************************
c     *                                                              *
c     *       subroutine umat - neo-hookean hyperelasticity          *
c     *                                                              *
c     ****************************************************************
c
c      
      subroutine umat_4( stress, statev, ddsdde, sse, spd, scd, rpl,
     1 ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,
     2 dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, 
     3 props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, 
     4 noel, npt, layer, kspt, kstep, kinc, kiter, kout )
c
      implicit real*8(a-h,o-z)
      parameter (nprecd=2)
c
      character*8 cmname
      dimension stress(*), statev(nstatv), ddsdde(ntens,ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 predef(1), dpred(1), props(nprops), coords(3), drot(3,3),
     3 dfgrd0(3,3), dfgrd1(3,3)

c              local arrays:
c
c              eelas  - logarithmic elastic strains
c              eelasp - principal elastic strains
c              bbar   - deviatoric right cauchy-green tensor
c              bbarp  - principal values of bbar
c              bbarn  - principal direction of bbar (and eelas)
c              distgr - deviatoric deformation gradient (distortion tensor)
c
      dimension eelas(6), eelasp(3), bbar(6), bbarp(3), bbarn(3,3),
     &          distgr(3,3)
c
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     &          six=6.d0)
c
c ----------------------------------------------------------------
c    umat for 3-D compressible neo-hookean hyperelasticity
c
c      props(1) - e
c      props(2) - nu
c ----------------------------------------------------------------
c
c              WARP3D handles thermal strains for this umat.
c              No adjustment needed for dstran. Just return.
c
      if( kiter .eq. 0 ) return
c
c 
c            elastic properties
c
      emod = props(1)
      enu  = props(2)
      c10  = emod/(four*(one+enu))
      d1   = six*(one-two*enu)/emod
c 
c            jacobian and distortion tensor
c
c      if( noel .eq. 5 .and. npt .eq. 3 ) then
c        write(kout,*) '... [F]n:'
c        write(kout,9000) (dfgrd0(i,1:3),i=1,3)
c        write(kout,*) '... [F]n+1:'
c        write(kout,9000) (dfgrd1(i,1:3),i=1,3)
c      end if
c 9000 format(3x,3e14.6)

      det = dfgrd1(1,1)*dfgrd1(2,2)*dfgrd1(3,3)
     &     -dfgrd1(1,2)*dfgrd1(2,1)*dfgrd1(3,3)
      det = det+dfgrd1(1,2)*dfgrd1(2,3)*dfgrd1(3,1)
     &         +dfgrd1(1,3)*dfgrd1(3,2)*dfgrd1(2,1)
     &         -dfgrd1(1,3)*dfgrd1(3,1)*dfgrd1(2,2)
     &         -dfgrd1(2,3)*dfgrd1(3,2)*dfgrd1(1,1)
      scale = det**(-one/three)
      distgr(1:3,1:3) = scale * dfgrd1(1:3,1:3)
c
c            calculate deviatoric left cauchy-green deformation tensor
c
      bbar(1) = distgr(1,1)**2+distgr(1,2)**2+distgr(1,3)**2
      bbar(2) = distgr(2,1)**2+distgr(2,2)**2+distgr(2,3)**2
      bbar(3) = distgr(3,3)**2+distgr(3,1)**2+distgr(3,2)**2
      bbar(4) = distgr(1,1)*distgr(2,1)+distgr(1,2)*distgr(2,2)
     &             + distgr(1,3)*distgr(2,3)
      bbar(5) = distgr(1,1)*distgr(3,1)+distgr(1,2)*distgr(3,2)
     &             + distgr(1,3)*distgr(3,3)
      bbar(6) = distgr(2,1)*distgr(3,1)+distgr(2,2)*distgr(3,2)
     &             + distgr(2,3)*distgr(3,3)
c
c                calculate the stress
c
      trbbar = (bbar(1)+bbar(2)+bbar(3))/three
      eg     = two*c10/det
      ek     = two/d1*(two*det-one)
      pr     = two/d1*(det-one)
      stress(1) = eg*(bbar(1)-trbbar)+pr
      stress(2) = eg*(bbar(2)-trbbar)+pr
      stress(3) = eg*(bbar(3)-trbbar)+pr
      stress(4) = eg*bbar(4)
      stress(5) = eg*bbar(5)
      stress(6) = eg*bbar(6)
c
c               calculate upper triangle of symmetric tangent stiffness.
c               fill in lower-triangle.
c
      eg23 = eg*two/three
      ddsdde(1,1) =  eg23*(bbar(1)+trbbar)+ek
      ddsdde(2,2) =  eg23*(bbar(2)+trbbar)+ek
      ddsdde(3,3) =  eg23*(bbar(3)+trbbar)+ek
      ddsdde(1,2) = -eg23*(bbar(1)+bbar(2)-trbbar)+ek
      ddsdde(1,3) = -eg23*(bbar(1)+bbar(3)-trbbar)+ek
      ddsdde(2,3) = -eg23*(bbar(2)+bbar(3)-trbbar)+ek
      ddsdde(1,4) =  eg23*bbar(4)/two
      ddsdde(2,4) =  eg23*bbar(4)/two
      ddsdde(3,4) = -eg23*bbar(4)
      ddsdde(4,4) =  eg*(bbar(1)+bbar(2))/two
      ddsdde(1,5) =  eg23*bbar(5)/two
      ddsdde(2,5) = -eg23*bbar(5)
      ddsdde(3,5) =  eg23*bbar(5)/two
      ddsdde(1,6) = -eg23*bbar(6)
      ddsdde(2,6) =  eg23*bbar(6)/two
      ddsdde(3,6) =  eg23*bbar(6)/two
      ddsdde(5,5) =  eg*(bbar(1)+bbar(3))/two
      ddsdde(6,6) =  eg*(bbar(2)+bbar(3))/two
      ddsdde(4,5) =  eg*bbar(6)/two
      ddsdde(4,6) =  eg*bbar(5)/two
      ddsdde(5,6) =  eg*bbar(4)/two
c
      do k1 = 1, 6
        do k2 = 1, k1-1
          ddsdde(k1,k2) = ddsdde(k2,k1)
        end do 
      end do
c
c            calculate logarithmic elastic strains (optional)
c
      call sprind(bbar, bbarp, bbarn, 1, ndi, nshr)
      eelasp(1) = log(sqrt(bbarp(1))/scale)
      eelasp(2) = log(sqrt(bbarp(2))/scale)
      eelasp(3) = log(sqrt(bbarp(3))/scale)
      eelas(1)  = eelasp(1)*bbarn(1,1)**2+eelasp(2)*bbarn(2, 1)**2
     &            + eelasp(3)*bbarn(3, 1)**2
      eelas(2)  = eelasp(1)*bbarn(1, 2)**2+eelasp(2)*bbarn(2, 2)**2
     &            + eelasp(3)*bbarn(3, 2)**2
      eelas(3)  = eelasp(1)*bbarn(1, 3)**2+eelasp(2)*bbarn(2, 3)**2
     &            + eelasp(3)*bbarn(3, 3)**2
      eelas(4)  = two*(eelasp(1)*bbarn(1, 1)*bbarn(1, 2)
     &            + eelasp(2)*bbarn(2, 1)*bbarn(2, 2)
     &            + eelasp(3)*bbarn(3, 1)*bbarn(3, 2))
      eelas(5) = two*(eelasp(1)*bbarn(1, 1)*bbarn(1, 3)
     &            + eelasp(2)*bbarn(2, 1)*bbarn(2, 3)
     &            + eelasp(3)*bbarn(3, 1)*bbarn(3, 3))
      eelas(6) = two*(eelasp(1)*bbarn(1, 2)*bbarn(1, 3)
     &            + eelasp(2)*bbarn(2, 2)*bbarn(2, 3)
     &            + eelasp(3)*bbarn(3, 2)*bbarn(3, 3))
c
c                   store elastic strains in state variable array
c
      statev(1:6) = eelas(1:6)
c
      return
      end



c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine xit  (called by UMATs)            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 3/22/12                     *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine xit
      implicit integer (a-z)
c
      write(*,*) '>>> UMAT called to abort execution'
      call die_abort
      stop
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine rotsig                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 03/29/12                       *
c     *                                                              *
c     *     tensor rotation routine for Abaqus compatible UMAT       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rotsig( in_vec, drot, out_vec, type, nrow, ncol )
      implicit none
c
c                      parameter declarations
c
      double precision in_vec(6), drot(3,3), out_vec(6)
      integer type, nrow, ncol
c
c                     locally defined arrays-variables
c
      double precision factor, one, half, two, a(3,3), 
     & t(3,3), c(3,3)
      logical local_debug
      data one, half, two, local_debug / 1.0, 0.5, 2.0, .true. /
      data one, half, two, local_debug / 1.0d00, 0.5d00, 
     &                                       2.0d00, .true. /
c
c
c                     out = drot * in * trans(drot)
c
c                     put input tensor into 3x3 matrix form.
c                     type = 1 input tensor has engr strain
c                     type = 2 input tensor has stress
c                     put into 3 x 3 form
c                     Abaqus umat ordering: x,y,z,xy,xz,yz
c
      factor = one
      if( type .eq. 1 ) factor = half
c
      a(1,1) = in_vec(1)
      a(2,1) = in_vec(4) * factor
      a(3,1) = in_vec(5) * factor
      a(1,2) = a(2,1)
      a(2,2) = in_vec(2)
      a(3,2) = in_vec(6) * factor
      a(1,3) = a(3,1)
      a(2,3) = a(3,2)
      a(3,3) = in_vec(3)
c
c                     t = a * trans(drot)
c
      t(1,1) = a(1,1)*drot(1,1) + a(1,2)*drot(1,2)  + a(1,3)*drot(1,3)
      t(2,1) = a(2,1)*drot(1,1) + a(2,2)*drot(1,2)  + a(2,3)*drot(1,3)
      t(3,1) = a(3,1)*drot(1,1) + a(3,2)*drot(1,2)  + a(3,3)*drot(1,3)
c
      t(1,2) = a(1,1)*drot(2,1) + a(1,2)*drot(2,2)  + a(1,3)*drot(2,3)
      t(2,2) = a(2,1)*drot(2,1) + a(2,2)*drot(2,2)  + a(2,3)*drot(2,3)
      t(3,2) = a(3,1)*drot(2,1) + a(3,2)*drot(2,2)  + a(3,3)*drot(2,3)
c
      t(1,3) = a(1,1)*drot(3,1) + a(1,2)*drot(3,2)  + a(1,3)*drot(3,3)
      t(2,3) = a(2,1)*drot(3,1) + a(2,2)*drot(3,2)  + a(2,3)*drot(3,3)
      t(3,3) = a(3,1)*drot(3,1) + a(3,2)*drot(3,2)  + a(3,3)*drot(3,3)
c
c                     c = drot * t
c
      c(1,1) = drot(1,1)*t(1,1) + drot(1,2)*t(2,1) + drot(1,3)*t(3,1)
      c(2,1) = drot(2,1)*t(1,1) + drot(2,2)*t(2,1) + drot(2,3)*t(3,1)
      c(3,1) = drot(3,1)*t(1,1) + drot(3,2)*t(2,1) + drot(3,3)*t(3,1)
c
      c(1,2) = drot(1,1)*t(1,2) + drot(1,2)*t(2,2) + drot(1,3)*t(3,2)
      c(2,2) = drot(2,1)*t(1,2) + drot(2,2)*t(2,2) + drot(2,3)*t(3,2)
      c(3,2) = drot(3,1)*t(1,2) + drot(3,2)*t(2,2) + drot(3,3)*t(3,2)
c
      c(1,3) = drot(1,1)*t(1,3) + drot(1,2)*t(2,3) + drot(1,3)*t(3,3)
      c(2,3) = drot(2,1)*t(1,3) + drot(2,2)*t(2,3) + drot(2,3)*t(3,3)
      c(3,3) = drot(3,1)*t(1,3) + drot(3,2)*t(2,3) + drot(3,3)*t(3,3)
c
      factor = one
      if( type .eq. 1 ) factor = two
      out_vec(1) = c(1,1)
      out_vec(2) = c(2,2)
      out_vec(3) = c(3,3)
      out_vec(4) = c(2,1) * two
      out_vec(5) = c(3,1) * two
      out_vec(6) = c(3,2) * two
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine sinv                         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 05/14/12                       *
c     *                                                              *
c     *              invariants 1 & 2 of 3D stress tensor            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine sinv( sig, sinv1, sinv2, ndi, nshr )
      implicit none
c
c                      parameter declarations
c
      double precision sig(*), sinv1, sinv2
      integer ndi, nshr
c
c                     locally defined arrays-variables
c
      double precision sig_dev(6), t1, t2, one_third, two, oneptfive
      data one_third, two, oneptfive
c     &    /  0.3333333, 2.0, 1.5 /
     &    /  0.3333333333333333d00, 2.0d00, 1.5d00 /
c
      sinv1 = one_third * ( sig(1) + sig(2) + sig(3) )
c
      sig_dev(1) = sig(1) - sinv1
      sig_dev(2) = sig(2) - sinv1
      sig_dev(3) = sig(3) - sinv1
      sig_dev(4) = sig(4)
      sig_dev(5) = sig(5)
      sig_dev(6) = sig(6)
c
      t1 = sig_dev(1)**2 + sig_dev(2)**2 + sig_dev(3)**2
      t2 = sig_dev(4)**2 + sig_dev(5)**2 + sig_dev(6)**2
      sinv2 = sqrt( oneptfive * (t1 + two * t2 ) )
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine sprinc                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/14/12                    *
c     *                                                              *
c     *     compute principal strain or stresses compatible with     *
c     *     UMAT specifications                                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine sprinc( s, ps, lstr, ndi, nshr )
      implicit integer (a-z)
      double precision s(*), ps(*)
c
c                    locally allocated
c
      double precision temp(6), wk(3), evec(3,3), half, one, factor
c
      data  half, one 
c     &   / 0.5, 1.0 / 
     &   / 0.5d00, 1.0d00 /
c
c        calculate the principal strains or stresses for UMAT support.
c        Abaqus stress/strain ordering: xx, yy, zz, xy, xz, yz
c        lstr = 1 (stresses), lstr = 2 (strains)
c
      factor = one
      if( lstr .eq. 2 ) factor = half
c
      temp(1) = s(1)
      temp(2) = s(4) * factor
      temp(3) = s(2)
      temp(4) = s(5) * factor
      temp(5) = s(6) * factor
      temp(6) = s(3)
      call ou3dpr( temp, 3, 0, ps, evec, 3, wk, ier )   ! warp3d routine
c
      return
      end  
c     ****************************************************************
c     *                                                              *
c     *                      subroutine sprind                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/14/12                    *
c     *                                                              *
c     *     compute principal strain or stresses compatible with     *
c     *     UMAT specifications                                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine sprind( s, ps, an, lstr, ndi, nshr )
      implicit integer (a-z)
      double precision s(*), ps(*), an(3,3)
c
c                    locally allocated
c
      double precision temp(6), wk(3), half, one, factor
c
      data  half, one 
c     &   / 0.5, 1.0 / 
     &   / 0.5d00, 1.0d00 /
c
c        calculate the principal strains or stresses and eigenvectors
c        for UMAT support.
c        Abaqus stress/strain ordering: xx, yy, zz, xy, xz, yz
c        lstr = 1 (stresses), lstr = 2 (strains)
c
      factor = one
      if( lstr .eq. 2 ) factor = half
c
      temp(1) = s(1)
      temp(2) = s(4) * factor
      temp(3) = s(2)
      temp(4) = s(5) * factor
      temp(5) = s(6) * factor
      temp(6) = s(3)
      call ou3dpr( temp, 3, 1, ps, an, 3, wk, ier )   ! warp3d routine
c
      return
      end  



      

c
c *******************************************************************
c *                                                                 *
c *        optional UMAT output rotuine                             *  
c *                                                                 *
c *   set upt to 3 material model dependent output values           *
c *                                                                 *
c *******************************************************************
c
c
      subroutine ou_umat( gpn, mxvl, span, iout, elestr,
     &                      stress, history )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer gpn, mxvl, span, iout
c
c     
      double precision stress(mxvl,*), elestr(mxvl,*), history(mxvl,*)
c
c               description of parameters
c               -------------------------
c
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     span              : number of elements in current block
c     iout              : write messages to this device number
c     stress            : current stresses for all
c                         elements in block for this gauss point
c (*) elestr            : stresses to be output for elements
c     history           : current history values for all
c                         elements in block for this gauss point
c     
c    (*)  values to be updated by this material model
c
c
c   stress ordering                elestr ordering
c     (1) sig-xx                   (1) sig-xx
c     (2) sig-yy                   (2) sig-yy
c     (3) sig-zz                   (3) sig-zz
c     (4) tau-xy                   (4) tau-xy
c     (5) tau-yz                   (5) tau-yz        
c     (6) tau-xz                   (6) tau-xz
c     (7) total work density       (7) total work density   
c                                  (8) mises equiv. stress
c                             (*)  (9) mat_val1
c                             (*) (10) mat_val2
c                             (*) (11) mat_val3
c
c  NOTE:  do NOT modify "stress" array or columns 1-8 of the "elestr"
c         array. only modify columns 9-11 of "elestr". These are the
c         3 "material model" dependent values output in the stress
c         values for a gauss point. See Section 2.12 of manual and
c         description of each material model. The output labels for 
c         columns 9-11 are "c1", "c2", "c3".
c
c
c
c
      integer i 
c
c           c1, current scalar plastic strain
c           c2, current yielding state. = 0 (never yielded),
c               = 1 actively yielding, = 2 (previously yielded
c                 but currently linear-elastic)
       do i = 1, span
         elestr(i,9)   = history(i,19) ! c1
         elestr(i,10)  = history(i,20) ! c2
       end do

c     
c
       return
       end
c c
c c
c c     Exmaple test file for about umats.
c c
c c
c c     square cross-section, straight bar along X-axis.
c c     fixed at both ends.
c c
c c     impose bending displacement increments at mid-span
c c     and temperature change over cross-section at mid-span.
c c
c c     umat_2 in the umats.f file. mises plasticity with
c c     kinematic, bilinear hardening.umat takes
c c     care of thermal strains. um_5 is CTE and is constant.
c c     WARP3D makes alpha (CTE) = 0.0 and this cannot adjust for
c c     thermal strain increment over a load step.
c c
c c     umat_3 in the umats.f file. mises plasticity with
c c     isotropic hardening. segmental hardening curve. 
c c     um_3 is total number of points on curve. (um_4, um_5)
c c     are first point. um_4 is yield stress, um_5 must be 0.0
c c     (plastic strain). (um_6, um_7) are next point on curve.
c c     (um_8, um_9) next point. etc.
c c     values for curve here define same plastic hardening slope
c c     H' as for model umat_2 just above.
c c     umat lets WARP3D handle thermal strains thus alpha is specified
c c     not as a umat property.
c c
c structure beam
c c
c material test_umat_2
c     properties umat  rho 0.0 ,
c        um_1 30000 um_2 0.3 um_3 60 um_4 1034.48276,
c        um_5 1.0e-04 
c c
c material test_umat_3
c     properties umat  rho 0.0  alpha 1.0e-04 ,
c        um_1 30000 um_2 0.3 um_3 2,
c       um_4 60 um_5 0.0,
c       um_6 1094.48276  um_7 1.0
c c
c number of nodes 99 elements 40
c elements
c   1-40  type  l3disop linear  material test_umat_2,
c               order 2x2x2  short  bbar
c c 1-40  type  l3disop linear  material test_umat_3, 
c c           order 2x2x2  short  bbar
c c 
c coordinates
c *echo off
c       1  0.000000000E+00  0.000000000E+00  0.000000000E+00
c       2  0.100000000E+01  0.000000000E+00  0.000000000E+00
c       3  0.200000000E+01  0.000000000E+00  0.000000000E+00
c       4  0.300000000E+01  0.000000000E+00  0.000000000E+00
c       5  0.400000000E+01  0.000000000E+00  0.000000000E+00
c       6  0.500000000E+01  0.000000000E+00  0.000000000E+00
c       7  0.600000000E+01  0.000000000E+00  0.000000000E+00
c       8  0.700000048E+01  0.000000000E+00  0.000000000E+00
c       9  0.800000095E+01  0.000000000E+00  0.000000000E+00
c      10  0.900000095E+01  0.000000000E+00  0.000000000E+00
c      11  0.100000000E+02  0.000000000E+00  0.000000000E+00
c      12  0.000000000E+00  0.500000000E+00  0.000000000E+00
c      13  0.100000000E+01  0.500000000E+00  0.000000000E+00
c      14  0.200000000E+01  0.500000000E+00  0.000000000E+00
c      15  0.300000000E+01  0.500000000E+00  0.000000000E+00
c      16  0.400000000E+01  0.500000000E+00  0.000000000E+00
c      17  0.500000000E+01  0.500000000E+00  0.000000000E+00
c      18  0.600000000E+01  0.500000000E+00  0.000000000E+00
c      19  0.700000048E+01  0.500000000E+00  0.000000000E+00
c      20  0.800000095E+01  0.500000000E+00  0.000000000E+00
c      21  0.900000095E+01  0.500000000E+00  0.000000000E+00
c      22  0.100000000E+02  0.500000000E+00  0.000000000E+00
c      23  0.000000000E+00  0.100000000E+01  0.000000000E+00
c      24  0.100000000E+01  0.100000000E+01  0.000000000E+00
c      25  0.200000000E+01  0.100000000E+01  0.000000000E+00
c      26  0.300000000E+01  0.100000000E+01  0.000000000E+00
c      27  0.400000000E+01  0.100000000E+01  0.000000000E+00
c      28  0.500000000E+01  0.100000000E+01  0.000000000E+00
c      29  0.600000000E+01  0.100000000E+01  0.000000000E+00
c      30  0.700000048E+01  0.100000000E+01  0.000000000E+00
c      31  0.800000095E+01  0.100000000E+01  0.000000000E+00
c      32  0.900000095E+01  0.100000000E+01  0.000000000E+00
c      33  0.100000000E+02  0.100000000E+01  0.000000000E+00
c      34  0.000000000E+00  0.000000000E+00  0.500000000E+00
c      35  0.100000000E+01  0.000000000E+00  0.500000000E+00
c      36  0.200000000E+01  0.000000000E+00  0.500000000E+00
c      37  0.300000000E+01  0.000000000E+00  0.500000000E+00
c      38  0.400000000E+01  0.000000000E+00  0.500000000E+00
c      39  0.500000000E+01  0.000000000E+00  0.500000000E+00
c      40  0.600000000E+01  0.000000000E+00  0.500000000E+00
c      41  0.700000048E+01  0.000000000E+00  0.500000000E+00
c      42  0.800000095E+01  0.000000000E+00  0.500000000E+00
c      43  0.900000095E+01  0.000000000E+00  0.500000000E+00
c      44  0.100000000E+02  0.000000000E+00  0.500000000E+00
c      45  0.000000000E+00  0.500000000E+00  0.500000000E+00
c      46  0.100000000E+01  0.500000000E+00  0.500000000E+00
c      47  0.200000000E+01  0.500000000E+00  0.500000000E+00
c      48  0.300000000E+01  0.500000000E+00  0.500000000E+00
c      49  0.400000000E+01  0.500000000E+00  0.500000000E+00
c      50  0.500000000E+01  0.500000000E+00  0.500000000E+00
c      51  0.600000000E+01  0.500000000E+00  0.500000000E+00
c      52  0.700000048E+01  0.500000000E+00  0.500000000E+00
c      53  0.800000095E+01  0.500000000E+00  0.500000000E+00
c      54  0.900000095E+01  0.500000000E+00  0.500000000E+00
c      55  0.100000000E+02  0.500000000E+00  0.500000000E+00
c      56  0.000000000E+00  0.100000000E+01  0.500000000E+00
c      57  0.100000000E+01  0.100000000E+01  0.500000000E+00
c      58  0.200000000E+01  0.100000000E+01  0.500000000E+00
c      59  0.300000000E+01  0.100000000E+01  0.500000000E+00
c      60  0.400000000E+01  0.100000000E+01  0.500000000E+00
c      61  0.500000000E+01  0.100000000E+01  0.500000000E+00
c      62  0.600000000E+01  0.100000000E+01  0.500000000E+00
c      63  0.700000048E+01  0.100000000E+01  0.500000000E+00
c      64  0.800000095E+01  0.100000000E+01  0.500000000E+00
c      65  0.900000095E+01  0.100000000E+01  0.500000000E+00
c      66  0.100000000E+02  0.100000000E+01  0.500000000E+00
c      67  0.000000000E+00  0.000000000E+00  0.100000000E+01
c      68  0.100000000E+01  0.000000000E+00  0.100000000E+01
c      69  0.200000000E+01  0.000000000E+00  0.100000000E+01
c      70  0.300000000E+01  0.000000000E+00  0.100000000E+01
c      71  0.400000000E+01  0.000000000E+00  0.100000000E+01
c      72  0.500000000E+01  0.000000000E+00  0.100000000E+01
c      73  0.600000000E+01  0.000000000E+00  0.100000000E+01
c      74  0.700000048E+01  0.000000000E+00  0.100000000E+01
c      75  0.800000095E+01  0.000000000E+00  0.100000000E+01
c      76  0.900000095E+01  0.000000000E+00  0.100000000E+01
c      77  0.100000000E+02  0.000000000E+00  0.100000000E+01
c      78  0.000000000E+00  0.500000000E+00  0.100000000E+01
c      79  0.100000000E+01  0.500000000E+00  0.100000000E+01
c      80  0.200000000E+01  0.500000000E+00  0.100000000E+01
c      81  0.300000000E+01  0.500000000E+00  0.100000000E+01
c      82  0.400000000E+01  0.500000000E+00  0.100000000E+01
c      83  0.500000000E+01  0.500000000E+00  0.100000000E+01
c      84  0.600000000E+01  0.500000000E+00  0.100000000E+01
c      85  0.700000048E+01  0.500000000E+00  0.100000000E+01
c      86  0.800000095E+01  0.500000000E+00  0.100000000E+01
c      87  0.900000095E+01  0.500000000E+00  0.100000000E+01
c      88  0.100000000E+02  0.500000000E+00  0.100000000E+01
c      89  0.000000000E+00  0.100000000E+01  0.100000000E+01
c      90  0.100000000E+01  0.100000000E+01  0.100000000E+01
c      91  0.200000000E+01  0.100000000E+01  0.100000000E+01
c      92  0.300000000E+01  0.100000000E+01  0.100000000E+01
c      93  0.400000000E+01  0.100000000E+01  0.100000000E+01
c      94  0.500000000E+01  0.100000000E+01  0.100000000E+01
c      95  0.600000000E+01  0.100000000E+01  0.100000000E+01
c      96  0.700000048E+01  0.100000000E+01  0.100000000E+01
c      97  0.800000095E+01  0.100000000E+01  0.100000000E+01
c      98  0.900000095E+01  0.100000000E+01  0.100000000E+01
c      99  0.100000000E+02  0.100000000E+01  0.100000000E+01
c incidences 
c        1       1       2      13      12      34      35      46      45
c        2       2       3      14      13      35      36      47      46
c        3       3       4      15      14      36      37      48      47
c        4       4       5      16      15      37      38      49      48
c        5       5       6      17      16      38      39      50      49
c        6       6       7      18      17      39      40      51      50
c        7       7       8      19      18      40      41      52      51
c        8       8       9      20      19      41      42      53      52
c        9       9      10      21      20      42      43      54      53
c       10      10      11      22      21      43      44      55      54
c       11      12      13      24      23      45      46      57      56
c       12      13      14      25      24      46      47      58      57
c       13      14      15      26      25      47      48      59      58
c       14      15      16      27      26      48      49      60      59
c       15      16      17      28      27      49      50      61      60
c       16      17      18      29      28      50      51      62      61
c       17      18      19      30      29      51      52      63      62
c       18      19      20      31      30      52      53      64      63
c       19      20      21      32      31      53      54      65      64
c       20      21      22      33      32      54      55      66      65
c       21      34      35      46      45      67      68      79      78
c       22      35      36      47      46      68      69      80      79
c       23      36      37      48      47      69      70      81      80
c       24      37      38      49      48      70      71      82      81
c       25      38      39      50      49      71      72      83      82
c       26      39      40      51      50      72      73      84      83
c       27      40      41      52      51      73      74      85      84
c       28      41      42      53      52      74      75      86      85
c       29      42      43      54      53      75      76      87      86
c       30      43      44      55      54      76      77      88      87
c       31      45      46      57      56      78      79      90      89
c       32      46      47      58      57      79      80      91      90
c       33      47      48      59      58      80      81      92      91
c       34      48      49      60      59      81      82      93      92
c       35      49      50      61      60      82      83      94      93
c       36      50      51      62      61      83      84      95      94
c       37      51      52      63      62      84      85      96      95
c       38      52      53      64      63      85      86      97      96
c       39      53      54      65      64      86      87      98      97
c       40      54      55      66      65      87      88      99      98
c blocking $ scalar: only sparse solver allowed
c c          small block sizes used to test parallel over blocks
c        1       4       1
c        2       4       5
c        3       4       9
c        4       4      13
c        5       4      17
c        6       4      21
c        7       4      25
c        8       4      29
c        9       4      33
c       10       4      37 
c c
c  initial conditions
c   temperature
c     nodes 1-99 temperature 100.0
c c
c  list 'nodes-center'  x=5
c c
c loading unit_temp
c   nodal loads
c     'nodes-center' temperature 1.0
c c
c constraints
c  plane x=0 fixed
c  plane x=10 fixed
c  'nodes-center' v -1.0
c c 
c *echo on
c c
c c      for displacement loading, first yield between constraints x
c c      0.02 and 0.03
c c      temperature loading first yields between 40-50 delta T
c c
c  loading test
c   nonlinear
c     step 1 unit_temp 50 constraints 0.03
c     step 2 unit_temp 50 constraints 0.03
c     step 3 unit_temp 50 constraints 0.03
c c
c  nonlinear analysis parameters
c    solution technique direct sparse
c    maximum iterations 4
c    minimum iterations 1
c    convergence test norm res tol 0.001 
c    time step  1.0e10
c    trace solution on lpcg_solution off 
c    extrapolate off
c    adaptive on
c c
c list 'left-end'  x=0
c list 'right-end' x=10
c    compute displacements for loading test step 3
c    output displacements nodes 'nodes-center'
c    output reactions 'left-end'
c    output reactions 'right-end'
c stop
