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
