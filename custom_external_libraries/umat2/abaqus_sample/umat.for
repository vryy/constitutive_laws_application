      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     * RPL,DDSDDT,DRPLDE,DRPLDT,
     * STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     * NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     * CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     * DDSDDE(NTENS,NTENS),
     * DDSDDT(NTENS),DRPLDE(NTENS),
     * STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     * PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DIMENSION DSTRES(6),D(3,3)
!
! EVALUATE NEW STRESS TENSOR
!
      EV = 0.
      DEV = 0.
      DO K1=1,NDI
        EV = EV + STRAN(K1)
        DEV = DEV + DSTRAN(K1)
      END DO
!
      TERM1 = .5*DTIME + PROPS(5)
      TERM1I = 1./TERM1
      TERM2 = (.5*DTIME*PROPS(1)+PROPS(3))*TERM1I*DEV
      TERM3 = (DTIME*PROPS(2)+2.*PROPS(4))*TERM1I
!
      DO K1=1,NDI
        DSTRES(K1) = TERM2+TERM3*DSTRAN(K1)
     *   +DTIME*TERM1I*(PROPS(1)*EV
     *   +2.*PROPS(2)*STRAN(K1)-STRESS(K1))
        STRESS(K1) = STRESS(K1) + DSTRES(K1)
      END DO
!
      TERM2 = (.5*DTIME*PROPS(2) + PROPS(4))*TERM1I
      I1 = NDI
      DO K1=1,NSHR
        I1 = I1+1
        DSTRES(I1) = TERM2*DSTRAN(I1)+
     *   DTIME*TERM1I*(PROPS(2)*STRAN(I1)-STRESS(I1))
        STRESS(I1) = STRESS(I1)+DSTRES(I1)
      END DO
!
! CREATE NEW JACOBIAN
!
      TERM2 = (DTIME*(.5*PROPS(1)+PROPS(2))+PROPS(3)+
     *   2.*PROPS(4))*TERM1I
      TERM3 = (.5*DTIME*PROPS(1)+PROPS(3))*TERM1I
      DO K1=1,NTENS
        DO K2=1,NTENS
          DDSDDE(K2,K1) = 0.
        END DO
      END DO
!
      DO K1=1,NDI
        DDSDDE(K1,K1) = TERM2
      END DO
!
      DO K1=2,NDI
        N2 = K1-1
        DO K2=1,N2
          DDSDDE(K2,K1) = TERM3
          DDSDDE(K1,K2) = TERM3
        END DO
      END DO
      TERM2 = (.5*DTIME*PROPS(2)+PROPS(4))*TERM1I
      I1 = NDI
      DO K1=1,NSHR
        I1 = I1+1
        DDSDDE(I1,I1) = TERM2
      END DO
!
! TOTAL CHANGE IN SPECIFIC ENERGY
!
      TDE = 0.
      DO K1=1,NTENS
        TDE = TDE + (STRESS(K1)-.5*DSTRES(K1))*DSTRAN(K1)
      END DO
!
! CHANGE IN SPECIFIC ELASTIC STRAIN ENERGY
!
      TERM1 = PROPS(1) + 2.*PROPS(2)
      DO K1=1,NDI
        D(K1,K1) = TERM1
      END DO
      DO K1=2,NDI
        N2 = K1-1
        DO K2=1,N2
          D(K1,K2) = PROPS(1)
          D(K2,K1) = PROPS(1)
        END DO
      END DO
      DEE = 0.
      DO K1=1,NDI
        TERM1 = 0.
        TERM2 = 0.
        DO K2=1,NDI
          TERM1 = TERM1 + D(K1,K2)*STRAN(K2)
          TERM2 = TERM2 + D(K1,K2)*DSTRAN(K2)
        END DO
        DEE = DEE + (TERM1+.5*TERM2)*DSTRAN(K1)
      END DO
      I1 = NDI
      DO K1=1,NSHR
        I1 = I1+1
        DEE = DEE + PROPS(2)*(STRAN(I1)+.5*DSTRAN(I1))*DSTRAN(I1)
      END DO
      SSE = SSE + DEE
      SCD = SCD + TDE - DEE
      RETURN
      END
