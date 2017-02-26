      Subroutine UdsmMC ( IDTask, iMod, IsUndr,
     *                    iStep, iTer, iEl, Int,
     *                    X, Y, Z,
     *                    Time0, dTime,
     *                    Props, Sig0, Swp0, StVar0,
     *                    dEps, D, BulkW,
     *                    Sig, Swp, StVar, ipl,
     *                    nStat,
     *                    NonSym, iStrsDep, iTimeDep, iTang,
     *                    iAbort )
!
! Purpose: User supplied soil model for Mohr-Coulomb
! Based on HYPLAS of Souza de Neto
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : return number of state variables
!                       5 : inquire matrix properties
!                           return switch for non-symmetric D-matrix
!                           stress/time dependent matrix
!                       6 : calculate elastic material stiffness matrix
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iMod     I   I    : model number (1..10)
!  IsUndr   I   I    : =1 for undrained, 0 otherwise
!  iStep    I   I    : Global step number
!  iter     I   I    : Global iteration number
!  iel      I   I    : Global element number
!  Int      I   I    : Global integration point number
!  X        I   R    : X-Position of integration point
!  Y        I   R    : Y-Position of integration point
!  Z        I   R    : Z-Position of integration point
!  Time0    I   R    : Time at start of step
!  dTime    I   R    : Time increment
!  Props    I   R()  : List with model parameters
!  Sig0     I   R()  : Stresses at start of step
!  Swp0     I   R    : Excess pore pressure start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulkmodulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : Non-Symmetric D-matrix ?
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iTang    O   I    : =1 for tangent D-matrix
!  iAbort   O   I    : =1 to force stopping of calculation
!
      Implicit Double Precision (A-H, O-Z)
!
      Dimension Props(*), Sig0(*), StVar0(*), dEps(*), D(6,6),
     *          Sig(*),   StVar(*),
     *          Strial(6), EEtD(6), UNIDEV(6), SOID(6), FOID(6,6),
     *          EIGVEC(3,3), EIGPRJ(3,6), Stra(6), PSTRS(3),
     *          DPSTRS(3,3)

      DATA TOL   /1.D-08/
      DATA MAXRT / 50 /
      DATA FOID / 1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *            0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     *            0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,
     *            0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,
     *            0.0d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,
     *            0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.5d0 /
      DATA SOID / 1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0 /

      If (IDTask .Eq. 1) Then
!       {Initialize state variables StVar0 }
        StVar(1)=0d0
        StVar(2)=0d0
        StVar(3)=0d0
        StVar(4)=0d0
        StVar(5)=0d0
        StVar(6)=0d0
        StVar(7)=0d0
        StVar(8)=0d0
        StVar(9)=0d0
      End If  ! IDTask = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 2) Then
!       { Calculate constitutive stresses Sig (and Swp) }
!       Initialize some algorithmic and internal variables
        DGAMA=R0
        DGAMB=R0
        IFPLAS=.FALSE.
        SUFAIL=.FALSE.
        EDGE=.FALSE.
        APEX=.FALSE.
        TWORET=.FALSE.
        APEXRET=.FALSE.
        EPBARN=StVar(7)
        EPBAR=EPBARN
!       Set some material properties
        YOUNG=PROPS(1)
        POISS=PROPS(2)
        SINPHI=PROPS(3)
        COSPHI=PROPS(4)
        SINPSI=PROPS(5)
        C=PROPS(6)
        HARD=PROPS(7)
!       Set some constants
        GMODU=YOUNG/(R2*(R1+POISS))
        BULK=YOUNG/(R3*(R1-R2*POISS))
        R2G=R2*GMODU
        R4G=R4*GMODU
        R2BULK=R2*BULK
        R2CPHI=R2*COSPHI
        R1D3=R1/R3
!       Compute elastic trial state
!       ---------------------------
!       Elastic trial volumetric strain and pressure stress
        EEtV=dEps(1)+dEps(2)+dEps(3)+dEps(7)+dEps(8)+dEps(9)
        PT=BULK*EEtV
!       Spectral decomposition of the elastic trial stress
        EEtvD3=EEtV*R1D3
        Strial(1)=R2G*(dEps(1)+dEps(7)-EEtvD3)+PT
        Strial(2)=R2G*(dEps(2)+dEps(8)-EEtvD3)+PT
        Strial(3)=R2G*(dEps(3)+dEps(9)-EEtvD3)+PT
!       Shear component
        Strial(4)=GMODU*(dEps(4)+dEps(10))
        Strial(5)=GMODU*(dEps(5)+dEps(11))
        Strial(6)=GMODU*(dEps(6)+dEps(12))
        iOpt = 1
!       Compute the principal stress and eigenprojection tensors
        Call PrnSig(iOpt,Strial,EIGVEC(3,:),EIGVEC(2,:),EIGVEC(1,:),
     *          PSTRS3,PSTRS2,PSTRS1,P,Q)
        DO I=1,3
          EIGPRJ(I,1)=EIGVEC(I,1)*EIGVEC(I,1)
          EIGPRJ(I,2)=EIGVEC(I,2)*EIGVEC(I,2)
          EIGPRJ(I,3)=EIGVEC(I,3)*EIGVEC(I,3)
          EIGPRJ(I,4)=EIGVEC(I,1)*EIGVEC(I,2)
          EIGPRJ(I,5)=EIGVEC(I,2)*EIGVEC(I,3)
          EIGPRJ(I,6)=EIGVEC(I,3)*EIGVEC(I,1)
        ENDDO
!       Compute trial yield function and check for plastic consistency
        COHE=C+HARD*EPBARN
        SMCT=PSTRS1-PSTRS3+(PSTRS1+PSTRS3)*SINPHI
        PHIA=SMCT-R2CPHI*COHE
        RES=PHIA
        IF(COHE.NE.R0)RES=RES/ABS(COHE)
        IF(RES.GT.TOL)THEN
!         Plastic step: Apply return mapping
!         identify possible edge return: either right or left of main plane
          SCAPRD=PSTRS1*(R1-SINPSI)+PSTRS2*(-R2)+PSTRS3*(R1+SINPSI)
          IF(SCAPRD.GE.R0)THEN
            RIGHT=.TRUE.
          ELSE
            RIGHT=.FALSE.
          ENDIF
!         Apply one-vector return mapping first (return to MAIN PLANE)
          SPHSPS=SINPHI*SINPSI
          CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS
          R4C2PH=R2CPHI*R2CPHI
!         Start Newton-Raphson iterations for DGAMA
          SUFAIL=.TRUE.
          DO NRITER=1,MXITER
!           Compute residual derivative
            DENOM=-CONSTA-R4C2PH*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
!           Compute Newton-Raphson increment and update variable DGAMA
            DDGAMA=-PHIA/DENOM
            DGAMA=DGAMA+DDGAMA
!           Compute new residual
            EPBAR=EPBARN+R2CPHI*DGAMA
            COHE=C+HARD*EPBAR
            PHIA=SMCT-CONSTA*DGAMA-R2CPHI*COHE
!           Check convergence
            RESNOR=ABS(PHIA)
            IF(SMCT.NE.R0)RESNOR=RESNOR/ABS(SMCT)
            IF(RESNOR.LE.TOL)THEN
              SUFAIL=.FALSE.
!             Check validity of 1-vector return (check sextant of converged stress)
              S1=PSTRS1-(R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI)*DGAMA
              S2=PSTRS2+(R4G*R1D3-R2BULK)*SINPSI*DGAMA
              S3=PSTRS3+(R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI)*DGAMA
              DELTA=DMAX1(ABS(S1),ABS(S2),ABS(S3))*SMALL
              IF(S1+DELTA.GE.S2.AND.S2+DELTA.GE.S3)THEN
!               converged stress is in the same sextant as trial stress -> 1-vector
!               return is valid.
                P=(S1+S2+S3)*R1D3
                TWORET=.FALSE.
                EXIT
              ELSE
!               converged stress is not in the same sextant -> 1-vector result is
!               not valid. Go to two-vector return map to edge
                TWORET=.TRUE.
                EXIT
              ENDIF
            ENDIF
          END DO
          IF(SUFAIL.EQ..TRUE.)THEN
!           TODO: throw some errors
          ELSE
            IF(TWORET.EQ..TRUE.)THEN
!             Apply two-vector return mapping to appropriate EDGE
              DGAMA=R0
              EPBAR=EPBARN
              COHE=C+HARD*EPBARN
              SMCTA=PSTRS1-PSTRS3+(PSTRS1+PSTRS3)*SINPHI
              IF(RIGHT.EQ..TRUE.)THEN
                SMCTB=PSTRS1-PSTRS2+(PSTRS1+PSTRS2)*SINPHI
              ELSE
                SMCTB=PSTRS2-PSTRS3+(PSTRS2+PSTRS3)*SINPHI
              ENDIF
              PHIA=SMCTA-R2CPHI*COHE
              PHIB=SMCTB-R2CPHI*COHE
              IF(RIGHT.EQ..TRUE.)THEN
                CONSTB=R2G*(R1+SINPHI+SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS
              ELSE
                CONSTB=R2G*(R1-SINPHI-SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS
              ENDIF
!             Start Newton-Raphson iterations for DGAMA and DGAMB
              SUFAIL=.TRUE.
              DO NRITER=1,MXITER
!               Compute residual derivative matrix
                FACTA=R4C2PH*HARD
                DRVAA=-CONSTA-FACTA
                DRVAB=-CONSTB-FACTA
                DRVBA=-CONSTB-FACTA
                DRVBB=-CONSTA-FACTA
!               Compute Newton-Raphson increment and update variables DGAMA and DGAMB
                R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA)
                DDGAMA=(-DRVBB*PHIA+DRVAB*PHIB)*R1DDET
                DDGAMB=(DRVBA*PHIA-DRVAA*PHIB)*R1DDET
                DGAMA=DGAMA+DDGAMA
                DGAMB=DGAMB+DDGAMB
!               Compute new residual
                EPBAR=EPBARN+R2CPHI*(DGAMA+DGAMB)
                COHE=C+HARD*EPBAR
                PHIA=SMCTA-CONSTA*DGAMA-CONSTB*DGAMB-R2CPHI*COHE
                PHIB=SMCTB-CONSTB*DGAMA-CONSTA*DGAMB-R2CPHI*COHE
!               Check convergence
                RESNOR=(ABS(PHIA)+ABS(PHIB))
                FACTOR=(ABS(SMCTA)+ABS(SMCTB))
                IF(FACTOR.NE.R0)RESNOR=RESNOR/FACTOR
                IF(RESNOR.LE.TOL)THEN
                  SUFAIL=.FALSE.
!                 Check validity of 2-vector return to edge
                  AUX1=R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI
                  AUX2=(R4G*R1D3-R2BULK)*SINPSI
                  AUX3=R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI
                  IF(RIGHT.EQ..TRUE.)THEN
                    S1=PSTRS1-AUX1*(DGAMA+DGAMB)
                    S2=PSTRS2+AUX2*DGAMA+AUX3*DGAMB
                    S3=PSTRS3+AUX3*DGAMA+AUX2*DGAMB
                  ELSE
                    S1=PSTRS1-AUX1*DGAMA+AUX2*DGAMB
                    S2=PSTRS2+AUX2*DGAMA-AUX1*DGAMB
                    S3=PSTRS3+AUX3*(DGAMA+DGAMB)
                  ENDIF
                  DELTA=DMAX1(ABS(S1),ABS(S2),ABS(S3))*SMALL
                  IF(S1+DELTA.GE.S2.AND.S2+DELTA.GE.S3)THEN
!                   converged stress is in the same sextant as trial stress -> 2-vector
!                   return to edge is valid.
                    EDGE=.TRUE.
                    APEXRET=.FALSE.
                    P=(S1+S2+S3)*R1D3
                    EXIT
                  ELSE
!                   converged stress is not in the same sextant -> 2-vector return to edge
!                   is not valid. Go to two-vector return map to APEX
                    APEXRET=.TRUE.
                    EXIT
                  ENDIF
                ENDIF
              END DO
              IF(SUFAIL.EQ..TRUE.)THEN
!               throw some errors
              ELSE
                IF(APEXRET.EQ..TRUE.)THEN
!                 Apply multi-vector return mapping to APEX
!                 Check conditions for which return to apex does not make sense
!                 Set initial guess for volumetric plastic strain increment DEPV
                  DEPV=R0
                  EPBAR=EPBARN
                  COHE=C+HARD*EPBAR
                  COTPHI=COSPHI/SINPHI
                  RES=COTPHI*COHE-PT
!                 Newton-Raphson iterations for DEPV
                  SUFAIL=.TRUE.
                  DO NRITER=1,MXITER
                    DENOM=COSPHI*COTPHI/SINPSI*HARD+BULK
                    DDEPV=-RES/DENOM
                    DEPV=DEPV+DDEPV
                    EPBAR=EPBARN+COSPHI/SINPSI*DEPV
                    COHE=C+HARD*EPBAR
                    P=PT-BULK*DEPV
                    RES=COTPHI*COHE-P
!                   check for convergence
                    RESNOR=ABS(RES)
                    IF(PT.NE.R0)RESNOR=RESNOR/ABS(PT)
                    IF(RESNOR.LE.TOL)THEN
                      SUFAIL=.FALSE.
                      APEX=.TRUE.
                      DGAMA=DEPV
                      DGAMB=R0
!                     update principal stresses
                      S1=P
                      S2=P
                      S3=P
                      EXIT
                    ENDIF
                  ENDDO
                  IF(SUFAIL.EQ..TRUE.)THEN
!                   throw some errors
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
!         update internal variable EPBAR  and stress components
!         -----------------------------------------------------
          DO I=1,6
            Sig(I)=S1*EIGPRJ(1,I)+S2*EIGPRJ(2,I)+S3*EIGPRJ(3,I)
          ENDDO
!         and elastic engineering strain
          EEvD3=P/BULK*R1D3
          StVar(1)=(Sig(1)-P)/R2G+EEvD3
          StVar(2)=(Sig(2)-P)/R2G+EEvD3
          StVar(3)=(Sig(3)-P)/R2G+EEvD3
          StVar(4)=Sig(4)/GMODU
          StVar(5)=Sig(5)/GMODU
          StVar(6)=Sig(6)/GMODU
          StVar(7)=EPBAR
          StVar(8)=DGAMA
          StVar(9)=DGAMB
!         update plastic state
          IF(APEX.EQ..TRUE.)THEN
            ipl=2
            StVar(10)=-1d0
          ELSE
            ipl=1
            IF(EDGE.EQ..FALSE.)THEN
              StVar(10)=1d0
            ELSE
              IF(RIGHT.EQ..TRUE.)THEN
                StVar(10)=2d0
              ELSE
                StVar(10)=3d0
              ENDIF
            ENDIF
          ENDIF
        ELSE
!         Elastic step: update stress using linear elastic law
          Sig(1)=Strial(1)
          Sig(2)=Strial(2)
          Sig(3)=Strial(3)
          Sig(4)=Strial(4)
          Sig(5)=Strial(5)
          Sig(6)=Strial(6)
!         elastic engineering strain
          StVar(1)=dEps(1)+dEps(7)
          StVar(2)=dEps(2)+dEps(8)
          StVar(3)=dEps(3)+dEps(9)
          StVar(4)=dEps(4)+dEps(10)
          StVar(5)=dEps(5)+dEps(11)
          StVar(6)=dEps(6)+dEps(12)
          StVar(7)=0d0
          StVar(8)=0d0
          StVar(9)=0d0
          StVar(10)=0d0
!         update plastic state
          ipl = 0
        END IF

      End If  ! IDTask = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 3) Then
!       { Create effective material stiffness matrix D }
!       Current accumulated plastic strain
        EPBAR=StVar(7)
!       Set material properties
        YOUNG=PROPS(1)
        POISS=PROPS(2)
        SINPHI=PROPS(3)
        COSPHI=PROPS(4)
        SINPSI=PROPS(5)
        C=PROPS(6)
        HARD=PROPS(7)
!       Set some constants
        GMODU=YOUNG/(R2*(R1+POISS))
        BULK=YOUNG/(R3*(R1-R2*POISS))
        R2G=R2*GMODU
        R4G=R4*GMODU
        R2BULK=R2*BULK
        R2CPHI=R2*COSPHI
        R4C2PH=R2CPHI*R2CPHI
        R1D3=R1/R3
        R2D3=R2*R1D3
        R2GD3=R2G*R1D3
        R4GD3=R4G*R1D3
        IF(ipl.GT.0)THEN
!         Compute elastoplastic consistent tangent
!         Spectral decomposition of the elastic trial strain
          Stra(1)=dEps(1)+dEps(7)
          Stra(2)=dEps(2)+dEps(8)
          Stra(3)=dEps(3)+dEps(9)
          Stra(4)=(dEps(4)+dEps(10))*RP5
          Stra(5)=(dEps(5)+dEps(11))*RP5
          Stra(6)=(dEps(6)+dEps(12))*RP5
          iOpt = 1
!         Compute the principal strain and eigenprojection tensors
          Call PrnSig(iOpt,Stra,EIGVEC(3,:),EIGVEC(2,:),EIGVEC(1,:),
     *          PSTRA3,PSTRA2,PSTRA1,P,Q)
          DO I=1,3
            EIGPRJ(I,1)=EIGVEC(I,1)*EIGVEC(I,1)
            EIGPRJ(I,2)=EIGVEC(I,2)*EIGVEC(I,2)
            EIGPRJ(I,3)=EIGVEC(I,3)*EIGVEC(I,3)
            EIGPRJ(I,4)=EIGVEC(I,1)*EIGVEC(I,2)
            EIGPRJ(I,5)=EIGVEC(I,2)*EIGVEC(I,3)
            EIGPRJ(I,6)=EIGVEC(I,3)*EIGVEC(I,1)
          ENDDO
!         and current total stress
          DO I=1,3
            PSTRS(I)=Sig(1)*EIGPRJ(1,1)+Sig(2)*EIGPRJ(1,2)+
     *             Sig(3)*EIGPRJ(1,3)+RP5*(Sig(4)*EIGPRJ(1,4)+
     *             Sig(5)*EIGPRJ(1,5)+Sig(6)*EIGPRJ(1,6))
          ENDDO
          IF(ipl.EQ.1)THEN
            IF(StVar(10).GT.1d0)THEN
!             Tangent consistent with 2-vector return to edge
              SPHSPS=SINPHI*SINPSI
              CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS
              IF(StVar(10).EQ.2d0)THEN
                CONSTB=R2G*(R1+SINPHI+SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS
              ELSE
                CONSTB=R2G*(R1-SINPHI-SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS
              ENDIF
              FACTA=R4C2PH*HARD
              DRVAA=-CONSTA-FACTA
              DRVAB=-CONSTB-FACTA
              DRVBA=-CONSTB-FACTA
              DRVBB=-CONSTA-FACTA
              AUX1=R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI
              AUX2=(R4GD3-R2BULK)*SINPSI
              AUX3=R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI
              R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA)
              IF(StVar(10).EQ.2d0)THEN
!               ...returned to right edge
                DPSTRS(1,1)=BULK+R4GD3+AUX1*(-DRVAB+DRVBB+DRVAA-DRVBA)*
     *                  (R2G+(R2BULK+R2GD3)*SINPHI)*R1DDET
                DPSTRS(1,2)=BULK-R2GD3+AUX1*(R2G*(DRVAB-DRVAA)+
     *                  ((-DRVAB+DRVBB+DRVAA-DRVBA)*(R2BULK+R2GD3)+
     *                  (DRVBA-DRVBB)*R2G)*SINPHI)*R1DDET
                DPSTRS(1,3)=BULK-R2GD3+AUX1*(R2G*(DRVBA-DRVBB)+
     *                  ((-DRVAB+DRVBB+DRVAA-DRVBA)*(R2BULK+R2GD3)+
     *                  (DRVAB-DRVAA)*R2G)*SINPHI)*R1DDET
                DPSTRS(2,1)=BULK-R2GD3+(AUX2*(DRVAB-DRVBB)+AUX3*(DRVBA-
     *                  DRVAA))*(R2G+(R2BULK+R2GD3)*SINPHI)*R1DDET
                DPSTRS(2,2)=BULK+R4GD3+(AUX2*((R2BULK*(DRVAB-DRVBB)+
     *                  (DRVAB*R2GD3+DRVBB*R4GD3))*SINPHI-DRVAB*R2G)+
     *                  AUX3*(DRVAA*R2G+(R2BULK*(DRVBA-DRVAA)-
     *                  (DRVAA*R2GD3+DRVBA*R4GD3))*SINPHI))*R1DDET
                DPSTRS(2,3)=BULK-R2GD3+(AUX2*((R2BULK*(DRVAB-DRVBB)-
     *                  (DRVBB*R2GD3+DRVAB*R4GD3))*SINPHI+DRVBB*R2G)+
     *                  AUX3*((R2BULK*(DRVBA-DRVAA)+(DRVAA*R4GD3+
     *                  DRVBA*R2GD3))*SINPHI-DRVBA*R2G))*R1DDET
                DPSTRS(3,1)=BULK-R2GD3+((AUX2*(DRVBA-DRVAA)+AUX3*(DRVAB-
     *                  DRVBB))*((R2BULK+R2GD3)*SINPHI+R2G))*R1DDET
                DPSTRS(3,2)=BULK-R2GD3+(AUX2*(((R2BULK*(DRVBA-DRVAA)-
     *                  (DRVBA*R4GD3+DRVAA*R2GD3))*SINPHI)+DRVAA*R2G)+
     *                  AUX3*(((R2BULK*(DRVAB-DRVBB)+(DRVAB*R2GD3+
     *                  DRVBB*R4GD3))*SINPHI)-DRVAB*R2G))*R1DDET
                DPSTRS(3,3)=BULK+R4GD3+(AUX2*(((R2BULK*(DRVBA-DRVAA)+
     *                  (DRVAA*R4GD3+DRVBA*R2GD3))*SINPHI)-DRVBA*R2G)+
     *                  AUX3*(((R2BULK*(DRVAB-DRVBB)-(DRVAB*R4GD3+
     *                  DRVBB*R2GD3))*SINPHI)+DRVBB*R2G))*R1DDET
              ELSE
!               ...returned to left edge
                DPSTRS(1,1)=BULK+R4GD3+(AUX1*(((R2BULK*(DRVBB-DRVAB)+
     *                  (DRVAB*R4GD3+DRVBB*R2GD3))*SINPHI)+DRVBB*R2G)+
     *                  AUX2*(((R2BULK*(DRVBA-DRVAA)+(DRVAA*R4GD3+
     *                  DRVBA*R2GD3))*SINPHI)+DRVBA*R2G))*R1DDET
                DPSTRS(1,2)=BULK-R2GD3+(AUX1*(((R2BULK*(DRVBB-DRVAB)-
     *                  (DRVAB*R2GD3+DRVBB*R4GD3))*SINPHI)-DRVAB*R2G)+
     *                  AUX2*(((R2BULK*(DRVBA-DRVAA)-(DRVAA*R2GD3+
     *                  DRVBA*R4GD3))*SINPHI)-DRVAA*R2G))*R1DDET
                DPSTRS(1,3)=BULK-R2GD3+((AUX1*(DRVBB-DRVAB)+AUX2*(DRVBA-
     *                  DRVAA))*(((R2BULK+R2GD3)*SINPHI)-R2G))*R1DDET
                DPSTRS(2,1)=BULK-R2GD3+(AUX1*(((R2BULK*(DRVAA-DRVBA)-
     *                  (DRVAA*R4GD3+DRVBA*R2GD3))*SINPHI)-DRVBA*R2G)+
     *                  AUX2*(((R2BULK*(DRVAB-DRVBB)-(DRVAB*R4GD3+
     *                  DRVBB*R2GD3))*SINPHI)-DRVBB*R2G))*R1DDET
                DPSTRS(2,2)=BULK+R4GD3+(AUX1*(((R2BULK*(DRVAA-DRVBA)+
     *                  (DRVAA*R2GD3+DRVBA*R4GD3))*SINPHI)+DRVAA*R2G)+
     *                  AUX2*(((R2BULK*(DRVAB-DRVBB)+(DRVAB*R2GD3+
     *                  DRVBB*R4GD3))*SINPHI)+DRVAB*R2G))*R1DDET
                DPSTRS(2,3)=BULK-R2GD3+((AUX1*(DRVAA-DRVBA)+AUX2*(DRVAB-
     *                  DRVBB))*(((R2BULK+R2GD3)*SINPHI)-R2G))*R1DDET
                DPSTRS(3,1)=BULK-R2GD3+(AUX3*(((R2BULK*(DRVAB-DRVBB-DRVAA+
     *                  DRVBA)+(DRVAA-DRVAB)*R4GD3+(DRVBA-DRVBB)*
     *                  R2GD3)*SINPHI)+(DRVBA-DRVBB)*R2G))*R1DDET
                DPSTRS(3,2)=BULK-R2GD3+(AUX3*(((R2BULK*(DRVAB-DRVBB-DRVAA+
     *                  DRVBA)+(DRVAB-DRVAA)*R2GD3+(DRVBB-DRVBA)*
     *                  R4GD3)*SINPHI)+(DRVAB-DRVAA)*R2G))*R1DDET
                DPSTRS(3,3)=BULK+R4GD3+(AUX3*(DRVAB-DRVBB-DRVAA+DRVBA)*
     *                  (((R2BULK+R2GD3)*SINPHI)-R2G))*R1DDET
              ENDIF
            ELSE
!             Tangent consistent with 1-vector return to main active plane
              SPHSPS=SINPHI*SINPSI
              CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS
              DENOM=-CONSTA-R4C2PH*HARD
              B1=(R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI)/DENOM
              B2=(R4G*R1D3-R2BULK)*SINPSI/DENOM
              B3=(R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI)/DENOM
              DPSTRS(1,1)=R2G*(R2D3+B1*(R1+R1D3*SINPHI))+
     *                BULK*(R1+R2*B1*SINPHI)
              DPSTRS(1,2)=R1D3*(R3*BULK-R2G)*(R1+R2*B1*SINPHI)
              DPSTRS(1,3)=R2G*(-R1D3-B1*(R1-R1D3*SINPHI))+
     *                BULK*(R1+R2*B1*SINPHI)
              DPSTRS(2,1)=R2G*(-R1D3-B2*(R1+R1D3*SINPHI))+
     *                BULK*(R1-R2*B2*SINPHI)
              DPSTRS(2,2)=R4G*R1D3*(R1+B2*SINPHI)+BULK*(R1-R2*B2*SINPHI)
              DPSTRS(2,3)=R2G*(-R1D3+B2*(R1-R1D3*SINPHI))+
     *                BULK*(R1-R2*B2*SINPHI)
              DPSTRS(3,1)=R2G*(-R1D3-B3*(R1+R1D3*SINPHI))+
     *                BULK*(R1-R2*B3*SINPHI)
              DPSTRS(3,2)=R1D3*(R3*BULK-R2G)*(R1-R2*B3*SINPHI)
              DPSTRS(3,3)=R2G*(R2D3+B3*(R1-R1D3*SINPHI))+
     *                BULK*(R1-R2*B3*SINPHI)
            ENDIF
          ELSEIF(ipl.EQ.2)THEN
!           Tangent consistent with multi-vector return to apex
            COTPHI=COSPHI/SINPHI
            DSIDEJ=BULK*(R1-(BULK/(BULK+HARD*COTPHI*COSPHI/SINPSI)))
            DPSTRS(1,1)=DSIDEJ
            DPSTRS(1,2)=DSIDEJ
            DPSTRS(1,3)=DSIDEJ
            DPSTRS(2,1)=DSIDEJ
            DPSTRS(2,2)=DSIDEJ
            DPSTRS(2,3)=DSIDEJ
            DPSTRS(3,1)=DSIDEJ
            DPSTRS(3,2)=DSIDEJ
            DPSTRS(3,3)=DSIDEJ
          ENDIF
!         TODO: I dont't know how to handle DGISO2 here
          
        ELSE
!         Compute elastic matrix
          FACTOR=BULK-R2G*R1D3
          DO I=1,6
            DO J=I,6
              D(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J)
            END DO
          END DO
          DO J=1,5
            DO I=J+1,6
              D(I,J)=D(J,I)
            END DO
          END DO
        ENDIF
      End If  ! IDTask = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 4) Then
!       { Return the number of state variables nStat }
        nStat=10
      End If  ! IDTask = 4
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 5) Then
!       { Return matrix attributes NonSym, iStrsDep, iTimeDep }
        NonSym=1
        iTang=1
        iStrsDep=0
        iTimeDep=0
      End If  ! IDTask = 5
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 6) Then
!       { Create elastic material stiffness matrix De }
        G       =   Props(1)       ! G
        xNu     =   Props(2)       ! nu
        F1  = 2*G*(1-xNu)/(1-2*xNu)
        F2  = 2*G*( xNu )/(1-2*xNu)
        Call MZeroR(D,36)
        Do i=1,3
          Do j=1,3
            D(i,j) = F2
          End Do
          D(i,i) = F1
          D(i+3,i+3) = G
        End Do
      End If  ! IDTask = 6
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 6!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Return
      End ! UdsmMC


