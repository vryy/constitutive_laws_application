      Subroutine UdsmDP ( IDTask, iMod, IsUndr,
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
! Purpose: User supplied soil model for Drucker-Prager
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
      Logical SMOOTH,APEX,SUFAIL

      Dimension Props(*), Sig0(*), StVar0(*), dEps(*), D(6,6),
     *          Sig(*),   StVar(*),
     *          Strial(6), EEtD(6), UNIDEV(6), SOID(6), FOID(6,6)

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
      End If  ! IDTask = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 2) Then
!       { Calculate constitutive stresses Sig (and Swp) }
!       Initialize some algorithmic and internal variables
        DGAMA=0d0
        ipl=0
        EPBARN=StVar0(7)
        EPBAR=EPBARN
!       Set some material properties
        YOUNG=PROPS(1)
        POISS=PROPS(2)
        ETA=PROPS(3)
        XI=PROPS(4)
        ETABAR=PROPS(5)
        C=PROPS(6)
        HARD=PROPS(7)
!       and some constants
        GMODU=YOUNG/(2d0*(1d0+POISS))
        BULK=YOUNG/(3d0*(1d0-2d0*POISS))
        R2G=2d0*GMODU
        R1D3=1d0/3d0
!       Compute elastic trial state
!       ---------------------------
!       Elastic trial volumetric strain and pressure stress
        EEtV=dEps(1)+dEps(2)+dEps(3)+dEps(7)+dEps(8)+dEps(9)
        PT=BULK*EEtV
!       Elastic trial deviatoric stress
        EEvD3=EEtV*R1D3
        Strial(1)=R2G*(dEps(1)+dEps(7)-EEVD3)
        Strial(2)=R2G*(dEps(2)+dEps(8)-EEVD3)
        Strial(3)=R2G*(dEps(3)+dEps(9)-EEVD3)
!       Shear component
        Strial(4)=GMODU*(dEps(4)+dEps(10))
        Strial(5)=GMODU*(dEps(5)+dEps(11))
        Strial(6)=GMODU*(dEps(6)+dEps(12))
!       Compute elastic trial stress J2 invariant and cohesion
        VARJ2t=Strial(4)*Strial(4)+Strial(5)*Strial(5)+
     *         Strial(6)*Strial(6)+
     *         0.5d0*(Strial(1)*Strial(1)+Strial(2)*Strial(2)+
     *                Strial(3)*Strial(3))
        COHE=C+HARD*EPBARN
!       Check for plastic consistency
!       -----------------------------
        SQRJ2t=SQRT(VARJ2t)
        PHI=SQRJ2t+ETA*PT-XI*COHE
        RES=PHI
        IF(COHE.NE.0d0)RES=RES/ABS(COHE)
        IF(RES.LT.TOL)THEN
!         Elastic step: update stress using linear elastic law
          Sig(1)=Strial(1)+PT
          Sig(2)=Strial(2)+PT
          Sig(3)=Strial(3)+PT
          Sig(4)=Strial(4)
          Sig(5)=Strial(5)
          Sig(6)=Strial(6)
!         elastic engineering strain; TODO: check
!          StVar(1)=StVar0(1)+dEps(1)
!          StVar(2)=StVar0(2)+dEps(2)
!          StVar(3)=StVar0(3)+dEps(3)
!          StVar(4)=StVar0(4)+dEps(4)
!          StVar(5)=StVar0(5)+dEps(5)
!          StVar(6)=StVar0(6)+dEps(6)
          StVar(1)=dEps(1)+dEps(7)
          StVar(2)=dEps(2)+dEps(8)
          StVar(3)=dEps(3)+dEps(9)
          StVar(4)=dEps(4)+dEps(10)
          StVar(5)=dEps(5)+dEps(11)
          StVar(6)=dEps(6)+dEps(12)
!         update the plastic state
          ipl = 0
!         update the plastic state
          StVar(7)=0d0
          StVar(8)=0d0
        ELSE
!         Plastic step: Use return mapping
!         Apply return mapping to smooth portion of cone
          SMOOTH=.FALSE.
          APEX=.FALSE.
!         -------------------------------------------------------------------
          SUFAIL=.TRUE.
          DO IPTER1=1,MAXRT
!           Compute residual derivative
            DENOM=-GMODU-BULK*ETABAR*ETA-XI*XI*HARD
!           Compute Newton-Raphson increment and update variable DGAMA
            DDGAMA=-PHI/DENOM
            DGAMA=DGAMA+DDGAMA
!           Compute new residual
            EPBAR=EPBARN+XI*DGAMA
            COHE=C+HARD*EPBAR
            SQRJ2=SQRJ2t-GMODU*DGAMA
            P=PT-BULK*ETABAR*DGAMA
            PHI=SQRJ2+ETA*P-XI*COHE
!           Check convergence
            RESNOR=ABS(PHI)
            IF(COHE.NE.0d0)RESNOR=RESNOR/ABS(COHE)
            IF(RESNOR.LE.TOL)THEN
              SUFAIL=.FALSE.
!             Check validity of return to smooth portion
              IF(SQRJ2.GE.0d0)THEN
!               results are valid, update stress components and other variables
                IF(SQRJ2t.EQ.0d0)THEN
                  FACTOR=0d0
                ELSE
                  FACTOR=1d0-GMODU*DGAMA/SQRJ2t
                END IF
                SMOOTH=.TRUE.
                APEX=.FALSE.
                EXIT
              ELSE
!               smooth wall return not valid - go to apex return procedure
                SMOOTH=.FALSE.
                APEX=.TRUE.
                EXIT
              END IF
            END IF
          END DO
          IF(SUFAIL.EQV..TRUE.)THEN
!           TODO: throw some errors
          END IF
          IF(APEX.EQV..TRUE.)THEN
!           Apply return mapping to APEX
!           --------------------------------------------------
!           perform checks and set some variables
            ALPHA=XI/ETABAR
            BETA=XI/ETA
!           Set initial guess for unknown DEPV and start iterations
            DEPV=0d0
            EPBAR=EPBARN
            COHE=C+HARD*EPBAR
            RES=BETA*COHE-PT
            SUFAIL=.TRUE.
            DO IPTER2=1,MAXRT
              DENOM=ALPHA*BETA*HARD+BULK
!             Compute Newton-Raphson increment and update variable DEPV
              DDEPV=-RES/DENOM
              DEPV=DEPV+DDEPV
!             Compute new residual
              EPBAR=EPBARN+ALPHA*DEPV
              COHE=C+HARD*EPBAR
              P=PT-BULK*DEPV
              RES=BETA*COHE-P
!             Check convergence
              RESNOR=ABS(RES)
              IF(COHE.NE.0d0)RESNOR=RESNOR/ABS(COHE)
              IF(RESNOR.LE.TOL)THEN
                SUFAIL=.FALSE.
!               update stress components and other variables
                DGAMA=DEPV/ETABAR
                FACTOR=0d0
                EXIT
              END IF
            END DO
            IF(SUFAIL.EQV..TRUE.)THEN
!             TODO: throw some errors
            END IF
          END IF
!          WRITE(*,*),"SMOOTH:",SMOOTH,"APEX:",APEX
!         Stress update successful; store converged stress components and other state variables
          Sig(1)=FACTOR*Strial(1)+P
          Sig(2)=FACTOR*Strial(2)+P
          Sig(3)=FACTOR*Strial(3)+P
          Sig(4)=FACTOR*Strial(4)
          Sig(5)=FACTOR*Strial(5)
          Sig(6)=FACTOR*Strial(6)
!         update EPBAR
          StVar(7)=EPBAR
!         update DGAM
          StVar(8)=DGAMA
!         compute converged elastic (engineering) strain components
          FACTOR=FACTOR/R2G
          EEvD3=P/(BULK*3d0)
          StVar(1)=FACTOR*Strial(1)+EEvD3
          StVar(2)=FACTOR*Strial(2)+EEvD3
          StVar(3)=FACTOR*Strial(3)+EEvD3
          StVar(4)=FACTOR*Strial(4)*2d0
          StVar(5)=FACTOR*Strial(5)*2d0
          StVar(6)=FACTOR*Strial(6)*2d0
!         update the plastic state
          IF(SMOOTH.EQV..TRUE.)THEN
            ipl=1
          END IF
          IF(APEX.EQV..TRUE.)THEN
            ipl=2
          END IF
        END IF ! RES.GT.TOL
      End If  ! IDTask = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 3) Then
!       { Create effective material stiffness matrix D }
!       Retrieve accumulated plastic strain, DGAMA and APEX algorithm flag
        EPBAR=StVar(7)
        DGAMA=StVar(8)
!       Set some material properties
        YOUNG=PROPS(1)
        POISS=PROPS(2)
        ETA=PROPS(3)
        XI=PROPS(4)
        ETABAR=PROPS(5)
        C=PROPS(6)
        HARD=PROPS(7)
!       and some constants
        GMODU=YOUNG/(2d0*(1d0+POISS))
        BULK=YOUNG/(3d0*(1d0-2d0*POISS))
        R2G=2d0*GMODU
        R1D3=1d0/3d0
        ROOT2=SQRT(2d0)
        IF(ipl.EQ.0)THEN
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
        ELSE
!         Compute elastoplastic consistent tangent
!         Hardening slope
          HSLOPE=HARD
          IF(ipl.EQ.2)THEN
!           Elastoplastic tangent consistent with apex return
            ALPHA=XI/ETABAR
            BETA=XI/ETA
            AFACT=BULK*(1d0-BULK/(BULK+ALPHA*BETA*HSLOPE))
            DO I=1,6
              DO J=1,6
                D(I,J)=AFACT*SOID(I)*SOID(J)
              END DO
            END DO
          ELSEIF(ipl.EQ.1)THEN
!           Elastoplastic tangent consistent with smooth cone wall return
!           Elastic trial deviatoric (physical) strain
            EEvD3=(dEps(1)+dEps(2)+dEps(3)+dEps(7)+dEps(8)+dEps(9))*R1D3
            EEtD(1)=dEps(1)+dEps(7)-EEvD3
            EEtD(2)=dEps(2)+dEps(8)-EEvD3
            EEtD(3)=dEps(3)+dEps(9)-EEvD3
            EEtD(4)=(dEps(4)+dEps(10))*0.5d0
            EEtD(5)=(dEps(5)+dEps(11))*0.5d0
            EEtD(6)=(dEps(6)+dEps(12))*0.5d0
            ETDNOR=SQRT(EEtD(1)*EEtD(1)+EEtD(2)*EEtD(2)+
     *             EEtD(3)*EEtD(3)+2d0*(EEtD(4)*EEtD(4)+
     *             EEtD(5)*EEtD(5)+EEtD(6)*EEtD(6)))
!           Unit deviatoric flow vector
            IF(ETDNOR.NE.0d0)THEN
              EDNINV=1d0/ETDNOR
            ELSE
              EDNINV=0d0
            ENDIF
            DO I=1,6
              UNIDEV(I)=EETD(I)*EDNINV
            END DO
!           Assemble tangent
            AUX=1d0/(GMODU+BULK*ETA*ETABAR+XI*XI*HSLOPE)
            AFACT=R2G*(1d0-DGAMA/(ROOT2*ETDNOR))
            AFACD3=AFACT*R1D3
            BFACT=R2G*(DGAMA/(ROOT2*ETDNOR)-GMODU*AUX)
            CFACT=-ROOT2*GMODU*BULK*AUX
            DFACT=BULK*(1d0-BULK*ETA*ETABAR*AUX)
            DO I=1,6
              DO J=1,6
                D(I,J)=AFACT*FOID(I,J)+BFACT*UNIDEV(I)*UNIDEV(J)+
     *                     CFACT*(ETA*UNIDEV(I)*SOID(J)+
     *                           ETABAR*SOID(I)*UNIDEV(J))+
     *                   (DFACT-AFACD3)*SOID(I)*SOID(J)
              END DO
            END DO
          ELSE
!           TODO: throw some errors
          ENDIF
        END IF
      End If  ! IDTask = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!END IDTask = 3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (IDTask .Eq. 4) Then
!       { Return the number of state variables nStat }
        nStat=8
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
      End ! UdsmDP


