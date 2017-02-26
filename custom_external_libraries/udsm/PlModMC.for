      Subroutine PlModMC ( IDTask, iMod, IsUndr,
     *                     iStep, iTer, iEl, Int,
     *                     X, Y, Z,
     *                     Time0, dTime,
     *                     Props, Sig0, Swp0, StVar0,
     *                     dEps, D, BulkW,
     *                     Sig, Swp, StVar, ipl,
     *                     nStat,
     *                     NonSym, iStrsDep, iTimeDep, iTang,
     *                     iAbort )
!
! Purpose: User supplied soil model
!          (Example: Mohr-Coulomb with tension cut-off)
! Remark: be careful, this is an demonstration example; it is not the full Mohr-Coulomb with tension cut-off model
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
     *          Sig(*),   StVar(*)
!
!---  Local variables
!
      Dimension DE(6,6), dSig(6), Prs_E(3), Prs(3),
     *          xN1(3), xN2(3), xN3(3)

      Data Pi/3.14159 26535 89793 23846 26433 83279 50288 41971 69399d0/

      io=0
      If (iEl.Eq.87 .And. Int.Eq.12) io=1
      nStatV = 0

! Contents of Props() (iMod=4, MC)
!  1 : G       shear modulus
!  2 : xNu     Poisson's ratio
!  3 : C       Cohesion
!  4 : Phi     Friction angle (degrees)
!  5 : Psi     Dilation angle (degrees)
!  6 : Tens    Allowable tensile stress
!  7 : cCosPhi C*Cos(Phi)  will be filled during IDTask = 1
!  8 : sPhi    Sin(Phi)    will be filled during IDTask = 1
!  9 : sPsi    Sin(Psi)    will be filled during IDTask = 1

      If (IDTask .Eq. 1) Then ! Initialize state variables
        ! Nothing to do here but also derive some properties
        Call MZeroR( StVar0, nStatV )
        Call MZeroR( StVar , nStatV )
        Rad  = 180d0 / Pi
!       G       = Props(1)       ! G
!       xNu     = Props(2)       ! nu
        C       = Props(3)       ! C
        Phi     = Props(4) / Rad ! Phi in radians
        Psi     = Props(5) / Rad ! Psi in radians
!       sTens   =-Props(6)       ! allowable tensile stress
        sPhi    = Sin(Phi)
        sPsi    = Sin(Psi)
        cCosPhi = C*Cos(Phi)

        Props(7) = sPhi
        Props(8) = sPsi
        Props(9) = cCosPhi
!        Call WriVal( io, 'Phi',phi)
!        Call WriVal( io, 'Psi',psi)
!        Call WriVal( io, 'sPhi',sphi)
!        Call WriVal( io, 'sPsi',spsi)
!        Call WriVec( io, 'Props', Props, 10)
      End If  ! IDTask = 1

      If (IDTask .Eq. 2) Then ! Calculate stresses
        Call CopyRVec( StVar0, StVar, nStatV )
        ipl     =   0
        G       =   Props(1)       ! G
        xNu     =   Props(2)       ! nu
        sTens   =   Props(6)       ! tensile strength (change sign)
        sPhi    =   Props(7)
        sPsi    =   Props(8)
        cCosPhi =   Props(9)
        If (sPhi.Gt.0) Then
          If (sTens.Gt.cCosPhi/sPhi) sTens = cCosPhi/sPhi
        End If
        sTens = - sTens

        If (IsUndr.Eq.1) Then
          xNu_U = 0.495d0 ! Undrained Poissons' ratio
          Fac=(1+xNu_U)/(1-2*xNu_U) - (1+xNu)/(1-2*xNu)
          Fac=2D0*G/3D0  * Fac
          BulkW = Fac
          dEpsV = dEps(1) + dEps(2) + dEps(3)
          dSwp  = BulkW * dEpsV
          Swp   = Swp0 + dSwp
        Else
          Swp = Swp0
        End If

        ! Fill elastic material matrix
        F1  = 2*G*(1-xNu)/(1-2*xNu)
        F2  = 2*G*( xNu )/(1-2*xNu)
        Call MZeroR(DE,36)
        Do i=1,3
          Do j=1,3
            DE(i,j) = F2
          End Do
          DE(i,i) = F1
          DE(i+3,i+3) = G
        End Do

!        If (iEl+Int+iter.Eq.3 .And. iStep.Lt.10) Then
!          Call WriMat( io, 'DE66', DE, 6, 6, 6 )
!        End If
        ! elastic stress increment
        Call MatVec( DE, 6, dEps, 6, dSig)
        ! elastic stress
        Call AddVec( Sig0, dSig, 1d0, 1d0, 6, Sig )

        ! calculate principal stresses and directions
        iOpt = 1
        Call PrnSig(iOpt, Sig, xN1, xN2, xN3, S1, S2, S3, P, Q)
        ! Sig    : tension     positive
        ! Prs(E) : compression positive
        Prs_E(1) = - S1 ! minus sign
        Prs_E(2) = - S2
        Prs_E(3) = - S3
        iArea = 2
        Call MC_Tens( iArea, G, xNu, sPhi, sPsi, cCosPhi, sTens,
     *                Prs_E, Prs, ipl )
        If (ipl.Ne.0) Then ! some plasticity
          ! Check Sig1 > Sig2 > Sig3
          If (                 Prs(2).Lt.Prs(3)) iarea = 1   ! Tr. compression
          If (IArea.Eq.2 .And. Prs(1).Lt.Prs(2)) iarea = 3   ! Tr. extension
          If (iArea.Ne.2) Then
            Call MC_Tens( iArea, G, xNu, sPhi, sPsi, cCosPhi, sTens,
     *                    Prs_E, Prs, ipl )
          End If
          ! Prs : compression positive
          S1 = - Prs(1) ! minus sign
          S2 = - Prs(2)
          S3  =- Prs(3)
          ! back to Cartesian stresses
          Call CarSig(S1,S2,S3,xN1,xN2,xN3,Sig)
          ! Sig    : tension positive
        End If
      End If ! IDTask = 2; get stresses

      If ( IDTask .Eq. 3 .Or.
     *     IDTask .Eq. 6     ) Then ! Calculate D-Matrix

        ! Always Elastic D-matrix | Why??
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
        BulkW = 0
        If (IsUndr.Eq.1) Then
          ! BulkW = ...
          xNu_U = 0.495d0
          Fac=(1+xNu_U)/(1-2*xNu_U) - (1+xNu)/(1-2*xNu)
          Fac=2D0*G/3D0  * Fac
          BulkW = Fac
        End If
      End If  ! IDTask = 3, 6

      If (IDTask .Eq. 4) Then ! Number of state parameters
        nStat    = nStatV
      End If  ! IDTask = 4

      If (IDTask .Eq. 5) Then ! matrix type
        NonSym   = 0  ! 1 for non-symmetric D-matrix
        iStrsDep = 0  ! 1 for stress dependent D-matrix
        iTang    = 0  ! 1 for tangent D-matrix
        iTimeDep = 0  ! 1 for time dependent D-matrix
      End If  ! IDTask = 5

      Return
      End ! MyMod_MC
