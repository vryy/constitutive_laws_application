C**********************************************************************
C234567890123456789012345678901234567890123456789012345678901234567890
C**********************************************************************

      Subroutine PlDP ( IDTask, iMod, IsUndr,
     *                      iStep, iTer, iEl, Int,
     *                      X, Y, Z,
     *                      Time0, dTime,
     *                      Props, Sig0, Swp0, StVar0,
     *                      dEps, D, BulkW,
     *                      Sig, Swp, StVar, ipl,
     *                      nStat, NonSym, iStrsDep, iTimeDep,iTang,
     *                      iPrjDir, iPrjLen, iAbort )

      Implicit Double Precision (A-H, O-Z)

      Dimension Props(*), Sig0(*), StVar0(*), dEps(*), D(6,6),
     *          Sig(*),   StVar(*)

!     Local variables

      Dimension dSig(6), xN1(3), xN2(3), xN3(3)

!     Expected contents of Props(1..?)
!     1 : G    Shear modulus
!     2 : xNu  Poissons' ratio
!     3 : Alf  Steepness of yield surface in p-q-plane
!     4 : Bet  Steepness of plastic potential function in p-q-plane
!     5 : C    'cohesion' term f = q + Alf*p - C
!    10 : BFac Ratio K/G to calculate Bulk modulus

      nStatV = 1
      If (IDTask .Eq. 1) Then ! Initialize state variables
        ! When there are state variables one could check a maximum
        ! For instance:
        !   StVar(1) contains maximum isotropic stress
        !   StVar(2) contains maximum shear mobilisation
        ! p = -(Sig0(1)+Sig0(2)+Sig0(3))/3
        ! or Call PrnSig(Sig0, ...., P0,Q0)
        !    P = -P0 ! compression positive
        ! ratio = Q0/Max(P,1d0) ! mobilisation degree
        ! StVar0(1) = Max (StVar0(1), p )
        ! StVar0(2) = Max (StVar0(2), ratio )
        G   = Props(1) ! G
        xNu = Props(2) ! nu
        BFac= (2+2*xNu)/(3-6*xNu) ! Bulk modulus/G
        Props(10) = BFac
        Bulk= BFac * G
      End If  ! IDTask = 1

      If (IDTask .Eq. 2) Then ! Calculate stresses
        ! Initialize StateVar
!        Write(*,*) 'Initialize StateVar'
        If (StVar0(1) .EQ. 0) StVar0(1) = 0.5
        ! Calculate new pore volume
        dStran_v = (dEps(1)+dEps(2)+dEps(3)+
     2              dEps(1)*dEps(2)+dEps(1)*
     3              dEps(3)+dEps(2)*dEps(3)-
     4              dEps(1)*dEps(2)*dEps(3))
        StVar(1) = StVar0(1)+(1d0+StVar0(1))*dStran_v

!        Write(*,*) 'after StVar'

        !Call CopyRVec( StVar0, StVar, nStatV )
        If (IsUndr.Eq.1) Then
          dEpsV = dEps(1) + dEps(2) + dEps(3)
          dSwp  = BulkW * dEpsV
          Swp   = Swp0 + dSwp
        Else
          Swp = Swp0
        End If

!        Write(*,*) 'before MatVec'
!        Write(*,*) 'dEps: ', dEps(1), dEps(2), dEps(3), dEps(4), dEps(5), dEps(6)
!        Write(*,*) 'dEps: ', dSig(1), dSig(2), dSig(3), dSig(4), dSig(5), dSig(6)

        Call MatVec( D, 6, dEps, 6, dSig)

!        Write(*,*) 'after MatVec'

        Call AddVec( Sig0, dSig, 1d0, 1d0, 6, Sig )

!        Write(*,*) 'after AddVec'

        iOpt = 1
        Call PrnSig(iOpt, Sig, xN1, xN2, xN3, S1, S2, S3, P, Q)

        ipl = 0
        G    = Props(1)
        xNu  = Props(2) ! nu
        Alf  = Props(3)
        Bet  = Props(4)
        C    = Props(5)

        P0   = - ( Sig0(1) + Sig0(2) + Sig0(3) ) /3
        P0   = 100
        G    = P0 / 100 * Props(1)
        Bulk = Props(10) * G
        f = q + Alf * p - C

        If (f .Gt. 1d-6) Then ! plastic stress correction
          !  f = q + Alf * p - C
          !  p = pe - xLam * Bulk * Bet
          !  q = qe - xLam * 3*G
          pe = p
          qe = q
          dfdp = Alf
          dfdq = 1
          dpdl = -Bulk*Bet
          dqdl = -3*G
          dfdl = (dfdp*dpdl + dfdq*dqdl)
          xLam = -f / dfdl
          p = pe - xLam * Bulk * Bet
          q = qe - xLam * 3 * G
          f = q + Alf * p - C
          ipl = 1
          If (p.Gt.C/Alf) Then ! apex
            p = C/Alf
            q = 0
            ipl = 2
          End If
          If (qe.Ne.0) Then
            s1 =  p + q/qe*(s1-pe)
            s2 =  p + q/qe*(s2-pe)
            s3 =  p + q/qe*(s3-pe)
          Else
            s1 =  p
            s2 =  p
            s3 =  p
          End If
        End If
        ! as example: update state variables
        ! Pt = -P ! compression positive
        ! ratio = Q/Max(Pt,1d0) ! mobilisation degree
        ! StVar(2) = Max (StVar0(1), pt )
        ! StVar(3) = Max (StVar0(2), ratio )

        ! back to Cartesian stresses
        Call CarSig(S1,S2,S3,xN1,xN2,xN3,Sig)

!        If (iEl.EQ.1 .AND. Int.EQ.1) Then
!          Write(1,*)'My Output for S1',S1, Time0+dTime
!          Flush(1)
!        End If

!        If (ipl.Gt.0) Then ! just checking f~0
!          iOpt = 0
!          Call PrnSig(iOpt,Sig,xN1,xN2,xN3,S1,S2,S3,P,Q)
!          f = q + Alf * p - C
!        End If
      End If ! IDTask = 2

      If ( IDTask .Eq. 3 .Or.
     *     IDTask .Eq. 6     ) Then ! Calculate D-Matrix
        G   = Props(1) ! G
        P0  = - ( Sig0(1) + Sig0(2) + Sig0(3) ) /3
        P0  = 100
        G   = P0 / 100 * Props(1)
        xNu = Props(2) ! nu

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
      End ! Subroutine DruPra

