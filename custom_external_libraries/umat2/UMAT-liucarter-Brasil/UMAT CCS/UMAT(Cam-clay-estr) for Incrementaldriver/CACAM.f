C     UMAT MADE BY CRISTHIAN MENDOZA, UNIVERSIDADE DE BRASILIA, 2011
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      USE tensors
      INCLUDE 'ABA_PARAM.INC'
C
C --------------------------------------------------------------- 
C     Declarating UMAT vaiables and constants 
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFBULKGRD1(3,3),
     4 MAT(2,2)
C ---------------------------------------------------------------      
      DOUBLE PRECISION T(3,3), DEPS(3,3), JAC(3,3,3,3), evoid, epsp(3,3)
     1 , ps, lambda, kappa, psinitial, M, p, q, nu, b, w, cgn, a
     2 , Kbulk, E, ELMOD(3,3,3,3), Idev(3,3,3,3), TrialSTRESS(3,3)
     3 , Df(3,3), con, D(3,3), F, dotTrialSTRESS(3,3), DG(3,3)
     4 , Trialp, Trialq, TrialF, DeltaPhi, caab, eta, caag, psnew
     5 , normaldev(3,3), Depsp(3,3), NormDepsp, deini, de, ab, ad
     6 , En, Tdev(3,3), dotepsp(3,3), dotT(3,3), CEP(3,3,3,3)    
      integer maxiter, kcounter
C     
      CALL D1(DSTRAN, D, dtime, NDI, NSHR, NTENS) 
      CALL Initial(STRESS,T, DSTRAN, DEPS, NTENS,NDI, NSHR)
C ---------------------------------------------------------------      
C     Material Parameters
      lambda=PROPS(1) !Compresion slope
      kappa=PROPS(2) !swelling slope
      psinitial=PROPS(3) !initial yield stress
      M=PROPS(4) !Critic state slope
      nu=PROPS(5) !Poisson modulus
      b=PROPS(6) !destructuring index b
      w=PROPS(7) !the effect of additional voids ratio on flow rule
      deini=PROPS(8) !additional voids ratio sustained by soil structure
C ---------------------------------------------------------------
C     Save State variables
C     Save Void ratio
      evoid=statev(1)
C     Plastic strains as state variable      
      CALL Vectortomatrix(STATEV(2:7),epsp) 
C     Overconsolidated mean stress state variable
      ps=statev(8)
C     additional voids ratio sustained by structure as state variable
      de=statev(9)
	  psnew=statev(10)
C ---------------------------------------------------------------       
C Trial Elastic step
C ---------------------------------------------------------------      
      CALL pq(T,p,q)  ! Subrotine p and q 
C ---------------------------------------------------------------      
      if (ps==0.) then
      ps=abs(psinitial) !overconsolidated mean stress             
      endif 
      if (abs(p)<1.0d0) then 
      p=1.0d0
      endif
C     update of yield stress 
      if (abs(p)>abs(ps)) then 
      ps=p
      endif
      eta=q/p
C ---------------------------------------------------------------      
      psnew=((q**2.0d0)/((M**2.0d0)*p))+p
      if (psnew>psinitial) then
      ps=psnew
      else
      ps=psinitial
      end if
C ---------------------------------------------------------------      
C     OJO, PREGUNTAR ESTA	  
C     update of additional voids ratio
c      if (tr(D)==0.0d0) then
c      de=deini
c      else
      if (psinitial<psnew) then
      de=deini*((psinitial/psnew)**b)                       
      else
      de=deini 
      endif
c      endif
C --------------------------------------------------------------- 
C     Calculate the material modules (Material constants)
C     bulk modulus       
      kbulk=abs((1.0d0+evoid)*p)/(kappa)   
      if (abs(kbulk)<100.0d0) then 
      kbulk=100.0d0
      endif 
C     young modulus
	  E=3.0d0*kbulk*(1.0d0-2.0d0*nu)       
C     Subrotine of elastic tensor
      CALL ELMOD1(E, nu, ELMOD)
C ---------------------------------------------------------------
C     update of initial yield stress
      if ((ps>psinitial).or.ps==0.0d0) then
      ps=ps                             
      else
      ps=psinitial 
      endif
C ---------------------------------------------------------------
C     Step trial(elasticity)
C ---------------------------------------------------------------
C     Trail stress 
      dotTrialSTRESS=ELMOD.double.D
      TrialSTRESS=T+dotTrialSTRESS*dtime      
C ---------------------------------------------------------------     
C     Yield function F 
      F=((q**2.0d0)/(M**2.0d0))+(abs(p)*(abs(p)-abs(ps)))     
C ---------------------------------------------------------------     
      call pq(TrialSTRESS,Trialp,Trialq)    
C     Trial yield function F 
      TrialF=((Trialq**2.0d0)/M**2.0d0)
     1 +(abs(Trialp)*(abs(Trialp)-abs(ps)))
C --------------------------------------------------------------- 
C     Check yield function F=0
      if (TrialF<0.0d0) then
C	  Trial step is ok (now saving state variables)
      T=TrialSTRESS
      JAC=ELMOD
C     Change of void ratio
      evoid=evoid+Tr(D)*(1.0d0+evoid)*dtime
C --------------------------------------------------------------- 
C     Saving State variables 
C     Alpha isotropic as state variable (scalar)
      STATEV(1)=evoid
C     Plastic strains as state variable (tensor R3xR3)
      CALL Matrixtovector(STATEV(2:7),epsp)
      STATEV(8)=ps
C     Salving additional voids ratio
      statev(9)=de
	  statev(10)=psnew
C --------------------------------------------------------------- 
C     Correction for yield
C ---------------------------------------------------------------       
      else ! Plastic corrector step
C ---------------------------------------------------------------
C     Calculation of the plastic multiplier     
C ---------------------------------------------------------------
C     Tensor deviator of fourth order 
      CALL Idev1(Idev)
      Tdev=Idev.double.T
      normaldev=Tdev/norm(Tdev)
      normaldev=-normaldev
C ---------------------------------------------------------------         
C     Constant a
      ab=(1.0d0-(w*de))
      ad=(1.0d0-(2.0d0*w*de))
	  a=ab/ad
C --------------------------------------------------------------- 
C     w correction
      if (ab>1.0d0) then
      w=0.0d0                              
      else
      w=w 
      endif
      if (ab<0.0d0) then
      w=1.0d0/deini                              
      else
      w=w 
      endif	  
C ---------------------------------------------------------------
C     The plastic multiplier 
      DeltaPhi=0.0d0 !initiation delta phi
      con=(M**2.0d0)*((2.0d0*p)-ps) !constante
      cgn=(M**2.0d0)*((2.0d0*p)-((2.0d0-(1.0d0/a))*(ps**(1.0d0/a))
     1 *(p**((a-1.0d0)/a)))) !constante
      caab=sqrt(3.0d0/2.0d0)*(2.0d0*q)!constant
      caag=sqrt(3.0d0/2.0d0)*(2.0d0*(1.0d0-(2.0d0*w*de))*q)!constant
      En=-(((1.0d0+evoid)*(p)*(M**2.0d0)*ps)/((lambda-kappa)
     1 +(b*de*(1.0d0+(eta/(M-eta))))))*((M**2.0d0)*((2.0d0*p)
     2 -((2.0d0-(1.0d0/a))*(ps**(1.0d0/a))*(p**((a-1.0d0)/a)))))
      DG=((1.0d0/3.0d0)*cgn*delta)+(caag*normaldev)
      Df=((1.0d0/3.0d0)*con*delta)+(caab*normaldev)
      DeltaPhi=(Df.double.ELMOD).double.D
      DeltaPhi=DeltaPhi/(((Df.double.ELMOD).double.DG)-En)
c      DeltaPhi=abs(DeltaPhi)
C ---------------------------------------------------------------        
C     Plastic strain
      if  ((norm(Tdev))/=0.0d0) then
      dotepsp=DeltaPhi*DG
      Depsp=dotepsp*dtime
      epsP=epsP+Depsp
      NormDepsp=norm(Depsp)
      else
      Depsp=0.0d0
      endif
C ---------------------------------------------------------------
C     Calculate the stress with plastic strains
      dotT=ELMOD.double.(D-dotepsp)
      T=T+dotT*dtime
C --------------------------------------------------------------- 
C     Saving State variables 
C     Void ratio
      evoid=evoid+Tr(D)*(1.0d0+evoid)*dtime
C     salving void ratio (scalar)
      STATEV(1)=evoid
C     salving plastic strain (tensor R3xR3)
      CALL Matrixtovector(STATEV(2:7),epsp)
      STATEV(8)=ps
C     Salving additional voids ratio
      STATEV(9)=de
	  STATEV(10)=psnew
C ---------------------------------------------------------------  
C     Calculate the elasto-plastic module
      if  (NormDepsp/=0) then
      cep=(ELMOD.double.DG).dyad.(Df.double.ELMOD)
      cep=cep/(((Df.double.ELMOD).double.DG)-En)
      CEP=ELMOD-cep
      Jac=CEP
       else
      JAC=ELMOD
       endif
       endif
C ---------------------------------------------------------------      
      Call Solution(NTENS, NDI, NSHR, T, STRESS, JAC
     1 , DDSDDE) 
      END SUBROUTINE UMAT