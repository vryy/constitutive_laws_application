CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                           MADE BY CRISTHIAN MENDOZA, BRASILÍA UNIVERSITY, 2010                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                             UMAT NOT MADE FOR UNDRAINED ANALYSIS ​                                    C​
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      
CXXXX       1. TENSOR LIBRARY BY WILLIAM FUENTES, UNIVERSIDAD DE LOS ANDES, 2009         XXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      
      MODULE tensors
      DOUBLE PRECISION delta(3,3), sq2, sq3, sq4, PI
      INTEGER I,J,K,L
      Public delta, sq2, sq3, sq4, PI
      data delta/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/ 
      parameter (PI=3.1415926535897932384626433832795d0)
      parameter (sq2=1.4142135623730950488016887242097d0)
      parameter (sq3=1.7320508075688772935274463415059d0)
      parameter (sq6=2.4494897427831780981972840747059d0)
                  
c------------------------------------------------------ 
c------------------------------------------------------        
      INTERFACE OPERATOR (.dot.)     
      MODULE PROCEDURE dot1, dot2
      END INTERFACE            
c------------------------------------------------------ 
c------------------------------------------------------        
      INTERFACE OPERATOR (.dyad.)     
      MODULE PROCEDURE dyad33x33
      END INTERFACE 
c------------------------------------------------------ 
c------------------------------------------------------        
      INTERFACE OPERATOR (.double.)     
      MODULE PROCEDURE double3333x33, double33x33, 
     1 double3333x3333, double33x3333
      END INTERFACE 
c------------------------------------------------------ 
c------------------------------------------------------    
      INTERFACE norm     
      MODULE PROCEDURE norm33
      END INTERFACE 
c------------------------------------------------------ 
c------------------------------------------------------ 
      INTERFACE Tr    
      MODULE PROCEDURE Tr33
      END INTERFACE 
c------------------------------------------------------ 
c------------------------------------------------------ 
      INTERFACE det    
      MODULE PROCEDURE det33
      END INTERFACE     
c------------------------------------------------------ 
c------------------------------------------------------ 
      INTERFACE gMC    
      MODULE PROCEDURE gMC
      END INTERFACE          
c------------------------------------------------------ 
c------------------------------------------------------ 
      INTERFACE Inverse    
      MODULE PROCEDURE Inverse3333, inverse33
      END INTERFACE       
c------------------------------------------------------ 
c------------------------------------------------------         
c------------------------------------------------------ 
c------------------------------------------------------                   
      CONTAINS
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION gMC(theta, phic) RESULT (g) 
      DOUBLE PRECISION, INTENT(IN):: theta, phic  
      double precision g    
c     Returns the function g for the Mohr Coloumb shape
c     Phic in radians        
      g=sin(phic)/(cos(theta)+sin(theta)*sin(phic)/sq3)
      return
      END FUNCTION gMC        
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION Tr33(a) RESULT (b) 
      DOUBLE PRECISION, INTENT(IN):: a(3,3)   
      double precision b      
      b=a(1,1)+a(2,2)+a(3,3)
      return
      END FUNCTION Tr33  
c------------------------------------------------------ 
c------------------------------------------------------            
      FUNCTION det33(T) RESULT (det) 
      DOUBLE PRECISION, INTENT(IN):: T(3,3)   
      double precision det      
      det=T(1,1)*T(2,2)*T(3,3)+T(1,2)*T(2,3)*T(3,1)
     1 +T(1,3)*T(2,1)*T(3,2)-T(1,3)*T(2,2)*T(3,1)
     2 -T(1,2)*T(2,2)*T(3,3)-T(1,1)*T(2,3)*T(3,2)
      return
      END FUNCTION det33        
c------------------------------------------------------ 
c------------------------------------------------------       
      FUNCTION dot1(a1,a2) RESULT (a3)
      DOUBLE PRECISION, INTENT(IN),dimension(3,3) :: a1, a2
      DOUBLE PRECISION,dimension(3,3) :: a3
      INTEGER i,j, k
      a3=0.0d0
      Do i=1, 3
      Do j=1, 3
      Do k=1, 3
      a3(i,j)=a3(i,j)+a1(i,k)*a2(k,j)
      Enddo
      Enddo
      Enddo
      END FUNCTION dot1
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION dot2(a1,a2) RESULT (a3)
      DOUBLE PRECISION, INTENT(IN),dimension(3) :: a1, a2
      DOUBLE PRECISION :: a3
      INTEGER i
      a3=0.0d0
      Do i=1, 3
      a3=a3+a1(i)*a2(i)
      Enddo
      END FUNCTION dot2  
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE ELMOD1(E, nu, ELMOD) 
      DOUBLE PRECISION ELMOD(3,3,3,3), E, nu,
     1 Idev(3,3,3,3), Ivol(3,3,3,3), Isym(3,3,3,3)
      INTEGER I,J,K,L
      Idev=0.0d0
      Ivol=0.0d0
      ELMOD=0.0D0
      Isym=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Isym(i,j,k,l)=1.0d0/2.0d0*(delta(i,k)*delta(j,l)+
     1 delta(i,l)*delta(j,k))              
      Enddo
      Enddo
      Enddo
      Enddo
      Ivol=0.0D0
      Ivol=1.0d0/3.0d0*delta.dyad.delta    
      Idev=Isym-Ivol         
      ELMOD=E/(1.0D0-2.0D0*nu)*Ivol+E/(1.0d0+nu)*Idev
      END SUBROUTINE ELMOD1  
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Isym1(Isym) 
      DOUBLE PRECISION Isym(3,3,3,3)
      INTEGER I,J,K,L
      Isym=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Isym(i,j,k,l)=1.0d0/2.0d0*(delta(i,k)*delta(j,l)+
     1 delta(i,l)*delta(j,k))              
      Enddo
      Enddo
      Enddo
      Enddo     
      END SUBROUTINE Isym1      
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Iunit1(Iunit) 
      DOUBLE PRECISION Iunit(3,3,3,3)
      INTEGER I,J,K,L
      Iunit=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Iunit(i,j,k,l)=delta(i,k)*delta(j,l)           
      Enddo
      Enddo
      Enddo
      Enddo     
      END SUBROUTINE Iunit1                                
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Idev1(Idev) 
      DOUBLE PRECISION Idev(3,3,3,3), 
     1  Ivol(3,3,3,3), Isym(3,3,3,3)
      INTEGER I,J,K,L
      Idev=0.0d0
      Ivol=0.0d0
      Isym=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Isym(i,j,k,l)=1.0d0/2.0d0*(delta(i,k)*delta(j,l)+
     1 delta(i,l)*delta(j,k))              
      Enddo
      Enddo
      Enddo
      Enddo
      Ivol=0.0D0
      Ivol=1.0d0/3.0d0*delta.dyad.delta    
      Idev=Isym-Ivol         
      END SUBROUTINE Idev1             
c------------------------------------------------------ 
c------------------------------------------------------                
      SUBROUTINE Matrixtovector(vec1,mat1)      
      double precision vec1(6),  mat1(3,3)   
      vec1(1)=mat1(1,1)  
      vec1(2)=mat1(2,2) 
      vec1(3)=mat1(3,3) 
      vec1(4)=mat1(1,2)  
      vec1(5)=mat1(1,3) 
      vec1(6)=mat1(2,3)  
      return          
      END SUBROUTINE Matrixtovector
c------------------------------------------------------ 
c------------------------------------------------------       
      SUBROUTINE Vectortomatrix(vec1,mat1)      
      double precision vec1(6),  mat1(3,3)   
      mat1(1,1)=vec1(1)  
      mat1(2,2)=vec1(2)  
      mat1(3,3)=vec1(3)  
      mat1(1,2)=vec1(4)  
      mat1(2,1)=vec1(4)  
      mat1(1,3)=vec1(5)  
      mat1(3,1)=vec1(5)  
      mat1(2,3)=vec1(6)  
      mat1(3,2)=vec1(6)  
      return          
      END SUBROUTINE Vectortomatrix               
c------------------------------------------------------ 
c------------------------------------------------------ 
      subroutine tensortomatrix(a3333,  b66)   ! returns b(6,6)
      double precision a3333(3,3,3,3),b66(6,6)
      integer i,j,i9(6),j9(6)
      data i9/1,2,3,1,1,2/
     .     j9/1,2,3,2,3,3/
      do  i=1,6   !  switch to matrix notation
      do  j=1,6
      b66(i,j)=a3333(i9(i),j9(i),i9(j),j9(j))  
      enddo
      enddo   
      return
      end subroutine tensortomatrix
c------------------------------------------------------ 
c------------------------------------------------------ 
      subroutine F1(T, Fm)
      double precision T(3,3),hatT(3,3), hTd(3,3), NormhatTd, trhtd3,
     1 cos3lode, Tanpsi, Tanpsi2, term1, Fm, Fm2, tensor(3,3)
C     Relative stress
      hatT=T/Tr(T)
C     Relative deviator stress 
      hTd=hatT-1.0d0/3.0d0*delta  
C     Norm relative deviator stress 
      NormhatTd=norm(hTd)
C     trace of 2 times the dot product of hTd  
      trhtd3=Tr((hTd.dot.hTd).dot.hTd)
C     Lode's angle 
      if (NormhatTd.le.1.d-10) then
        cos3lode=1.0d0
      else
        cos3lode=(-sq6*trhtd3/(NormhatTd)**(3.d0))
        if (cos3lode.gt.1.0d0)  cos3lode =  1.0d0
        if (cos3lode.lt.-1.0d0) cos3lode = -1.0d0 
      endif
C     Psi angle  
      Tanpsi=sq3*NormhatTd
C     F function
C      Tanpsi=Tan(psi)
      Tanpsi2=Tanpsi**2.d0
      term1=Tanpsi2/8.+(2.-Tanpsi2)/(2.+sq2*Tanpsi*
     1 cos3lode)
      if (term1.lt.0) then
      Fm=1.d-10
      else
      Fm=sqrt(term1)-Tanpsi/(2.0d0*sq2) 
      endif
      end subroutine F1
c------------------------------------------------------ 
c------------------------------------------------------ 
      subroutine pq(T,p,q)   ! p and q
      double precision T(3,3),p,q, DevSTRESS(3,3),
     1 NormSt
C      mean stress 
      p=-1.0d0/3.0d0*Tr(T)
C     deviator stress
      DevSTRESS=T+p*delta
C     Norm of the deviator stress 
      NormSt=norm(DevSTRESS)
C     Trial q 
      q=sqrt(3.0d0/2.0d0)*NormSt    
      return
      end subroutine pq   
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Initial(STRESS,T, DSTRAN, DEPS, NTENS,NDI, NSHR)      
      double precision STRESS(ntens), T(3,3)
     1 ,DSTRAN(ntens), DEPS(3,3) 
      Integer ntens, nshr, ndi   
      DEPS=0.0D0
      T=0.0D0
C     
      do i=1,ndi
      T(i,i)=stress(i)
      DEPS(i,i)=DSTRAN(i)
      enddo 
C     
      if (nshr.ge.1) then
      T(1,2)=stress(4)
      T(2,1)=stress(4)    
      DEPS(1,2)=0.5d0*DSTRAN(4)
      DEPS(2,1)=0.5d0*DSTRAN(4)           
      endif
      if (nshr.ge.2) then
      T(1,3)=stress(5)
      T(3,1)=stress(5)   
      DEPS(1,3)=0.5d0*DSTRAN(5)
      DEPS(3,1)=0.5d0*DSTRAN(5)          
      endif
      if (nshr.ge.3) then
      T(2,3)=stress(6)
      T(3,2)=stress(6)    
      DEPS(2,3)=0.5d0*DSTRAN(6)
      DEPS(3,2)=0.5d0*DSTRAN(6)         
      endif   
      return          
      END SUBROUTINE Initial                
c------------------------------------------------------ 
c------------------------------------------------------ 
       SUBROUTINE Solution(NTENS, NDI, NSHR, T, STRESS, JAC, DDSDDE) 
       integer NTENS, NDI, NSHR, i, j, k, l
C     Subroutine for filling the stress and Jacobian matrix        
       double precision T(3,3), JAC(3,3,3,3), STRESS(NTENS),
     1 DDSDDE(NTENS,NTENS), JAC66(6,6)
      k=1
      l=1
c------------------------------------------------------
      do i=1,ndi
      stress(i)=T(i,i)
      enddo 
C     
      if (nshr.ge.1) then
      stress(ndi+1)=T(1,2)         
      endif
      if (nshr.ge.2) then
      stress(ndi+2)=T(1,3)         
      endif
      if (nshr.ge.3) then
      stress(ndi+3)=T(2,3)         
      endif   
      call tensortomatrix(jac,  jac66)   
        do i=1,ndi
        do j=1,ndi
          ddsdde(i,j)=jac66(i,j)
        enddo
      enddo 
      do i=ndi+1,ndi+nshr
        do j=1,ndi
          ddsdde(i,j)=jac66(3+k,j)
        enddo
        k=k+1
      enddo  
      do i=1,ndi
      l=1
        do j=ndi+1,ndi+nshr
          ddsdde(i,j)=jac66(i,3+l)
          l=l+1
        enddo
      enddo   
      k=1
      do i=ndi+1,ndi+nshr
        l=1
        do j=ndi+1,ndi+nshr
          ddsdde(i,j)=jac66(3+k,3+l)
          l=l+1
        enddo
        k=k+1
      enddo   
       Return       
       END SUBROUTINE Solution  
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION dyad33x33(a1,a2) RESULT (a3)
      DOUBLE PRECISION, INTENT(IN),dimension(3,3) :: a1, a2
      DOUBLE PRECISION ,dimension(3,3,3,3) :: a3
      INTEGER i,j,k,l
      a3=0.0d0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      a3(i,j,k,l)=a3(i,j,k,l)+a1(i,j)*a2(k,l)
      Enddo
      Enddo
      Enddo
      Enddo
      END FUNCTION dyad33x33 
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION double3333x33(a,b) RESULT (c) 
      DOUBLE PRECISION, INTENT(IN):: a(3,3,3,3),b(3,3)    
      double precision c(3,3)
      do i=1,3
      do j=1,3
      c(i,j)=a(i,j,1,1)*b(1,1)+
     1 a(i,j,1,2)*b(1,2)+
     2 a(i,j,1,3)*b(1,3)+
     3 a(i,j,2,1)*b(2,1)+
     4 a(i,j,2,2)*b(2,2)+
     5 a(i,j,2,3)*b(2,3)+
     6 a(i,j,3,1)*b(3,1)+
     7 a(i,j,3,2)*b(3,2)+
     8 a(i,j,3,3)*b(3,3)
      enddo
      enddo
      return
      END FUNCTION double3333x33 
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION double33x3333(b,a) RESULT (c) 
      DOUBLE PRECISION, INTENT(IN):: b(3,3),a(3,3,3,3)    
      double precision c(3,3), d(3,3)
      integer i1,i2
      do 10 i1=1,3
      do 10 i2=1,3
   10 d(i1,i2)       =  a(1,1,i1,i2)*b(1,1)+
     1  a(1,2,i1,i2)*b(1,2)+
     2  a(1,3,i1,i2)*b(1,3)+
     3  a(2,1,i1,i2)*b(2,1)+
     4  a(2,2,i1,i2)*b(2,2)+
     5  a(2,3,i1,i2)*b(2,3)+
     6  a(3,1,i1,i2)*b(3,1)+
     7  a(3,2,i1,i2)*b(3,2)+
     8  a(3,3,i1,i2)*b(3,3)
      do 20 i1=1,3
      do 20 i2=1,3
   20   c(i1,i2) = d(i1,i2)  ! d() prevents errors if: call x99x91(a,b,b)
      return     
c      do i=1,3
c      do j=1,3
c      do k=1,3
c      do l=1,3
c      c(k,l)=b(i,j)*a(j,i,k,l)
c      enddo
c      enddo     
c      enddo
c      enddo
      return
      END FUNCTION double33x3333       
c------------------------------------------------------ 
c------------------------------------------------------ 
      FUNCTION double33x33(a,b) RESULT (c) 
      DOUBLE PRECISION, INTENT(IN):: a(3,3),b(3,3)    
      double precision c
      c=0.0d0
      do i=1,3
      do j=1,3
      c=c+a(i,j)*b(i,j)
      enddo
      enddo
      return
      END FUNCTION double33x33 
c------------------------------------------------------ 
c------------------------------------------------------
      FUNCTION double3333x3333(a,b) RESULT (c) 
C     Double contraction fourth order tensors  
      integer i1,i2,j1,j2
      DOUBLE PRECISION, INTENT(IN):: a(3,3,3,3),b(3,3,3,3)      
      DOUBLE PRECISION c(3,3,3,3)
      c=0.0d0
      do 10 i1=1,3
      do 10 i2=1,3
      do 10 j1=1,3
      do 10 j2=1,3
  10  c(i1,i2,j1,j2) =c(i1,i2,j1,j2)+
     1                  a(i1,i2,1,1)*b(1,1,j1,j2)+
     2                  a(i1,i2,1,2)*b(1,2,j1,j2)+
     3                  a(i1,i2,1,3)*b(1,3,j1,j2)+
     4                  a(i1,i2,2,1)*b(2,1,j1,j2)+
     5                  a(i1,i2,2,2)*b(2,2,j1,j2)+
     6                  a(i1,i2,2,3)*b(2,3,j1,j2)+
     7                  a(i1,i2,3,1)*b(3,1,j1,j2)+
     8                  a(i1,i2,3,2)*b(3,2,j1,j2)+
     9                  a(i1,i2,3,3)*b(3,3,j1,j2)
      return
      END FUNCTION double3333x3333     
c------------------------------------------------------ 
c------------------------------------------------------
      FUNCTION norm33(a) RESULT (b) 
      DOUBLE PRECISION, INTENT(IN):: a(3,3)   
      double precision b      
      b=sqrt(a(1,1)*a(1,1)+a(1,2)*a(1,2)+a(1,3)*a(1,3)+
     1  a(2,1)*a(2,1)+a(2,2)*a(2,2)+a(2,3)*a(2,3)+
     2  a(3,1)*a(3,1)+a(3,2)*a(3,2)+a(3,3)*a(3,3)) 
      return
      END FUNCTION norm33      
c------------------------------------------------------ 
c------------------------------------------------------                  
      SUBROUTINE D1(DSTRAN, D, dtime, NDI, NSHR, NTENS) 
C     Strain rate tensor D
      integer i,j, NDI, NSHR, NTENS
      double precision D(3,3), DSTRAN(6), dtime
      if (dtime==0.0d0) then
      D=0.0D0
      else
      D=0.0D0
      Do i=1,ndi
      D(i,i)=dstran(i)/dtime ! covariant components, matrix format
      Enddo
      if (nshr.ge.1) then
      D(1,2)=dstran(4)/(2.0d0*dtime)
      D(2,1)=D(1,2)
      endif
      if (nshr.ge.2) then
      D(1,3)=dstran(5)/(2.0d0*dtime)
      D(3,1)=D(1,3)   
      endif
      if (nshr.ge.3) then
      D(2,3)=dstran(6)/(2.0d0*dtime)
      D(3,2)=D(2,3)  
      endif  
      endif    
      END SUBROUTINE D1
c------------------------------------------------------ 
c------------------------------------------------------
      FUNCTION inverse3333(a3333) RESULT (b3333) 
      DOUBLE PRECISION, INTENT(IN):: a3333(3,3,3,3)   
      DOUBLE PRECISION b3333(3,3,3,3)          
      DOUBLE PRECISION a(9,9),b(9,9),
     1 c(9,18),cc,dd
      integer i,j,k,inv,i9(9),j9(9)
      data i9/1,2,3,1,2,1,3,2,3/,
     .     j9/1,2,3,2,1,3,1,3,2/
      do 10 i=1,9   !  switch to matrix notation
      do 10 j=1,9
  10  a(i,j)=a3333(i9(i),j9(i),i9(j),j9(j))
c (1) Preparam. c(9,18)
      do i=1,9
      do j=1,9
      c(i,j+9)=0.0d0
      c(i,j)=a(i,j)
      enddo
      c(i,9+i)=1.0d0
      enddo
c (2) invert c with result  going to the right half of it
      do i=1,9
      cc = c(i,i)
      if (abs(cc).lt.1d-4) then
      b3333=0.0d0
      inv =0              ! inversion failed
      Print*, "UMAT: Inversion failed"
c      pause
      CALL XIT
      return
      endif
      c(i,i)=cc-1.0d0
      do k=i+1,18
      dd=c(i,k)/cc
      do j=1,9
      c(j,k)=c(j,k)-dd*c(j,i)
      enddo
      enddo
      enddo
c (3) copy result to b()
      do i=1,9
      do j=1,9
      b(i,j)= c(i,j+9)
      enddo
      enddo
      inv=1                    ! inversion successfull
      do 20 i=1,9     ! switch to tensorial notation
      do 20 j=1,9
  20  b3333(i9(i),j9(i),i9(j),j9(j))=b(i,j)
      return
      END FUNCTION inverse3333
c------------------------------------------------------ 
c------------------------------------------------------ 
c     
      FUNCTION inverse33(a33) RESULT (b33)  
C
      integer m, n, i, j
      PARAMETER (M=3,N=3)
      Double precision a33(3,3),b33(3,3), p
C
      DO 5 I=1,M
      DO 5 J=1,N
      b33(I,J)=b33(I,J)
5     CONTINUE
      DO 10 K=1,M
      P=b33(K,K)
      b33(K,K)=1.0d0
      DO 20 J=1,N
      b33(K,J)=b33(K,J)/P
20    CONTINUE
      DO 10 I=1,M
      IF(I .EQ. K) GO TO 10
      P=b33(I,K)
      b33(I,K)=0.0d0
      DO 30 J=1,N
      b33(I,J)=b33(I,J)-b33(K,J)*P
30    CONTINUE
10    CONTINUE
      return
      END FUNCTION inverse33    
c------------------------------------------------------ 
c------------------------------------------------------                                                 
      SUBROUTINE dpdTdqdT(T, dpdT, dqdT)   
      Double precision T(3,3),dpdT(3,3), dqdT(3,3), p, q,
     1 Tdevunit(3,3), Tdev(3,3), Idev(3,3,3,3)
c     Computes the derivatives dpdT dqdT
c     Negative for compression always      
      CALL pq(T,p,q)  
      dpdT=-1.0d0/3.0d0*delta
      call Idev1(Idev)
      Tdev=Idev.double.T
      if (abs(norm(Tdev))<1.0d-20) then
      dqdT=0.0d0
      else
      dqdT=sqrt(3.0d0/2.0d0)*Tdev/norm(Tdev)    
      endif
      END SUBROUTINE dpdTdqdT  
c------------------------------------------------------ 
c------------------------------------------------------                                                 
      SUBROUTINE Lode(T, theta)   
      Double precision T(3,3), p, q,
     1  Tdev(3,3), Idev(3,3,3,3), theta
c     Computes the derivatives dpdT dqdT
c     Negative for compression always      
      CALL pq(T,p,q)  
      call Idev1(Idev)
      Tdev=Idev.double.T
      if (norm(Tdev)<=1.0e-20) then
      theta=0.0d0
      else
      theta=-1.0d0/3.0d0*asin(3.0d0*sq3/2.0d0*det(Tdev)
     1 /(q/sq3)**3.0d0)
      endif
      END SUBROUTINE Lode        
c------------------------------------------------------   
      SUBROUTINE CEP1(Elmod, dFdT, dgdT, H, CEP)   
      Double precision Elmod(3,3,3,3), dfdT(3,3), dgdT(3,3)
     1 ,H, CEP(3,3,3,3), chi, tensor1(3,3), tensor2(3,3)
c     Computes the elastoplastic tangent moduli
      chi=((dfdT.double.Elmod).double.dgdT)-H
      tensor1=ELMOD.double.dgdT
      tensor2=dfdT.double.ELMOD
      if (norm(tensor1)<1.0d-10) then
      CEP=elmod
      return
      endif
      CEP=ELMOD-(tensor1.dyad.tensor2)
      CEP=CEP/chi 
      return
      END SUBROUTINE CEP1                                                  
      END MODULE tensors  
c------------------------------------------------------ 
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      
CXXXX     2. UMAT CAM CLAY WITH STRUCTURE- (LIU & CARTER), 2005     XXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      
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
