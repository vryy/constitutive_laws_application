! Copyright (C)  2007  Andrzej Niemunis
!
! usovle is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! usolve is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
!  USA.

       MODULE unsymmetric_module                                        ! \com written by Andrzej Niemunis, may 2006
       IMPLICIT NONE

       interface outerprod
         module procedure outerproda,outerprodb
       end interface


       private
       public  matrix , set_to_identity, usolver, outerprod, 
     &        map2stran, map2D, map2stress, map2T,i6,j6
     
     
       integer, parameter,dimension(1:6) :: i6=(/ 1,2,3,1,1,2/),         !  \com aux. tables for transformations tensor <--> matrix
     &                                     j6=(/ 1,2,3,2,3,3/)        
     
       CONTAINS
       
       
        function map2stran(a,ntens)                                       ! \com ..........{\large MAP2STRAN}
        implicit none                                                   ! \com converts D(3,3)  to stran(6) with $\gamma_{12} = 2 \epsilon_{12}$ etc.
        real(8), intent(in), dimension(1:3,1:3) :: a
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens) :: map2stran
        integer :: i
        map2stran(1)=a(1,1)
        map2stran(2)=a(2,2)
        map2stran(3)=a(3,3)
        do i=4,ntens
        map2stran(i) = a(i6(i),j6(i)) * 2.0d0                           ! \com  abaqus needs gammas  as strains
        enddo                                                           ! \com   i6=(/ 1,2,3, 1,1,2/), j6=(/ 1,2,3, 2,3,3/)
      end function map2stran

      function map2D(a,ntens)                                           ! \com ..........{\large MAP2D}
        implicit none                                                   ! \com convert strain rate from vector dstran(1:ntens) to  D(3,3)
        real(8),  dimension(1:3,1:3) :: map2D
        integer, intent(in) :: ntens
        real(8), intent(in), dimension(:) :: a
        integer :: i
        map2D=0
        map2D(1,1) = a(1)
        map2D(2,2) = a(2)
        map2D(3,3) = a(3)
        do i=4,ntens
         map2D(i6(i),j6(i))=a(i)/2.0d0                                   ! \com  abaqus uses  gammas  for shear strains
         map2D(j6(i),i6(i))=a(i)/2.0d0                                   ! \com   i6=(/ 1,2,3, 1,1,2/), j6=(/ 1,2,3, 2,3,3/)
        enddo
      end function map2D

      function map2stress(a,ntens)                                      ! \com ..........{\large MAP2STRESS}
        implicit none                                                   ! \com convert tensor T(3,3)  to matrix stress(ntens)
        real(8), intent(in), dimension(1:3,1:3) :: a
        integer, intent(in) :: ntens
        real(8),  dimension(1:ntens) :: map2stress
        integer :: i
        do i=1,ntens
          map2stress(i) = a(i6(i),j6(i))                                ! \com  abaqus stress
        enddo                                                           ! \com   i6=(/ 1,2,3, 1,1,2/), j6=(/ 1,2,3, 2,3,3/)
      end function map2stress

      function map2T(a,ntens)                                           ! \com ..........{\large MAP2T}
        implicit none                                                   ! \com  convert  matrix stress(1:ntens)  to tensor T(3,3)
        real(8),  dimension(1:3,1:3) :: map2T
        integer, intent(in) :: ntens
        real(8), intent(in), dimension(:) :: a
        integer :: i
        map2T=0
        map2T(1,1) = a(1)
        map2T(2,2) = a(2)
        map2T(3,3) = a(3)
        do i=4,ntens
        map2T(i6(i),j6(i))=a(i)
        map2T(j6(i),i6(i))=a(i)
        enddo
      end function map2T


       SUBROUTINE getHausholder( ab, rho, gamma)                        !  \com  Kielbasinski s. 119
!  \com  being given a vector ab(1:n) finds  vb(1:n) such that
!   \com       $ H . ab  = \{rho, 0, 0, 0, ...\}$
!  \com  where: $H = Identity - vb .out. vb^T / gamma $  is an orthogonal matrix  n#n
!  \com        $rho = - sign(ab(1)) * sqrt(ab^T . ab)$
!  \com      $ vb =  {1,0,0...0 } - ab / rho  $
!   \com     $  gamma =  (vb^T . vb )/2 $            (= vb(1) but may be too tricky)
!  \com   Notice that ab is overwritten by the values of vb
         implicit none
         INTEGER   :: n  !  \com , INTENT(IN)        :: n
         real(8), dimension(:), INTENT(INOUT)   :: ab
         real(8), INTENT(OUT)  :: rho, gamma
         real(8)   ::  aa, a1a1
         n = size(ab)
         IF(n <= 0) RETURN
         IF(n == 1) then                                                !  \com onedimensional input vector a
           gamma = 1.0d0
           rho = ab(1)
           ab(1) = 0.0d0                                                !  \com = vb
           RETURN
         endif
         a1a1 = dot_product(ab(2:n),ab(2:n))
         IF(a1a1 <= tiny(a1a1)) then                                    !   \com input vector a in the desired form
           rho = ab(1)                                                  !  \com possibly $\rho = 0 $
           ab = 0.0d0                                                   !  \com = vb
           gamma = 1.0d0
           RETURN
         endif
         aa= ab(1)*ab(1)  + a1a1
         rho = -sign(sqrt(aa) , ab(1))
         ab(1) = 1.0d0 - ab(1)/rho                                      !  \com  = vb(1)
         ab(2:n) = -ab(2:n)/rho                                         !  \com = vb(2:n)
         gamma=0.5d0*dot_product(ab,ab)
      END SUBROUTINE  getHausholder


      SUBROUTINE app_HausholderLeft(vb, gamma, A)                 !  \com Kielbasinski s. 119
         implicit none
         INTEGER :: m,n
         real(8), dimension(:), INTENT(IN)   :: vb
         real(8),  INTENT(IN)   :: gamma
         real(8), dimension(:,:), INTENT(INOUT)   :: A
         integer :: icol
         real(8) :: beta
         m = size(A,1)
         n = size(A,2)
         if(size(vb) /= m)  stop '  size(vb) /= m  '
         IF(n <= 0  .or. m <= 0) RETURN
         IF(m == 1)   RETURN                                            !  \com nothing to do  because vb=0  for  n=1
         do icol = 1,n
           beta = dot_product(vb(1:m), A(1:m,icol))/gamma
           A(1:m,icol)  =   A(1:m,icol) - vb(1:m) * beta
         enddo
       END SUBROUTINE  app_HausholderLeft

      SUBROUTINE app_HausholderRight( vb, gamma, A)                !  \com Kielbasinski s. 119
         implicit none
         INTEGER  :: m,n
         real(8), dimension(:), INTENT(IN)   :: vb
         real(8),  INTENT(IN)   :: gamma
         real(8), dimension(:,:), INTENT(INOUT)   :: A
         integer :: irow
         real(8) :: beta
         m = size(A,1)
         n = size(A,2)
         if(size(vb) /= n)  stop '  size(vb) /= n '
         IF(n <= 0  .or. m <= 0) RETURN
         IF(n == 1)  RETURN                                             !  \com nothing to do  because vb=0  for  n=1
         do irow = 1,m
           beta = dot_product(A(irow,1:n),vb(1:n))/gamma
            A(irow, 1:n ) =   A(irow,1:n) - vb(1:n) * beta
         enddo
       END SUBROUTINE  app_HausholderRight


       function  getH( v, gamma)
        implicit none
       real(8) , dimension(:) ,intent(in) :: v
       real(8)  ,intent(in) :: gamma
       real(8) , dimension(size(v), size(v)) :: getH
       integer :: i
       getH = 0.0d0
       do i= 1,size(v)
       getH(i,i) = 1.0d0
       enddo
       getH = getH - outerprod(v,v)/gamma
       end function  getH

      FUNCTION outerproda(a,b)                                           !  \com outer product of two vectors
       implicit none
       REAL(8), DIMENSION(:), INTENT(IN) :: a,b
       REAL(8), DIMENSION(size(a),size(b)) :: outerproda
       outerproda  =  spread(a,2,size(b))*spread(b,1,size(a))
      END FUNCTION outerproda

      FUNCTION outerprodb(a,b)                                          !  \com outer product of two vectors
       implicit none
       complex(8), DIMENSION(:), INTENT(IN) :: a,b
       complex(8), DIMENSION(size(a),size(b)) :: outerprodb
       outerprodb  =  spread(a,2,size(b))*spread(b,1,size(a))
      END FUNCTION outerprodb



      SUBROUTINE  app_jacobian_similarity(A, p,q, c, s)                 !    \com jacobian similarity tranformation of a square symmetric matrix A
       implicit none
      INTEGER, INTENT(IN)        :: p,q                                        !  \com ( $ A : =  G^T .A . G $ with    Givens   rotation  G = {{c,s},{-s,c}}  )
      real(8), INTENT(IN)      :: c ,s                                         !  \com  Kielbasinski  385
      real(8), dimension(:,:),intent(inout) :: A
      real(8), dimension(:), allocatable :: prow ,qrow
      real(8) :: App, Apq, Aqq
      integer :: m, n
      m = size(A,1)
      n = size(A,2)
      if(min(m,n) < 2) stop 'error: jacobian_similarity  A is too small'
      if(p == q)  stop 'error: jacobian_similarity  p == q'
      if(m /= n ) stop 'error: jacobian_similarity  A is not square'
      if(p<1 .or. p>n) stop 'error: jacobian_similarity p out of range'
      if(q<1 .or. q>n) stop 'error: jacobian_similarity q out of range'
      allocate(prow(n) ,qrow(n) )
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
      deallocate(prow ,qrow)
      END SUBROUTINE  app_jacobian_similarity


      SUBROUTINE  get_jacobian_rot(A, p,q, c, s)             !   \com returns jacobian similarity  tranformation param.
       implicit none
      INTEGER, INTENT(OUT)        :: p,q                      !  \com  for iterative diagonalization of  a square symm.  A
      real(8), INTENT(OUT)      :: c ,s                      !  \com Kielbasinski 385-386
      real(8), dimension(:,:),intent(in) :: A
      real(8) :: App, Apq, Aqq, d, t,maxoff
      integer :: m,n,i,j
      m = size(A,1)
      n = size(A,2)
      if(m /= n ) stop 'error: get_jacobian_tridiag  A is not square'
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

      subroutine set_to_identity(A)
       implicit none
      real(8) :: A(:,:)
      integer :: i,m,n
      m=size(A,1)
      n=size(A,2)
      if(m /= n) stop 'error UNIT:  A is not square'
      A = 0.0d0
      do i=1,n
      A(i,i) = 1.0d0
      enddo
      end  subroutine  set_to_identity

      SUBROUTINE  app_GivensLeft(A, p,q, c, s)              !  \com left  Givens rotation  of   square symmetric matrix A
       implicit none
      INTEGER, INTENT(IN)        :: p,q                     !  \com $A := G. A $ ,$ G = {{c,s},{-s,c}}$  Kielbasinski  p. 121
      real(8), INTENT(IN)      :: c ,s
      real(8), dimension(:,:),intent(inout) :: A
      real(8), dimension(:), allocatable :: prow ,qrow
      integer ::  n,m
      m = size(A,1)
      n = size(A,2)
      if(min(m,n) < 2) stop 'error:  app_Givens  A is too small'
      if(p == q)  stop 'error:  app_Givens  p == q'
      if(size(A,1) /= n ) stop 'error:  app_Givens  A is not square'
      if(p<1 .or. p>n) stop 'error:  app_Givens p out of range'
      if(q<1 .or. q>n) stop 'error:  app_Givens q out of range'
      allocate(prow(n) ,qrow(n) )
      prow(1:n) = c*A(p,1:n) -  s*A(q,1:n)
      qrow(1:n) =  s*A(p,1:n) + c*A(q,1:n)
      A(p,1:n) =   prow(1:n)
      A(q,1:n) =   qrow(1:n)
      deallocate(prow ,qrow)
      END SUBROUTINE  app_GivensLeft

      function inv_bidiagonal(bid, ok)               !  \com inverts a bidiagonal matrix  bid(n,2)
       implicit none
      logical, intent(out) :: ok                                !  \com returns square n#n  inverse matrix
!  \com this algorithm is based on the analytical solution of
!   \com Inverse[  {{a1,b1,0,0,0,0}, {0,a2,b2,0,0,0},{0,0,a3,b3,0,0},{0,0,0,a4,b4,0},{0,0,0,0,a5,b5}, {0,0,0,0,0,a6}}  ]
      real(8),dimension(:,:), intent(in) :: bid
      real(8),dimension(1:size(bid,1),1:size(bid,1) )::Inv_bidiagonal
      real(8),dimension(:), allocatable :: eta
      integer(4) :: n,i,j
      ok = .false.
      n = size(bid,1)
      allocate(eta(1:n))
      if( minval(  abs(bid(:,1))    ) <= tiny(1.0d0)   )  return
      ok = .true.
      if (n==1) then
       inv_bidiagonal =  1.0d0/bid
       return
      endif
      eta = bid(1:n,2)/ bid(1:n,1)
      inv_bidiagonal = 0.0d0
      do i=1,n
       inv_bidiagonal(i,i) = 1.0d0/bid(i,1)
       do j=i+1,n
         inv_bidiagonal(i,j) =  product( eta(i:j-1))/bid(j,1) *
     &                        sign(1.0d0,-2.0d0*mod(i+j,2)+1.0d0)     !  \com corrected Niemunis mai 2006
       enddo
      enddo
      end function inv_bidiagonal


       SUBROUTINE  matrix(task, x, n, ok )                              !  \com improve conditioning number or inverse (or both)
  !  use nag_mat_inv                                                  !   \com professional  inversion
       implicit none
       integer, intent(in) :: n
       character(7), intent(in) :: task                                  !  \com ='improve' or 'inverse'  or 'singula'
       real(8), INTENT(IN OUT)  :: x(n,n)
       LOGICAL, intent(out) ::   ok

       INTEGER ::  i, j,m, L,iter
       real(8) ::   cosine, sine
       real(8),dimension(:,:),allocatable:: U, Vt, HL, HR,GT,
     &                                       bidiagonal,gamma
      real(8),dimension( : ) , allocatable::  vb
      integer ,dimension( : ) , allocatable::  jarray
      real(8) ::   max_diag, min_diag, condition_number, thres ,
     &                max_off , av_diag
       m = n                                                            !  \com square matrix

  !     if(task == 'inverse')  then        !  \com iiiiiiiiiiiiiiiii inverse iiiiiiiiiiiiiiiiiiiiii
  !       call nag_gen_mat_inv(x(1:n,1:n))                               ! \com  professional inversion
  !       ok = .true.                                                    ! \com  professional inversion
  !        return                                                        ! \com professional inversion
  !    endif


      allocate( bidiagonal(1:n, 1:2), gamma(1:n, 1:2),
     &    U (1:m,1:m),  Vt (1:n,1:n), HL(1:m,1:m), HR(1:n,1:n) ,
     &     vb(m) ,GT(n,n) ,  jarray(1:n) )

      bidiagonal = 0.0d0
      gamma = 1.0d0

      do L = 1,n-1
      call getHausholder( x(L:m, L), bidiagonal(L,1) , gamma(L,1))      ! \com  poddiag czesc L-tej kolumny  zamieniona na  v dla 1- v .out. v /gamma
      call app_HausholderLeft( x(L:n, L), gamma(L,1), x(L:m, L+1:n))
      call getHausholder( x(L,L+1:n),bidiagonal(L,2),gamma(L,2))
      call app_HausholderRight( x(L,L+1:n),gamma(L,2), x(L+1:m,L+1:n))
      enddo




       bidiagonal(n,1:2)  = (/ x(n,n), 0.0d0 /)                         !  \com last line
       gamma(n,1:2)  = 1.0d0

      call set_to_identity(U)
      call set_to_identity(Vt)

        do L = n-1 , 1, -1                                              !  \com collect    left   reflections  into  Uout
           vb = 0
           vb(L:m) = x(L:m,L)
           HL = getH( vb, gamma(L,1))
           U  = matmul(HL,U )
        enddo

        do L = n-2 , 1, -1                                              ! \com collect   right   reflections  into  Vtout
            vb = 0
            vb(L+1:m) = x(L,L+1:m)
            HR = getH( vb, gamma(L,2))
            Vt= matmul(Vt,HR)
        enddo


      if(task == 'inverse')  then        !  \com iiiiiiiiiiiiiiiii inverse iiiiiiiiiiiiiiiiiiiiii
         x = inv_bidiagonal(bidiagonal, ok)
         x = matmul(matmul(transpose(Vt) , x), transpose(U))
         return
      endif                           ! \com iiiiiiiiiiiiiiiii end  inverse iiiiiiiiiiiiiiiiiiiiii

      if(task == 'improve')  then    ! \com mmmmmmmmmmm  improve mmmmmmm
          ok = .true.
          if( maxval(abs( bidiagonal(1:n,1)) ) <= tiny(1.0d0) )
     &    stop 'error: improve: all diagonal elements are zero'
          do i=1,n
              min_diag  =  minval(abs( bidiagonal(1:n,1)) )
              max_diag  =  maxval(abs( bidiagonal(1:n,1)) )
              condition_number = min_diag/max_diag
              thres = 1.0d-5 * max_diag
              if (condition_number > 1.0d-5) then
                ok=.true.
                exit
              else
                jarray = minloc(abs( bidiagonal(1:n,1)) )
                j = jarray(1)                                           !  \com position where
                bidiagonal(j,1) = sign( 1.0d0, bidiagonal(j,1) ) *thres
                cycle
              endif
          enddo
          x = 0.0d0                                                     ! \com restore the original martix
          do i=1,n-1
               x(i,i:i+1) = bidiagonal(i,1:2)
          enddo
          x(n,n) =  bidiagonal(n,1)
          x = matmul(matmul(U, x),Vt)
          return
      endif           ! \com mmmmmmmmmmmmmmmmmmmmm  end improve mmmmmmmmmmmmmmmm
      if (task=='singula') then   ! \com sssssssssssssssssssssss  singular values sssssssssssssss
          x = 0.0d0
          do i=1,n-1
               x(i,i:i+1) = bidiagonal(i,1:2)
          enddo
          x(n,n) =  bidiagonal(n,1)
          x = matmul(transpose(x) , x)  ! \com construct  a tridiagonal matrix $B^T . B $
          do  iter = 1,30
            call  get_jacobian_rot(x, i ,j, cosine, sine)
            call  app_jacobian_similarity(x, i,j, cosine, sine)
          enddo
          max_off = 0.0d0
          av_diag = 0.0d0
          do i=1,n-1
              max_off = max( max_off, maxval( abs( x(i,i+1:n) )  ))
              av_diag = av_diag + abs( x(i,i) )
          enddo
          av_diag = (av_diag + abs(  x(n,n) ) )/n
          ok = .true.
          if(max_off > av_diag * 1.0d-10 ) ok = .false.
          do i=1,n
                x(i,i) = sqrt( x(i,i) )
          enddo
          return
      endif                 ! \com  sssssssssssssssssssssss  end singular values sssssssssssssss

      END SUBROUTINE matrix


       subroutine USOLVER(KK,u,rhs,is,ntens)
!  \com    KK - stiffness  is spoiled  within the subroutine
!  \com    u - strain   rhs - stress
!  \com    is(i)= 1 means rhs(i) is prescribed,
!  \com    is(i)= 0  means u(i) is prescribed

      implicit none
      integer :: ntens , is(ntens),i,j,ii,jj, nis,nes
      logical :: ok
      real(8), dimension(1:ntens,1:ntens), intent(inout):: KK
      real(8), dimension(1:ntens), intent(inout)::  u,rhs
      real(8), dimension(1:ntens):: rhs1
      real(8), allocatable :: rhsPrim(:), KKprim(:,:), uprim(:)

      nis = 0                                                           ! \com count nonzero is
      do i=1,ntens
      if(is(i) == 1) nis=nis+1                                          ! \com number of prescribed stress components
      enddo
      nes = ntens - nis                                                 ! \com  number of prescribed strain components
      if (nis == 0) then                                                ! \com a special case with full strain control
      rhs =  matmul(KK,u)
      return
      endif
      if (nes == 0) then                                                ! \com a special case with full stress control
      allocate( KKprim(nis,nis) )
      KKprim = KK
      call matrix('inverse', KKprim, ntens, ok )
      u =  matmul(KKprim,rhs)
      deallocate(KKprim)
      return
      endif
      ! \com modify the rhs  to rhs1
      rhs1 = rhs
      do i=1,ntens
      if (is(i) == 0) then                                              ! \com  modify rhs wherever strain control
        rhs1 = rhs1 - u(i)*KK(:,i)
      endif
      enddo
      allocate(KKprim(nis,nis), rhsprim(nis), uprim(nis) )              ! \com re-dimension  stiffness and rhs
      ii = 0
      do i=1,ntens                                                      ! \com  over rows of the original KK
        if (is(i) == 1) then                                            ! \com  the i-th row is qualified    (must be inverted)
           ii=ii+1
           rhsPrim(ii) = rhs1(i)
           jj=0
           do j=1,ntens                                                 !  \com over elements of each qualified row
              if (is(j) == 1) then
                 jj = jj+1
                 KKprim(ii,jj) = KK(i,j)
              endif
           enddo     ! \com  j
        endif        ! \com  is(i)==0
      enddo          ! \com  i
      if (nis == 1) then
        uprim = rhsprim / KKprim(1,1)
      else
        call matrix('inverse', KKprim, nis, ok )
        uprim = matmul(KKprim,rhsprim)
      endif
      ii = 0
      do i=1,ntens
        if ( is(i) == 1 ) then                                          ! \com put uprim into the u - vector
           ii=ii+1
           u(i) = uprim(ii)
        endif
      enddo
      do i=1,ntens
        if ( is(i) == 0 ) then
           rhs(i) = dot_product( KK(i,:), u)                            ! \com   calculate rhs where u prescribed
        endif
      enddo
      !  rhs1 = rhs -  matmul(KK,u)                                     !  \com check the iversion ( except trivial cases)
      deallocate(KKprim,rhsprim,uprim)
      return
      end subroutine USOLVER

      END MODULE unsymmetric_module




