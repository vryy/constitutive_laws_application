c----------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      real*8	statev(10),dummy
      integer	nstatv,noel,npt,recnum


      statev(1)  = 1.1d0    ! ... e void ratio
      statev(2)  = 0.0d0    ! ... plastic strains component 1 
      statev(3)  = 0.0d0    ! ... plastic strains component 2
      statev(4)  = 0.0d0    ! ... 
      statev(5)  = 0.0d0    ! ... 
      statev(6)  = 0.0d0    ! ... plastic strains component 5
      statev(7)  = 0.0d0    ! ... plastic strains component 6 
      statev(8)  = 0.0d0    ! ... overconsolidated mean stress pc 
      statev(9)  = 0.0d0    ! ... additional voids ratio by soil structure de  
      statev(10)  = 0.0d0    ! ... yielding stress py 
      end
