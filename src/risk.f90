!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sum L' 
! output: L'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine risk(nobs,rs,nrs,x,xb,fjj)
      IMPLICIT NONE
      integer :: nobs
      INTEGER::nrs
      INTEGER::rs(nrs)
      integer :: j
      integer :: j1
      integer :: j2
      integer :: k
      double precision :: x(nobs)
      double precision :: xb(nobs)                                       
      double precision :: fjj                                            
      double precision :: se                                           
      double precision :: sxe
      fjj=0.0
      se = 0.0
      sxe = 0.0                                            
      do k=nrs,1,-1                                              
          j2=nobs
          if(k<nrs)j2=rs(k+1)-1                                                             
          j1=rs(k)                                                  
          do j=j2,j1,-1                                             
              sxe=sxe+x(j)*xb(j)
              se=se+xb(j)                                                  
          enddo                                                        
          fjj=fjj+(x(rs(k))-sxe/se)/nobs
      enddo
      return                                                          
      end                                                  





