!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sum L' 
! output: L'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine drv(nobs,nvars,rs,nrs,x,xb,vl)
      IMPLICIT NONE
      integer :: nobs
      integer :: nvars
      integer :: nrs
      integer :: rs(nrs)
      integer :: j
      integer :: j1
      integer :: j2
      integer :: k
      double precision :: x(nobs,nvars)
      double precision :: xb(nobs)                                       
      double precision :: vl(nvars)                                              
      double precision :: se                                           
      double precision :: sxe(nvars)                                                                                                                                                                     
      vl=0.0
      se = 0.0
      sxe = 0.0
    do k=nrs,1,-1                                              
        j2=nobs
        if(k<nrs)j2=rs(k+1)-1                                                             
        j1=rs(k)                                                  
        do j=j2,j1,-1                                             
            sxe=sxe+x(j,:)*xb(j)
            se=se+xb(j)                                                  
        enddo                                                        
        vl=vl+(x(rs(k),:)-sxe/se)/nobs
    enddo                                                     
      end                                                  





