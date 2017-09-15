!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sum L' 
! output: L'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine derivative(no,ni,n_s,i_s,idx,x,e,irs,r)
      IMPLICIT NONE
      integer :: no
      integer :: ni
      integer :: n_s
      integer :: j
      integer :: j1
      integer :: j2
      integer :: k
      double precision :: x(no,ni)
      double precision :: e(no)                                       
      double precision :: r(ni)                                              
      double precision :: se                                           
      double precision :: sxe(ni)                                                                                      
      integer :: i_s(n_s)
      integer :: idx(no)                                         
      integer :: irs(n_s)                                                                                     
      r=0.0
      se = 0.0
      sxe = 0.0
      do 15421 k=n_s,1,-1                                              
      j2=i_s(k)                                                        
      j1=1                                                             
      if(k>1) j1=i_s(k-1)+1                                         
      do 15431 j=j2,j1,-1                                             
      sxe=sxe+x(idx(j),:)*e(idx(j))
      se=se+e(idx(j))                                                  
15431 continue                                                        
      r=r+(x(irs(k),:)-sxe/se)/no
15421 continue                                                        
      continue
      return                                                          
      end                                                  





