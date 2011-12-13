!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sum loglik
! output: loglik
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine objective(no,ni,n_s,i_s,idx,x,f,e,irs,r)
      IMPLICIT NONE
      integer :: no
      integer :: ni
      integer :: l
      integer :: n_s
      integer :: j
      integer :: j1
      integer :: j2
      integer :: k
      double precision :: x(no,ni)
      double precision :: f(no)                                       
      double precision :: e(no)                                       
      double precision :: r                                          
      double precision :: se                                           
      integer :: i_s(n_s)
      integer :: idx(no)                                         
      integer :: irs(n_s)                                                                                     
	  r=0.0
	  se = 0.0
15420 do 15421 k=n_s,1,-1                                              
      j2=i_s(k)                                                        
      j1=1                                                             
      if(k>1) j1=i_s(k-1)+1                                         
15430 do 15431 j=j2,j1,-1                                             
      se=se+e(idx(j))                                                  
15431 continue                                                        
      r=r+(-f(irs(k))+log(se))/no
15421 continue                                                        
15422 continue
      return                                                          
      end                                                  





