!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute d for tied failure time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                    
      subroutine failure(no,n_s,i_s,idx,irs)
      IMPLICIT NONE
      integer :: no
      integer :: n_s
      integer :: k                               
      integer :: i_s(n_s)
      integer :: idx(no)
      integer :: irs(n_s)                                         
      irs(1)=idx(1)
      do k=2,n_s                                                 
      irs(k)=idx(i_s(k-1)+1)
      enddo                                                        
      return                                                          
      end