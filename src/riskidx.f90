!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pre-process survival data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine riskidx(no,y,d,n_s,i_s,idx,t0,jerr)
      IMPLICIT NONE
      integer :: no
      integer :: n_s
      integer :: nidx
      integer :: j0
      integer :: j
      integer :: jerr
      integer :: idx(no)                                    
      integer :: i_s(no)
      double precision :: t0                        
      double precision :: yk                        
      double precision :: y(no)                          
      double precision :: d(no)                                        
15210 do 15211 j=1,no                                                 
      idx(j)=j                                                         
15211 continue                                                        
15212 continue
!     sort vector y and store the order in idx                                                        
      call quicksort(y,idx,1,no)        
!     CHECK if all d==0; j is the first index where d = 1
      j=1                                                             
15250 continue                                                        
15251 if(d(idx(j))>0.0)goto 15252                                   
      j=j+1                                                           
      if(j > no)goto 15252                                           
      goto 15251                                                      
15252 continue
      if(j < no-1)goto 15271                                       
      jerr=-5                                                         
      return                                                          
15271 continue                                                        
      j0=j-1       ! # of d == 0 before d[j]                                                 
      nidx=no-j0                                                        
15280 do 15281 j=1,nidx                                                 
      idx(j)=idx(j+j0)   ! exclude first j - 1 elements (whose d == 1 ) from the vector.                                                 
15281 continue                                                        
15282 continue
      jerr=0                                                          
      n_s=0        ! initailize n_s                                                    
      t0=y(idx(1))                                                     
      yk=t0                                                           
      j=2                                                             
15290 continue                                                     
15291 continue                                                        
15300 continue
15301 if(d(idx(j))>0.0.and.y(idx(j))>yk)goto 15302                
      j=j+1                      ! find the index j uch that d[2]>0 and y[2] > y[1] (no tie)                                         
      if(j>nidx)goto 15302                                           
      goto 15301                                                      
15302 continue       
      n_s=n_s+1                                                         
      i_s(n_s)=j-1                                                      
      if(j>nidx)goto 15292                                           
      if(j  /=  nidx)goto 15321                                         
      n_s=n_s+1                                                         
      i_s(n_s)=nidx                                                       
      goto 15292                                                      
15321 continue                                                        
      yk=y(idx(j))                                                     
      j=j+1                                                           
      goto 15291                 
15292 continue                                         
      return
      end