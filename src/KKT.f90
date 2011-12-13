      subroutine KKT(no,ni,x,y,d,a,nlam,flog)
      IMPLICIT NONE    
      integer :: no
      integer :: ni
      integer :: nlam
      integer :: jerr 
      double precision :: x(no,ni)
      double precision :: y(no)
      double precision :: d(no)
      double precision :: a(ni,nlam)
      double precision :: flog(ni,nlam)
!     working variable
      integer :: i
      integer :: j
      integer :: n_s
      integer :: lam
      double precision :: fmax
      double precision :: t0   
      double precision :: f(no,nlam)
      double precision :: e(no,nlam)
      double precision :: r(ni)
      double precision :: xm
      integer :: idx(no)
      integer :: i_s(no)
      integer :: irs(no)
      r = 0.0
      do j=1,ni                                                 
          xm=sum(x(:,j))/no                                        
          x(:,j)=x(:,j)-xm                                           
      enddo                                            
      f=0.0
      fmax=log(huge(f(1,1))*0.1)
      call riskidx(no,y,d,n_s,i_s,idx,t0,jerr)
      call failure(no,n_s,i_s,idx,irs)
      if(jerr /= 0) return 
      f=matmul(x,a)
      e=exp(sign(min(abs(f),fmax),f))                    
      do lam=1,nlam
          call derivative(no,ni,n_s,i_s,idx,x,e(:,lam),irs,r)
          flog(:,lam) = r                                                         
      enddo
      return
      end                                                        
       