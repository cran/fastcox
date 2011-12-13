! --------------------------------------------------
SUBROUTINE standardcox(nobs,nvars,x,ju,isd,xmean,xnorm)     
! --------------------------------------------------
	IMPLICIT NONE
! - - - arg types - - -
	INTEGER ::  nobs,nvars,isd,ju(nvars)
    DOUBLE PRECISION ::  x(nobs,nvars),xmean(nvars),xnorm(nvars)
! - - - local declarations - - -
	INTEGER :: j                                  
! - - - begin - - -                                                        
	DO j=1,nvars                                                       
    	IF(ju(j)==1) THEN                                        
    		xmean(j)=sum(x(:,j))/nobs     !mean                                       
    		x(:,j)=x(:,j)-xmean(j)                                                                                      
    		xnorm(j)=sqrt(dot_product(x(:,j),x(:,j))/nobs)    !standard deviation                        
      		IF(isd==1) x(:,j)=x(:,j)/xnorm(j)                                                                                               
		ENDIF                                                             
	ENDDO                  
END SUBROUTINE standardcox
