! --------------------------------------------------
SUBROUTINE chkvars(nobs,nvars,x,ju)                
! --------------------------------------------------
	IMPLICIT NONE
! - - - arg types - - -
	INTEGER ::  nobs,nvars,ju(nvars)
    DOUBLE PRECISION ::  x(nobs,nvars)       
! - - - local declarations - - -
	INTEGER :: i,j
	DOUBLE PRECISION :: t
! - - - begin - - -
    DO j=1,nvars                                                       
    	ju(j)=0                                                               
    	t=x(1,j)                                                              
 		DO i=2,nobs                                                       
      		IF (x(i,j)/=t) THEN
				ju(j)=1
				EXIT
			ENDIF
        ENDDO
	ENDDO                                                              
END SUBROUTINE chkvars
