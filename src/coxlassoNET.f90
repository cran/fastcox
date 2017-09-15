! ------------------------------------------------------------------------------
! coxlassoNET.f90: the cocktail algorithm for elastic net penalized Cox's model.
! ------------------------------------------------------------------------------
! 
! USAGE:
! 
! coxlassoNET(rs,nrs,alpha,nobs,nvars,x,jd,pf,dfmax,pmax,nlam,flmin,ulam,&
!                eps,isd,maxit,nalam,beta,ibeta,nbeta,alam,npass,jerr)
! 
! INPUT ARGUMENTS:
!    
!    rs = the index of failure times, whose length is nrs
!    nrs = the number of risk sets
!    alpha = the elasticnet mixing parameter, with (0 <= alpha <= 1)
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = ordered matrix of predictors, of dimension N * p; each row is 
!                     an observation vector.
!    jd(jd(1)+1) = predictor variable deletion flag
!                  jd(1) = 0  => use all variables
!                  jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
!    pf(nvars) = relative penalties for each predictor variable
!                pf(j) = 0 => jth variable unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => minimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    isd = standarization flag:
!          isd = 0 => fit on original predictor variables
!          isd = 1 => fit on standardized predictor variables
!          Note: output solutions always reference original
!                variables locations and scales.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    beta(pmax, nalam) = compressed coefficient values for each solution
!    ibeta(pmax) = pointers to compressed coefficients
!    nbeta(nalam) = number of compressed coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 7777 => all used predictors have zero variance
!                    jerr = 10000 => maxval(vp) <= 0.0
!           jerr < 0 => non fatal error - partial output:
!                    Solutions for larger lambdas (1:(k-1)) returned.
!                    jerr = -k => convergence for kth lambda value not reached
!                           after maxit (see above) iterations.
!                    jerr = -10000-k => number of non zero coefficients along path
!                           exceeds pmax (see above) at kth lambda value.
! 
! LICENSE: GNU GPL (version 2)
! 
! AUTHORS:
!    Yi Yang (yi.yang6@mcgill.ca) and Hui Zou (hzou@stat.umn.edu), 
!    School of Statistics, University of Minnesota.
! 
! REFERENCES:
!    Yang, Y. and Zou, H. (2012). A Cocktail Algorithm for Solving 
!    The Elastic Net Penalized Cox's Regression in High Dimensions
!    Statistics and Its Interface. 6:2, 167-173. 


! --------------------------------------------------
SUBROUTINE coxlassoNET(rs,nrs,alpha,nobs,nvars,x,jd,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,isd,maxit,nalam,beta,ibeta,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: dfmax
    INTEGER :: pmax
    INTEGER :: nlam
    INTEGER :: isd
    INTEGER :: nalam
    INTEGER :: npass
    INTEGER :: jerr
    INTEGER :: maxit
    INTEGER :: jd(*)
    INTEGER :: ibeta(pmax)
    INTEGER :: nbeta(nlam)
    DOUBLE PRECISION :: alpha
    DOUBLE PRECISION :: flmin
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: x(nobs,nvars)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: ulam(nlam)
    DOUBLE PRECISION :: beta(pmax,nlam)
    DOUBLE PRECISION :: alam(nlam)
    ! - - - local declarations - - -
    INTEGER :: j
    INTEGER :: l
    INTEGER :: nk
    INTEGER :: nrs
    INTEGER :: rs(nrs)
    INTEGER, DIMENSION (:), ALLOCATABLE :: ju
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION ::xnorm(nvars)
! - - - begin - - -    
! - - - allocate variables - - -      
    ALLOCATE(ju(1:nvars))                        
    xnorm=0.0
    xmean=0.0
    CALL chkvars(nobs,nvars,x,ju)       
    IF(jd(1)>0) ju(jd(2:(jd(1)+1)))=0             
    IF(maxval(ju) <= 0) THEN                        
        jerr=7777                                     
        RETURN          
    ENDIF
    IF(maxval(pf) <= 0.0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)                                  
    pf=pf*nvars/sum(pf)
    CALL standardcox(nobs,nvars,x,ju,isd,xmean,xnorm) 
    CALL coxlassoNETpath(rs,nrs,alpha,nobs,nvars,x,ju,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,nalam,beta,ibeta,nbeta,alam,npass,jerr)                                 
    IF(jerr>0) RETURN    ! check error after calling function    
! - - - organize beta afterward - - -                      
    DO l=1,nalam                                 
          nk=nbeta(l)                                                                   
          IF(isd==1) THEN                           
             DO j=1,nk                                  
                  beta(j,l)=beta(j,l)/xnorm(ibeta(j))                        
            ENDDO
         ENDIF                                     
    ENDDO    
    DEALLOCATE(ju)                             
    RETURN                                     
END SUBROUTINE coxlassoNET
! --------------------------------------------------
SUBROUTINE coxlassoNETpath(rs,nrs,alpha,nobs,nvars,x,ju,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,nalam,beta,m,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
        ! - - - arg types - - -    
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl=1.0E-6
    INTEGER, PARAMETER :: mnlam=6
    INTEGER :: mnl  
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: dfmax
    INTEGER :: pmax
    INTEGER :: nlam
    INTEGER :: maxit
    INTEGER :: nalam
    INTEGER :: npass
    INTEGER :: jerr
    INTEGER :: ju(nvars)
    INTEGER :: m(pmax)
    INTEGER :: nbeta(nlam)
    INTEGER :: nrs
    INTEGER :: rs(nrs)
    DOUBLE PRECISION :: alpha
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: x(nobs,nvars)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: beta(pmax,nlam)
    DOUBLE PRECISION :: ulam(nlam)
    DOUBLE PRECISION ::  alam(nlam)
    ! - - - local declarations - - -              
    DOUBLE PRECISION :: d
    DOUBLE PRECISION :: dif
    DOUBLE PRECISION :: oldb
    DOUBLE PRECISION :: u
    DOUBLE PRECISION :: v
    DOUBLE PRECISION :: al
    DOUBLE PRECISION :: alf
    DOUBLE PRECISION :: flmin
    DOUBLE PRECISION :: fj
    double precision :: oma
    double precision :: sb
    double precision :: sa
    DOUBLE PRECISION :: tlam
    DOUBLE PRECISION :: al0
!     DOUBLE PRECISION::fjj(nvars)
    DOUBLE PRECISION :: xb(nobs)
    DOUBLE PRECISION :: b(nvars)
    DOUBLE PRECISION :: oldbeta(nvars)
    DOUBLE PRECISION :: maj(nvars)
    double precision :: ga(nvars)
    double precision :: vl(nvars)
    INTEGER :: s
    INTEGER :: k
    INTEGER :: j
    INTEGER :: l
    INTEGER :: ix
    INTEGER :: ni
    INTEGER :: me
    INTEGER :: mm(nvars)
    INTEGER :: one(nobs)
    INTEGER :: ixx(nvars)
! - - - begin - - -
! - - - some initial setup - - -
    ixx=0      
    oma=1.0-alpha
    one=0
    b=0.0                                   
    oldbeta=0.0
    m=0                                       
    mm=0                                      
    npass=0                                     
    ni=npass
    mnl=min(mnlam,nlam)   
    xb = 1.0     
    maj = 0.0
    al = 0.0D0
    alf = 0.0D0
    DO j = 1,nvars
        DO s = 1, nrs
            maj(j) = maj(j) + (maxval(x(rs(s):nobs,j))-minval(x(rs(s):nobs,j)))**2
        ENDDO
    ENDDO
    maj = 0.25*maj/nobs
    IF(flmin < 1.0) THEN
        flmin=max(mfl,flmin)          
        alf=flmin**(1.0/(nlam-1.0))                                                       
    ENDIF
    vl = 0.0
    call drv(nobs,nvars,rs,nrs,x,xb,vl)
    ga = abs(vl)
    DO l=1,nlam
        al0 = al   
        IF(flmin>=1.0) THEN                 
            al=ulam(l)                                 
        ELSE 
            IF(l > 2) THEN                          
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN  
                al0=0.0
                DO j = 1,nvars
                    IF(ju(j)==0) cycle
                    IF(pf(j)>0.0) al0 = max(al0, ga(j)/pf(j)) 
                END DO
                al0 = al0 / alpha
                al = al0 * alf
            ENDIF
        ENDIF
        sa=alpha*al
        sb=oma*al
        tlam=alpha*(2.0*al-al0)                                    
        do k=1,nvars
            if(ixx(k) == 1) cycle
            if(ju(k) == 0) cycle
            if(ga(k) > tlam*pf(k)) ixx(k)=1
        enddo
        ! --------- outer loop ----------------------------                                                              
        DO  
            IF(ni>0) oldbeta(m(1:ni))=b(m(1:ni))    
        ! --middle loop-------------------------------------
            DO
                npass=npass+1              
                dif=0.0      
                DO k=1,nvars
                    if(ixx(k) == 0) cycle                                                                                  
                    oldb=b(k)     
                    call  risk(nobs,rs,nrs,x(:,k),xb,fj)       
                    u = maj(k)*b(k)+fj
                    v = pf(k)*sa   
                    if(abs(u) <= v) then
                        b(k) = 0                              
                    else
                        b(k) = sign(abs(u) - v,u)/(maj(k)+pf(k)*sb)                                     
                    endif
                    d=b(k)-oldb                                    
                    IF(abs(d)>0.0) THEN                          
                        dif=max(dif,maj(k)*d**2)    
                        xb=xb*exp(x(:,k)*d)
                        IF(mm(k)==0) THEN                         
                            ni=ni+1      
                            IF(ni>pmax) EXIT                           
                            mm(k)=ni                                    
                            m(ni)=k      !indicate which one is non-zero
                        ENDIF
                    ENDIF
                ENDDO  
                IF(ni>pmax) EXIT                              
                IF(dif<eps) EXIT
                if(npass > maxit) then                               
                      jerr=-l                                                  
                      return
                endif    
        ! --inner loop----------------------                                      
                DO                    
                    npass=npass+1
                    dif=0.0     
                    DO j=1,ni                                
                        k=m(j)                                    
                        oldb=b(k)
                        call  risk(nobs,rs,nrs,x(:,k),xb,fj)                                 
                        u = maj(k)*b(k)+fj
                        v = pf(k)*sa   
                        if(abs(u) <= v) then
                            b(k) = 0                              
                        else
                            b(k) = sign(abs(u) - v,u)/(maj(k)+pf(k)*sb)                                     
                        endif         
                        d=b(k)-oldb                                    
                        IF(abs(d)>0.0) THEN                        
                            dif=max(dif,maj(k)*d**2)    
                            xb=xb*exp(x(:,k)*d)
                        ENDIF                         
                    ENDDO                  
                    IF(dif<eps) EXIT 
                      if(npass > maxit) then                               
                            jerr=-l                                                  
                            return
                      endif  
                ENDDO
            ENDDO                                 
            IF(ni>pmax) EXIT  
        !--- this is the final check ------------------------
            ix=0
            DO j=1,ni                                
                IF(maj(m(j))*(b(m(j))-oldbeta(m(j)))**2 < eps) cycle                 
                ix=1                                 
                EXIT
            ENDDO        
            IF(ix/=0) cycle     ! change this line when strong rule is added
            call drv(nobs,nvars,rs,nrs,x,xb,vl)
            do k=1,nvars                                            
                if(ixx(k)==1)cycle                                  
                if(ju(k)==0)cycle                                   
                ga(k)=abs(vl(k))
                if(ga(k) > sa*pf(k))then                          
                      ixx(k)=1                                                   
                    ix=1                                                       
                endif                                                   
            enddo 
            if(ix == 1) cycle
            exit
        ENDDO     
    ! final update variable save results------------                                                
        IF(ni>pmax) THEN                          
            jerr=-10000-l                                 
            EXIT
        ENDIF
        IF(ni>0) beta(1:ni,l)=b(m(1:ni))                    
        nbeta(l)=ni                                 
        alam(l)=al                                  
        nalam=l
        IF(l<mnl) CYCLE
        IF(flmin>=1.0) CYCLE       
        me = count(beta(1:ni,l)/=0.0)                                                        
        IF(me>dfmax) EXIT                                                              
    ENDDO    
    RETURN                                    
END SUBROUTINE coxlassoNETpath
