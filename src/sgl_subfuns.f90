MODULE sgl_subfuns

   USE spmatmul
   IMPLICIT NONE

   CONTAINS

   SUBROUTINE strong_rule (is_in_E_set, ga, pf, tlam)
      IMPLICIT NONE
      INTEGER :: g, k
      INTEGER, DIMENSION (:), INTENT(inout) :: is_in_E_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: ga
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: pf
      DOUBLE PRECISION, INTENT(in) :: tlam
      DOUBLE PRECISION :: z
      k = SIZE(is_in_E_set)
      !z = tlam * (1 - alsparse)
      z = tlam
      DO g = 1, k
         IF (is_in_E_set(g) .EQ. 1) CYCLE
         IF (ga(g) > pf(g) * z) is_in_E_set(g) = 1
      ENDDO
      RETURN
   END SUBROUTINE strong_rule

   SUBROUTINE kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf,pfl1, bs, lama1,&
   lama2,lama3, ga, nvars,drgix,drgiy,cn,cn_s,cn_e,b)
      IMPLICIT NONE
      INTEGER :: g, startix, endix, nvars
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) :: bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: vl
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lama1,lama2,lama3
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: theta_tilt
      INTEGER:: startdrg, enddrg, drginx, drgstartinx, drgendinx, theta_tilt_start, theta_tilt_end, theta_drg_norm
      INTEGER:: cn
      INTEGER:: drgix(cn)
      INTEGER:: drgiy(cn)
      INTEGER:: cn_s(bn)
      INTEGER:: cn_e(bn),drginxsub
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(in) :: b
            
      DO g = 1, bn
         IF (is_in_E_set(g) .EQ. 1) CYCLE
         startix = ix(g)
         endix = iy(g)
         
         ALLOCATE(theta_tilt (bs(g)))



          startdrg=cn_s(g)
        	enddrg=cn_e(g)
        	theta_tilt=b(startix:endix)
          !theta_tilt(1)= SIGN(b(startix)/b(startix), b(startix))-lama*pfl1(1)
          theta_tilt(1)=0
        	DO drginx =  startdrg, enddrg
        		drgstartinx=drgix(drginx)
        		drgendinx=drgiy(drginx)
            theta_drg_norm= SQRT(DOT_PRODUCT(b(drgstartinx: drgendinx), b(drgstartinx: drgendinx)))
            theta_tilt_start= drgstartinx- startix+1
            theta_tilt_end= drgendinx- startix+1
            
            
            DO drginxsub = theta_tilt_start, theta_tilt_end
            theta_tilt(drginxsub)= b(drginxsub+startix-1)*lama2*pfl1(drginxsub+startix-1)/ theta_drg_norm
            ENDDO
            !theta_tilt(theta_tilt_start: theta_tilt_end)= b(drgstartinx: drgendinx)/ theta_drg_norm
        	ENDDO

         
         ALLOCATE(s(bs(g)))
         s = vl(startix:endix)
         !CALL softthresh(s, lama*pfl1(startix:endix), bs(g))
         !CALL softthresh(s, lama*pfl1(startix:endix)+ lama*pfl1* theta_tilt, bs(g))
         
         CALL softthresh(s, lama3*pfl1(startix:endix)+ theta_tilt, bs(g))
         snorm = SQRT(DOT_PRODUCT(s,s))
         ga(g) = snorm
         IF(ga(g) > pf(g)*lama1 ) THEN
            is_in_E_set(g) = 1
            violation = 1
         ENDIF
         DEALLOCATE(s)
         DEALLOCATE(theta_tilt)
      ENDDO
      RETURN
   END SUBROUTINE kkt_check


   SUBROUTINE update_step(bsg, startix, endix, b, lama1,lama2,lama3, t_for_sg, pfg, pfl1, x,&
         isDifZero, nobs, r, gamg, maxDif,nvars, lb, ub, startsubgp, endsubgp, drgix, drgiy,drgdim)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars
      INTEGER, INTENT(in) :: startix, endix
      DOUBLE PRECISION :: gamg
      DOUBLE PRECISION, INTENT(inout) :: maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb, s, dd
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama1,lama2,lama3, t_for_sg, pfg, lb, ub
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: x(nobs,nvars)
      INTEGER, INTENT(inout) :: isDifZero
      INTEGER :: k
      INTEGER :: i
      INTEGER :: subgp_ix, subgp_iy
      INTEGER, INTENT(in) :: startsubgp, endsubgp, drgdim
      INTEGER, INTENT(in) :: drgix(drgdim),drgiy(drgdim)
      INTEGER :: s_endinx


      ALLOCATE(s(bsg))
      ALLOCATE(oldb(bsg))
      isDifZero = 0
      oldb = b(startix:endix)
      s = MATMUL(r, x(:, startix:endix))/nobs
      s = s*t_for_sg + b(startix:endix)
      
      s_endinx = endix-startix+1
      CALL softthresh(s(2:s_endinx), lama3*t_for_sg*pfl1, bsg-1)
      
      DO i = startsubgp, endsubgp
      subgp_ix = drgix(i) - drgix(startsubgp) +2
      subgp_iy = drgiy(i) - drgix(startsubgp) +2
      snorm = SQRT(DOT_PRODUCT(s(subgp_ix: subgp_iy), s(subgp_ix: subgp_iy)))
      tea = snorm - t_for_sg * lama2 * pfg 
      IF (tea > 0.0D0) THEN
         s(subgp_ix: subgp_iy) = s(subgp_ix: subgp_iy)* tea / snorm
      ELSE
         s(subgp_ix: subgp_iy) = 0.0D0
      ENDIF
      ENDDO

      
      
      snorm = SQRT(DOT_PRODUCT(s,s))
      tea = snorm - t_for_sg * lama1 * pfg
      IF (tea > 0.0D0) THEN
         b(startix:endix) = s * tea / snorm
         DO k = startix, endix
            b(k) = MIN(MAX(lb, b(k)), ub)
         ENDDO
         
         IF(startix.EQ.1) THEN
            b(1)=s(1)
         ENDIF
         
      ELSE
         b(startix:endix) = 0.0D0
         
         IF(startix.EQ.1) THEN
            b(1)=s(1)
         ENDIF
         
      ENDIF
      
      
      
      
      ALLOCATE(dd(bsg))
      dd = b(startix:endix) - oldb
      IF (ANY(dd .ne. 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         r = r - MATMUL(x(:,startix:endix), dd)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE update_step


   SUBROUTINE strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,pfl1,bs,&
         lama1,lama2,lama3,ga,is_in_S_set,x,r,nobs,nvars,vl,drgix,drgiy,cn,cn_s,cn_e,b)
      IMPLICIT NONE
      INTEGER, INTENT(in)::nobs
      INTEGER, INTENT(in)::nvars
      DOUBLE PRECISION,INTENT(in):: x(nobs, nvars)
      DOUBLE PRECISION, INTENT(in):: r(nobs)
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: vl
      INTEGER :: g, startix, endix
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) ::bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      INTEGER, DIMENSION(:), INTENT(in) :: is_in_S_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lama1,lama2,lama3
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: theta_tilt
      INTEGER:: startdrg, enddrg, drginx, drgstartinx, drgendinx, theta_tilt_start, theta_tilt_end, theta_drg_norm
      INTEGER:: cn,drginxsub
      INTEGER:: drgix(cn)
      INTEGER:: drgiy(cn)
      INTEGER:: cn_s(bn)
      INTEGER:: cn_e(bn)
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      
      violation = 0
      DO g = 1, bn
         IF (is_in_S_set(g) .EQ. 1) THEN
            startix = ix(g)
            endix = iy(g)
            
            ALLOCATE(theta_tilt (bs(g)))


            
            startdrg=cn_s(g)
            enddrg=cn_e(g)
            theta_tilt=b(startix:endix)
            !theta_tilt(1)= SIGN(b(startix)/b(startix), b(startix))-lama*pfl1(1)
            theta_tilt(1)=0
          	DO drginx =  startdrg, enddrg
          		drgstartinx=drgix(drginx)
          		drgendinx=drgiy(drginx)
              theta_drg_norm= SQRT(DOT_PRODUCT(b(drgstartinx: drgendinx), b(drgstartinx: drgendinx)))
              theta_tilt_start= drgstartinx- startix+1
              theta_tilt_end= drgendinx- startix+1
              
              
               DO drginxsub = theta_tilt_start, theta_tilt_end
                theta_tilt(drginxsub)= b(drginxsub+startix-1)*lama2*pfl1(drginxsub+startix-1)/ theta_drg_norm
                ENDDO
                !theta_tilt(theta_tilt_start: theta_tilt_end)= b(drgstartinx: drgendinx)/ theta_drg_norm
              
          	ENDDO

            
            
            
            ALLOCATE(s(bs(g)))
            s = MATMUL(r, x(:,startix:endix)) / nobs
            vl(startix:endix) = s
            CALL softthresh(s, lama3*pfl1(startix:endix)+ theta_tilt, bs(g))
            
            !CALL softthresh(s, lama*pfl1(startix:endix)+ lama*pfl1* theta_tilt, bs(g))
            !CALL softthresh(s, lama*pfl1(startix:endix), bs(g))
            
            snorm = SQRT(dot_PRODUCT(s,s))
            ga(g) = snorm
            DEALLOCATE(s)
            DEALLOCATE(theta_tilt)
                     
            IF (is_in_E_set(g) .EQ. 1) CYCLE
            IF (ga(g) > pf(g)*lama1) THEN
               is_in_E_set(g) = 1
               violation = 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
   END SUBROUTINE strong_kkt_check

   SUBROUTINE sp_update_step(bsg, startix, endix, b, lama1,lama2,lama3, t_for_sg, pfg, pfl1, x,&
         xidx, xcptr, nnz, isDifZero, nobs, r, gamg, maxDif, nvars, lb, ub, startsubgp, endsubgp, drgix, drgiy,drgdim)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars, nnz
      INTEGER, INTENT(in) :: startix, endix
      DOUBLE PRECISION :: gamg
      DOUBLE PRECISION, INTENT(inout) :: maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb, s, dd
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama1,lama2,lama3, t_for_sg, pfg, lb, ub
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      INTEGER, INTENT(in) :: xidx(nnz)
      INTEGER, INTENT(in) :: xcptr(nvars + 1)
      INTEGER, INTENT(inout) :: isDifZero
      INTEGER :: k
      INTEGER :: i
      INTEGER :: subgp_ix, subgp_iy
      INTEGER, INTENT(in) :: startsubgp, endsubgp, drgdim
      INTEGER, INTENT(in) :: drgix(drgdim),drgiy(drgdim)
      INTEGER :: s_endinx






      ALLOCATE(s(bsg))
      ALLOCATE(oldb(bsg))
      isDifZero = 0
      s = 0.0D0
      oldb = b(startix:endix)
      ! print *, oldb

      CALL spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
      s = s * t_for_sg / nobs + b(startix:endix)
      
      
      s_endinx = endix-startix+1
      CALL softthresh(s(2:s_endinx), lama3*t_for_sg*pfl1, bsg-1)
      
      
      
      DO i = startsubgp, endsubgp
      subgp_ix = drgix(i) - drgix(startsubgp) +2
      subgp_iy = drgiy(i) - drgix(startsubgp) +2
      snorm = SQRT(DOT_PRODUCT(s(subgp_ix: subgp_iy), s(subgp_ix: subgp_iy)))
      tea = snorm - t_for_sg * lama2 * pfg 
      IF (tea > 0.0D0) THEN
         s(subgp_ix: subgp_iy) = s(subgp_ix: subgp_iy)* tea / snorm
      ELSE
         s(subgp_ix: subgp_iy) = 0.0D0
      ENDIF
      ENDDO

      
      
      
      
      
      
      snorm = SQRT(DOT_PRODUCT(s,s))
      tea = snorm - t_for_sg * lama1 * pfg
      IF (tea > 0.0D0) THEN
         b(startix:endix) = s * tea / snorm
         DO k = startix, endix
            b(k) = MIN(MAX(b(k), lb), ub)
         ENDDO
      ELSE
         b(startix:endix) = 0.0D0
      ENDIF
      
      IF(startix.EQ.1) THEN
         b(1)=s(1)
      ENDIF
         
         
      ALLOCATE(dd(bsg))
      dd = b(startix:endix) - oldb
      IF(ANY(dd .ne. 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         CALL ymspax(x, xidx, xcptr, nobs, nvars, nnz, dd, r, startix, endix, bsg)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE sp_update_step


   SUBROUTINE sp_strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,pfl1,bs,&
         lama1,lama2,lama3,ga,is_in_S_set,x,xidx,xcptr,nnz,r,nobs,nvars,vl,drgix,drgiy,cn,cn_s,cn_e,b)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: nobs, nvars, nnz
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      INTEGER, INTENT(in) :: xidx(nnz)
      INTEGER, INTENT(in) :: xcptr(nvars + 1)
      DOUBLE PRECISION, INTENT(in):: r(nobs)
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: vl
      INTEGER :: g, startix, endix
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) ::bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      INTEGER, DIMENSION(:), INTENT(in) :: is_in_S_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lama1,lama2,lama3
      
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: theta_tilt
      INTEGER:: startdrg, enddrg, drginx, drgstartinx, drgendinx, theta_tilt_start, theta_tilt_end, theta_drg_norm
      INTEGER:: cn,drginxsub
      INTEGER:: drgix(cn)
      INTEGER:: drgiy(cn)
      INTEGER:: cn_s(bn)
      INTEGER:: cn_e(bn)      
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      

      violation = 0
      DO g = 1, bn
         IF(is_in_S_set(g) .EQ. 1) THEN
            startix = ix(g)
            endix = iy(g)
            ALLOCATE(s(bs(g)))
            
            startdrg=cn_s(g)
            enddrg=cn_e(g)
            theta_tilt=b(startix:endix)
            !theta_tilt(1)= SIGN(b(startix)/b(startix), b(startix))-lama*pfl1(1)
            theta_tilt(1)=0
          	DO drginx =  startdrg, enddrg
          		drgstartinx=drgix(drginx)
          		drgendinx=drgiy(drginx)
              theta_drg_norm= SQRT(DOT_PRODUCT(b(drgstartinx: drgendinx), b(drgstartinx: drgendinx)))
              theta_tilt_start= drgstartinx- startix+1
              theta_tilt_end= drgendinx- startix+1
              
              
              DO drginxsub = theta_tilt_start, theta_tilt_end
                theta_tilt(drginxsub)= b(drginxsub+startix-1)*lama2*pfl1(drginxsub+startix-1)/ theta_drg_norm
              ENDDO
          	ENDDO




            s = 0.0D0
            CALL spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
            s= s / nobs
            
            !vl(startix:endix) = s / nobs
            vl(startix:endix) = s
            ! CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
            CALL softthresh(s, lama3*pfl1(startix:endix)+ theta_tilt, bs(g))
            snorm = SQRT(dot_PRODUCT(s,s))
            ! print *, "kkt snorm = ", snorm
            ga(g) = snorm
            DEALLOCATE(s)
            DEALLOCATE(theta_tilt)
            
            
            IF(is_in_E_set(g) .EQ. 1) CYCLE
            IF(ga(g) > pf(g) * lama1) THEN
               is_in_E_set(g) = 1
               violation = 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
   END SUBROUTINE sp_strong_kkt_check

END MODULE sgl_subfuns


