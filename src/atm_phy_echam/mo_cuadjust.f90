MODULE mo_cuadjust
#ifndef __xlC__
#define SWDIV_NOCHK(a,b) (a)/(b)
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#endif

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: vtmpc1
  USE mo_echam_convect_tables,     ONLY: prepare_ua_index_spline, lookup_ua_spline, &
                                   lookup_ua_list_spline, lookup_ubc, lookup_ubc_list
#ifdef _PROFILE
  USE mo_profile,        ONLY: trace_start, trace_stop
#endif
                 
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: cuadjtq

CONTAINS

SUBROUTINE cuadjtq(  kproma, kbdim, klev, kk,                 &
           pp,       pt,       pq,       ldidx, ldcnt,  kcall)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          D.SALMOND         CRAY(UK))      12/8/91
!
!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
!               KCALL=0    ENV. T AND QS IN*CUINI*
!               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!
!          EXTERNALS
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SETPHYS*.
!

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kcall, kk, klev, kproma, kbdim

  !  Array arguments with intent(In):
  REAL(wp), INTENT (IN) :: pp(kbdim)
  INTEGER, INTENT (IN) :: ldidx(kbdim)
  INTEGER, INTENT (IN) :: ldcnt

  !  Array arguments with intent(InOut):
  REAL(wp), INTENT (INOUT) :: pq(kbdim,klev), pt(kbdim,klev)

  !  Local scalars:
  REAL(wp):: zcond1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  INTEGER :: jl, nl, nsum

  !  Local arrays:
  REAL(wp) :: zcond(kbdim),zppi(kbdim),za(kbdim)
  REAL(wp) :: ua(kbdim),dua(kbdim),ub(kbdim),uc(kbdim)
  INTEGER  :: idx(kbdim),ncond(kbdim)


  !  Executable statements

#ifdef _PROFILE
  CALL trace_start ('cuadjtq', 32)
#endif
!
!----------------------------------------------------------------------
!
!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------
!

  IF (kcall >= 0.AND. kcall <= 2 ) THEN

     CALL lookup_ubc_list('cuadjtq (1)',kproma,ldcnt,ldidx(1),pt(1,kk),ub(1),uc(1))
     CALL lookup_ua_list_spline('cuadjtq (1)',kproma,ldcnt,ldidx(1),pt(1,kk),ua(1),dua(1))

!DIR$ IVDEP
!OCL NOVREC
!IBM* ASSERT(NODEPS)
     DO 111 nl=1,ldcnt
        jl = ldidx(nl)
        zppi(jl)=1._wp/pp(jl)
111  END DO

     IF (kcall.EQ.0) THEN

        ! mpuetz: 40% of cuadtjq (lookup of tlucua und tlucub !!)

!DIR$ IVDEP
!OCL NOVREC
!IBM ASSERT(NODEPS)
        DO 112 nl=1,ldcnt
           jl = ldidx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_wp,zes)
           zcor = SWDIV_NOCHK(1._wp,(1._wp-vtmpc1*zes))
           zqsat= zes*zcor

           zdqsdt = zppi(jl)*zcor**2*dua(nl)
           zlcdqsdt  = FSEL(zes-0.4_wp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))

           zcond(jl) = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._wp+zlcdqsdt))

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond(jl)
           pq(jl,kk) = pq(jl,kk) - zcond(jl)

           ! mpuetz: zcond is almost always > 0
           ncond(nl) = INT(FSEL(-ABS(zcond(jl)),0._wp,1._wp))
112     END DO

     ELSE IF (kcall.EQ.1) THEN

!DIR$ IVDEP
!OCL NOVREC
!IBM ASSERT(NODEPS)
        DO 212 nl=1,ldcnt
           jl = ldidx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_wp,zes)
           zcor = SWDIV_NOCHK(1._wp,(1._wp-vtmpc1*zes))
           zqsat= zes*zcor

           zdqsdt = zppi(jl)*zcor**2*dua(nl)
           zlcdqsdt  = FSEL(zes-0.4_wp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))

           zcond(jl) = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._wp+zlcdqsdt))
           zcond(jl) = MAX(zcond(jl),0._wp)

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond(jl)
           pq(jl,kk) = pq(jl,kk) - zcond(jl)
           ! mpuetz: zcond is almost always > 0
           ncond(nl) = INT(FSEL(-ABS(zcond(jl)),0._wp,1._wp))
212     END DO

     ELSE

!DIR$ IVDEP
!OCL NOVREC
!IBM* ASSERT(NODEPS)
        DO 312 nl=1,ldcnt
           jl = ldidx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_wp,zes)
           zcor = SWDIV_NOCHK(1._wp,(1._wp-vtmpc1*zes))
           zqsat= zes*zcor

           zdqsdt = zppi(jl)*zcor**2*dua(nl)
           zlcdqsdt  = FSEL(zes-0.4_wp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))
           zcond(jl) = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._wp+zlcdqsdt))

           zcond(jl) = MIN(zcond(jl),0._wp)

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond(jl)
           pq(jl,kk) = pq(jl,kk) - zcond(jl)
           ! mpuetz: zcond is almost always > 0
           ncond(nl) = INT(FSEL(-ABS(zcond(jl)),0._wp,1._wp))
312     END DO

     END IF

     nsum = 1
     DO nl=1,ldcnt
        idx(nsum) = ldidx(nl)
        nsum = nsum + ncond(nl)
     END DO
     nsum = nsum - 1

#ifdef __ibmdbg__
     print *,'cuad(',kcall,')',ldcnt,nsum,ldcnt
#endif

     IF (nsum > 0) THEN

        CALL lookup_ubc_list('cuadjtq (2)',kproma,nsum,idx(1),pt(1,kk),ub(1),uc(1))
        CALL lookup_ua_list_spline('cuadjtq (2)',kproma,nsum,idx(1),pt(1,kk),ua(1),dua(1))

!DIR$ IVDEP
!OCL NOVREC
!IBM* ASSERT(NODEPS)
!IBM* UNROLL(3)
        DO 116 nl=1,nsum
           jl = idx(nl)

           zes  = ua(nl)*zppi(jl)
           zes  = MIN(0.5_wp,zes)
           zcor = SWDIV_NOCHK(1._wp,(1._wp-vtmpc1*zes))
           zqsat= zes*zcor

           zdqsdt = zppi(jl)*zcor**2*dua(nl)
           zlcdqsdt = FSEL(zes-0.4_wp,zqsat*zcor*ub(nl),zdqsdt*uc(nl))
           zcond1   = SWDIV_NOCHK((pq(jl,kk)-zqsat),(1._wp+zlcdqsdt))

           pt(jl,kk) = pt(jl,kk) + uc(nl)*zcond1
           pq(jl,kk) = pq(jl,kk) - zcond1
116     END DO

     END IF

  END IF

#ifdef _PROFILE
  CALL trace_stop ('cuadjtq', 32)
#endif

  RETURN
END SUBROUTINE cuadjtq

END MODULE mo_cuadjust