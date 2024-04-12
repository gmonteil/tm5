!*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 1999-2009 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
!
!*******************************************************************************
!  Content:
!      Intel(R) Math Kernel Library (MKL) Fortran interface for LAPACK routines
!*******************************************************************************

!DEC$IF .NOT. DEFINED (MKL_LAPACK_HEADER)

!DEC$DEFINE MKL_LAPACK_HEADER

! LAPACK 3.2 routines

      INTERFACE
      SUBROUTINE CBDSQR(UPLO,N,NCVT,NRU,NCC,D,E,VT,LDVT,U,LDU,C,LDC, &
                        RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDC,LDU,LDVT,N,NCC,NCVT,NRU
      REAL               D(*),E(*),RWORK(*)
      COMPLEX            C(LDC,*),U(LDU,*),VT(LDVT,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBBRD(VECT,M,N,NCC,KL,KU,AB,LDAB,D,E,Q,LDQ,PT,LDPT,C, &
                        LDC,WORK,RWORK,INFO)
      CHARACTER          VECT
      INTEGER            INFO,KL,KU,LDAB,LDC,LDPT,LDQ,M,N,NCC
      REAL               D(*),E(*),RWORK(*)
      COMPLEX            AB(LDAB,*),C(LDC,*),PT(LDPT,*),Q(LDQ,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBCON(NORM,N,KL,KU,AB,LDAB,IPIV,ANORM,RCOND,WORK, &
                        RWORK,INFO)
      CHARACTER          NORM
      INTEGER            INFO,KL,KU,LDAB,N
      REAL               ANORM,RCOND
      INTEGER            IPIV(*)
      REAL               RWORK(*)
      COMPLEX            AB(LDAB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBEQUB(M,N,KL,KU,AB,LDAB,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      REAL               AMAX,COLCND,ROWCND
      REAL               C(*),R(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBEQU(M,N,KL,KU,AB,LDAB,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      REAL               AMAX,COLCND,ROWCND
      REAL               C(*),R(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBRFS(TRANS,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB,IPIV,B,LDB, &
                        X,LDX,FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NRHS
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),WORK(*),X(LDX, &
                         *)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBRFSX(TRANS,EQUED,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB, &
                         IPIV,R,C,B,LDB,X,LDX,RCOND,BERR,N_ERR_BNDS, &
                         ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS,PARAMS, &
                         WORK,RWORK,INFO)
      CHARACTER          EQUED,TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NPARAMS,NRHS, &
                         N_ERR_BNDS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),C(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),R(*),RWORK(*)
      COMPLEX            AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),WORK(*),X(LDX, &
                         *)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBSV(N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)
      INTEGER            INFO,KL,KU,LDAB,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            AB(LDAB,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBSVX(FACT,TRANS,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB,IPIV, &
                        EQUED,R,C,B,LDB,X,LDX,RCOND,FERR,BERR,WORK, &
                        RWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NRHS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),C(*),FERR(*),R(*),RWORK(*)
      COMPLEX            AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),WORK(*),X(LDX, &
                         *)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBSVXX(FACT,TRANS,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB,IPIV, &
                         EQUED,R,C,B,LDB,X,LDX,RCOND,RPVGRW,BERR, &
                         N_ERR_BNDS,ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS, &
                         PARAMS,WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NPARAMS,NRHS, &
                         N_ERR_BNDS
      REAL               RCOND,RPVGRW
      INTEGER            IPIV(*)
      REAL               BERR(*),C(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),R(*),RWORK(*)
      COMPLEX            AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),WORK(*),X(LDX, &
                         *)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      INTEGER            IPIV(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      INTEGER            IPIV(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGBTRS(TRANS,N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,KL,KU,LDAB,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            AB(LDAB,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEBAK(JOB,SIDE,N,ILO,IHI,SCALE,M,V,LDV,INFO)
      CHARACTER          JOB,SIDE
      INTEGER            IHI,ILO,INFO,LDV,M,N
      REAL               SCALE(*)
      COMPLEX            V(LDV,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEBAL(JOB,N,A,LDA,ILO,IHI,SCALE,INFO)
      CHARACTER          JOB
      INTEGER            IHI,ILO,INFO,LDA,N
      REAL               SCALE(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEBD2(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      REAL               D(*),E(*)
      COMPLEX            A(LDA,*),TAUP(*),TAUQ(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEBRD(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      REAL               D(*),E(*)
      COMPLEX            A(LDA,*),TAUP(*),TAUQ(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGECON(NORM,N,A,LDA,ANORM,RCOND,WORK,RWORK,INFO)
      CHARACTER          NORM
      INTEGER            INFO,LDA,N
      REAL               ANORM,RCOND
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEEQUB(M,N,A,LDA,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,LDA,M,N
      REAL               AMAX,COLCND,ROWCND
      REAL               C(*),R(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEEQU(M,N,A,LDA,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,LDA,M,N
      REAL               AMAX,COLCND,ROWCND
      REAL               C(*),R(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEES(JOBVS,SORT,SELECT,N,A,LDA,SDIM,W,VS,LDVS,WORK, &
                       LWORK,RWORK,BWORK,INFO)
      CHARACTER          JOBVS,SORT
      INTEGER            INFO,LDA,LDVS,LWORK,N,SDIM
      LOGICAL            BWORK(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),VS(LDVS,*),W(*),WORK(*)
      LOGICAL            SELECT
      EXTERNAL           SELECT
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEESX(JOBVS,SORT,SELECT,SENSE,N,A,LDA,SDIM,W,VS,LDVS, &
                        RCONDE,RCONDV,WORK,LWORK,RWORK,BWORK,INFO)
      CHARACTER          JOBVS,SENSE,SORT
      INTEGER            INFO,LDA,LDVS,LWORK,N,SDIM
      REAL               RCONDE,RCONDV
      LOGICAL            BWORK(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),VS(LDVS,*),W(*),WORK(*)
      LOGICAL            SELECT
      EXTERNAL           SELECT
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEEV(JOBVL,JOBVR,N,A,LDA,W,VL,LDVL,VR,LDVR,WORK,LWORK, &
                       RWORK,INFO)
      CHARACTER          JOBVL,JOBVR
      INTEGER            INFO,LDA,LDVL,LDVR,LWORK,N
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),VL(LDVL,*),VR(LDVR,*),W(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEEVX(BALANC,JOBVL,JOBVR,SENSE,N,A,LDA,W,VL,LDVL,VR, &
                        LDVR,ILO,IHI,SCALE,ABNRM,RCONDE,RCONDV,WORK, &
                        LWORK,RWORK,INFO)
      CHARACTER          BALANC,JOBVL,JOBVR,SENSE
      INTEGER            IHI,ILO,INFO,LDA,LDVL,LDVR,LWORK,N
      REAL               ABNRM
      REAL               RCONDE(*),RCONDV(*),RWORK(*),SCALE(*)
      COMPLEX            A(LDA,*),VL(LDVL,*),VR(LDVR,*),W(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEGS(JOBVSL,JOBVSR,N,A,LDA,B,LDB,ALPHA,BETA,VSL,LDVSL, &
                       VSR,LDVSR,WORK,LWORK,RWORK,INFO)
      CHARACTER          JOBVSL,JOBVSR
      INTEGER            INFO,LDA,LDB,LDVSL,LDVSR,LWORK,N
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),VSL(LDVSL, &
                         *),VSR(LDVSR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEGV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHA,BETA,VL,LDVL,VR, &
                       LDVR,WORK,LWORK,RWORK,INFO)
      CHARACTER          JOBVL,JOBVR
      INTEGER            INFO,LDA,LDB,LDVL,LDVR,LWORK,N
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),VL(LDVL,*), &
                         VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEHD2(N,ILO,IHI,A,LDA,TAU,WORK,INFO)
      INTEGER            IHI,ILO,INFO,LDA,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEHRD(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            IHI,ILO,INFO,LDA,LWORK,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGELQ2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGELQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGELSD(M,N,NRHS,A,LDA,B,LDB,S,RCOND,RANK,WORK,LWORK, &
                        RWORK,IWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS,RANK
      REAL               RCOND
      INTEGER            IWORK(*)
      REAL               RWORK(*),S(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGELS(TRANS,M,N,NRHS,A,LDA,B,LDB,WORK,LWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGELSS(M,N,NRHS,A,LDA,B,LDB,S,RCOND,RANK,WORK,LWORK, &
                        RWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS,RANK
      REAL               RCOND
      REAL               RWORK(*),S(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGELSX(M,N,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,RWORK, &
                        INFO)
      INTEGER            INFO,LDA,LDB,M,N,NRHS,RANK
      REAL               RCOND
      INTEGER            JPVT(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGELSY(M,N,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK, &
                        RWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS,RANK
      REAL               RCOND
      INTEGER            JPVT(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEQL2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEQLF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEQP3(M,N,A,LDA,JPVT,TAU,WORK,LWORK,RWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      INTEGER            JPVT(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEQPF(M,N,A,LDA,JPVT,TAU,WORK,RWORK,INFO)
      INTEGER            INFO,LDA,M,N
      INTEGER            JPVT(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEQR2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGERFS(TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX, &
                        FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGERFSX(TRANS,EQUED,N,NRHS,A,LDA,AF,LDAF,IPIV,R,C,B, &
                         LDB,X,LDX,RCOND,BERR,N_ERR_BNDS,ERR_BNDS_NORM, &
                         ERR_BNDS_COMP,NPARAMS,PARAMS,WORK,RWORK,INFO)
      CHARACTER          EQUED,TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),C(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),R(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGERQ2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGERQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGESC2(N,A,LDA,RHS,IPIV,JPIV,SCALE)
      INTEGER            LDA,N
      REAL               SCALE
      INTEGER            IPIV(*),JPIV(*)
      COMPLEX            A(LDA,*),RHS(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGESDD(JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,RWORK, &
                        IWORK,INFO)
      CHARACTER          JOBZ
      INTEGER            INFO,LDA,LDU,LDVT,LWORK,M,N
      INTEGER            IWORK(*)
      REAL               RWORK(*),S(*)
      COMPLEX            A(LDA,*),U(LDU,*),VT(LDVT,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGESVD(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK, &
                        RWORK,INFO)
      CHARACTER          JOBU,JOBVT
      INTEGER            INFO,LDA,LDU,LDVT,LWORK,M,N
      REAL               RWORK(*),S(*)
      COMPLEX            A(LDA,*),U(LDU,*),VT(LDVT,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      INTEGER            INFO,LDA,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGESVX(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,R,C, &
                        B,LDB,X,LDX,RCOND,FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),C(*),FERR(*),R(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGESVXX(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,R,C, &
                         B,LDB,X,LDX,RCOND,RPVGRW,BERR,N_ERR_BNDS, &
                         ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS,PARAMS, &
                         WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND,RPVGRW
      INTEGER            IPIV(*)
      REAL               BERR(*),C(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),R(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGETC2(N,A,LDA,IPIV,JPIV,INFO)
      INTEGER            INFO,LDA,N
      INTEGER            IPIV(*),JPIV(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGETF2(M,N,A,LDA,IPIV,INFO)
      INTEGER            INFO,LDA,M,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGETRF(M,N,A,LDA,IPIV,INFO)
      INTEGER            INFO,LDA,M,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDA,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGBAK(JOB,SIDE,N,ILO,IHI,LSCALE,RSCALE,M,V,LDV,INFO)
      CHARACTER          JOB,SIDE
      INTEGER            IHI,ILO,INFO,LDV,M,N
      REAL               LSCALE(*),RSCALE(*)
      COMPLEX            V(LDV,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGBAL(JOB,N,A,LDA,B,LDB,ILO,IHI,LSCALE,RSCALE,WORK, &
                        INFO)
      CHARACTER          JOB
      INTEGER            IHI,ILO,INFO,LDA,LDB,N
      REAL               LSCALE(*),RSCALE(*),WORK(*)
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGES(JOBVSL,JOBVSR,SORT,SELCTG,N,A,LDA,B,LDB,SDIM, &
                       ALPHA,BETA,VSL,LDVSL,VSR,LDVSR,WORK,LWORK,RWORK, &
                       BWORK,INFO)
      CHARACTER          JOBVSL,JOBVSR,SORT
      INTEGER            INFO,LDA,LDB,LDVSL,LDVSR,LWORK,N,SDIM
      LOGICAL            BWORK(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),VSL(LDVSL, &
                         *),VSR(LDVSR,*),WORK(*)
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGESX(JOBVSL,JOBVSR,SORT,SELCTG,SENSE,N,A,LDA,B,LDB, &
                        SDIM,ALPHA,BETA,VSL,LDVSL,VSR,LDVSR,RCONDE, &
                        RCONDV,WORK,LWORK,RWORK,IWORK,LIWORK,BWORK,INFO)
      CHARACTER          JOBVSL,JOBVSR,SENSE,SORT
      INTEGER            INFO,LDA,LDB,LDVSL,LDVSR,LIWORK,LWORK,N,SDIM
      LOGICAL            BWORK(*)
      INTEGER            IWORK(*)
      REAL               RCONDE(2),RCONDV(2),RWORK(*)
      COMPLEX            A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),VSL(LDVSL, &
                         *),VSR(LDVSR,*),WORK(*)
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHA,BETA,VL,LDVL,VR, &
                       LDVR,WORK,LWORK,RWORK,INFO)
      CHARACTER          JOBVL,JOBVR
      INTEGER            INFO,LDA,LDB,LDVL,LDVR,LWORK,N
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),VL(LDVL,*), &
                         VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGEVX(BALANC,JOBVL,JOBVR,SENSE,N,A,LDA,B,LDB,ALPHA, &
                        BETA,VL,LDVL,VR,LDVR,ILO,IHI,LSCALE,RSCALE, &
                        ABNRM,BBNRM,RCONDE,RCONDV,WORK,LWORK,RWORK, &
                        IWORK,BWORK,INFO)
      CHARACTER          BALANC,JOBVL,JOBVR,SENSE
      INTEGER            IHI,ILO,INFO,LDA,LDB,LDVL,LDVR,LWORK,N
      REAL               ABNRM,BBNRM
      LOGICAL            BWORK(*)
      INTEGER            IWORK(*)
      REAL               LSCALE(*),RCONDE(*),RCONDV(*),RSCALE(*), &
                         RWORK(*)
      COMPLEX            A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),VL(LDVL,*), &
                         VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGGLM(N,M,P,A,LDA,B,LDB,D,X,Y,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      COMPLEX            A(LDA,*),B(LDB,*),D(*),WORK(*),X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGHRD(COMPQ,COMPZ,N,ILO,IHI,A,LDA,B,LDB,Q,LDQ,Z,LDZ, &
                        INFO)
      CHARACTER          COMPQ,COMPZ
      INTEGER            IHI,ILO,INFO,LDA,LDB,LDQ,LDZ,N
      COMPLEX            A(LDA,*),B(LDB,*),Q(LDQ,*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGLSE(M,N,P,A,LDA,B,LDB,C,D,X,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      COMPLEX            A(LDA,*),B(LDB,*),C(*),D(*),WORK(*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGQRF(N,M,P,A,LDA,TAUA,B,LDB,TAUB,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      COMPLEX            A(LDA,*),B(LDB,*),TAUA(*),TAUB(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGRQF(M,P,N,A,LDA,TAUA,B,LDB,TAUB,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      COMPLEX            A(LDA,*),B(LDB,*),TAUA(*),TAUB(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGSVD(JOBU,JOBV,JOBQ,M,N,P,K,L,A,LDA,B,LDB,ALPHA,BETA, &
                        U,LDU,V,LDV,Q,LDQ,WORK,RWORK,IWORK,INFO)
      CHARACTER          JOBQ,JOBU,JOBV
      INTEGER            INFO,K,L,LDA,LDB,LDQ,LDU,LDV,M,N,P
      INTEGER            IWORK(*)
      REAL               ALPHA(*),BETA(*),RWORK(*)
      COMPLEX            A(LDA,*),B(LDB,*),Q(LDQ,*),U(LDU,*),V(LDV,*), &
                         WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGGSVP(JOBU,JOBV,JOBQ,M,P,N,A,LDA,B,LDB,TOLA,TOLB,K,L, &
                        U,LDU,V,LDV,Q,LDQ,IWORK,RWORK,TAU,WORK,INFO)
      CHARACTER          JOBQ,JOBU,JOBV
      INTEGER            INFO,K,L,LDA,LDB,LDQ,LDU,LDV,M,N,P
      REAL               TOLA,TOLB
      INTEGER            IWORK(*)
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),B(LDB,*),Q(LDQ,*),TAU(*),U(LDU,*), &
                         V(LDV,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGTCON(NORM,N,DL,D,DU,DU2,IPIV,ANORM,RCOND,WORK,INFO)
      CHARACTER          NORM
      INTEGER            INFO,N
      REAL               ANORM,RCOND
      INTEGER            IPIV(*)
      COMPLEX            D(*),DL(*),DU(*),DU2(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGTRFS(TRANS,N,NRHS,DL,D,DU,DLF,DF,DUF,DU2,IPIV,B,LDB, &
                        X,LDX,FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDB,LDX,N,NRHS
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            B(LDB,*),D(*),DF(*),DL(*),DLF(*),DU(*),DU2(*), &
                         DUF(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGTSV(N,NRHS,DL,D,DU,B,LDB,INFO)
      INTEGER            INFO,LDB,N,NRHS
      COMPLEX            B(LDB,*),D(*),DL(*),DU(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGTSVX(FACT,TRANS,N,NRHS,DL,D,DU,DLF,DF,DUF,DU2,IPIV,B, &
                        LDB,X,LDX,RCOND,FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          FACT,TRANS
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            B(LDB,*),D(*),DF(*),DL(*),DLF(*),DU(*),DU2(*), &
                         DUF(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGTTRF(N,DL,D,DU,DU2,IPIV,INFO)
      INTEGER            INFO,N
      INTEGER            IPIV(*)
      COMPLEX            D(*),DL(*),DU(*),DU2(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGTTRS(TRANS,N,NRHS,DL,D,DU,DU2,IPIV,B,LDB,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            B(LDB,*),D(*),DL(*),DU(*),DU2(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CGTTS2(ITRANS,N,NRHS,DL,D,DU,DU2,IPIV,B,LDB)
      INTEGER            ITRANS,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            B(LDB,*),D(*),DL(*),DU(*),DU2(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBEVD(JOBZ,UPLO,N,KD,AB,LDAB,W,Z,LDZ,WORK,LWORK,RWORK, &
                        LRWORK,IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,KD,LDAB,LDZ,LIWORK,LRWORK,LWORK,N
      INTEGER            IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AB(LDAB,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBEV(JOBZ,UPLO,N,KD,AB,LDAB,W,Z,LDZ,WORK,RWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,KD,LDAB,LDZ,N
      REAL               RWORK(*),W(*)
      COMPLEX            AB(LDAB,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBEVX(JOBZ,RANGE,UPLO,N,KD,AB,LDAB,Q,LDQ,VL,VU,IL,IU, &
                        ABSTOL,M,W,Z,LDZ,WORK,RWORK,IWORK,IFAIL,INFO)
      CHARACTER          JOBZ,RANGE,UPLO
      INTEGER            IL,INFO,IU,KD,LDAB,LDQ,LDZ,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            IFAIL(*),IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AB(LDAB,*),Q(LDQ,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBGST(VECT,UPLO,N,KA,KB,AB,LDAB,BB,LDBB,X,LDX,WORK, &
                        RWORK,INFO)
      CHARACTER          UPLO,VECT
      INTEGER            INFO,KA,KB,LDAB,LDBB,LDX,N
      REAL               RWORK(*)
      COMPLEX            AB(LDAB,*),BB(LDBB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBGVD(JOBZ,UPLO,N,KA,KB,AB,LDAB,BB,LDBB,W,Z,LDZ,WORK, &
                        LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,KA,KB,LDAB,LDBB,LDZ,LIWORK,LRWORK,LWORK,N
      INTEGER            IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AB(LDAB,*),BB(LDBB,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBGV(JOBZ,UPLO,N,KA,KB,AB,LDAB,BB,LDBB,W,Z,LDZ,WORK, &
                       RWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,KA,KB,LDAB,LDBB,LDZ,N
      REAL               RWORK(*),W(*)
      COMPLEX            AB(LDAB,*),BB(LDBB,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBGVX(JOBZ,RANGE,UPLO,N,KA,KB,AB,LDAB,BB,LDBB,Q,LDQ, &
                        VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,RWORK,IWORK, &
                        IFAIL,INFO)
      CHARACTER          JOBZ,RANGE,UPLO
      INTEGER            IL,INFO,IU,KA,KB,LDAB,LDBB,LDQ,LDZ,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            IFAIL(*),IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AB(LDAB,*),BB(LDBB,*),Q(LDQ,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHBTRD(VECT,UPLO,N,KD,AB,LDAB,D,E,Q,LDQ,WORK,INFO)
      CHARACTER          UPLO,VECT
      INTEGER            INFO,KD,LDAB,LDQ,N
      REAL               D(*),E(*)
      COMPLEX            AB(LDAB,*),Q(LDQ,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHECON(UPLO,N,A,LDA,IPIV,ANORM,RCOND,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      REAL               ANORM,RCOND
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEEQUB(UPLO,N,A,LDA,S,SCOND,AMAX,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEEVD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,LRWORK, &
                        IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,LDA,LIWORK,LRWORK,LWORK,N
      INTEGER            IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,LDA,LWORK,N
      REAL               RWORK(*),W(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEEVR(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W, &
                        Z,LDZ,ISUPPZ,WORK,LWORK,RWORK,LRWORK,IWORK, &
                        LIWORK,INFO)
      CHARACTER          JOBZ,RANGE,UPLO
      INTEGER            IL,INFO,IU,LDA,LDZ,LIWORK,LRWORK,LWORK,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            ISUPPZ(*),IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            A(LDA,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W, &
                        Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
      CHARACTER          JOBZ,RANGE,UPLO
      INTEGER            IL,INFO,IU,LDA,LDZ,LWORK,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            IFAIL(*),IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            A(LDA,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEGS2(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,ITYPE,LDA,LDB,N
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,ITYPE,LDA,LDB,N
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEGVD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK, &
                        RWORK,LRWORK,IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,ITYPE,LDA,LDB,LIWORK,LRWORK,LWORK,N
      INTEGER            IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK, &
                       INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,ITYPE,LDA,LDB,LWORK,N
      REAL               RWORK(*),W(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHEGVX(ITYPE,JOBZ,RANGE,UPLO,N,A,LDA,B,LDB,VL,VU,IL,IU, &
                        ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL, &
                        INFO)
      CHARACTER          JOBZ,RANGE,UPLO
      INTEGER            IL,INFO,ITYPE,IU,LDA,LDB,LDZ,LWORK,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            IFAIL(*),IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHERFS(UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX,FERR, &
                        BERR,WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHERFSX(UPLO,EQUED,N,NRHS,A,LDA,AF,LDAF,IPIV,S,B,LDB,X, &
                         LDX,RCOND,BERR,N_ERR_BNDS,ERR_BNDS_NORM, &
                         ERR_BNDS_COMP,NPARAMS,PARAMS,WORK,RWORK,INFO)
      CHARACTER          EQUED,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),RWORK(*),S(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHESV(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,LWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDB,LWORK,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHESVX(FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX, &
                        RCOND,FERR,BERR,WORK,LWORK,RWORK,INFO)
      CHARACTER          FACT,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,LWORK,N,NRHS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHESVXX(FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,S,B, &
                         LDB,X,LDX,RCOND,RPVGRW,BERR,N_ERR_BNDS, &
                         ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS,PARAMS, &
                         WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND,RPVGRW
      INTEGER            IPIV(*)
      REAL               BERR(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),RWORK(*),S(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHETD2(UPLO,N,A,LDA,D,E,TAU,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      REAL               D(*),E(*)
      COMPLEX            A(LDA,*),TAU(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHETF2(UPLO,N,A,LDA,IPIV,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHETRD(UPLO,N,A,LDA,D,E,TAU,WORK,LWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LWORK,N
      REAL               D(*),E(*)
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHETRF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LWORK,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHETRI(UPLO,N,A,LDA,IPIV,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHETRS(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHFRK(TRANSR,UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C)
      CHARACTER          TRANS,TRANSR,UPLO
      INTEGER            K,LDA,N
      REAL               ALPHA,BETA
      COMPLEX            A(LDA,*),C(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHGEQZ(JOB,COMPQ,COMPZ,N,ILO,IHI,H,LDH,T,LDT,ALPHA, &
                        BETA,Q,LDQ,Z,LDZ,WORK,LWORK,RWORK,INFO)
      CHARACTER          COMPQ,COMPZ,JOB
      INTEGER            IHI,ILO,INFO,LDH,LDQ,LDT,LDZ,LWORK,N
      REAL               RWORK(*)
      COMPLEX            ALPHA(*),BETA(*),H(LDH,*),Q(LDQ,*),T(LDT,*), &
                         WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CHLA_TRANSTYPE(TRANS)
      CHARACTER          CHLA_TRANSTYPE
      INTEGER            TRANS
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPCON(UPLO,N,AP,IPIV,ANORM,RCOND,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      REAL               ANORM,RCOND
      INTEGER            IPIV(*)
      COMPLEX            AP(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPEVD(JOBZ,UPLO,N,AP,W,Z,LDZ,WORK,LWORK,RWORK,LRWORK, &
                        IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,LDZ,LIWORK,LRWORK,LWORK,N
      INTEGER            IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AP(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPEV(JOBZ,UPLO,N,AP,W,Z,LDZ,WORK,RWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,LDZ,N
      REAL               RWORK(*),W(*)
      COMPLEX            AP(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPEVX(JOBZ,RANGE,UPLO,N,AP,VL,VU,IL,IU,ABSTOL,M,W,Z, &
                        LDZ,WORK,RWORK,IWORK,IFAIL,INFO)
      CHARACTER          JOBZ,RANGE,UPLO
      INTEGER            IL,INFO,IU,LDZ,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            IFAIL(*),IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AP(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPGST(ITYPE,UPLO,N,AP,BP,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,ITYPE,N
      COMPLEX            AP(*),BP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPGVD(ITYPE,JOBZ,UPLO,N,AP,BP,W,Z,LDZ,WORK,LWORK, &
                        RWORK,LRWORK,IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,ITYPE,LDZ,LIWORK,LRWORK,LWORK,N
      INTEGER            IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AP(*),BP(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPGV(ITYPE,JOBZ,UPLO,N,AP,BP,W,Z,LDZ,WORK,RWORK,INFO)
      CHARACTER          JOBZ,UPLO
      INTEGER            INFO,ITYPE,LDZ,N
      REAL               RWORK(*),W(*)
      COMPLEX            AP(*),BP(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPGVX(ITYPE,JOBZ,RANGE,UPLO,N,AP,BP,VL,VU,IL,IU, &
                        ABSTOL,M,W,Z,LDZ,WORK,RWORK,IWORK,IFAIL,INFO)
      CHARACTER          JOBZ,RANGE,UPLO
      INTEGER            IL,INFO,ITYPE,IU,LDZ,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            IFAIL(*),IWORK(*)
      REAL               RWORK(*),W(*)
      COMPLEX            AP(*),BP(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPRFS(UPLO,N,NRHS,AP,AFP,IPIV,B,LDB,X,LDX,FERR,BERR, &
                        WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AFP(*),AP(*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPSV(UPLO,N,NRHS,AP,IPIV,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            AP(*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPSVX(FACT,UPLO,N,NRHS,AP,AFP,IPIV,B,LDB,X,LDX,RCOND, &
                        FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          FACT,UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AFP(*),AP(*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPTRD(UPLO,N,AP,D,E,TAU,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      REAL               D(*),E(*)
      COMPLEX            AP(*),TAU(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPTRF(UPLO,N,AP,IPIV,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      INTEGER            IPIV(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPTRI(UPLO,N,AP,IPIV,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      INTEGER            IPIV(*)
      COMPLEX            AP(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHPTRS(UPLO,N,NRHS,AP,IPIV,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            AP(*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHSEIN(SIDE,EIGSRC,INITV,SELECT,N,H,LDH,W,VL,LDVL,VR, &
                        LDVR,MM,M,WORK,RWORK,IFAILL,IFAILR,INFO)
      CHARACTER          EIGSRC,INITV,SIDE
      INTEGER            INFO,LDH,LDVL,LDVR,M,MM,N
      LOGICAL            SELECT(*)
      INTEGER            IFAILL(*),IFAILR(*)
      REAL               RWORK(*)
      COMPLEX            H(LDH,*),VL(LDVL,*),VR(LDVR,*),W(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CHSEQR(JOB,COMPZ,N,ILO,IHI,H,LDH,W,Z,LDZ,WORK,LWORK, &
                        INFO)
      CHARACTER          COMPZ,JOB
      INTEGER            IHI,ILO,INFO,LDH,LDZ,LWORK,N
      COMPLEX            H(LDH,*),W(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLABRD(M,N,NB,A,LDA,D,E,TAUQ,TAUP,X,LDX,Y,LDY)
      INTEGER            LDA,LDX,LDY,M,N,NB
      REAL               D(*),E(*)
      COMPLEX            A(LDA,*),TAUP(*),TAUQ(*),X(LDX,*),Y(LDY,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLACGV(N,X,INCX)
      INTEGER            INCX,N
      COMPLEX            X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLACN2(N,V,X,EST,KASE,ISAVE)
      INTEGER            KASE,N
      REAL               EST
      INTEGER            ISAVE(3)
      COMPLEX            V(*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLACON(N,V,X,EST,KASE)
      INTEGER            KASE,N
      REAL               EST
      COMPLEX            V(N),X(N)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLACP2(UPLO,M,N,A,LDA,B,LDB)
      CHARACTER          UPLO
      INTEGER            LDA,LDB,M,N
      REAL               A(LDA,*)
      COMPLEX            B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLACPY(UPLO,M,N,A,LDA,B,LDB)
      CHARACTER          UPLO
      INTEGER            LDA,LDB,M,N
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLACRM(M,N,A,LDA,B,LDB,C,LDC,RWORK)
      INTEGER            LDA,LDB,LDC,M,N
      REAL               B(LDB,*),RWORK(*)
      COMPLEX            A(LDA,*),C(LDC,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLACRT(N,CX,INCX,CY,INCY,C,S)
      INTEGER            INCX,INCY,N
      COMPLEX            C,S
      COMPLEX            CX(*),CY(*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLADIV(X,Y)
      COMPLEX            CLADIV
      COMPLEX            X,Y
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAED0(QSIZ,N,D,E,Q,LDQ,QSTORE,LDQS,RWORK,IWORK,INFO)
      INTEGER            INFO,LDQ,LDQS,N,QSIZ
      INTEGER            IWORK(*)
      REAL               D(*),E(*),RWORK(*)
      COMPLEX            Q(LDQ,*),QSTORE(LDQS,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAED7(N,CUTPNT,QSIZ,TLVLS,CURLVL,CURPBM,D,Q,LDQ,RHO, &
                        INDXQ,QSTORE,QPTR,PRMPTR,PERM,GIVPTR,GIVCOL, &
                        GIVNUM,WORK,RWORK,IWORK,INFO)
      INTEGER            CURLVL,CURPBM,CUTPNT,INFO,LDQ,N,QSIZ,TLVLS
      REAL               RHO
      INTEGER            GIVCOL(2,*),GIVPTR(*),INDXQ(*),IWORK(*), &
                         PERM(*),PRMPTR(*),QPTR(*)
      REAL               D(*),GIVNUM(2,*),QSTORE(*),RWORK(*)
      COMPLEX            Q(LDQ,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAED8(K,N,QSIZ,Q,LDQ,D,RHO,CUTPNT,Z,DLAMDA,Q2,LDQ2,W, &
                        INDXP,INDX,INDXQ,PERM,GIVPTR,GIVCOL,GIVNUM,INFO)
      INTEGER            CUTPNT,GIVPTR,INFO,K,LDQ,LDQ2,N,QSIZ
      REAL               RHO
      INTEGER            GIVCOL(2,*),INDX(*),INDXP(*),INDXQ(*),PERM(*)
      REAL               D(*),DLAMDA(*),GIVNUM(2,*),W(*),Z(*)
      COMPLEX            Q(LDQ,*),Q2(LDQ2,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAEIN(RIGHTV,NOINIT,N,H,LDH,W,V,B,LDB,RWORK,EPS3, &
                        SMLNUM,INFO)
      LOGICAL            NOINIT,RIGHTV
      INTEGER            INFO,LDB,LDH,N
      REAL               EPS3,SMLNUM
      COMPLEX            W
      REAL               RWORK(*)
      COMPLEX            B(LDB,*),H(LDH,*),V(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAESY(A,B,C,RT1,RT2,EVSCAL,CS1,SN1)
      COMPLEX            A,B,C,CS1,EVSCAL,RT1,RT2,SN1
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAEV2(A,B,C,RT1,RT2,CS1,SN1)
      REAL               CS1,RT1,RT2
      COMPLEX            A,B,C,SN1
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAG2Z(M,N,SA,LDSA,A,LDA,INFO)
      INTEGER            INFO,LDA,LDSA,M,N
      COMPLEX            SA(LDSA,*)
      COMPLEX*16         A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAGS2(UPPER,A1,A2,A3,B1,B2,B3,CSU,SNU,CSV,SNV,CSQ,SNQ)
      LOGICAL            UPPER
      REAL               A1,A3,B1,B3,CSQ,CSU,CSV
      COMPLEX            A2,B2,SNQ,SNU,SNV
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAGTM(TRANS,N,NRHS,ALPHA,DL,D,DU,X,LDX,BETA,B,LDB)
      CHARACTER          TRANS
      INTEGER            LDB,LDX,N,NRHS
      REAL               ALPHA,BETA
      COMPLEX            B(LDB,*),D(*),DL(*),DU(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAHEF(UPLO,N,NB,KB,A,LDA,IPIV,W,LDW,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KB,LDA,LDW,N,NB
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),W(LDW,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAHQR(WANTT,WANTZ,N,ILO,IHI,H,LDH,W,ILOZ,IHIZ,Z,LDZ, &
                        INFO)
      LOGICAL            WANTT,WANTZ
      INTEGER            IHI,IHIZ,ILO,ILOZ,INFO,LDH,LDZ,N
      COMPLEX            H(LDH,*),W(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAHR2(N,K,NB,A,LDA,TAU,T,LDT,Y,LDY)
      INTEGER            K,LDA,LDT,LDY,N,NB
      COMPLEX            A(LDA,*),T(LDT,NB),TAU(NB),Y(LDY,NB)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAHRD(N,K,NB,A,LDA,TAU,T,LDT,Y,LDY)
      INTEGER            K,LDA,LDT,LDY,N,NB
      COMPLEX            A(LDA,*),T(LDT,NB),TAU(NB),Y(LDY,NB)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAIC1(JOB,J,X,SEST,W,GAMMA,SESTPR,S,C)
      INTEGER            J,JOB
      REAL               SEST,SESTPR
      COMPLEX            C,GAMMA,S
      COMPLEX            W(J),X(J)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLALS0(ICOMPQ,NL,NR,SQRE,NRHS,B,LDB,BX,LDBX,PERM, &
                        GIVPTR,GIVCOL,LDGCOL,GIVNUM,LDGNUM,POLES,DIFL, &
                        DIFR,Z,K,C,S,RWORK,INFO)
      INTEGER            GIVPTR,ICOMPQ,INFO,K,LDB,LDBX,LDGCOL,LDGNUM,NL, &
                         NR,NRHS,SQRE
      REAL               C,S
      INTEGER            GIVCOL(LDGCOL,*),PERM(*)
      REAL               DIFL(*),DIFR(LDGNUM,*),GIVNUM(LDGNUM,*), &
                         POLES(LDGNUM,*),RWORK(*),Z(*)
      COMPLEX            B(LDB,*),BX(LDBX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLALSA(ICOMPQ,SMLSIZ,N,NRHS,B,LDB,BX,LDBX,U,LDU,VT,K, &
                        DIFL,DIFR,Z,POLES,GIVPTR,GIVCOL,LDGCOL,PERM, &
                        GIVNUM,C,S,RWORK,IWORK,INFO)
      INTEGER            ICOMPQ,INFO,LDB,LDBX,LDGCOL,LDU,N,NRHS,SMLSIZ
      INTEGER            GIVCOL(LDGCOL,*),GIVPTR(*),IWORK(*),K(*), &
                         PERM(LDGCOL,*)
      REAL               C(*),DIFL(LDU,*),DIFR(LDU,*),GIVNUM(LDU,*), &
                         POLES(LDU,*),RWORK(*),S(*),U(LDU,*),VT(LDU,*), &
                         Z(LDU,*)
      COMPLEX            B(LDB,*),BX(LDBX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLALSD(UPLO,SMLSIZ,N,NRHS,D,E,B,LDB,RCOND,RANK,WORK, &
                        RWORK,IWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS,RANK,SMLSIZ
      REAL               RCOND
      INTEGER            IWORK(*)
      REAL               D(*),E(*),RWORK(*)
      COMPLEX            B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANGB(NORM,N,KL,KU,AB,LDAB,WORK)
      REAL               CLANGB
      CHARACTER          NORM
      INTEGER            KL,KU,LDAB,N
      REAL               WORK(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANGE(NORM,M,N,A,LDA,WORK)
      REAL               CLANGE
      CHARACTER          NORM
      INTEGER            LDA,M,N
      REAL               WORK(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANGT(NORM,N,DL,D,DU)
      REAL               CLANGT
      CHARACTER          NORM
      INTEGER            N
      COMPLEX            D(*),DL(*),DU(*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANHB(NORM,UPLO,N,K,AB,LDAB,WORK)
      REAL               CLANHB
      CHARACTER          NORM,UPLO
      INTEGER            K,LDAB,N
      REAL               WORK(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANHE(NORM,UPLO,N,A,LDA,WORK)
      REAL               CLANHE
      CHARACTER          NORM,UPLO
      INTEGER            LDA,N
      REAL               WORK(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANHF(NORM,TRANSR,UPLO,N,A,WORK)
      REAL               CLANHF
      CHARACTER          NORM,TRANSR,UPLO
      INTEGER            N
      REAL               WORK(0:*)
      COMPLEX            A(0:*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANHP(NORM,UPLO,N,AP,WORK)
      REAL               CLANHP
      CHARACTER          NORM,UPLO
      INTEGER            N
      REAL               WORK(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANHS(NORM,N,A,LDA,WORK)
      REAL               CLANHS
      CHARACTER          NORM
      INTEGER            LDA,N
      REAL               WORK(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANHT(NORM,N,D,E)
      REAL               CLANHT
      CHARACTER          NORM
      INTEGER            N
      REAL               D(*)
      COMPLEX            E(*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANSB(NORM,UPLO,N,K,AB,LDAB,WORK)
      REAL               CLANSB
      CHARACTER          NORM,UPLO
      INTEGER            K,LDAB,N
      REAL               WORK(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANSP(NORM,UPLO,N,AP,WORK)
      REAL               CLANSP
      CHARACTER          NORM,UPLO
      INTEGER            N
      REAL               WORK(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANSY(NORM,UPLO,N,A,LDA,WORK)
      REAL               CLANSY
      CHARACTER          NORM,UPLO
      INTEGER            LDA,N
      REAL               WORK(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANTB(NORM,UPLO,DIAG,N,K,AB,LDAB,WORK)
      REAL               CLANTB
      CHARACTER          DIAG,NORM,UPLO
      INTEGER            K,LDAB,N
      REAL               WORK(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANTP(NORM,UPLO,DIAG,N,AP,WORK)
      REAL               CLANTP
      CHARACTER          DIAG,NORM,UPLO
      INTEGER            N
      REAL               WORK(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      FUNCTION CLANTR(NORM,UPLO,DIAG,M,N,A,LDA,WORK)
      REAL               CLANTR
      CHARACTER          DIAG,NORM,UPLO
      INTEGER            LDA,M,N
      REAL               WORK(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAPLL(N,X,INCX,Y,INCY,SSMIN)
      INTEGER            INCX,INCY,N
      REAL               SSMIN
      COMPLEX            X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAPMT(FORWRD,M,N,X,LDX,K)
      LOGICAL            FORWRD
      INTEGER            LDX,M,N
      INTEGER            K(*)
      COMPLEX            X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQGB(M,N,KL,KU,AB,LDAB,R,C,ROWCND,COLCND,AMAX,EQUED)
      CHARACTER          EQUED
      INTEGER            KL,KU,LDAB,M,N
      REAL               AMAX,COLCND,ROWCND
      REAL               C(*),R(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQGE(M,N,A,LDA,R,C,ROWCND,COLCND,AMAX,EQUED)
      CHARACTER          EQUED
      INTEGER            LDA,M,N
      REAL               AMAX,COLCND,ROWCND
      REAL               C(*),R(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQHB(UPLO,N,KD,AB,LDAB,S,SCOND,AMAX,EQUED)
      CHARACTER          EQUED,UPLO
      INTEGER            KD,LDAB,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQHE(UPLO,N,A,LDA,S,SCOND,AMAX,EQUED)
      CHARACTER          EQUED,UPLO
      INTEGER            LDA,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQHP(UPLO,N,AP,S,SCOND,AMAX,EQUED)
      CHARACTER          EQUED,UPLO
      INTEGER            N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQP2(M,N,OFFSET,A,LDA,JPVT,TAU,VN1,VN2,WORK)
      INTEGER            LDA,M,N,OFFSET
      INTEGER            JPVT(*)
      REAL               VN1(*),VN2(*)
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQPS(M,N,OFFSET,NB,KB,A,LDA,JPVT,TAU,VN1,VN2,AUXV,F, &
                        LDF)
      INTEGER            KB,LDA,LDF,M,N,NB,OFFSET
      INTEGER            JPVT(*)
      REAL               VN1(*),VN2(*)
      COMPLEX            A(LDA,*),AUXV(*),F(LDF,*),TAU(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQR0(WANTT,WANTZ,N,ILO,IHI,H,LDH,W,ILOZ,IHIZ,Z,LDZ, &
                        WORK,LWORK,INFO)
      LOGICAL            WANTT,WANTZ
      INTEGER            IHI,IHIZ,ILO,ILOZ,INFO,LDH,LDZ,LWORK,N
      COMPLEX            H(LDH,*),W(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQR1(N,H,LDH,S1,S2,V)
      INTEGER            LDH,N
      COMPLEX            S1,S2
      COMPLEX            H(LDH,*),V(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQR2(WANTT,WANTZ,N,KTOP,KBOT,NW,H,LDH,ILOZ,IHIZ,Z, &
                        LDZ,NS,ND,SH,V,LDV,NH,T,LDT,NV,WV,LDWV,WORK, &
                        LWORK)
      LOGICAL            WANTT,WANTZ
      INTEGER            IHIZ,ILOZ,KBOT,KTOP,LDH,LDT,LDV,LDWV,LDZ,LWORK, &
                         N,ND,NH,NS,NV,NW
      COMPLEX            H(LDH,*),SH(*),T(LDT,*),V(LDV,*),WORK(*), &
                         WV(LDWV,*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQR3(WANTT,WANTZ,N,KTOP,KBOT,NW,H,LDH,ILOZ,IHIZ,Z, &
                        LDZ,NS,ND,SH,V,LDV,NH,T,LDT,NV,WV,LDWV,WORK, &
                        LWORK)
      LOGICAL            WANTT,WANTZ
      INTEGER            IHIZ,ILOZ,KBOT,KTOP,LDH,LDT,LDV,LDWV,LDZ,LWORK, &
                         N,ND,NH,NS,NV,NW
      COMPLEX            H(LDH,*),SH(*),T(LDT,*),V(LDV,*),WORK(*), &
                         WV(LDWV,*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQR4(WANTT,WANTZ,N,ILO,IHI,H,LDH,W,ILOZ,IHIZ,Z,LDZ, &
                        WORK,LWORK,INFO)
      LOGICAL            WANTT,WANTZ
      INTEGER            IHI,IHIZ,ILO,ILOZ,INFO,LDH,LDZ,LWORK,N
      COMPLEX            H(LDH,*),W(*),WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQR5(WANTT,WANTZ,KACC22,N,KTOP,KBOT,NSHFTS,S,H,LDH, &
                        ILOZ,IHIZ,Z,LDZ,V,LDV,U,LDU,NV,WV,LDWV,NH,WH, &
                        LDWH)
      LOGICAL            WANTT,WANTZ
      INTEGER            IHIZ,ILOZ,KACC22,KBOT,KTOP,LDH,LDU,LDV,LDWH, &
                         LDWV,LDZ,N,NH,NSHFTS,NV
      COMPLEX            H(LDH,*),S(*),U(LDU,*),V(LDV,*),WH(LDWH,*), &
                         WV(LDWV,*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQSB(UPLO,N,KD,AB,LDAB,S,SCOND,AMAX,EQUED)
      CHARACTER          EQUED,UPLO
      INTEGER            KD,LDAB,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQSP(UPLO,N,AP,S,SCOND,AMAX,EQUED)
      CHARACTER          EQUED,UPLO
      INTEGER            N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAQSY(UPLO,N,A,LDA,S,SCOND,AMAX,EQUED)
      CHARACTER          EQUED,UPLO
      INTEGER            LDA,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAR1V(N,B1,BN,LAMBDA,D,L,LD,LLD,PIVMIN,GAPTOL,Z, &
                        WANTNC,NEGCNT,ZTZ,MINGMA,R,ISUPPZ,NRMINV,RESID, &
                        RQCORR,WORK)
      LOGICAL            WANTNC
      INTEGER            B1,BN,N,NEGCNT,R
      REAL               GAPTOL,LAMBDA,MINGMA,NRMINV,PIVMIN,RESID, &
                         RQCORR,ZTZ
      INTEGER            ISUPPZ(*)
      REAL               D(*),L(*),LD(*),LLD(*),WORK(*)
      COMPLEX            Z(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAR2V(N,X,Y,Z,INCX,C,S,INCC)
      INTEGER            INCC,INCX,N
      REAL               C(*)
      COMPLEX            S(*),X(*),Y(*),Z(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARCM(M,N,A,LDA,B,LDB,C,LDC,RWORK)
      INTEGER            LDA,LDB,LDC,M,N
      REAL               A(LDA,*),RWORK(*)
      COMPLEX            B(LDB,*),C(LDC,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARFB(SIDE,TRANS,DIRECT,STOREV,M,N,K,V,LDV,T,LDT,C, &
                        LDC,WORK,LDWORK)
      CHARACTER          DIRECT,SIDE,STOREV,TRANS
      INTEGER            K,LDC,LDT,LDV,LDWORK,M,N
      COMPLEX            C(LDC,*),T(LDT,*),V(LDV,*),WORK(LDWORK,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARF(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
      CHARACTER          SIDE
      INTEGER            INCV,LDC,M,N
      COMPLEX            TAU
      COMPLEX            C(LDC,*),V(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARFG(N,ALPHA,X,INCX,TAU)
      INTEGER            INCX,N
      COMPLEX            ALPHA,TAU
      COMPLEX            X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARFP(N,ALPHA,X,INCX,TAU)
      INTEGER            INCX,N
      COMPLEX            ALPHA,TAU
      COMPLEX            X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARFT(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
      CHARACTER          DIRECT,STOREV
      INTEGER            K,LDT,LDV,N
      COMPLEX            T(LDT,*),TAU(*),V(LDV,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARFX(SIDE,M,N,V,TAU,C,LDC,WORK)
      CHARACTER          SIDE
      INTEGER            LDC,M,N
      COMPLEX            TAU
      COMPLEX            C(LDC,*),V(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARGV(N,X,INCX,Y,INCY,C,INCC)
      INTEGER            INCC,INCX,INCY,N
      REAL               C(*)
      COMPLEX            X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARNV(IDIST,ISEED,N,X)
      INTEGER            IDIST,N
      INTEGER            ISEED(4)
      COMPLEX            X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARRV(N,VL,VU,D,L,PIVMIN,ISPLIT,M,DOL,DOU,MINRGP, &
                        RTOL1,RTOL2,W,WERR,WGAP,IBLOCK,INDEXW,GERS,Z, &
                        LDZ,ISUPPZ,WORK,IWORK,INFO)
      INTEGER            DOL,DOU,INFO,LDZ,M,N
      REAL               MINRGP,PIVMIN,RTOL1,RTOL2,VL,VU
      INTEGER            IBLOCK(*),INDEXW(*),ISPLIT(*),ISUPPZ(*), &
                         IWORK(*)
      REAL               D(*),GERS(*),L(*),W(*),WERR(*),WGAP(*),WORK(*)
      COMPLEX            Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARTG(F,G,CS,SN,R)
      REAL               CS
      COMPLEX            F,G,R,SN
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARTV(N,X,INCX,Y,INCY,C,S,INCC)
      INTEGER            INCC,INCX,INCY,N
      REAL               C(*)
      COMPLEX            S(*),X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARZB(SIDE,TRANS,DIRECT,STOREV,M,N,K,L,V,LDV,T,LDT,C, &
                        LDC,WORK,LDWORK)
      CHARACTER          DIRECT,SIDE,STOREV,TRANS
      INTEGER            K,L,LDC,LDT,LDV,LDWORK,M,N
      COMPLEX            C(LDC,*),T(LDT,*),V(LDV,*),WORK(LDWORK,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARZ(SIDE,M,N,L,V,INCV,TAU,C,LDC,WORK)
      CHARACTER          SIDE
      INTEGER            INCV,L,LDC,M,N
      COMPLEX            TAU
      COMPLEX            C(LDC,*),V(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLARZT(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
      CHARACTER          DIRECT,STOREV
      INTEGER            K,LDT,LDV,N
      COMPLEX            T(LDT,*),TAU(*),V(LDV,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLASCL(TYPE,KL,KU,CFROM,CTO,M,N,A,LDA,INFO)
      CHARACTER          TYPE
      INTEGER            INFO,KL,KU,LDA,M,N
      REAL               CFROM,CTO
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLASET(UPLO,M,N,ALPHA,BETA,A,LDA)
      CHARACTER          UPLO
      INTEGER            LDA,M,N
      COMPLEX            ALPHA,BETA
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLASR(SIDE,PIVOT,DIRECT,M,N,C,S,A,LDA)
      CHARACTER          DIRECT,PIVOT,SIDE
      INTEGER            LDA,M,N
      REAL               C(*),S(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLASSQ(N,X,INCX,SCALE,SUMSQ)
      INTEGER            INCX,N
      REAL               SCALE,SUMSQ
      COMPLEX            X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLASWP(N,A,LDA,K1,K2,IPIV,INCX)
      INTEGER            INCX,K1,K2,LDA,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLASYF(UPLO,N,NB,KB,A,LDA,IPIV,W,LDW,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KB,LDA,LDW,N,NB
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),W(LDW,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLATBS(UPLO,TRANS,DIAG,NORMIN,N,KD,AB,LDAB,X,SCALE, &
                        CNORM,INFO)
      CHARACTER          DIAG,NORMIN,TRANS,UPLO
      INTEGER            INFO,KD,LDAB,N
      REAL               SCALE
      REAL               CNORM(*)
      COMPLEX            AB(LDAB,*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLATDF(IJOB,N,Z,LDZ,RHS,RDSUM,RDSCAL,IPIV,JPIV)
      INTEGER            IJOB,LDZ,N
      REAL               RDSCAL,RDSUM
      INTEGER            IPIV(*),JPIV(*)
      COMPLEX            RHS(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLATPS(UPLO,TRANS,DIAG,NORMIN,N,AP,X,SCALE,CNORM,INFO)
      CHARACTER          DIAG,NORMIN,TRANS,UPLO
      INTEGER            INFO,N
      REAL               SCALE
      REAL               CNORM(*)
      COMPLEX            AP(*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLATRD(UPLO,N,NB,A,LDA,E,TAU,W,LDW)
      CHARACTER          UPLO
      INTEGER            LDA,LDW,N,NB
      REAL               E(*)
      COMPLEX            A(LDA,*),TAU(*),W(LDW,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLATRS(UPLO,TRANS,DIAG,NORMIN,N,A,LDA,X,SCALE,CNORM, &
                        INFO)
      CHARACTER          DIAG,NORMIN,TRANS,UPLO
      INTEGER            INFO,LDA,N
      REAL               SCALE
      REAL               CNORM(*)
      COMPLEX            A(LDA,*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLATRZ(M,N,L,A,LDA,TAU,WORK)
      INTEGER            L,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLATZM(SIDE,M,N,V,INCV,TAU,C1,C2,LDC,WORK)
      CHARACTER          SIDE
      INTEGER            INCV,LDC,M,N
      COMPLEX            TAU
      COMPLEX            C1(LDC,*),C2(LDC,*),V(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAUU2(UPLO,N,A,LDA,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CLAUUM(UPLO,N,A,LDA,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBCON(UPLO,N,KD,AB,LDAB,ANORM,RCOND,WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,N
      REAL               ANORM,RCOND
      REAL               RWORK(*)
      COMPLEX            AB(LDAB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBEQU(UPLO,N,KD,AB,LDAB,S,SCOND,AMAX,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBRFS(UPLO,N,KD,NRHS,AB,LDAB,AFB,LDAFB,B,LDB,X,LDX, &
                        FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,LDAFB,LDB,LDX,N,NRHS
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),WORK(*),X(LDX, &
                         *)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBSTF(UPLO,N,KD,AB,LDAB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,N
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBSV(UPLO,N,KD,NRHS,AB,LDAB,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,LDB,N,NRHS
      COMPLEX            AB(LDAB,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBSVX(FACT,UPLO,N,KD,NRHS,AB,LDAB,AFB,LDAFB,EQUED,S,B, &
                        LDB,X,LDX,RCOND,FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,UPLO
      INTEGER            INFO,KD,LDAB,LDAFB,LDB,LDX,N,NRHS
      REAL               RCOND
      REAL               BERR(*),FERR(*),RWORK(*),S(*)
      COMPLEX            AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),WORK(*),X(LDX, &
                         *)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBTF2(UPLO,N,KD,AB,LDAB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,N
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBTRF(UPLO,N,KD,AB,LDAB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,N
      COMPLEX            AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPBTRS(UPLO,N,KD,NRHS,AB,LDAB,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,KD,LDAB,LDB,N,NRHS
      COMPLEX            AB(LDAB,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPFTRF(TRANSR,UPLO,N,A,INFO)
      CHARACTER          TRANSR,UPLO
      INTEGER            INFO,N
      COMPLEX            A(0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPFTRI(TRANSR,UPLO,N,A,INFO)
      CHARACTER          TRANSR,UPLO
      INTEGER            INFO,N
      COMPLEX            A(0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPFTRS(TRANSR,UPLO,N,NRHS,A,B,LDB,INFO)
      CHARACTER          TRANSR,UPLO
      INTEGER            INFO,LDB,N,NRHS
      COMPLEX            A(0:*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOCON(UPLO,N,A,LDA,ANORM,RCOND,WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      REAL               ANORM,RCOND
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOEQUB(N,A,LDA,S,SCOND,AMAX,INFO)
      INTEGER            INFO,LDA,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOEQU(N,A,LDA,S,SCOND,AMAX,INFO)
      INTEGER            INFO,LDA,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPORFS(UPLO,N,NRHS,A,LDA,AF,LDAF,B,LDB,X,LDX,FERR,BERR, &
                        WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPORFSX(UPLO,EQUED,N,NRHS,A,LDA,AF,LDAF,S,B,LDB,X,LDX, &
                         RCOND,BERR,N_ERR_BNDS,ERR_BNDS_NORM, &
                         ERR_BNDS_COMP,NPARAMS,PARAMS,WORK,RWORK,INFO)
      CHARACTER          EQUED,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND
      REAL               BERR(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),RWORK(*),S(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDB,N,NRHS
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOSVX(FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,EQUED,S,B,LDB,X, &
                        LDX,RCOND,FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      REAL               RCOND
      REAL               BERR(*),FERR(*),RWORK(*),S(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOSVXX(FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,EQUED,S,B,LDB,X, &
                         LDX,RCOND,RPVGRW,BERR,N_ERR_BNDS,ERR_BNDS_NORM, &
                         ERR_BNDS_COMP,NPARAMS,PARAMS,WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND,RPVGRW
      REAL               BERR(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),RWORK(*),S(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOTF2(UPLO,N,A,LDA,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOTRF(UPLO,N,A,LDA,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOTRI(UPLO,N,A,LDA,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDB,N,NRHS
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPCON(UPLO,N,AP,ANORM,RCOND,WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      REAL               ANORM,RCOND
      REAL               RWORK(*)
      COMPLEX            AP(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPEQU(UPLO,N,AP,S,SCOND,AMAX,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPRFS(UPLO,N,NRHS,AP,AFP,B,LDB,X,LDX,FERR,BERR,WORK, &
                        RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AFP(*),AP(*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPSV(UPLO,N,NRHS,AP,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS
      COMPLEX            AP(*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPSVX(FACT,UPLO,N,NRHS,AP,AFP,EQUED,S,B,LDB,X,LDX, &
                        RCOND,FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               RCOND
      REAL               BERR(*),FERR(*),RWORK(*),S(*)
      COMPLEX            AFP(*),AP(*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPTRF(UPLO,N,AP,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPTRI(UPLO,N,AP,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPPTRS(UPLO,N,NRHS,AP,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS
      COMPLEX            AP(*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPSTF2(UPLO,N,A,LDA,PIV,RANK,TOL,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N,RANK
      REAL               TOL
      INTEGER            PIV(N)
      REAL               WORK(2*N)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPSTRF(UPLO,N,A,LDA,PIV,RANK,TOL,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N,RANK
      REAL               TOL
      INTEGER            PIV(N)
      REAL               WORK(2*N)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTCON(N,D,E,ANORM,RCOND,RWORK,INFO)
      INTEGER            INFO,N
      REAL               ANORM,RCOND
      REAL               D(*),RWORK(*)
      COMPLEX            E(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTEQR(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
      CHARACTER          COMPZ
      INTEGER            INFO,LDZ,N
      REAL               D(*),E(*),WORK(*)
      COMPLEX            Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTRFS(UPLO,N,NRHS,D,E,DF,EF,B,LDB,X,LDX,FERR,BERR, &
                        WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               BERR(*),D(*),DF(*),FERR(*),RWORK(*)
      COMPLEX            B(LDB,*),E(*),EF(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTSV(N,NRHS,D,E,B,LDB,INFO)
      INTEGER            INFO,LDB,N,NRHS
      REAL               D(*)
      COMPLEX            B(LDB,*),E(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTSVX(FACT,N,NRHS,D,E,DF,EF,B,LDB,X,LDX,RCOND,FERR, &
                        BERR,WORK,RWORK,INFO)
      CHARACTER          FACT
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               RCOND
      REAL               BERR(*),D(*),DF(*),FERR(*),RWORK(*)
      COMPLEX            B(LDB,*),E(*),EF(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTTRF(N,D,E,INFO)
      INTEGER            INFO,N
      REAL               D(*)
      COMPLEX            E(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTTRS(UPLO,N,NRHS,D,E,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS
      REAL               D(*)
      COMPLEX            B(LDB,*),E(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CPTTS2(IUPLO,N,NRHS,D,E,B,LDB)
      INTEGER            IUPLO,LDB,N,NRHS
      REAL               D(*)
      COMPLEX            B(LDB,*),E(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CROT(N,CX,INCX,CY,INCY,C,S)
      INTEGER            INCX,INCY,N
      REAL               C
      COMPLEX            S
      COMPLEX            CX(*),CY(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPCON(UPLO,N,AP,IPIV,ANORM,RCOND,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      REAL               ANORM,RCOND
      INTEGER            IPIV(*)
      COMPLEX            AP(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      CHARACTER          UPLO
      INTEGER            INCX,INCY,N
      COMPLEX            ALPHA,BETA
      COMPLEX            AP(*),X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPR(UPLO,N,ALPHA,X,INCX,AP)
      CHARACTER          UPLO
      INTEGER            INCX,N
      COMPLEX            ALPHA
      COMPLEX            AP(*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPRFS(UPLO,N,NRHS,AP,AFP,IPIV,B,LDB,X,LDX,FERR,BERR, &
                        WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AFP(*),AP(*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPSV(UPLO,N,NRHS,AP,IPIV,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            AP(*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPSVX(FACT,UPLO,N,NRHS,AP,AFP,IPIV,B,LDB,X,LDX,RCOND, &
                        FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          FACT,UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AFP(*),AP(*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPTRF(UPLO,N,AP,IPIV,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      INTEGER            IPIV(*)
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPTRI(UPLO,N,AP,IPIV,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,N
      INTEGER            IPIV(*)
      COMPLEX            AP(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSPTRS(UPLO,N,NRHS,AP,IPIV,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            AP(*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSRSCL(N,SA,SX,INCX)
      INTEGER            INCX,N
      REAL               SA
      COMPLEX            SX(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSTEDC(COMPZ,N,D,E,Z,LDZ,WORK,LWORK,RWORK,LRWORK,IWORK, &
                        LIWORK,INFO)
      CHARACTER          COMPZ
      INTEGER            INFO,LDZ,LIWORK,LRWORK,LWORK,N
      INTEGER            IWORK(*)
      REAL               D(*),E(*),RWORK(*)
      COMPLEX            WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSTEGR(JOBZ,RANGE,N,D,E,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ, &
                        ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,RANGE
      INTEGER            IL,INFO,IU,LDZ,LIWORK,LWORK,M,N
      REAL               ABSTOL,VL,VU
      INTEGER            ISUPPZ(*),IWORK(*)
      REAL               D(*),E(*),W(*),WORK(*)
      COMPLEX            Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAIL, &
                        INFO)
      INTEGER            INFO,LDZ,M,N
      INTEGER            IBLOCK(*),IFAIL(*),ISPLIT(*),IWORK(*)
      REAL               D(*),E(*),W(*),WORK(*)
      COMPLEX            Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSTEMR(JOBZ,RANGE,N,D,E,VL,VU,IL,IU,M,W,Z,LDZ,NZC, &
                        ISUPPZ,TRYRAC,WORK,LWORK,IWORK,LIWORK,INFO)
      CHARACTER          JOBZ,RANGE
      LOGICAL            TRYRAC
      INTEGER            IL,INFO,IU,LDZ,LIWORK,LWORK,M,N,NZC
      REAL               VL,VU
      INTEGER            ISUPPZ(*),IWORK(*)
      REAL               D(*),E(*),W(*),WORK(*)
      COMPLEX            Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSTEQR(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
      CHARACTER          COMPZ
      INTEGER            INFO,LDZ,N
      REAL               D(*),E(*),WORK(*)
      COMPLEX            Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYCON(UPLO,N,A,LDA,IPIV,ANORM,RCOND,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      REAL               ANORM,RCOND
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYEQUB(UPLO,N,A,LDA,S,SCOND,AMAX,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      REAL               AMAX,SCOND
      REAL               S(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER          UPLO
      INTEGER            INCX,INCY,LDA,N
      COMPLEX            ALPHA,BETA
      COMPLEX            A(LDA,*),X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
      CHARACTER          UPLO
      INTEGER            INCX,LDA,N
      COMPLEX            ALPHA
      COMPLEX            A(LDA,*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYRFS(UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX,FERR, &
                        BERR,WORK,RWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYRFSX(UPLO,EQUED,N,NRHS,A,LDA,AF,LDAF,IPIV,S,B,LDB,X, &
                         LDX,RCOND,BERR,N_ERR_BNDS,ERR_BNDS_NORM, &
                         ERR_BNDS_COMP,NPARAMS,PARAMS,WORK,RWORK,INFO)
      CHARACTER          EQUED,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),RWORK(*),S(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYSV(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,LWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDB,LWORK,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYSVX(FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX, &
                        RCOND,FERR,BERR,WORK,LWORK,RWORK,INFO)
      CHARACTER          FACT,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,LWORK,N,NRHS
      REAL               RCOND
      INTEGER            IPIV(*)
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYSVXX(FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,S,B, &
                         LDB,X,LDX,RCOND,RPVGRW,BERR,N_ERR_BNDS, &
                         ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS,PARAMS, &
                         WORK,RWORK,INFO)
      CHARACTER          EQUED,FACT,UPLO
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      REAL               RCOND,RPVGRW
      INTEGER            IPIV(*)
      REAL               BERR(*),ERR_BNDS_COMP(NRHS,*), &
                         ERR_BNDS_NORM(NRHS,*),PARAMS(*),RWORK(*),S(*)
      COMPLEX            A(LDA,*),AF(LDAF,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYTF2(UPLO,N,A,LDA,IPIV,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYTRF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LWORK,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYTRI(UPLO,N,A,LDA,IPIV,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CSYTRS(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LDB,N,NRHS
      INTEGER            IPIV(*)
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTBCON(NORM,UPLO,DIAG,N,KD,AB,LDAB,RCOND,WORK,RWORK, &
                        INFO)
      CHARACTER          DIAG,NORM,UPLO
      INTEGER            INFO,KD,LDAB,N
      REAL               RCOND
      REAL               RWORK(*)
      COMPLEX            AB(LDAB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTBRFS(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,B,LDB,X,LDX, &
                        FERR,BERR,WORK,RWORK,INFO)
      CHARACTER          DIAG,TRANS,UPLO
      INTEGER            INFO,KD,LDAB,LDB,LDX,N,NRHS
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AB(LDAB,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTBTRS(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,B,LDB,INFO)
      CHARACTER          DIAG,TRANS,UPLO
      INTEGER            INFO,KD,LDAB,LDB,N,NRHS
      COMPLEX            AB(LDAB,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTFSM(TRANSR,SIDE,UPLO,TRANS,DIAG,M,N,ALPHA,A,B,LDB)
      CHARACTER          DIAG,SIDE,TRANS,TRANSR,UPLO
      INTEGER            LDB,M,N
      COMPLEX            ALPHA
      COMPLEX            A(0:*),B(0:LDB-1,0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTFTRI(TRANSR,UPLO,DIAG,N,A,INFO)
      CHARACTER          DIAG,TRANSR,UPLO
      INTEGER            INFO,N
      COMPLEX            A(0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTFTTP(TRANSR,UPLO,N,ARF,AP,INFO)
      CHARACTER          TRANSR,UPLO
      INTEGER            INFO,N
      COMPLEX            AP(0:*),ARF(0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTFTTR(TRANSR,UPLO,N,ARF,A,LDA,INFO)
      CHARACTER          TRANSR,UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(0:LDA-1,0:*),ARF(0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGEVC(SIDE,HOWMNY,SELECT,N,S,LDS,P,LDP,VL,LDVL,VR, &
                        LDVR,MM,M,WORK,RWORK,INFO)
      CHARACTER          HOWMNY,SIDE
      INTEGER            INFO,LDP,LDS,LDVL,LDVR,M,MM,N
      LOGICAL            SELECT(*)
      REAL               RWORK(*)
      COMPLEX            P(LDP,*),S(LDS,*),VL(LDVL,*),VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGEX2(WANTQ,WANTZ,N,A,LDA,B,LDB,Q,LDQ,Z,LDZ,J1,INFO)
      LOGICAL            WANTQ,WANTZ
      INTEGER            INFO,J1,LDA,LDB,LDQ,LDZ,N
      COMPLEX            A(LDA,*),B(LDB,*),Q(LDQ,*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGEXC(WANTQ,WANTZ,N,A,LDA,B,LDB,Q,LDQ,Z,LDZ,IFST,ILST, &
                        INFO)
      LOGICAL            WANTQ,WANTZ
      INTEGER            IFST,ILST,INFO,LDA,LDB,LDQ,LDZ,N
      COMPLEX            A(LDA,*),B(LDB,*),Q(LDQ,*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGSEN(IJOB,WANTQ,WANTZ,SELECT,N,A,LDA,B,LDB,ALPHA, &
                        BETA,Q,LDQ,Z,LDZ,M,PL,PR,DIF,WORK,LWORK,IWORK, &
                        LIWORK,INFO)
      LOGICAL            WANTQ,WANTZ
      INTEGER            IJOB,INFO,LDA,LDB,LDQ,LDZ,LIWORK,LWORK,M,N
      REAL               PL,PR
      LOGICAL            SELECT(*)
      INTEGER            IWORK(*)
      REAL               DIF(*)
      COMPLEX            A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),Q(LDQ,*), &
                         WORK(*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGSJA(JOBU,JOBV,JOBQ,M,P,N,K,L,A,LDA,B,LDB,TOLA,TOLB, &
                        ALPHA,BETA,U,LDU,V,LDV,Q,LDQ,WORK,NCYCLE,INFO)
      CHARACTER          JOBQ,JOBU,JOBV
      INTEGER            INFO,K,L,LDA,LDB,LDQ,LDU,LDV,M,N,NCYCLE,P
      REAL               TOLA,TOLB
      REAL               ALPHA(*),BETA(*)
      COMPLEX            A(LDA,*),B(LDB,*),Q(LDQ,*),U(LDU,*),V(LDV,*), &
                         WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGSNA(JOB,HOWMNY,SELECT,N,A,LDA,B,LDB,VL,LDVL,VR,LDVR, &
                        S,DIF,MM,M,WORK,LWORK,IWORK,INFO)
      CHARACTER          HOWMNY,JOB
      INTEGER            INFO,LDA,LDB,LDVL,LDVR,LWORK,M,MM,N
      LOGICAL            SELECT(*)
      INTEGER            IWORK(*)
      REAL               DIF(*),S(*)
      COMPLEX            A(LDA,*),B(LDB,*),VL(LDVL,*),VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGSY2(TRANS,IJOB,M,N,A,LDA,B,LDB,C,LDC,D,LDD,E,LDE,F, &
                        LDF,SCALE,RDSUM,RDSCAL,INFO)
      CHARACTER          TRANS
      INTEGER            IJOB,INFO,LDA,LDB,LDC,LDD,LDE,LDF,M,N
      REAL               RDSCAL,RDSUM,SCALE
      COMPLEX            A(LDA,*),B(LDB,*),C(LDC,*),D(LDD,*),E(LDE,*), &
                         F(LDF,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTGSYL(TRANS,IJOB,M,N,A,LDA,B,LDB,C,LDC,D,LDD,E,LDE,F, &
                        LDF,SCALE,DIF,WORK,LWORK,IWORK,INFO)
      CHARACTER          TRANS
      INTEGER            IJOB,INFO,LDA,LDB,LDC,LDD,LDE,LDF,LWORK,M,N
      REAL               DIF,SCALE
      INTEGER            IWORK(*)
      COMPLEX            A(LDA,*),B(LDB,*),C(LDC,*),D(LDD,*),E(LDE,*), &
                         F(LDF,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTPCON(NORM,UPLO,DIAG,N,AP,RCOND,WORK,RWORK,INFO)
      CHARACTER          DIAG,NORM,UPLO
      INTEGER            INFO,N
      REAL               RCOND
      REAL               RWORK(*)
      COMPLEX            AP(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTPRFS(UPLO,TRANS,DIAG,N,NRHS,AP,B,LDB,X,LDX,FERR,BERR, &
                        WORK,RWORK,INFO)
      CHARACTER          DIAG,TRANS,UPLO
      INTEGER            INFO,LDB,LDX,N,NRHS
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            AP(*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTPTRI(UPLO,DIAG,N,AP,INFO)
      CHARACTER          DIAG,UPLO
      INTEGER            INFO,N
      COMPLEX            AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTPTRS(UPLO,TRANS,DIAG,N,NRHS,AP,B,LDB,INFO)
      CHARACTER          DIAG,TRANS,UPLO
      INTEGER            INFO,LDB,N,NRHS
      COMPLEX            AP(*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTPTTF(TRANSR,UPLO,N,AP,ARF,INFO)
      CHARACTER          TRANSR,UPLO
      INTEGER            INFO,N
      COMPLEX            AP(0:*),ARF(0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTPTTR(UPLO,N,AP,A,LDA,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*),AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRCON(NORM,UPLO,DIAG,N,A,LDA,RCOND,WORK,RWORK,INFO)
      CHARACTER          DIAG,NORM,UPLO
      INTEGER            INFO,LDA,N
      REAL               RCOND
      REAL               RWORK(*)
      COMPLEX            A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTREVC(SIDE,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,MM,M, &
                        WORK,RWORK,INFO)
      CHARACTER          HOWMNY,SIDE
      INTEGER            INFO,LDT,LDVL,LDVR,M,MM,N
      LOGICAL            SELECT(*)
      REAL               RWORK(*)
      COMPLEX            T(LDT,*),VL(LDVL,*),VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTREXC(COMPQ,N,T,LDT,Q,LDQ,IFST,ILST,INFO)
      CHARACTER          COMPQ
      INTEGER            IFST,ILST,INFO,LDQ,LDT,N
      COMPLEX            Q(LDQ,*),T(LDT,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRRFS(UPLO,TRANS,DIAG,N,NRHS,A,LDA,B,LDB,X,LDX,FERR, &
                        BERR,WORK,RWORK,INFO)
      CHARACTER          DIAG,TRANS,UPLO
      INTEGER            INFO,LDA,LDB,LDX,N,NRHS
      REAL               BERR(*),FERR(*),RWORK(*)
      COMPLEX            A(LDA,*),B(LDB,*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRSEN(JOB,COMPQ,SELECT,N,T,LDT,Q,LDQ,W,M,S,SEP,WORK, &
                        LWORK,INFO)
      CHARACTER          COMPQ,JOB
      INTEGER            INFO,LDQ,LDT,LWORK,M,N
      REAL               S,SEP
      LOGICAL            SELECT(*)
      COMPLEX            Q(LDQ,*),T(LDT,*),W(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRSNA(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,S,SEP, &
                        MM,M,WORK,LDWORK,RWORK,INFO)
      CHARACTER          HOWMNY,JOB
      INTEGER            INFO,LDT,LDVL,LDVR,LDWORK,M,MM,N
      LOGICAL            SELECT(*)
      REAL               RWORK(*),S(*),SEP(*)
      COMPLEX            T(LDT,*),VL(LDVL,*),VR(LDVR,*),WORK(LDWORK,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRSYL(TRANA,TRANB,ISGN,M,N,A,LDA,B,LDB,C,LDC,SCALE, &
                        INFO)
      CHARACTER          TRANA,TRANB
      INTEGER            INFO,ISGN,LDA,LDB,LDC,M,N
      REAL               SCALE
      COMPLEX            A(LDA,*),B(LDB,*),C(LDC,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRTI2(UPLO,DIAG,N,A,LDA,INFO)
      CHARACTER          DIAG,UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRTRI(UPLO,DIAG,N,A,LDA,INFO)
      CHARACTER          DIAG,UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRTRS(UPLO,TRANS,DIAG,N,NRHS,A,LDA,B,LDB,INFO)
      CHARACTER          DIAG,TRANS,UPLO
      INTEGER            INFO,LDA,LDB,N,NRHS
      COMPLEX            A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRTTF(TRANSR,UPLO,N,A,LDA,ARF,INFO)
      CHARACTER          TRANSR,UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(0:LDA-1,0:*),ARF(0:*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTRTTP(UPLO,N,A,LDA,AP,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,N
      COMPLEX            A(LDA,*),AP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTZRQF(M,N,A,LDA,TAU,INFO)
      INTEGER            INFO,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CTZRZF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNG2L(M,N,K,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,K,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNG2R(M,N,K,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,K,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGBR(VECT,M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
      CHARACTER          VECT
      INTEGER            INFO,K,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGHR(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            IHI,ILO,INFO,LDA,LWORK,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGL2(M,N,K,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,K,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGLQ(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,K,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGQL(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,K,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGQR(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,K,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGR2(M,N,K,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,K,LDA,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGRQ(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,K,LDA,LWORK,M,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNGTR(UPLO,N,A,LDA,TAU,WORK,LWORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDA,LWORK,N
      COMPLEX            A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNM2L(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNM2R(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMBR(VECT,SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK, &
                        LWORK,INFO)
      CHARACTER          SIDE,TRANS,VECT
      INTEGER            INFO,K,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMHR(SIDE,TRANS,M,N,ILO,IHI,A,LDA,TAU,C,LDC,WORK, &
                        LWORK,INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            IHI,ILO,INFO,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNML2(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMLQ(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK, &
                        INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMQL(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK, &
                        INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMQR(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK, &
                        INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMR2(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMR3(SIDE,TRANS,M,N,K,L,A,LDA,TAU,C,LDC,WORK,INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,L,LDA,LDC,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMRQ(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK, &
                        INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMRZ(SIDE,TRANS,M,N,K,L,A,LDA,TAU,C,LDC,WORK,LWORK, &
                        INFO)
      CHARACTER          SIDE,TRANS
      INTEGER            INFO,K,L,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUNMTR(SIDE,UPLO,TRANS,M,N,A,LDA,TAU,C,LDC,WORK,LWORK, &
                        INFO)
      CHARACTER          SIDE,TRANS,UPLO
      INTEGER            INFO,LDA,LDC,LWORK,M,N
      COMPLEX            A(LDA,*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUPGTR(UPLO,N,AP,TAU,Q,LDQ,WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDQ,N
      COMPLEX            AP(*),Q(LDQ,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE CUPMTR(SIDE,UPLO,TRANS,M,N,AP,TAU,C,LDC,WORK,INFO)
      CHARACTER          SIDE,TRANS,UPLO
      INTEGER            INFO,LDC,M,N
      COMPLEX            AP(*),C(LDC,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DBDSDC(UPLO,COMPQ,N,D,E,U,LDU,VT,LDVT,Q,IQ,WORK,IWORK, &
                        INFO)
      CHARACTER          COMPQ,UPLO
      INTEGER            INFO,LDU,LDVT,N
      INTEGER            IQ(*),IWORK(*)
      DOUBLE PRECISION   D(*),E(*),Q(*),U(LDU,*),VT(LDVT,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DBDSQR(UPLO,N,NCVT,NRU,NCC,D,E,VT,LDVT,U,LDU,C,LDC, &
                        WORK,INFO)
      CHARACTER          UPLO
      INTEGER            INFO,LDC,LDU,LDVT,N,NCC,NCVT,NRU
      DOUBLE PRECISION   C(LDC,*),D(*),E(*),U(LDU,*),VT(LDVT,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DDISNA(JOB,M,N,D,SEP,INFO)
      CHARACTER          JOB
      INTEGER            INFO,M,N
      DOUBLE PRECISION   D(*),SEP(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBBRD(VECT,M,N,NCC,KL,KU,AB,LDAB,D,E,Q,LDQ,PT,LDPT,C, &
                        LDC,WORK,INFO)
      CHARACTER          VECT
      INTEGER            INFO,KL,KU,LDAB,LDC,LDPT,LDQ,M,N,NCC
      DOUBLE PRECISION   AB(LDAB,*),C(LDC,*),D(*),E(*),PT(LDPT,*),Q(LDQ, &
                         *),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBCON(NORM,N,KL,KU,AB,LDAB,IPIV,ANORM,RCOND,WORK, &
                        IWORK,INFO)
      CHARACTER          NORM
      INTEGER            INFO,KL,KU,LDAB,N
      DOUBLE PRECISION   ANORM,RCOND
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   AB(LDAB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBEQUB(M,N,KL,KU,AB,LDAB,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      DOUBLE PRECISION   AMAX,COLCND,ROWCND
      DOUBLE PRECISION   AB(LDAB,*),C(*),R(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBEQU(M,N,KL,KU,AB,LDAB,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      DOUBLE PRECISION   AMAX,COLCND,ROWCND
      DOUBLE PRECISION   AB(LDAB,*),C(*),R(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBRFS(TRANS,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB,IPIV,B,LDB, &
                        X,LDX,FERR,BERR,WORK,IWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NRHS
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),BERR(*), &
                         FERR(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBRFSX(TRANS,EQUED,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB, &
                         IPIV,R,C,B,LDB,X,LDX,RCOND,BERR,N_ERR_BNDS, &
                         ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS,PARAMS, &
                         WORK,IWORK,INFO)
      CHARACTER          EQUED,TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NPARAMS,NRHS, &
                         N_ERR_BNDS
      DOUBLE PRECISION   RCOND
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),BERR(*),C(*), &
                         ERR_BNDS_COMP(NRHS,*),ERR_BNDS_NORM(NRHS,*), &
                         PARAMS(*),R(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBSV(N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)
      INTEGER            INFO,KL,KU,LDAB,LDB,N,NRHS
      INTEGER            IPIV(*)
      DOUBLE PRECISION   AB(LDAB,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBSVX(FACT,TRANS,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB,IPIV, &
                        EQUED,R,C,B,LDB,X,LDX,RCOND,FERR,BERR,WORK, &
                        IWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NRHS
      DOUBLE PRECISION   RCOND
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),BERR(*),C(*), &
                         FERR(*),R(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBSVXX(FACT,TRANS,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB,IPIV, &
                         EQUED,R,C,B,LDB,X,LDX,RCOND,RPVGRW,BERR, &
                         N_ERR_BNDS,ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS, &
                         PARAMS,WORK,IWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,KL,KU,LDAB,LDAFB,LDB,LDX,N,NPARAMS,NRHS, &
                         N_ERR_BNDS
      DOUBLE PRECISION   RCOND,RPVGRW
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   AB(LDAB,*),AFB(LDAFB,*),B(LDB,*),BERR(*),C(*), &
                         ERR_BNDS_COMP(NRHS,*),ERR_BNDS_NORM(NRHS,*), &
                         PARAMS(*),R(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      INTEGER            IPIV(*)
      DOUBLE PRECISION   AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
      INTEGER            INFO,KL,KU,LDAB,M,N
      INTEGER            IPIV(*)
      DOUBLE PRECISION   AB(LDAB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGBTRS(TRANS,N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,KL,KU,LDAB,LDB,N,NRHS
      INTEGER            IPIV(*)
      DOUBLE PRECISION   AB(LDAB,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEBAK(JOB,SIDE,N,ILO,IHI,SCALE,M,V,LDV,INFO)
      CHARACTER          JOB,SIDE
      INTEGER            IHI,ILO,INFO,LDV,M,N
      DOUBLE PRECISION   SCALE(*),V(LDV,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEBAL(JOB,N,A,LDA,ILO,IHI,SCALE,INFO)
      CHARACTER          JOB
      INTEGER            IHI,ILO,INFO,LDA,N
      DOUBLE PRECISION   A(LDA,*),SCALE(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEBD2(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      DOUBLE PRECISION   A(LDA,*),D(*),E(*),TAUP(*),TAUQ(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEBRD(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      DOUBLE PRECISION   A(LDA,*),D(*),E(*),TAUP(*),TAUQ(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGECON(NORM,N,A,LDA,ANORM,RCOND,WORK,IWORK,INFO)
      CHARACTER          NORM
      INTEGER            INFO,LDA,N
      DOUBLE PRECISION   ANORM,RCOND
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEEQUB(M,N,A,LDA,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,LDA,M,N
      DOUBLE PRECISION   AMAX,COLCND,ROWCND
      DOUBLE PRECISION   A(LDA,*),C(*),R(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEEQU(M,N,A,LDA,R,C,ROWCND,COLCND,AMAX,INFO)
      INTEGER            INFO,LDA,M,N
      DOUBLE PRECISION   AMAX,COLCND,ROWCND
      DOUBLE PRECISION   A(LDA,*),C(*),R(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEES(JOBVS,SORT,SELECT,N,A,LDA,SDIM,WR,WI,VS,LDVS, &
                       WORK,LWORK,BWORK,INFO)
      CHARACTER          JOBVS,SORT
      INTEGER            INFO,LDA,LDVS,LWORK,N,SDIM
      LOGICAL            BWORK(*)
      DOUBLE PRECISION   A(LDA,*),VS(LDVS,*),WI(*),WORK(*),WR(*)
      LOGICAL            SELECT
      EXTERNAL           SELECT
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEESX(JOBVS,SORT,SELECT,SENSE,N,A,LDA,SDIM,WR,WI,VS, &
                        LDVS,RCONDE,RCONDV,WORK,LWORK,IWORK,LIWORK, &
                        BWORK,INFO)
      CHARACTER          JOBVS,SENSE,SORT
      INTEGER            INFO,LDA,LDVS,LIWORK,LWORK,N,SDIM
      DOUBLE PRECISION   RCONDE,RCONDV
      LOGICAL            BWORK(*)
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),VS(LDVS,*),WI(*),WORK(*),WR(*)
      LOGICAL            SELECT
      EXTERNAL           SELECT
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEEV(JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK, &
                       LWORK,INFO)
      CHARACTER          JOBVL,JOBVR
      INTEGER            INFO,LDA,LDVL,LDVR,LWORK,N
      DOUBLE PRECISION   A(LDA,*),VL(LDVL,*),VR(LDVR,*),WI(*),WORK(*), &
                         WR(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEEVX(BALANC,JOBVL,JOBVR,SENSE,N,A,LDA,WR,WI,VL,LDVL, &
                        VR,LDVR,ILO,IHI,SCALE,ABNRM,RCONDE,RCONDV,WORK, &
                        LWORK,IWORK,INFO)
      CHARACTER          BALANC,JOBVL,JOBVR,SENSE
      INTEGER            IHI,ILO,INFO,LDA,LDVL,LDVR,LWORK,N
      DOUBLE PRECISION   ABNRM
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),RCONDE(*),RCONDV(*),SCALE(*),VL(LDVL, &
                         *),VR(LDVR,*),WI(*),WORK(*),WR(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEGS(JOBVSL,JOBVSR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA, &
                       VSL,LDVSL,VSR,LDVSR,WORK,LWORK,INFO)
      CHARACTER          JOBVSL,JOBVSR
      INTEGER            INFO,LDA,LDB,LDVSL,LDVSR,LWORK,N
      DOUBLE PRECISION   A(LDA,*),ALPHAI(*),ALPHAR(*),B(LDB,*),BETA(*), &
                         VSL(LDVSL,*),VSR(LDVSR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEGV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL, &
                       LDVL,VR,LDVR,WORK,LWORK,INFO)
      CHARACTER          JOBVL,JOBVR
      INTEGER            INFO,LDA,LDB,LDVL,LDVR,LWORK,N
      DOUBLE PRECISION   A(LDA,*),ALPHAI(*),ALPHAR(*),B(LDB,*),BETA(*), &
                         VL(LDVL,*),VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEHD2(N,ILO,IHI,A,LDA,TAU,WORK,INFO)
      INTEGER            IHI,ILO,INFO,LDA,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEHRD(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            IHI,ILO,INFO,LDA,LWORK,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEJSV(JOBA,JOBU,JOBV,JOBR,JOBT,JOBP,M,N,A,LDA,SVA,U, &
                        LDU,V,LDV,WORK,LWORK,IWORK,INFO)
      CHARACTER          JOBA,JOBP,JOBR,JOBT,JOBU,JOBV
      INTEGER            INFO,LDA,LDU,LDV,LWORK,M,N
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),SVA(N),U(LDU,*),V(LDV,*),WORK(LWORK)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGELQ2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGELQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGELSD(M,N,NRHS,A,LDA,B,LDB,S,RCOND,RANK,WORK,LWORK, &
                        IWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS,RANK
      DOUBLE PRECISION   RCOND
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),S(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGELS(TRANS,M,N,NRHS,A,LDA,B,LDB,WORK,LWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGELSS(M,N,NRHS,A,LDA,B,LDB,S,RCOND,RANK,WORK,LWORK, &
                        INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS,RANK
      DOUBLE PRECISION   RCOND
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),S(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGELSX(M,N,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,INFO)
      INTEGER            INFO,LDA,LDB,M,N,NRHS,RANK
      DOUBLE PRECISION   RCOND
      INTEGER            JPVT(*)
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGELSY(M,N,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK, &
                        INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,NRHS,RANK
      DOUBLE PRECISION   RCOND
      INTEGER            JPVT(*)
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEQL2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEQLF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEQP3(M,N,A,LDA,JPVT,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      INTEGER            JPVT(*)
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEQPF(M,N,A,LDA,JPVT,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      INTEGER            JPVT(*)
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEQR2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGERFS(TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX, &
                        FERR,BERR,WORK,IWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   A(LDA,*),AF(LDAF,*),B(LDB,*),BERR(*),FERR(*), &
                         WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGERFSX(TRANS,EQUED,N,NRHS,A,LDA,AF,LDAF,IPIV,R,C,B, &
                         LDB,X,LDX,RCOND,BERR,N_ERR_BNDS,ERR_BNDS_NORM, &
                         ERR_BNDS_COMP,NPARAMS,PARAMS,WORK,IWORK,INFO)
      CHARACTER          EQUED,TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      DOUBLE PRECISION   RCOND
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   A(LDA,*),AF(LDAF,*),B(LDB,*),BERR(*),C(*), &
                         ERR_BNDS_COMP(NRHS,*),ERR_BNDS_NORM(NRHS,*), &
                         PARAMS(*),R(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGERQ2(M,N,A,LDA,TAU,WORK,INFO)
      INTEGER            INFO,LDA,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGERQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,M,N
      DOUBLE PRECISION   A(LDA,*),TAU(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGESC2(N,A,LDA,RHS,IPIV,JPIV,SCALE)
      INTEGER            LDA,N
      DOUBLE PRECISION   SCALE
      INTEGER            IPIV(*),JPIV(*)
      DOUBLE PRECISION   A(LDA,*),RHS(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGESDD(JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,IWORK, &
                        INFO)
      CHARACTER          JOBZ
      INTEGER            INFO,LDA,LDU,LDVT,LWORK,M,N
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),S(*),U(LDU,*),VT(LDVT,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGESVD(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK, &
                        INFO)
      CHARACTER          JOBU,JOBVT
      INTEGER            INFO,LDA,LDU,LDVT,LWORK,M,N
      DOUBLE PRECISION   A(LDA,*),S(*),U(LDU,*),VT(LDVT,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      INTEGER            INFO,LDA,LDB,N,NRHS
      INTEGER            IPIV(*)
      DOUBLE PRECISION   A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGESVJ(JOBA,JOBU,JOBV,M,N,A,LDA,SVA,MV,V,LDV,WORK, &
                        LWORK,INFO)
      CHARACTER          JOBA,JOBU,JOBV
      INTEGER            INFO,LDA,LDV,LWORK,M,MV,N
      DOUBLE PRECISION   A(LDA,*),SVA(N),V(LDV,*),WORK(LWORK)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGESVX(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,R,C, &
                        B,LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NRHS
      DOUBLE PRECISION   RCOND
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   A(LDA,*),AF(LDAF,*),B(LDB,*),BERR(*),C(*), &
                         FERR(*),R(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGESVXX(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,R,C, &
                         B,LDB,X,LDX,RCOND,RPVGRW,BERR,N_ERR_BNDS, &
                         ERR_BNDS_NORM,ERR_BNDS_COMP,NPARAMS,PARAMS, &
                         WORK,IWORK,INFO)
      CHARACTER          EQUED,FACT,TRANS
      INTEGER            INFO,LDA,LDAF,LDB,LDX,N,NPARAMS,NRHS,N_ERR_BNDS
      DOUBLE PRECISION   RCOND,RPVGRW
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   A(LDA,*),AF(LDAF,*),B(LDB,*),BERR(*),C(*), &
                         ERR_BNDS_COMP(NRHS,*),ERR_BNDS_NORM(NRHS,*), &
                         PARAMS(*),R(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGETC2(N,A,LDA,IPIV,JPIV,INFO)
      INTEGER            INFO,LDA,N
      INTEGER            IPIV(*),JPIV(*)
      DOUBLE PRECISION   A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGETF2(M,N,A,LDA,IPIV,INFO)
      INTEGER            INFO,LDA,M,N
      INTEGER            IPIV(*)
      DOUBLE PRECISION   A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGETRF(M,N,A,LDA,IPIV,INFO)
      INTEGER            INFO,LDA,M,N
      INTEGER            IPIV(*)
      DOUBLE PRECISION   A(LDA,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LWORK,N
      INTEGER            IPIV(*)
      DOUBLE PRECISION   A(LDA,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDA,LDB,N,NRHS
      INTEGER            IPIV(*)
      DOUBLE PRECISION   A(LDA,*),B(LDB,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGBAK(JOB,SIDE,N,ILO,IHI,LSCALE,RSCALE,M,V,LDV,INFO)
      CHARACTER          JOB,SIDE
      INTEGER            IHI,ILO,INFO,LDV,M,N
      DOUBLE PRECISION   LSCALE(*),RSCALE(*),V(LDV,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGBAL(JOB,N,A,LDA,B,LDB,ILO,IHI,LSCALE,RSCALE,WORK, &
                        INFO)
      CHARACTER          JOB
      INTEGER            IHI,ILO,INFO,LDA,LDB,N
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),LSCALE(*),RSCALE(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGES(JOBVSL,JOBVSR,SORT,SELCTG,N,A,LDA,B,LDB,SDIM, &
                       ALPHAR,ALPHAI,BETA,VSL,LDVSL,VSR,LDVSR,WORK, &
                       LWORK,BWORK,INFO)
      CHARACTER          JOBVSL,JOBVSR,SORT
      INTEGER            INFO,LDA,LDB,LDVSL,LDVSR,LWORK,N,SDIM
      LOGICAL            BWORK(*)
      DOUBLE PRECISION   A(LDA,*),ALPHAI(*),ALPHAR(*),B(LDB,*),BETA(*), &
                         VSL(LDVSL,*),VSR(LDVSR,*),WORK(*)
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGESX(JOBVSL,JOBVSR,SORT,SELCTG,SENSE,N,A,LDA,B,LDB, &
                        SDIM,ALPHAR,ALPHAI,BETA,VSL,LDVSL,VSR,LDVSR, &
                        RCONDE,RCONDV,WORK,LWORK,IWORK,LIWORK,BWORK, &
                        INFO)
      CHARACTER          JOBVSL,JOBVSR,SENSE,SORT
      INTEGER            INFO,LDA,LDB,LDVSL,LDVSR,LIWORK,LWORK,N,SDIM
      LOGICAL            BWORK(*)
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),ALPHAI(*),ALPHAR(*),B(LDB,*),BETA(*), &
                         RCONDE(2),RCONDV(2),VSL(LDVSL,*),VSR(LDVSR,*), &
                         WORK(*)
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL, &
                       LDVL,VR,LDVR,WORK,LWORK,INFO)
      CHARACTER          JOBVL,JOBVR
      INTEGER            INFO,LDA,LDB,LDVL,LDVR,LWORK,N
      DOUBLE PRECISION   A(LDA,*),ALPHAI(*),ALPHAR(*),B(LDB,*),BETA(*), &
                         VL(LDVL,*),VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGEVX(BALANC,JOBVL,JOBVR,SENSE,N,A,LDA,B,LDB,ALPHAR, &
                        ALPHAI,BETA,VL,LDVL,VR,LDVR,ILO,IHI,LSCALE, &
                        RSCALE,ABNRM,BBNRM,RCONDE,RCONDV,WORK,LWORK, &
                        IWORK,BWORK,INFO)
      CHARACTER          BALANC,JOBVL,JOBVR,SENSE
      INTEGER            IHI,ILO,INFO,LDA,LDB,LDVL,LDVR,LWORK,N
      DOUBLE PRECISION   ABNRM,BBNRM
      LOGICAL            BWORK(*)
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),ALPHAI(*),ALPHAR(*),B(LDB,*),BETA(*), &
                         LSCALE(*),RCONDE(*),RCONDV(*),RSCALE(*), &
                         VL(LDVL,*),VR(LDVR,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGGLM(N,M,P,A,LDA,B,LDB,D,X,Y,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),D(*),WORK(*),X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGHRD(COMPQ,COMPZ,N,ILO,IHI,A,LDA,B,LDB,Q,LDQ,Z,LDZ, &
                        INFO)
      CHARACTER          COMPQ,COMPZ
      INTEGER            IHI,ILO,INFO,LDA,LDB,LDQ,LDZ,N
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),Q(LDQ,*),Z(LDZ,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGLSE(M,N,P,A,LDA,B,LDB,C,D,X,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),C(*),D(*),WORK(*),X(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGQRF(N,M,P,A,LDA,TAUA,B,LDB,TAUB,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),TAUA(*),TAUB(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGRQF(M,P,N,A,LDA,TAUA,B,LDB,TAUB,WORK,LWORK,INFO)
      INTEGER            INFO,LDA,LDB,LWORK,M,N,P
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),TAUA(*),TAUB(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGSVD(JOBU,JOBV,JOBQ,M,N,P,K,L,A,LDA,B,LDB,ALPHA,BETA, &
                        U,LDU,V,LDV,Q,LDQ,WORK,IWORK,INFO)
      CHARACTER          JOBQ,JOBU,JOBV
      INTEGER            INFO,K,L,LDA,LDB,LDQ,LDU,LDV,M,N,P
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),Q(LDQ,*), &
                         U(LDU,*),V(LDV,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGGSVP(JOBU,JOBV,JOBQ,M,P,N,A,LDA,B,LDB,TOLA,TOLB,K,L, &
                        U,LDU,V,LDV,Q,LDQ,IWORK,TAU,WORK,INFO)
      CHARACTER          JOBQ,JOBU,JOBV
      INTEGER            INFO,K,L,LDA,LDB,LDQ,LDU,LDV,M,N,P
      DOUBLE PRECISION   TOLA,TOLB
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),Q(LDQ,*),TAU(*),U(LDU,*), &
                         V(LDV,*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGSVJ0(JOBV,M,N,A,LDA,D,SVA,MV,V,LDV,EPS,SFMIN,TOL, &
                        NSWEEP,WORK,LWORK,INFO)
      CHARACTER          JOBV
      INTEGER            INFO,LDA,LDV,LWORK,M,MV,N,NSWEEP
      DOUBLE PRECISION   EPS,SFMIN,TOL
      DOUBLE PRECISION   A(LDA,*),D(N),SVA(N),V(LDV,*),WORK(LWORK)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGSVJ1(JOBV,M,N,N1,A,LDA,D,SVA,MV,V,LDV,EPS,SFMIN,TOL, &
                        NSWEEP,WORK,LWORK,INFO)
      CHARACTER          JOBV
      INTEGER            INFO,LDA,LDV,LWORK,M,MV,N,N1,NSWEEP
      DOUBLE PRECISION   EPS,SFMIN,TOL
      DOUBLE PRECISION   A(LDA,*),D(N),SVA(N),V(LDV,*),WORK(LWORK)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGTCON(NORM,N,DL,D,DU,DU2,IPIV,ANORM,RCOND,WORK,IWORK, &
                        INFO)
      CHARACTER          NORM
      INTEGER            INFO,N
      DOUBLE PRECISION   ANORM,RCOND
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   D(*),DL(*),DU(*),DU2(*),WORK(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGTRFS(TRANS,N,NRHS,DL,D,DU,DLF,DF,DUF,DU2,IPIV,B,LDB, &
                        X,LDX,FERR,BERR,WORK,IWORK,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDB,LDX,N,NRHS
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   B(LDB,*),BERR(*),D(*),DF(*),DL(*),DLF(*),DU(*), &
                         DU2(*),DUF(*),FERR(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGTSV(N,NRHS,DL,D,DU,B,LDB,INFO)
      INTEGER            INFO,LDB,N,NRHS
      DOUBLE PRECISION   B(LDB,*),D(*),DL(*),DU(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGTSVX(FACT,TRANS,N,NRHS,DL,D,DU,DLF,DF,DUF,DU2,IPIV,B, &
                        LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)
      CHARACTER          FACT,TRANS
      INTEGER            INFO,LDB,LDX,N,NRHS
      DOUBLE PRECISION   RCOND
      INTEGER            IPIV(*),IWORK(*)
      DOUBLE PRECISION   B(LDB,*),BERR(*),D(*),DF(*),DL(*),DLF(*),DU(*), &
                         DU2(*),DUF(*),FERR(*),WORK(*),X(LDX,*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGTTRF(N,DL,D,DU,DU2,IPIV,INFO)
      INTEGER            INFO,N
      INTEGER            IPIV(*)
      DOUBLE PRECISION   D(*),DL(*),DU(*),DU2(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE DGTTRS(TRANS,N,NRHS,DL,D,DU,DU2,IPIV,B,LDB,INFO)
      CHARACTER          TRANS
      INTEGER            INFO,LDB,N,NRHS
      INTEGER            IPIV(*)
      DOUBLE PRECISION   B(LDB,*),D(*),DL(*),DU(*),DU2(*)
      END
      END INTERFACE


!DEC$ENDIF
