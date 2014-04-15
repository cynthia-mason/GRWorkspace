! COPYRIGHT (c) 2007  Science and Technology Facilities Council
! Original date August 2007
! PACKAGE MF64A/AD
! AUTHORS Iain Duff (i.duff@rl.ac.uk) and
!         Jacko Koster (jacko.koster@uninett.no)
!
! Version 1.1.0
! See ChangeLog for version history.
!
      SUBROUTINE MF64ID(ICNTL)
      IMPLICIT NONE

C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MF64A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER ICNTL(10)
C
C  Local variables
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the defaut value, the incoming data is checked for
C     out-of-range indices and duplicates.  Setting ICNTL(4) to any
C     other will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) non-zero, if the data is
C     known to avoid these problems.
C
C    ICNTL(5) to ICNTL(10) are not used by MF64A/AD but are set to
C     zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MF64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MF64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
C  Purpose
C  =======
C
C This subroutine attempts to find a column permutation for an NxN
C sparse matrix A = {a_ij} that makes the permuted matrix have N
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a column permutation that maximizes the smallest element
C on the diagonal, maximizes the sum of the diagonal entries, or
C maximizes the product of the diagonal entries of the permuted matrix.
C For the latter option, the subroutine also finds scaling factors
C that may be used to scale the matrix so that the nonzero diagonal
C entries of the permuted matrix are one in absolute value and all the
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..N, for the
C rows and v(j), j=1..N, for the columns are returned so that the
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C
C  Parameters
C  ==========
C
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      COMPLEX*16 A(NE)
      DOUBLE PRECISION DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C     permuted matrix has as many entries on its diagonal as possible.
C     The values on the diagonal are of arbitrary size. HSL subroutine
C     MC21A/AD is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     The algorithm differs from the one used for JOB = 2 and may
C     have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C     of the diagonal entries of the permuted matrix is maximized.
C     See [3].
C   5 Compute a column permutation of the matrix so that the product
C     of the diagonal entries of the permuted matrix is maximized
C     and vectors to scale the matrix so that the nonzero diagonal
C     entries of the permuted matrix are one in absolute value and
C     all the off-diagonal entries are less than or equal to one in
C     absolute value. See [3].
C  Restriction: 1 <= JOB <= 5.
C
C N is an INTEGER variable which must be set by the user to the
C   order of the matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the
C   subroutine.
C   Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE.
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated
C   entries are not allowed. The array IRN is not altered by the
C   subroutine.
C
C A is a COMPLEX*16 array of length NE.
C   The user must set A(K), K=1..NE, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not used by the subroutine when JOB = 1.
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C CPERM is an INTEGER array of length N that need not be set by the
C   user. On successful exit, CPERM contains the column permutation.
C   Column ABS(CPERM(J)) of the original matrix is column J in the
C   permuted matrix, J=1..N. For the N-NUM  entries of CPERM that are
C   not matched the permutation array is set negative so that a full
C   permutation of the matrix is obtained even in the structurally
C   singular case.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LIW >= 5N
C     JOB = 2 :  LIW >= 4N
C     JOB = 3 :  LIW >= 10N + NE
C     JOB = 4 :  LIW >= 5N
C     JOB = 5 :  LIW >= 5N
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LDW is not used
C     JOB = 2 :  LDW >= N
C     JOB = 3 :  LDW >= NE
C     JOB = 4 :  LDW >= 2N + NE
C     JOB = 5 :  LDW >= 3N + NE
C
C DW is a DOUBLE PRECISION array of length LDW
C   that is used for workspace. If JOB = 5, on return,
C   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MF64A/AD and must be set by the user before calling
C   MF64A/AD. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC46I/ID is -1.
C
C   ICNTL(4) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC46I/ID is 0.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C        overflow when used to scale the matrix.
C        (For JOB = 5 entry only.)
C   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : NE < 1. Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -7 : repeated entries are found. INFO(2) contains the index of a
C        column in which such entries are found.
C  INFO(3) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C References:
C  [1]  I. S. Duff, (1981),
C       "Algorithm 575. Permutations for a zero-free diagonal",
C       ACM Trans. Math. Software 7(3), 387-390.
C  [2]  I. S. Duff and J. Koster, (1998),
C       "The design and use of algorithms for permuting large
C       entries to the diagonal of sparse matrices",
C       SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3]  I. S. Duff and J. Koster, (1999),
C       "On algorithms for permuting large entries to the diagonal
C       of sparse matrices",
C       Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England.

C Local variables and parameters
      INTEGER I,J,K
      DOUBLE PRECISION FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
C External routines and functions
      EXTERNAL MC21AD,MF64BD,MC64RD,MC64SD,MC64WD
C Intrinsic functions
      INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
      RINF = HUGE(RINF)

C Check value of JOB
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Check value of NE
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
C Check LIW
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
C If JOB = 1, do not check
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
C Check row indices. Use IW(1:N) as workspace
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF

C Print diagnostics on input
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF

C Set components of INFO to zero
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE

C Compute maximum matching with MC21A/AD
      IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:5N) is workspace
        CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.2) THEN
C IW(1:5N), DW(1:N) are workspaces
        CALL MF64BD(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value.
        CALL MC64RD(N,NE,IP,IW,DW)
C IW(NE+1:NE+10N) is workspace
        CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
C B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
C Check size of scaling factors
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE
C       GO TO 90
      ENDIF

   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
C Matrix is structurally singular, return with warning
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
C Scaling factors are large, return with warning
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF

C Print diagnostics on output
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MF64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MF64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MF64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MF64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MF64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MF64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MF64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MF64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MF64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MF64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MF64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      COMPLEX*16 A(NE)
      DOUBLE PRECISION D(N)

C N, NE, IP, IRN are described in MF64A/AD.
C A is a COMPLEX*16 array of length
C   NE. A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length N. On exit, it contains the
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C IW is an INTEGER work array of length 4N.
C DW is a DOUBLE PRECISION array of length N.

C Local variables
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
C Local parameters
      DOUBLE PRECISION RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
      INTRINSIC ABS,MIN
C External subroutines and/or functions
      EXTERNAL MC64DD,MC64ED,MC64FD


C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Initialization
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
C Scan columns of matrix;
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
C Update BV with smallest of all the largest maximum absolute values
C of the rows.
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Rescan unassigned columns; improve initial assignment
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
        CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
C Add row I to Q2
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
C Add row I to Q, and push it
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE

        DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64ED(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
C End of dummy loop; this point is never reached
          ENDIF
C Move row Q0
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
C Scan column that matches with row Q0
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Update D(I)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
C Delete row I from Q (if necessary); add row I to Q2
                IF (DI.NE.MINONE) THEN
                  LPOS = L(I)
                  CALL MC64FD(LPOS,QLEN,N,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
C Add row I to Q (if necessary); push row I up Q
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,1)
              ENDIF
C Update tree
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE

C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
        BV = MIN(BV,CSP)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is bottleneck value of final matching
      IF (NUM.EQ.N) GO TO 1000

C Matrix is structurally singular, complete IPERM.
C JPERM, PR are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = - I
  320 CONTINUE

 1000 RETURN
      END

