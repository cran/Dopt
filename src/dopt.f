C     (C) COPYRIGHT, 1994.  Royal Statistical Society, London.
C
C     Requests to incorporate any or all of this software into any
C     commercial product must be addressed to 
C
C     The Secretary,
C     The Royal Statistical Society,
C     25 Enford St,
C     LONDON,      UK,    W1H2BH

C-------------------------------------------------------------------------

      SUBROUTINE FEDROV(X, DIM1, NCAND, KIN, N, NBLOCK, IN, BLKSIZ, K,
     +    RSTART, NRBAR, D, RBAR, PICKED, LNDET, XX, TOL, ZPZ, WK,
     +    IFAULT)

C     [Formerly named DOPT, but renamed here to avoid loader table
C      conflicts.  W. Venables.]
C
C
C     Slightly changed from version distributed by W. Venables, for use
C     with R. All lines with changes marked #KJETIL.  Kjetil Halvorsen.
C     Subroutine RAND renamed to RRAND to avoid warning from gcc. 
C
C
C     Heuristic algorithm to pick N rows of X out of NCAND to maximize
C     the determinant of X'X, using the Fedorov exchange algorithm.
C
C     Code written by Alan Miller and Nam Nguyen from CSIRO, Melbourne,
C     Australia and to be published in Applied Statistics.
C
C     [Now appeared as AS295, Applied Statistics, vol 43, 669-677, 
C     (1994)  W. Venables.]
C
C     Auxiliary routines required: A random number generator (RRAND),
C     and routines CLEAR & REGCF from AS 274.
C
C     [All supplied herewith.  W. Venables.]
C
      implicit none
      INTEGER DIM1, NCAND, KIN, N, NBLOCK, IN(*), BLKSIZ(*), K, NRBAR,
     +    PICKED(N), IFAULT
      LOGICAL RSTART
C      REAL X(DIM1, KIN)    Original definition. We have to change to double;#KJETIL
      DOUBLE PRECISION X(DIM1, KIN)
      DOUBLE PRECISION D(K), RBAR(NRBAR), LNDET, XX(K), TOL(K), WK(K),
     +    ZPZ(NCAND, *)
C
C     Local variables
C
      INTEGER I, J, NIN, POINT, CASE, NB, BLOCK, L, POS, BEST, FIRST,
     +    LAST, CAND, LASTIN, LSTOUT, DROP, REMPOS, BL, RPOS, LAST1,
     +    LAST2, FIRST1, FIRST2, POS1, POS2, BLOCK1, BLOCK2, CASE1,
     +    CASE2, POSI, POSJ, BESTB1, BESTB2, BESTP1, BESTP2, RANK,
     +    MXRANK, INC
      DOUBLE PRECISION ONE, ZERO, MINUS1, TEMP, EPS, DMAX, ABOVE1,
     +    SUM, SMALL, HUNDRD
      DOUBLE PRECISION DELTA
      REAL RRAND
      LOGICAL CHANGE
      DATA ONE/1.D0/, ZERO/0.D0/, ABOVE1/1.0001D0/, EPS/1.D-15/,
     +    MINUS1/-1.D0/, SMALL/1.D-04/, HUNDRD/100.D0/
C
C     Argument checks
C
      IFAULT = 0
      IF (DIM1 .LT. NCAND) IFAULT = 1
      IF (K .GT. N) IFAULT = IFAULT + 2
      IF (NRBAR .LT. K* (K-1)/2) IFAULT = IFAULT + 4
      IF (K .NE. KIN+NBLOCK) IFAULT = IFAULT + 8
      IF (NBLOCK .GT. 1) THEN
        L = 0
        DO 10 BLOCK = 1, NBLOCK
   10   L = L + BLKSIZ(BLOCK)
        IF (N .NE. L) IFAULT = IFAULT + 16
      ELSE
        IF (N .NE. BLKSIZ(1)) IFAULT = IFAULT + 16
      END IF
C
C     NB = max(1, NBLOCK) so that we can force it to go through DO-loops
C     once.   NIN = no. of design points forced into the design.
C
      NB = MAX(1, NBLOCK)
      NIN = 0
      DO 20 I = 1, NB
        IF (IN(I) .LT. 0) GO TO 30
        IF (IN(I) .GT. 0) THEN
          IF (IN(I) .GT. BLKSIZ(I)) GO TO 30
          NIN = NIN + IN(I)
        END IF
   20 CONTINUE
      IF (NIN .LE. N) GO TO 40
   30 IFAULT = IFAULT + 32
   40 CONTINUE
      IF (IFAULT .NE. 0) RETURN
      CALL CLEAR(K, NRBAR, D, RBAR, IFAULT)
C
C     Set up an array of tolerances.
C
      DO 50 I = 1, K
   50 TOL(I) = ZERO
      BLOCK = 1
      DO 70 CASE = 1, NCAND
        CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, CASE)
        DO 60 I = 1, K
   60   TOL(I) = TOL(I) + ABS(XX(I))
   70 CONTINUE
      TEMP = FLOAT(N) * EPS / NCAND
      DO 80 I = 1, K
        IF (I .LE. NBLOCK) THEN
          TOL(I) = EPS
        ELSE
          TOL(I) = TOL(I) * TEMP
        END IF
   80 CONTINUE
C
C     Form initial Cholesky factorization.
C
      POS = 1
      DO 120 BLOCK = 1, NB
        IF (RSTART) THEN
          LAST1 = (IN(BLOCK) + BLKSIZ(BLOCK))/2
          INC = SQRT(FLOAT(NCAND) + SMALL)
        END IF
        DO 110 I = 1, BLKSIZ(BLOCK)
          IF (RSTART .AND. I .GT. IN(BLOCK)) THEN
            POINT = 1 + NCAND * RRAND()
C
C     If I <= LAST1, use a random point, otherwise find the candidate
C     which maximizes the rank, and then maximizes the subspace
C     determinant for that rank.
C
            IF (I .GT. LAST1) THEN
              MXRANK = 0
              LNDET = -HUNDRD
              DO 100 CAND = 1, NCAND, INC
                CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, POINT)
                CALL MODTR2(K, NRBAR, XX, D, RBAR, TOL, RANK, SUM)
                IF (RANK .LT. MXRANK) GO TO 90
                IF (RANK .EQ. MXRANK .AND. SUM .LT. LNDET) GO TO 90
                BEST = POINT
                MXRANK = RANK
                LNDET = SUM * ABOVE1
   90           POINT = POINT + INC
                IF (POINT .GT. NCAND) POINT = POINT - NCAND
  100         CONTINUE
              POINT = BEST
            END IF
            PICKED(POS) = POINT
          ELSE
C
C     Case in which a full design has been input, or points are to be
C     forced into the design.
C
            POINT = PICKED(POS)
          END IF
C
C     Augment the Cholesky factorization.
C
          CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, POINT)
          CALL MODTRI(K, NRBAR, ONE, XX, D, RBAR, TOL)
          POS = POS + 1
  110   CONTINUE
  120 CONTINUE
C
C     Adjust factorization in case of singular matrix.
C
      CALL SINGM(K, NRBAR, D, RBAR, TOL, WK, IFAULT)
C
C     If rank of input design < K, try replacing points.
C
      IF (IFAULT .EQ. 0) GO TO 280
C
C     Find first row of Cholesky factorization with a zero multiplier.
C
  180 DO 190 POS = 1, K
        IF (D(POS) .LT. TOL(POS)) GO TO 200
  190 CONTINUE
      GO TO 280
C
C     Find linear relationship between variable in position POS and the
C     previous variables.
C
  200 L = POS - 1
      DO 210 I = 1, POS-1
        WK(I) = RBAR(L)
        L = L + K - I - 1
  210 CONTINUE
      CALL REGCF(K, NRBAR, D, RBAR, WK, TOL, WK, POS-1, IFAULT)
C
C     Find a candidate point which does not satisfy this linear
C     relationship.   Use a random start.
C
      BL = 1
      CASE = 1 + NCAND * RRAND()
      DO 230 CAND = 1, NCAND
        CALL GETX(X, DIM1, KIN, NBLOCK, K, BL, XX, CASE)
        SUM = XX(POS)
        DO 220 I = 1, POS-1
  220   SUM = SUM - WK(I) * XX(I)
        IF (ABS(SUM) .GT. HUNDRD * TOL(POS)) GO TO 240
        CASE = CASE + 1
        IF (CASE .GT. NCAND) CASE = 1
  230 CONTINUE
C
C     Failed to find any candidate point which would make the design
C     of higher rank.
C
      IFAULT = -1
      RETURN
C
C     Before adding the point, find one which it can replace without
C     lowering the rank.
C
  240 BL = 0
      TEMP = ONE - SMALL
      POS = IN(1) + 1
      DO 270 BLOCK = 1, NB
        DO 260 J = IN(BLOCK)+1, BLKSIZ(BLOCK)
          L = PICKED(POS)
          CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, L)
          CALL BKSUB2(RBAR, NRBAR, K, XX, WK)
          SUM = ZERO
          DO 250 I = 1, K
  250     IF (D(I) .GT. TOL(I)) SUM = SUM + WK(I)**2 / D(I)
          IF (SUM .LT. TEMP) THEN
            TEMP = SUM
            REMPOS = POS
            BL = BLOCK
          END IF
          POS = POS + 1
  260   CONTINUE
        IF (BLOCK .LT. NBLOCK) POS = POS + IN(BLOCK+1)
  270 CONTINUE
C
C     If BL = 0 it means that any point removed from the existing design
C     would reduce the rank.
C
      IF (BL .EQ. 0) THEN
        IFAULT = -1
        RETURN
      END IF
C      
C     Add candidate CASE in block BL, then delete the design point
C     already in that position.
C
      CALL GETX(X, DIM1, KIN, NBLOCK, K, BL, XX, CASE)
      CALL MODTRI(K, NRBAR, ONE, XX, D, RBAR, TOL)
      L = PICKED(REMPOS)
      CALL GETX(X, DIM1, KIN, NBLOCK, K, BL, XX, L)
      CALL MODTRI(K, NRBAR, MINUS1, XX, D, RBAR, TOL)
      PICKED(REMPOS) = CASE

      GO TO 180
C----------------------------------------------------------------------
C
C     Design is now of full rank.
C 
C     Calculate z'z for all candidate points.
C     z is the solution of R'z = x, so that z'z = x'.inv(X'X).x
C     WK holds sqrt(D) times vector z on return from BKSUB2.
C
  280 DO 310 BLOCK = 1, NB
        DO 300 CASE = 1, NCAND
          CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, CASE)
          CALL BKSUB2(RBAR, NRBAR, K, XX, WK)
          TEMP = ZERO
          DO 290 I = 1, K
  290     TEMP = TEMP + WK(I)**2 / D(I)
          ZPZ(CASE, BLOCK) = TEMP
  300   CONTINUE
  310 CONTINUE

C----------------------------------------------------------------------
C
C     Start of Fedorov exchange algorithm
C
      LASTIN = 0
      LSTOUT = 0
  320 CHANGE = .FALSE.
      LAST = 0
      DO 420 BLOCK = 1, NB
        FIRST = LAST + 1 + IN(BLOCK)
        LAST = LAST + BLKSIZ(BLOCK)
        DMAX = SMALL
        BEST = 0
C
C     Start at a random position within the block.
C     I = no. of point being considered for deletion.
C
        POS = FIRST + (BLKSIZ(BLOCK) - IN(BLOCK)) * RRAND()
        DO 350 CASE = IN(BLOCK)+1, BLKSIZ(BLOCK)
          POS = POS + 1
          IF (POS .GT. LAST) POS = FIRST
          I = PICKED(POS)
          IF (I .EQ. LASTIN) GO TO 350
          CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, I)
          CALL BKSUB2(RBAR, NRBAR, K, XX, WK)
          CALL BKSUB1(RBAR, NRBAR, K, WK, WK, TOL, D)
C 
C     Cycle through the candidates for exchange, using a random start.
C     J = no. of point being considered for addition.
C
          J = 1 + NCAND * RRAND()
          DO 340 CAND = 1, NCAND
            J = J + 1
            IF (J .GT. NCAND) J = 1
            IF (J .EQ. I .OR. J .EQ. LSTOUT) GO TO 340
C... The Cauchy-Schwarz test.
            TEMP = ZPZ(J,BLOCK) - ZPZ(I,BLOCK)
            IF (TEMP .LT. DMAX) GO TO 340
C
            CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, J)
            SUM = ZERO
            DO 330 L = 1, K
  330       SUM = SUM + XX(L) * WK(L)
            TEMP = TEMP + SUM**2 - ZPZ(I,BLOCK)*ZPZ(J,BLOCK)
            IF (TEMP .GT. DMAX) THEN
              DMAX = TEMP * ABOVE1
              BEST = J
              REMPOS = POS
              DROP = I
            END IF
  340     CONTINUE
  350   CONTINUE
C
C     Exchange points BEST & DROP in position REMPOS, if the determinant
C     is increased.
C
        IF (BEST .NE. 0) THEN
          CHANGE = .TRUE.
          IF (NB .EQ. 1) THEN
            LASTIN = BEST
            LSTOUT = DROP
          END IF
C
C     Add the new point, BEST, first to avoid ill-conditioning.
C     Update z'z.
C
          CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, BEST)
          CALL BKSUB2(RBAR, NRBAR, K, XX, WK)
          CALL BKSUB1(RBAR, NRBAR, K, WK, WK, TOL, D)
          CALL MODTRI(K, NRBAR, ONE, XX, D, RBAR, TOL)
          TEMP = ONE + ZPZ(BEST,BLOCK)
          DO 360 BL = 1, NB
            DO 380 CASE = 1, NCAND
              CALL GETX(X, DIM1, KIN, NBLOCK, K, BL, XX, CASE)
              SUM = ZERO
              DO 370 L = 1, K
  370         SUM = SUM + XX(L)*WK(L)
              ZPZ(CASE,BL) = ZPZ(CASE,BL) - SUM**2 / TEMP
  380       CONTINUE
  360     CONTINUE
C
C     Remove the point DROP, and update z'z.
C
          CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, DROP)
          CALL BKSUB2(RBAR, NRBAR, K, XX, WK)
          CALL BKSUB1(RBAR, NRBAR, K, WK, WK, TOL, D)
          CALL MODTRI(K, NRBAR, MINUS1, XX, D, RBAR, TOL)
          TEMP = ONE - ZPZ(DROP,BLOCK)
          DO 390 BL = 1, NB
            DO 410 CASE = 1, NCAND                            
              CALL GETX(X, DIM1, KIN, NBLOCK, K, BL, XX, CASE)
              SUM = ZERO
              DO 400 L = 1, K
  400         SUM = SUM + XX(L)*WK(L)
              ZPZ(CASE,BL) = ZPZ(CASE,BL) + SUM**2 / TEMP
  410       CONTINUE
  390     CONTINUE
C
          PICKED(REMPOS) = BEST
        END IF
  420 CONTINUE
C
C     Repeat until there is no further improvement
C
      IF (CHANGE) GO TO 320
C---------------------------------------------------------------------
C
C     If there is more than one block, try swapping treatments between
C     blocks.   This is the Cook & Nachtsheim (1989) algorithm.
C
      IF (NBLOCK .LE. 1) GO TO 500
C
C
C...  RPOS is the position of the first element in RBAR after the rows
C     for the block constants.
      RPOS = NBLOCK * K - NBLOCK * (NBLOCK + 1)/2 + 1
C
  430 LAST1 = 0
C...  POS1 & POS2 will hold the positions of the start of the means of
C     the X-variables in the two blocks being considered in RBAR.
      POS1 = NBLOCK
      DMAX = SMALL
      CHANGE = .FALSE.
      DO 490 BLOCK1 = 1, NBLOCK-1
        FIRST1 = LAST1 + 1 + IN(BLOCK1)
        LAST1 = LAST1 + BLKSIZ(BLOCK1)
        LAST2 = LAST1
        POS2 = POS1 + K - 1 - BLOCK1
        DO 480 BLOCK2 = BLOCK1+1, NBLOCK
          FIRST2 = LAST2 + 1 + IN(BLOCK2)
          LAST2 = LAST2 + BLKSIZ(BLOCK2)
          DO 470 CASE1 = IN(BLOCK1)+1, BLKSIZ(BLOCK1)
            POSI = FIRST1 - 1 + CASE1
            I = PICKED(POSI)
            CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, I)
            DO 440 L = 1, KIN
  440       ZPZ(L,1) = XX(L+NBLOCK)
            DO 460 CASE2 = IN(BLOCK2)+1, BLKSIZ(BLOCK2)
              POSJ = FIRST2 - 1 + CASE2
              J = PICKED(POSJ)
              IF (I .EQ. J) GO TO 460
              CALL GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, J)
              DO 450 L = 1, KIN
  450         ZPZ(L,2) = XX(L+NBLOCK)
C...  Pass the orthogonal factorization to DELTA with the top NBLOCK
C     rows removed, i.e. without that part relating to the blocks.
              TEMP = DELTA(KIN, ZPZ(1,1), ZPZ(1,2), RBAR(POS1),
     +             RBAR(POS2), BLKSIZ(BLOCK1), BLKSIZ(BLOCK2), NRBAR,
     +             D(NBLOCK+1), RBAR(RPOS), ZPZ(1,3), WK, XX, ZPZ(1,2))
              IF (TEMP .GT. DMAX) THEN
                DMAX = TEMP * ABOVE1
                BESTB1 = BLOCK1
                BESTB2 = BLOCK2
                BESTP1 = POSI
                BESTP2 = POSJ
                CHANGE = .TRUE.
              END IF
  460       CONTINUE
  470     CONTINUE
          POS2 = POS2 + K - 1 - BLOCK2
  480   CONTINUE
        POS1 = POS1 + K - 1 - BLOCK1
  490 CONTINUE
C
C     If CHANGE = .TRUE. then make the swap, otherwise the search ends.
C
      IF (CHANGE) THEN
        I = PICKED(BESTP1)
        J = PICKED(BESTP2)
        CALL GETX(X, DIM1, KIN, NBLOCK, K, BESTB2, XX, I)
        CALL MODTRI(K, NRBAR, ONE, XX, D, RBAR, TOL)
        CALL GETX(X, DIM1, KIN, NBLOCK, K, BESTB1, XX, J)
        CALL MODTRI(K, NRBAR, ONE, XX, D, RBAR, TOL)
        CALL GETX(X, DIM1, KIN, NBLOCK, K, BESTB1, XX, I)
        CALL MODTRI(K, NRBAR, MINUS1, XX, D, RBAR, TOL)
        CALL GETX(X, DIM1, KIN, NBLOCK, K, BESTB2, XX, J)
        CALL MODTRI(K, NRBAR, MINUS1, XX, D, RBAR, TOL)
        PICKED(BESTP1) = J
        PICKED(BESTP2) = I
        GO TO 430
      END IF
C
C
C     Calculate log of determinant.
C
  500 LNDET = ZERO
      DO 510 I = 1, K
  510 LNDET = LNDET + LOG(D(I))
      RETURN
      END
C
C
C
C
C
C
C
      SUBROUTINE MODTRI(NP, NRBAR, WEIGHT, XROW, D, RBAR, TOL)
C
C     ALGORITHM AS75.1 with minor modifications 
C
C     Modify a triangular (Cholesky) decomposition.
C     Calling this routine updates d and rbar by adding another design
C     point with weight = WEIGHT, which may be negative.
C
C     *** WARNING  Array XROW is overwritten.
C
      implicit none 
      INTEGER NP, NRBAR
      DOUBLE PRECISION WEIGHT, XROW(NP), D(NP), RBAR(*), TOL(NP)
C
C     Local variables
C
      INTEGER I, K, NEXTR
      DOUBLE PRECISION ZERO, W, XI, DI, WXI, DPI, CBAR, SBAR, XK
C
      DATA ZERO/0.D0/
C
      W = WEIGHT
      NEXTR = 1
      DO 30 I = 1, NP
C
C     Skip unnecessary transformations.   Test on exact zeroes must be
C     used or stability can be destroyed.
C
        IF (W .EQ. ZERO) RETURN
C
        XI = XROW(I)
        IF (ABS(XI) .LT. TOL(I)) THEN
          NEXTR = NEXTR + NP - I
          GO TO 30
        END IF
C
        DI = D(I)
        WXI = W*XI
        DPI = DI + WXI*XI
C
C     Test for new singularity.
C
        IF (DPI .LT. TOL(I)) THEN
          DPI = ZERO
          CBAR = ZERO
          SBAR = ZERO
          W = ZERO
        ELSE
          CBAR = DI/DPI
          SBAR = WXI/DPI
          W = CBAR*W
        END IF
C
        D(I) = DPI
        DO 20 K = I + 1, NP
          XK = XROW(K)
          XROW(K) = XK - XI*RBAR(NEXTR)
          RBAR(NEXTR) = CBAR*RBAR(NEXTR) + SBAR*XK
          NEXTR = NEXTR + 1
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE MODTR2(NP, NRBAR, XROW, D, RBAR, TOL, RANK, LNDET)
C
C     ALGORITHM AS75.1 with modifications.
C     Calculate the effect of an update of a QR factorization upon the
C     rank and determinant, without changing D or RBAR.
C
C     *** WARNING  Array XROW is overwritten.
C
      implicit none
      INTEGER NP, NRBAR, RANK
      DOUBLE PRECISION XROW(NP), D(NP), RBAR(*), TOL(NP), LNDET
C
C     Local variables
C
      INTEGER I, K, NEXTR, J
      DOUBLE PRECISION ZERO, ONE, W, XI, DI, WXI, DPI, CBAR, XK
C
      DATA ZERO/0.D0/, ONE/1.D0/
C
      W = ONE
      RANK = 0
      LNDET = ZERO
      NEXTR = 1
      DO 30 I = 1, NP
C
C     Skip unnecessary transformations.   Test on exact zeroes must be
C     used or stability can be destroyed.
C
        IF (W .EQ. ZERO) THEN
          DO 10 J = I, NP
            IF (D(J) .GT. TOL(J)) THEN
              RANK = RANK + 1
              LNDET = LNDET + LOG(D(J))
            END IF
   10     CONTINUE
          RETURN
        END IF
C
        XI = XROW(I)
        IF (ABS(XI) .LT. TOL(I)) THEN
          IF (D(I) .GT. TOL(I)) THEN
            RANK = RANK + 1
            LNDET = LNDET + LOG(D(I))
          END IF
          NEXTR = NEXTR + NP - I
          GO TO 30
        END IF
C
        DI = D(I)
        WXI = W*XI
        DPI = DI + WXI*XI
C
C     Test for new singularity.
C
        IF (DPI .LT. TOL(I)) THEN
          DPI = ZERO
          CBAR = ZERO
          W = ZERO
        ELSE
          CBAR = DI/DPI
          W = CBAR*W
          LNDET = LNDET + LOG(DPI)
          RANK = RANK + 1
        END IF
C
        DO 20 K = I + 1, NP
          XK = XROW(K)
          XROW(K) = XK - XI*RBAR(NEXTR)
          NEXTR = NEXTR + 1
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE GETX(X, DIM1, KIN, NBLOCK, K, BLOCK, XX, CASE)
C
C     Copy one case from X to XX.
C
      implicit none
      INTEGER DIM1, KIN, NBLOCK, K, BLOCK, CASE
C      REAL X(DIM1, KIN)    #KJETIL   Must be changed to DOUBLE:
      DOUBLE PRECISION X(DIM1, KIN)
C
      DOUBLE PRECISION XX(K)
C
C     Local variables.
C
      INTEGER I, J
      DOUBLE PRECISION ZERO, ONE
      DATA ZERO/0.D0/, ONE/1.D0/
C
      DO 10 I = 1, NBLOCK
        IF (I .NE. BLOCK) THEN
          XX(I) = ZERO
        ELSE
          XX(I) = ONE
        END IF
   10 CONTINUE
      J = NBLOCK + 1
      DO 20 I = 1, KIN
        XX(J) = X(CASE, I)
        J = J + 1
   20 CONTINUE
C
      RETURN
      END
C      
      SUBROUTINE BKSUB1(RBAR, NRBAR, K, RHS, SOLN, TOL, D)
C 
C     Solves  D R y = z for y (SOLN), where z = RHS.
C     RBAR is an upper-triangular matrix with implicit 1's on it's
C     diagonal, stored by rows.
C
      implicit none
      INTEGER NRBAR, K
      DOUBLE PRECISION RBAR(NRBAR), RHS(K), SOLN(K), TOL(K), D(K)
C
C     Local variables
C
      INTEGER ROW, COL, POS
      DOUBLE PRECISION ZERO, TEMP
      DATA ZERO /0.D0/
C
      POS = K * (K - 1) / 2
      DO 20 ROW = K, 1, -1
        IF (D(ROW) .GT. TOL(ROW)) THEN
          TEMP = RHS(ROW) / D(ROW)
          DO 10 COL = K, ROW+1, -1
            TEMP = TEMP - RBAR(POS) * SOLN(COL)
            POS = POS - 1
   10     CONTINUE
          SOLN(ROW) = TEMP
        ELSE
          POS = POS - K + ROW
          SOLN(ROW) = ZERO
        END IF
   20 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE BKSUB2(RBAR, NRBAR, K, RHS, SOLN)
C
C     Solves  R'(sqrt(D).z) = x where (sqrt(D).z) = SOLN and x = RHS.
C     RBAR is an upper-triangular matrix with implicit 1's on it's
C     diagonal, stored by rows.
C 
      implicit none
      INTEGER NRBAR, K
      DOUBLE PRECISION RBAR(NRBAR), RHS(K), SOLN(K)
C
C     Local variables
C
      INTEGER ROW, COL, POS
      DOUBLE PRECISION TEMP
C
      SOLN(1) = RHS(1)
      DO 20 ROW = 2, K
        TEMP = RHS(ROW)
        POS = ROW - 1
        DO 10 COL = 1, ROW-1
          TEMP = TEMP - RBAR(POS) * SOLN(COL)
          POS = POS + K - COL - 1
   10   CONTINUE
        SOLN(ROW) = TEMP
   20 CONTINUE
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DELTA(K, XJ, XL, XBARI, XBARK, NI, NK,
     +    NRBAR, D, RBAR, Z, A, B, DIFF)
C
C     Calculate the delta function for the swap of case J in block I
C     with case L in block K.   Uses the method of Cook & Nachtsheim.
C
      implicit none
      INTEGER K, NI, NK, NRBAR
      DOUBLE PRECISION XJ(K), XL(K), XBARI(K), XBARK(K), D(K),
     +    RBAR(NRBAR), Z(K,3), A(K), B(K), DIFF(K)
      DOUBLE PRECISION DOTP
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION CONST, ONE, TWO, TEMP, E11, E12, E21, E22
      DATA ONE/1.D0/, TWO/2.D0/
C
C     Calculate vectors DIFF, A & B.
C
      CONST = TWO - ONE/NI - ONE/NK
      DO 10 I = 1, K
        TEMP = XJ(I) - XL(I)
        DIFF(I) = -TEMP
        A(I) = TEMP - XBARI(I) + XBARK(I)
        B(I) = A(I) - CONST * TEMP
   10 CONTINUE
C
C     Calculate the z-vectors by back-substitution.
C     Z1 for A, Z2 for B and Z3 for DIFF.
C     N.B. The solutions returned from BKSUB2 have the I-th element
C     multiplied by sqrt(D(I)).
C
      CALL BKSUB2(RBAR, NRBAR, K, A, Z(1,1))
      CALL BKSUB2(RBAR, NRBAR, K, B, Z(1,2))
      CALL BKSUB2(RBAR, NRBAR, K, DIFF, Z(1,3))
C
C     Calculate the elements E11, E12, E21 & E22 as dot-products of the
C     appropriate z-vectors.
C
      E11 = DOTP(K, Z(1,3), Z(1,1), D)
      E12 = DOTP(K, Z(1,3), Z(1,3), D)
      E21 = DOTP(K, Z(1,2), Z(1,1), D)
      E22 = DOTP(K, Z(1,2), Z(1,3), D)
C
C     Return the determinant of the matrix:   E11+1    E12
C     (minus 1)                                E21    E22+1
C
      DELTA = (E11+ONE) * (E22+ONE) - E12 * E21 - ONE
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION DOTP(K, X, Y, D)
C
C     [Formerly named DOTPRD, but renamed here to avoid loader table
C      conflicts.  W. Venables.]
C
C     Dot-product scaled by vector D.
C
      implicit none
      INTEGER K
      DOUBLE PRECISION X(K), Y(K), D(K)
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION ZERO
      DATA ZERO/0.D0/
C
      DOTP = ZERO
      DO 10 I = 1, K
   10 DOTP = DOTP + X(I) * Y(I) / D(I)
      RETURN
      END
C
C
      SUBROUTINE SINGM(NP, NRBAR, D, RBAR, TOL, WORK, IFAULT)
C
C     Modified from:
C     ALGORITHM AS274.5  APPL. STATIST. (1992) VOL.41, NO.2
C
C     Checks for singularities, and adjusts orthogonal
C     reductions produced by AS75.1.
C
C     Auxiliary routine called: MODTRI
C
      implicit none
      INTEGER NP, NRBAR, IFAULT
      DOUBLE PRECISION D(NP), RBAR(NRBAR), TOL(NP), WORK(NP)
C
C     Local variables
C
      DOUBLE PRECISION ZERO, TEMP
      INTEGER COL, POS, ROW, NP2, POS2, J
C
      DATA ZERO/0.D0/
C
C     Check input parameters
C
      IFAULT = 0
      IF (NP .LE. 0) IFAULT = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IFAULT = IFAULT + 2
      IF (IFAULT .NE. 0) RETURN
C
      DO 10 COL = 1, NP
   10 WORK(COL) = SQRT(D(COL))
C
      DO 40 COL = 1, NP
C
C     Set elements within RBAR to zero if they are less than TOL(COL) in
C     absolute value after being scaled by the square root of their row
C     multiplier.
C
        TEMP = TOL(COL)
        POS = COL - 1
        DO 20 ROW = 1, COL-1
          IF (ABS(RBAR(POS)) * WORK(ROW) .LT. TEMP) RBAR(POS) = ZERO
          POS = POS + NP - ROW - 1
   20   CONTINUE
C
C     If diagonal element is near zero, set it to zero, and use MODTRI
C     to augment the projections in the lower rows of the factorization.
C
        IF (WORK(COL) .LE. TEMP) THEN
          IFAULT = IFAULT - 1
          IF (COL .LT. NP) THEN
            NP2 = NP - COL
            POS2 = POS + NP - COL + 1
            IF (NP2 .GT. 1) THEN
              CALL MODTRI(NP2, NP2*(NP2-1)/2, D(COL), RBAR(POS+1),
     +             D(COL+1), RBAR(POS2), TOL) 
            ELSE
              CALL MODTRI(1, 0, D(COL), RBAR(POS+1), D(COL+1), RBAR(1),
     +            TOL) 
            END IF
            DO 30 J = POS+1, POS2-1
   30       RBAR(J) = ZERO
          END IF
          D(COL) = ZERO
        END IF
   40 CONTINUE
      RETURN
      END
C
      SUBROUTINE XXTR(NP, NRBAR, D, RBAR, NREQ, TRACE, RINV)

C
C     Calculate the trace of the inverse of X'X (= R'R).
C
      implicit none
      INTEGER NP, NRBAR, NREQ
      DOUBLE PRECISION D(NP), RBAR(NRBAR), TRACE, RINV(*)
C
C     Local variables
C
      INTEGER POS, ROW, COL
      DOUBLE PRECISION ONE, ZERO
      DATA ONE/1.D0/, ZERO/0.D0/
C
C     Get the inverse of R
C
      CALL INV(NP, NRBAR, RBAR, NREQ, RINV)
C
C     Trace = the sum of the diagonal elements of RINV * (1/D) * (RINV)'
C
      TRACE = ZERO
      POS = 1
      DO 20 ROW = 1, NREQ
        TRACE = TRACE + ONE / D(ROW)
        DO 10 COL = ROW+1, NREQ
          TRACE = TRACE + RINV(POS)**2 / D(COL)
          POS = POS + 1
   10   CONTINUE
   20 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE INV(NP, NRBAR, RBAR, NREQ, RINV)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Invert first NREQ rows and columns of Cholesky factorization
C     produced by AS 75.1.
C
      implicit none
      INTEGER NP, NRBAR, NREQ
C     We are using this with double precision, so must change. (gcc gives warnings) #KJETIL
C      REAL RBAR(*), RINV(*)
      DOUBLE PRECISION RBAR(*), RINV(*)
C
C     Local variables.
C
      INTEGER POS, ROW, COL, START, K, POS1, POS2
      DOUBLE PRECISION SUM, ZERO
C
      DATA ZERO/0.0/
C
C     Invert RBAR ignoring row multipliers, from the bottom up.
C
      POS = NREQ * (NREQ-1)/2
      DO 30 ROW = NREQ-1, 1, -1
        START = (ROW-1) * (NP+NP-ROW)/2 + 1
        DO 20 COL = NREQ, ROW+1, -1
          POS1 = START
          POS2 = POS
          SUM = ZERO
          DO 10 K = ROW+1, COL-1
            POS2 = POS2 + NREQ - K
            SUM = SUM - RBAR(POS1) * RINV(POS2)
            POS1 = POS1 + 1
   10     CONTINUE
          RINV(POS) = SUM - RBAR(POS1)
          POS = POS - 1
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END

C-------------------------------------------------------------------------

      SUBROUTINE CLEAR(NP, NRBAR, D, RBAR, IER)
C
C     Modified version of:
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Sets arrays to zero prior to calling INCLUD
C
      implicit none
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(*)
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      DO 10 I = 1, NP
	D(I) = ZERO
   10 CONTINUE
      DO 20 I = 1, NRBAR
   20 RBAR(I) = ZERO
      RETURN
      END
C
C
      SUBROUTINE REGCF(NP, NRBAR, D, RBAR, THETAB, TOL, BETA,
     +     NREQ, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Modified version of AS75.4 to calculate regression coefficients
C     for the first NREQ variables, given an orthogonal reduction from
C     AS75.1.
C
      implicit none 
      INTEGER NP, NRBAR, NREQ, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), TOL(NP),
     +     BETA(NP)
C
C     Local variables
C
      INTEGER I, J, NEXTR
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (NREQ .LT. 1 .OR. NREQ .GT. NP) IER = IER + 4
      IF (IER .NE. 0) RETURN
C
      DO 20 I = NREQ, 1, -1
	IF (SQRT(D(I)) .LT. TOL(I)) THEN
	  BETA(I) = ZERO
	  D(I) = ZERO
	  GO TO 20
	END IF
	BETA(I) = THETAB(I)
	NEXTR = (I-1) * (NP+NP-I)/2 + 1
	DO 10 J = I+1, NREQ
	  BETA(I) = BETA(I) - RBAR(NEXTR) * BETA(J)
	  NEXTR = NEXTR + 1
   10   CONTINUE
   20 CONTINUE
C
      RETURN
      END
C
C
C
C
      REAL FUNCTION RRAND()
C
C     The Wichmann & Hill random number generator
C     Algorithm AS183, Appl. Statist., 31, 188-190, 1982.
C     The cycle length is 6.95E+12.
C     This random number generator is very slow compared with most
C     others, but it is dependable, and the results are reproducible.
C
      implicit none
      INTEGER*4 IX, IY, IZ
      COMMON /RANDC/ IX, IY, IZ
C
      IX = MOD(171*IX, 30269)
      IY = MOD(172*IY, 30307)
      IZ = MOD(170*IZ, 30323)
      RRAND = MOD(FLOAT(IX)/30269. + FLOAT(IY)/30307. +
     1       FLOAT(IZ)/30323. , 1.0)
      RETURN
      END

      SUBROUTINE SETRAN(INIT)
C     
C     Simple initializer for use in S-PLUS routines - WNV April, 1995
C
      implicit none
      INTEGER*4 INIT(3), IX, IY, IZ
      COMMON /RANDC/ IX, IY, IZ
      IX = INIT(1)
      IY = INIT(2)
      IZ = INIT(3)
      RETURN
      END

C      subroutine tesran(n, x)
C      integer n
C      real rand
C      double precision x(*)
C      do i = 1,n
C         x(i) = rand()
C      end do
C      return
C      end



