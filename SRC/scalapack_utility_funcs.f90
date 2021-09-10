       SUBROUTINE pdlawrite( FILNAM, M, N, A, IA, JA, DESCA, IRWRIT, &
                             ICWRIT, WORK )
 !
 !  -- ScaLAPACK tools routine (version 1.8) --
 !     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
 !     and University of California, Berkeley.
 !
 !     written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
 !     adapted by Julie Langou, April 2007 (julie@cs.utk.edu)
 !
 !     .. Scalar Arguments ..
       INTEGER            IA, ICWRIT, IRWRIT, JA, M, N
 !     ..
 !     .. Array Arguments ..
       CHARACTER*(*)      FILNAM
       INTEGER            DESCA( * )
       DOUBLE PRECISION   A( * ), WORK( * )
 !     ..
 !
 !  Purpose
 !  =======
 !
 !  PDLAWRITE writes to a file named FILNAMa distributed matrix sub( A )
 !  denoting A(IA:IA+M-1,JA:JA+N-1). The local pieces are sent to and
 !  written by the process of coordinates (IRWWRITE, ICWRIT).
 !
 !  WORK must be of size >= MB_ = DESCA( MB_ ).
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       INTEGER            NOUT
       parameter( nout = 13 )
       INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,  &
                         lld_, mb_, m_, nb_, n_, rsrc_
       parameter( block_cyclic_2d = 1, dlen_ = 9, dt_ = 1, &
                  ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6, &
                  rsrc_ = 7, csrc_ = 8, lld_ = 9 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL, &
                          icurrow, ii, iia, in, j, jb, jj, jja, jn, k, &
                          lda, mycol, myrow, npcol, nprow
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           blacs_barrier, blacs_gridinfo, infog2l, &
                          dgerv2d, dgesd2d
 !     ..
 !     .. External Functions ..
       INTEGER            ICEIL
       EXTERNAL           iceil
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          min, mod
 !     ..
 !     .. Executable Statements ..
 !
 !     Get grid parameters
 !
       ictxt = desca( ctxt_ )
       CALL blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
 !
       IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
          OPEN( nout, file=filnam, status='UNKNOWN' )
          WRITE( nout, fmt = * ) m, n
       END IF
 !
       CALL infog2l( ia, ja, desca, nprow, npcol, myrow, mycol, &
                     iia, jja, iarow, iacol )
       icurrow = iarow
       icurcol = iacol
       ii = iia
       jj = jja
       lda = desca( lld_ )
 !
 !     Handle the first block of column separately
 !
       jn = min( iceil( ja, desca( nb_ ) ) * desca( nb_ ), ja+n-1 )
       jb = jn-ja+1
       DO 60 h = 0, jb-1
          in = min( iceil( ia, desca( mb_ ) ) * desca( mb_ ), ia+m-1 )
          ib = in-ia+1
          IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
             IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                DO 10 k = 0, ib-1
                   WRITE( nout, fmt = 9999 ) a( ii+k+(jj+h-1)*lda )
    10          CONTINUE
             END IF
          ELSE
             IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                CALL dgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda ), lda, &
                              irwrit, icwrit )
             ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                CALL dgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                              icurrow, icurcol )
                DO 20 k = 1, ib
                   WRITE( nout, fmt = 9999 ) work( k )
    20          CONTINUE
             END IF
          END IF
          IF( myrow.EQ.icurrow ) THEN
             ii = ii + ib
          END IF
          icurrow = mod( icurrow+1, nprow )
          CALL blacs_barrier( ictxt, 'All' )
 !
 !        Loop over remaining block of rows
 !
          DO 50 i = in+1, ia+m-1, desca( mb_ )
             ib = min( desca( mb_ ), ia+m-i )
             IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
                IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                   DO 30 k = 0, ib-1
                      WRITE( nout, fmt = 9999 ) a( ii+k+(jj+h-1)*lda )
    30             CONTINUE
                END IF
             ELSE
                IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                   CALL dgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda ), &
                                 lda, irwrit, icwrit )
                ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                   CALL dgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                                 icurrow, icurcol )
                   DO 40 k = 1, ib
                      WRITE( nout, fmt = 9999 ) work( k )
    40             CONTINUE
                END IF
             END IF
             IF( myrow.EQ.icurrow ) THEN
                ii = ii + ib
             END IF
             icurrow = mod( icurrow+1, nprow )
             CALL blacs_barrier( ictxt, 'All' )
    50    CONTINUE
 !
         ii = iia
         icurrow = iarow
    60 CONTINUE
 !
       IF( mycol.EQ.icurcol ) THEN
         jj = jj + jb
       END IF
       icurcol = mod( icurcol+1, npcol )
       CALL blacs_barrier( ictxt, 'All' )
 !
 !     Loop over remaining column blocks
 !
       DO 130 j = jn+1, ja+n-1, desca( nb_ )
          jb = min(  desca( nb_ ), ja+n-j )
          DO 120 h = 0, jb-1
             in = min( iceil( ia, desca( mb_ ) ) * desca( mb_ ), ia+m-1 )
             ib = in-ia+1
             IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
                IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                   DO 70 k = 0, ib-1
                      WRITE( nout, fmt = 9999 ) a( ii+k+(jj+h-1)*lda )
    70             CONTINUE
                END IF
             ELSE
                IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                   CALL dgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda ), &
                                lda, irwrit, icwrit )
                ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                   CALL dgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                                icurrow, icurcol )
                   DO 80 k = 1, ib
                      WRITE( nout, fmt = 9999 ) work( k )
    80             CONTINUE
                END IF
             END IF
             IF( myrow.EQ.icurrow ) THEN
                ii = ii + ib
             END IF
             icurrow = mod( icurrow+1, nprow )
             CALL blacs_barrier( ictxt, 'All' )
 !
 !           Loop over remaining block of rows
 !
             DO 110 i = in+1, ia+m-1, desca( mb_ )
                ib = min( desca( mb_ ), ia+m-i )
                IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
                   IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                      DO 90 k = 0, ib-1
                         WRITE( nout, fmt = 9999 ) a( ii+k+(jj+h-1)*lda)
    90                CONTINUE
                   END IF
                ELSE
                   IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                      CALL dgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda ),&
                                   lda, irwrit, icwrit )
                    ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                      CALL dgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                                    icurrow, icurcol )
                      DO 100 k = 1, ib
                         WRITE( nout, fmt = 9999 ) work( k )
   100                CONTINUE
                   END IF
                END IF
                IF( myrow.EQ.icurrow ) THEN
                   ii = ii + ib
                END IF
                icurrow = mod( icurrow+1, nprow )
                CALL blacs_barrier( ictxt, 'All' )
   110       CONTINUE
 !
             ii = iia
             icurrow = iarow
   120    CONTINUE
 !
          IF( mycol.EQ.icurcol ) THEN
             jj = jj + jb
          END IF
          icurcol = mod( icurcol+1, npcol )
          CALL blacs_barrier( ictxt, 'All' )
 !
   130 CONTINUE
 !
       IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
          CLOSE( nout )
       END IF
 !
  9999 FORMAT( d30.18 )
 !
       RETURN
 !
 !     End of PDLAWRITE
 !
       END









       SUBROUTINE pzlawrite( FILNAM, M, N, A, IA, JA, DESCA, IRWRIT, &
                             ICWRIT, WORK )
 !
 !  -- ScaLAPACK tools routine (version 1.8) --
 !     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
 !     and University of California, Berkeley.
 !
 !     written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
 !     adapted by Julie Langou, April 2007 (julie@cs.utk.edu)
 !
 !     .. Scalar Arguments ..
       INTEGER            IA, ICWRIT, IRWRIT, JA, M, N
 !     ..
 !     .. Array Arguments ..
       CHARACTER*(*)      FILNAM
       INTEGER            DESCA( * )
       COMPLEX*16         A( * ), WORK( * )
 !     ..
 !
 !  Purpose
 !  =======
 !
 !  PZLAWRITE writes to a file named FILNAMa distributed matrix sub( A )
 !  denoting A(IA:IA+M-1,JA:JA+N-1). The local pieces are sent to and
 !  written by the process of coordinates (IRWWRITE, ICWRIT).
 !
 !  WORK must be of size >= MB_ = DESCA( MB_ ).
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       INTEGER            NOUT
       parameter( nout = 13 )
       INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, &
                          lld_, mb_, m_, nb_, n_, rsrc_
       parameter( block_cyclic_2d = 1, dlen_ = 9, dt_ = 1,    &
                  ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ =6, &
                            rsrc_ = 7, csrc_ = 8, lld_ = 9 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL, &
                          icurrow, ii, iia, in, j, jb, jj, jja, jn, k, &
                          lda, mycol, myrow, npcol, nprow
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           blacs_barrier, blacs_gridinfo, infog2l, &
                          zgerv2d, zgesd2d
 !     ..
 !     .. External Functions ..
       INTEGER            ICEIL
       EXTERNAL           iceil
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          dble, dimag, min, mod
 !     ..
 !     .. Executable Statements ..
 !
 !     Get grid parameters
 !
       ictxt = desca( ctxt_ )
       CALL blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
 !
       IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
          OPEN( nout, file=filnam, status='UNKNOWN' )
          WRITE( nout, fmt = * ) m, n
       END IF
 !
       CALL infog2l( ia, ja, desca, nprow, npcol, myrow, mycol, &
                     iia, jja, iarow, iacol )
       icurrow = iarow
       icurcol = iacol
       ii = iia
       jj = jja
       lda = desca( lld_ )
 !
 !     Handle the first block of column separately
 !
       jn = min( iceil( ja, desca( nb_ ) ) * desca( nb_ ), ja+n-1 )
       jb = jn-ja+1
       DO 60 h = 0, jb-1
          in = min( iceil( ia, desca( mb_ ) ) * desca( mb_ ), ia+m-1 )
          ib = in-ia+1
          IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
             IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                DO 10 k = 0, ib-1
                   WRITE( nout, fmt = 9999 ) a( ii+k+(jj+h-1)*lda )
    10          CONTINUE
             END IF
          ELSE
             IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                CALL zgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda ), lda, &
                              irwrit, icwrit )
             ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                CALL zgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                              icurrow, icurcol )
                DO 20 k = 1, ib
                   WRITE( nout, fmt = 9999 ) dble(work( k )), &
                    dimag(work( k ))
    20          CONTINUE
             END IF
          END IF
          IF( myrow.EQ.icurrow ) THEN
             ii = ii + ib
          END IF
          icurrow = mod( icurrow+1, nprow )
          CALL blacs_barrier( ictxt, 'All' )
 !
 !        Loop over remaining block of rows
 !
          DO 50 i = in+1, ia+m-1, desca( mb_ )
             ib = min( desca( mb_ ), ia+m-i )
             IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
                IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                   DO 30 k = 0, ib-1
                      WRITE( nout, fmt = 9999 ) &
                       dble(a( ii+k+(jj+h-1)*lda )), &
                       dimag(a( ii+k+(jj+h-1)*lda ))      
    30             CONTINUE
                END IF
             ELSE
                IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                   CALL zgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda ), &
                                 lda, irwrit, icwrit )
                ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN 
                   CALL zgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                                 icurrow, icurcol )
                   DO 40 k = 1, ib
                      WRITE( nout, fmt = 9999 ) dble(work( k )), &
                                                dimag(work( k ))
    40             CONTINUE
                END IF
             END IF
             IF( myrow.EQ.icurrow ) THEN
                ii = ii + ib
             END IF
             icurrow = mod( icurrow+1, nprow )
             CALL blacs_barrier( ictxt, 'All' )
    50    CONTINUE
 !
         ii = iia
         icurrow = iarow
    60 CONTINUE
 !
       IF( mycol.EQ.icurcol ) THEN
          jj = jj + jb
       END IF
       icurcol = mod( icurcol+1, npcol )
       CALL blacs_barrier( ictxt, 'All' )
 !
 !     Loop over remaining column blocks
 !
       DO 130 j = jn+1, ja+n-1, desca( nb_ )
          jb = min(  desca( nb_ ), ja+n-j )
          DO 120 h = 0, jb-1
             in = min( iceil( ia, desca( mb_ ) ) * desca( mb_ ), ia+m-1)
             ib = in-ia+1
             IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
                IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                   DO 70 k = 0, ib-1
                      WRITE( nout, fmt = 9999 ) &
                     dble(a( ii+k+(jj+h-1)*lda )), &
                     dimag(a( ii+k+(jj+h-1)*lda ))
    70             CONTINUE
                END IF
             ELSE
                IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                   CALL zgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda ), &
                                 lda, irwrit, icwrit )
                ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                   CALL zgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                                 icurrow, icurcol )
                   DO 80 k = 1, ib
                      WRITE( nout, fmt = 9999 ) dble(work( k )), &
                                                dimag(work( k))
    80             CONTINUE
                END IF
             END IF
             IF( myrow.EQ.icurrow ) THEN
                ii = ii + ib
             END IF
             icurrow = mod( icurrow+1, nprow )
             CALL blacs_barrier( ictxt, 'All' )
 !
 !           Loop over remaining block of rows
 !
             DO 110 i = in+1, ia+m-1, desca( mb_ )
                ib = min( desca( mb_ ), ia+m-i )
                IF( icurrow.EQ.irwrit .AND. icurcol.EQ.icwrit ) THEN
                   IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                      DO 90 k = 0, ib-1
                         WRITE( nout, fmt = 9999 )  &
                          dble(a( ii+k+(jj+h-1)*lda )), &
                          dimag(a( ii+k+(jj+h-1)*lda ))
    90                CONTINUE
                   END IF
                ELSE
                   IF( myrow.EQ.icurrow .AND. mycol.EQ.icurcol ) THEN
                      CALL zgesd2d( ictxt, ib, 1, a( ii+(jj+h-1)*lda), &
                                    lda, irwrit, icwrit )
                    ELSE IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
                      CALL zgerv2d( ictxt, ib, 1, work, desca( mb_ ), &
                                    icurrow, icurcol )
                      DO 100 k = 1, ib
                         WRITE( nout, fmt = 9999 ) dble(work( k )), &
                                                   dimag(work( k ))
   100                CONTINUE
                   END IF
                END IF
                IF( myrow.EQ.icurrow ) THEN
                   ii = ii + ib
                END IF
                icurrow = mod( icurrow+1, nprow )
                CALL blacs_barrier( ictxt, 'All' )
   110       CONTINUE
 !
             ii = iia
             icurrow = iarow
   120    CONTINUE
 !
          IF( mycol.EQ.icurcol ) THEN
             jj = jj + jb
          END IF
          icurcol = mod( icurcol+1, npcol )
          CALL blacs_barrier( ictxt, 'All' )
 !
   130 CONTINUE
 !
       IF( myrow.EQ.irwrit .AND. mycol.EQ.icwrit ) THEN
          CLOSE( nout )
       END IF
 !
  9999 FORMAT( e15.8,e15.8 )
 !
       RETURN
 !
 !     End of PZLAWRITE
 !
       END       
