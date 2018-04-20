PROGRAM histo
  USE utilities
  IMPLICIT NONE

  INTEGER, PARAMETER :: sp =4, nmax=200, nsamp=2000000
  REAL(kind = sp):: len = 4.0
  REAL (KIND= sp), ALLOCATABLE :: X(:), dat(:), hist(:)
  REAL (KIND= sp):: checksum(4), my_sum(6), dx!, diff  ! the real is 4byte long but has been convert to double presicion type with kind =8
  INTEGER (KIND=sp):: i, j

  CALL random_SEED()
  ALLOCATE(X(nmax))
  ALLOCATE(dat(nsamp))
  ALLOCATE(hist(nmax))

  CALL grid(X)
  OPEN(unit=120,file='4_sum',status='unknown')
  OPEN(unit=180,file='6_sum',status='unknown')
  !$OMP PARALLEL DO
  DO i = 1, nsamp, 1
     DO j = 1, 6, 1
        CALL random_NUMBER(dx)
        dx = dx - 0.5
        my_sum(j) = dx
     END DO
     dat(i) = SUM(my_sum)
  END DO
  !$OMP END PARALLEL DO


  CALL histgrm(dat,hist,x,len)
  DO i = 1, nmax, 1
     WRITE(180,*) x(i), hist(i)
  END DO

  !$OMP PARALLEL DO
  DO i = 1, nsamp, 1
     DO j = 1, 4, 1
        CALL random_NUMBER(dx)
        dx = dx - 0.5
        checksum(j) = dx
     END DO
     dat(i) = SUM(checksum)
  END DO
  !$OMP END PARALLEL DO

  CALL histgrm(dat,hist,x,len)
  DO i = 1, nmax, 1
     WRITE(120,*) x(i), hist(i)
  END DO


  CLOSE(120)
  CLOSE(180)


  DEALLOCATE(X)
  DEALLOCATE(hist)
  DEALLOCATE(dat)

END PROGRAM histo
