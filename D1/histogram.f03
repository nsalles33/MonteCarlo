PROGRAM histo
  USE utilities

  IMPLICIT NONE

  INTEGER, PARAMETER :: sp =4, nmax=200, nsamp=2000000, n = 50
  REAL(kind = sp):: len = 4.0
  REAL (KIND= sp), ALLOCATABLE :: X(:), dat(:), hist(:), prob(:,:), g(:)
  REAL (KIND= sp):: my_sum(n), dx!, diff  ! the real is 4byte long but has been convert to double presicion type with kind =8
  INTEGER (KIND=sp):: i, j, k

  CALL random_SEED()
  ALLOCATE(X(nmax))
  ALLOCATE(g(nmax))
  ALLOCATE(dat(nsamp))
  ALLOCATE(hist(nmax))
  ALLOCATE(prob(nmax,n))

  CALL grid(X)

  OPEN(unit=120,file='out_put',status='unknown')

  DO i = 1, n, 1
     DO j = 1, nsamp, 1
        my_sum(:) = 0.
        DO k = 1, i, 1
           CALL random_NUMBER(dx)
           dx = dx - 0.5
           my_sum(j) = dx
        END DO
        dat(i) = SUM(my_sum)
     END DO
     CALL histgrm(dat,hist,x,len, prob)
  END DO

  DO i = 1, nmax, 1
     WRITE(180,*) x(i), hist(i)
  END DO

  !$OMP PARALLEL DO
  !$OMP END PARALLEL DO

  CALL histgrm(dat,hist,x,len, prob)
  DO i = 1, nmax, 1
     WRITE(120,*) x(i), hist(i)
  END DO


  CLOSE(120)

  DEALLOCATE(X)
  DEALLOCATE(hist)
  DEALLOCATE(prob)
  DEALLOCATE(dat)
  DEALLOCATE(g)

END PROGRAM histo
