PROGRAM histo
  USE utilities

  IMPLICIT NONE

  INTEGER, PARAMETER :: sp =4, nmax=200, nsamp=20000, n = 10
  REAL(kind = sp):: len = 4.0
  REAL (KIND= sp), ALLOCATABLE :: X(:), dat(:), hist(:,:), prob(:,:), g(:)
  REAL (KIND= sp):: my_sum(n), err(n), dx, sigma, mu
  INTEGER (KIND=sp):: i, j, k

  sigma = 1.0
  mu = 0.0

  CALL random_SEED()

  ALLOCATE(X(nmax))
  ALLOCATE(g(nmax))
  ALLOCATE(dat(nsamp))
  ALLOCATE(hist(nmax,n))
  ALLOCATE(prob(nmax,n))

  CALL grid(X)

  OPEN(unit=120,file='_put',status='unknown')
  OPEN(unit=180,file='hist_out',status='unknown')

  DO i = 1, n, 1
     DO j = 1, nsamp, 1
        my_sum(:) = 0.
        DO k = 1, i, 1
           CALL random_NUMBER(dx)
           dx = dx - 0.5
           my_sum(k) = dx
        END DO
        dat(j) = SUM(my_sum)
        dat(j) = SQRT(12.0/i)*dat(j)
     END DO
     CALL histgrm(dat,hist(:,i),x,len, prob(:,i))
  END DO

  CALL guss(x,mu,sigma,g)

  DO i = 1, nmax, 1
     WRITE(180,*) x(i), hist(i,6), g(i),prob(i,4), prob(i,6), prob(i,8), prob(i,10)
  END DO
  i = n
  j = nmax
  CALL error(err, prob, g, i, j)
  DO i = 1, n, 1
     WRITE(120,*) i, err(i)
  END DO
  ! CALL histgrm(dat,hist,x,len, prob)
  ! DO i = 1, nmax, 1
  !    WRITE(120,*) x(i), hist(i)
  ! END DO


  CLOSE(120)
  CLOSE(180)

  DEALLOCATE(X)
  DEALLOCATE(hist)
  DEALLOCATE(prob)
  DEALLOCATE(dat)
  DEALLOCATE(g)

END PROGRAM histo
