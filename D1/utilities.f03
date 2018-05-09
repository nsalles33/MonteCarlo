MODULE utilities
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: grid, histgrm, guss, error
CONTAINS

  SUBROUTINE grid(dat)
    IMPLICIT NONE
    REAL,INTENT(INOUT) :: dat(:)
    INTEGER :: num, i
    REAL :: dx


    num = SIZE(dat, dim=1)!(dat)
    dx = 4.0/num
    DO i = 1, num, 1
       dat(i) = -2 + i*dx
    END DO
  END SUBROUTINE grid

  SUBROUTINE histgrm(dat, hist, X, len, prob)
    IMPLICIT NONE
    REAL,INTENT(inout) :: dat(:), hist(:), X(:), prob(:), len
    INTEGER :: n, m, i, j
    REAL :: dx


    n = SIZE(dat, dim=1)!(dat)
    m = SIZE(hist, dim=1)!(dat)
    dx = len/m

    hist(:) = 0
    prob(:) = 0

    DO i = 2, m, 1
       DO j = 1, n, 1
          ! IF ( dat(j) .LE. X(1) ) THEN
          !    hist(1) = hist(1) + 1
          IF ( dat(j) .GT. x(i-1) .AND. dat(j) .LE. x(i) ) THEN
             hist(i) = hist(i) + 1
          END IF
       END DO
    END DO
    DO i = 1, n, 1
       IF ( dat(i) .LE. X(1) .AND. dat(i) .GT. (X(1) - dx)) THEN
          hist(1) = hist(1) + 1
       ENDIF
    END DO
    prob(:) = hist(:)/(n*dx)
  END SUBROUTINE histgrm

  SUBROUTINE guss(x,mu,sigma,g)
    IMPLICIT NONE
    REAL, INTENT(inout) :: mu, sigma, x(:), g(:)
    REAL:: pi, dx
    INTEGER:: n, i
    pi = 3.14159265359
    n = SIZE(x, dim=1)
    ! dx = 4.0/200
    ! x(:) = x(:) + dx/2
    DO i = 1, n, 1
       g(i)=1./SQRT(2*pi*sigma*sigma)*EXP(-(x(i)-mu)**2/(2.0*sigma*sigma))
    END DO
  END SUBROUTINE guss

  SUBROUTINE error(err, prob, g, na, max)
    IMPLICIT NONE
    REAL, INTENT(INOUT) ::err(:), prob(:,:), g(:)
    INTEGER, INTENT(INOUT) :: na, max
    REAL :: er(max)
    INTEGER:: i
    DO i = 1, na, 1
       er(:) = 100*ABS(prob(:,i) - g(:))/g(:)
       err(i) = MAXVAL(er)
    END DO
  END SUBROUTINE error

END MODULE utilities
