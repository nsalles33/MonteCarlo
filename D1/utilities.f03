MODULE utilities
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: grid, histgrm, corr, energy
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

  SUBROUTINE histgrm(dat, hist, X, len)
    IMPLICIT NONE
    REAL,INTENT(inout) :: dat(:), hist(:), X(:), len
    INTEGER :: n, m, i, j
    REAL :: dx


    n = SIZE(dat, dim=1)!(dat)
    m = SIZE(hist, dim=1)!(dat)
    dx = len/m
    hist(:) = 0
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
  END SUBROUTINE histgrm


  SUBROUTINE corr(num, J, T)
    IMPLICIT NONE

    INTEGER, PARAMETER :: sp = 4, dp = 8
    REAL(KIND= dp),INTENT(inout) ::  J(:,:)
    REAL(KIND= dp),INTENT(inout) :: T
    INTEGER(KIND= dp), INTENT(in) :: num
    INTEGER(KIND= dp) :: i1, i2
    REAL (KIND= dp):: x, pi, gamma, inv_pi

    pi  = 3.14159265359/(num*T)
    inv_pi = num/3.14159265359
    gamma = 1.25
    DO i1 = 1, num, 1
       DO i2 = i1+1, num, 1
          x = pi*(i1 - i2)
          x = inv_pi*ABS(SIN(x))
          J(i1,i2) = -x**(-gamma)
       END DO
    END DO

  END SUBROUTINE corr

  SUBROUTINE energy(J, A , N, M, H)
    IMPLICIT NONE

    INTEGER, PARAMETER :: sp = 4, dp = 8
    REAL(KIND= dp),INTENT(inout) ::  J(:,:), A(:,:), H(:)
    INTEGER(KIND= dp), INTENT(in) :: n, m
    INTEGER(KIND= dp) :: i1, i2, i3

    DO i1 = 1, m, 1
       DO i2 = 1, n, 1
          DO i3 = 1, n, 1
             H(i1) =+ J(i3,i2)*A(i1,i3)*A(i1,i2)
          END DO
       END DO
    END DO
    WRITE(*,*) H
  END SUBROUTINE energy

END MODULE utilities
