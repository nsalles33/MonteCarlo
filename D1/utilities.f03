MODULE utilities
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: grid, corr
CONTAINS

  SUBROUTINE grid(dat, l)
    IMPLICIT NONE
    REAL,DIMENSION(:),INTENT(inout) :: dat
    REAL, INTENT(in) :: l
    INTEGER :: num, i
    REAL :: dx


    num = SIZE(dat, dim=1)!(dat)
    dx = l/num
    DO i = 1, num, 1
       dat(i) = i*dx
    END DO
  END SUBROUTINE grid

  SUBROUTINE corr(num, J, T)
    IMPLICIT NONE

    INTEGER, PARAMETER :: sp = 4, dp = 8
    REAL(KIND= dp),INTENT(inout) ::  J(:,:), T
    INTEGER, INTENT(in) :: num
    INTEGER :: i1, i2
    REAL :: x, pi, gamma, inv_pi

    pi  = 3.14159265359/(num*T)
    inv_pi = num/3.14159265359
    gamma = 1.25

    DO i1 = 1, num, 1
       DO i2 = i1+1, num, 1
          x = pi*(i1 - i2)
          x = inv_pi*ABS(SIN(x))
          J(i1,i2) = x**(-gamma)
       END DO
    END DO
  END SUBROUTINE corr

END MODULE utilities
