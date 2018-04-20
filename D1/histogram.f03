PROGRAM readfloat
  USE utilities
  IMPLICIT NONE

  INTEGER, PARAMETER :: sp =4, nmax =100
  REAL(kind = sp), PARAMETER :: l = 4.0
  REAL (KIND= sp), ALLOCATABLE :: X(:)
  ! REAL (KIND= sp):: dx!checksum , my_sum, tolerance, dx  ! the real is 4byte long but has been convert to double presicion type with kind =8
  ! INTEGER (KIND=sp):: n

  ! dx = l/nmax

  ALLOCATE(X(nmax))

  CALL grid(X,l)



  DEALLOCATE(X)

END PROGRAM readfloat
