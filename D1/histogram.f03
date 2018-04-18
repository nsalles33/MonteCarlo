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




  ! READ(5,*) n         !reading the number of integers
  ! WRITE(*,*) "Total number of integers", n
  !
  ! READ(5,*)t          ! so, read them into the array t
  ! READ(5,*) checksum  ! the provided sum of the numbers
  !
  ! my_sum = SUM(t)
  ! x = ABS(my_sum - checksum)
  ! IF (x .LE. tolerance) THEN
  !    WRITE(*,*) "my_sum:",my_sum ,"  =  ","givensum:", checksum
  !    WRITE(*,*) "It's working :)"
  !    WRITE(*,*) "=============================================="
  ! ELSE
  !    WRITE(*,*) "the absolute difference is:", x
  ! ENDIF

  DEALLOCATE(X)

END PROGRAM readfloat
