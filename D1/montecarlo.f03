PROGRAM monte_carlo
  USE utilities
  IMPLICIT NONE

  INTEGER, PARAMETER :: sp = 4, dp = 8
  ! REAL(kind = sp), PARAMETER :: l = 4.0
  REAL (KIND= dp), ALLOCATABLE :: H(:,:), J(:,:), A(:,:)
  REAL (KIND= dp):: dx, x!checksum , my_sum, tolerance, dx
  REAL (KIND= dp):: time1, time2, time3, T
  INTEGER (KIND=sp):: i1, i2, k, N, M, nsweep

  PRINT *,'Sweeps N M Temp'
  READ *, nsweep, N, M, T
  !
  CALL CPU_TIME(time1)
  ! i = NINT(time1)
  ! CALL random_SEED(i)
  !

  ALLOCATE(H(n,n))
  ALLOCATE(J(m,n))
  ALLOCATE(A(n,m))

  OPEN(unit=130,file='A',status='unknown')
  DO i1=1,N,1
     DO i2 = 1, m, 1
        CALL random_NUMBER(x)
        IF ( (NINT(X)) .EQ. 1 ) THEN
           A(i2,i1)=-1
        END IF
     END DO
  END DO
  CLOSE(130)

  CALL corr(N,J, T)
  DO i1 = 1, nsweep, 1
     DO k = 1, N, 1

     END DO
  END DO



  DEALLOCATE(A)
  DEALLOCATE(H)
  DEALLOCATE(J)

  CALL CPU_TIME(TIME2)
  WRITE(*,*) TIME2-TIME1

END PROGRAM monte_carlo
